import re
import gzip
import shutil
import requests
import json
import numpy as np
import pandas as pd
import os
import collections

#main gdc api querry function
def get_projects_info(project_names):
    '''
    a method for retrieveing data about cases and their related samples for specified gdc projects.
    The method looks for which of the projects specified actually exist in the GDC and then for each
    of the project retrieves all the cases with data about them like demographic and data about
    samples related to each case.
    
    Input:
    List of project names for which data should be retrieved
    
    Output:
    dict:
        data dict: data in the form of a dictionary with the structure: {project_name: {case_id: {data about the case}}}
        image to sample: mapping from image file name to sample id for matching with labels
        case to images: mapping from case id to list of associated images for dataset creation
        labels: label dataframe
        mutational signatures: mutatiuonal signatures dataframe
        hugo symbols: hugo symbols dataframe
    '''
    
    #check if project_names is a list of strings
    if not(isinstance(project_names,list) and all(isinstance(i,str) for i in project_names)):
        raise TypeError("project_names expects a list of strings")
    
    #define api endpoints
    cases_endpt = 'https://api.gdc.cancer.gov/cases'
    projects_endpt = 'https://api.gdc.cancer.gov/projects'
    
    
    #check which of the specified project names are in gdc
    print("looking for project names in GDC database.")
    
    filters = {
        "op":"in",
        "content":{
            "field":"project_id",
            "value":project_names
        }
    }
    params = {
        "filters": json.dumps(filters),
        "fields":"project_id",
        "format":"json",
        "pretty":"true",
    }
    
    response = requests.get(projects_endpt, params = params)
    response = json.loads(response.content.decode("utf-8"))
    print(80*'-')
    
    if not(response['warnings']=={}):
        print("warnings:")
        print(response['warnings'])
        
    found_projects = list(map(lambda x: x['project_id'],response['data']['hits']))
    print("Projects found in GDC:",found_projects)
    out_samples = []
    projects_data={}
    image_to_sample={}
    case_to_images={}
    out_hugos = []
    
    #mutational signature file for the entire tcga dataset
    signatures = pd.read_csv(os.path.join("manifest", 'TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv'))
    signatures['case_id'] = signatures['Sample Names'].apply(lambda x : '-'.join(x.split('-')[:3]))
    signatures = signatures.iloc[:,3:]
    #search for cases for each project
    print(80*'-')
    for project in found_projects:
        
        print("looking for cases for project:",project)
        
        filters = {
            "op":"=",
            "content":{
                "field":"project.project_id",
                "value":project
            }
        }
        
        fields = ['case_id','project.project_id',"submitter_id","files.file_id","files.file_name"]
        fields = ','.join(fields)
        
        params = {
            "filters":json.dumps(filters),
            "fields" :fields,
            "format" :"json",
            "pretty" :"true",
            "size"   : "1100",
            "expand" : "demographic,samples,files,diagnoses"
            }

        response = requests.get(cases_endpt, params = params)
        cases = json.loads(response.content.decode("utf-8"))
        print("retrieved",cases['data']['pagination']['count'],"cases in the",project,"project")
                
        
        #add cases to dictionary
        found_cases = cases['data']['hits']
        out_cases = {}
        for case in found_cases:
            
            #section for creating a dataframe
            for sample in case['samples']:
                sample_dict = {}
                sample_dict['case_barcode'] = case['submitter_id']
                sample_dict['project'] = case['project']['project_id']
                sample_dict['case_id'] = case['case_id']

                #add demographic data for sample
                demographic_dict = case.get('demographic',{})
                for key in demographic_dict:
                    sample_dict["demographic."+key] = demographic_dict[key]

                #add diagnose data for sample
                diagnoses_dict = case.get('diagnoses',{})
                #the diagnostic value is a list with either length 0 or 1 if empty list then change it to empty
                #dict for convenience
                if len(diagnoses_dict)>0:
                    diagnoses_dict = diagnoses_dict[0]
                else:
                    diagnoses_dict = {}
                for key in diagnoses_dict:
                    sample_dict["diagnose."+key] = diagnoses_dict[key]
                sample_dict['sample.barcode']=sample['submitter_id']
                #print(sample['submitter_id'])
                for key in sample:
                    if key != 'submitter_id':
                        sample_dict['sample.'+key]=sample[key]
                #print(sample_dict['sample.barcode'])
                out_samples.append(sample_dict)
                
                
            #dict section
            demographic = case.get('demographic',{})
            samples = case.get('samples',[])
            mutational_signatures = signatures[signatures['case_id']==case['submitter_id']]
            
            if mutational_signatures.shape[0] != 0:
                mutational_signatures = mutational_signatures.iloc[0,:-1].to_list()
            else:
                mutational_signatures = []

            out_cases[case['submitter_id']]={"demographic":demographic,"samples":samples,
                                             "case_id":case['case_id'],"hugo_symbols":[],
                                             "mutational_signature":mutational_signatures,
                                             }
            
            #add slide image file names to cases for mapping from cases to images
            for sample in out_cases[case['submitter_id']]['samples']:
                #sample_mut_sig = 
                sample['image_files']=[]
            
            images_for_case = []
            for file in case['files']:
                if file['data_format']=='SVS':
                    sample_id = '-'.join(file['submitter_id'].split('-')[:4])
                    images_for_case.append(file['file_name'])
                    image_to_sample[file['file_name']] = sample_id
                    for sample in out_cases[case['submitter_id']]['samples']:
                        if  sample['submitter_id']==sample_id:
                            sample['image_files'].append(file['file_name'])
            #add list of image file names to mapping for that case             
            if images_for_case != []:
                case_to_images[case['submitter_id']]=images_for_case
            
            #delete samples with no image files
            for i,sample in enumerate(out_cases[case['submitter_id']]['samples']):
                if sample['image_files']==[]:
                    del out_cases[case['submitter_id']]['samples'][i]
            
            #delete cases with no samples
            if out_cases[case['submitter_id']]['samples'] == []:
                del out_cases[case['submitter_id']]
        
        #download maf file to add hugo symbols to each case
        print("downloading maf file for project",project)
        zip_file,maf_file = download_maf_for_proj(project)
        
        #if maf file for the given project was found and extracted then add it to 
        if maf_file != None:
            print("downloaded file:",zip_file)
            print("file extracted to:",maf_file)

            maf = pd.read_csv(maf_file,sep='\t',header=5)
            temp = maf.loc[:,['Hugo_Symbol','Tumor_Sample_Barcode']]
            temp['case_id'] = temp['Tumor_Sample_Barcode'].apply(lambda x : '-'.join(x.split('-')[:3]))

            #add hugo symbols to each case
            for i,r in temp.iterrows():
                out_cases[r['case_id']]["hugo_symbols"].append(r["Hugo_Symbol"])
            #add cases to hugo symbol dataframe
            for case in out_cases:
                symbols = out_cases[case]['hugo_symbols']
                symbols = dict(collections.Counter(symbols))
                symbols['case_barcode'] = case
                out_hugos.append(symbols)
        
        #add the cases doctionary for given project to the projects dictionary
        projects_data[project]=out_cases
        print(80*'-')
    
    samples = pd.DataFrame.from_records(out_samples)
    samples = samples.fillna(np.nan)
    print("done")   
    return {"data dict":projects_data,"image to sample":image_to_sample,
            "case to images":case_to_images,"labels":samples, "mutational signatures" : signatures,
            "hugo symbols" : pd.DataFrame(out_hugos)}

#file download functionality
def download_extract(file_id,project_name):
    download_endpt = "https://api.gdc.cancer.gov/data/{}".format(file_id)
    response = requests.get(download_endpt, headers = {"Content-Type": "application/json"})
    response_head_cd = response.headers["Content-Disposition"]
    file_name = os.path.join("manifest", re.findall("filename=(.+)", response_head_cd)[0])
    out_name = os.path.join("manifest", project_name+'-maf.tsv')
    
    with open(file_name, "wb") as output_file:
        output_file.write(response.content)
    
    with gzip.open(file_name, 'rb') as f_in:
        with open(out_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return file_name,out_name

def download_maf_for_proj(project_name):
    files_endpt = "https://api.gdc.cancer.gov/files"

    workflow = {"op":"=",
                    "content":{
                        "field":"analysis.workflow_type",
                        "value":"MuSE Variant Aggregation and Masking"
                        }
                    }

    category =  {"op":"=",
                    "content":{
                        "field":"data_category",
                        "value":"Simple Nucleotide Variation"
                        }
                    }

    dat_type = {"op":"=",
                    "content":{
                        "field":"data_type",
                        "value":"Masked Somatic Mutation"
                        }
                    }
    proj = {"op":"=",
                "content":{
                    "field":"cases.project.project_id",
                    "value":project_name
                }
            }
    filters = {
        "op":"and",
        "content":[workflow,dat_type,category,proj]
    }
    params = {
        "filters":json.dumps(filters),
        "fields" : "data_type,data_category,file_id,cases.project.project_id",
        "format" :"json",
        "pretty" :"true",
        "size"   : "5",
        "expand" : "analysis",
        }

    response = requests.get(files_endpt, params = params)
    file_data = json.loads(response.content.decode("utf-8"))
    
    #if no maf file found do nothing
    if(file_data['data']['pagination']['total']==0):
        print("no maf file found for project",project_name)
        return None,None
    
    file_id = file_data['data']['hits'][0]['file_id']
    return download_extract(file_id,project_name)

def download_image(file_name,path=""):
    file_path = os.path.join(path,file_name)
    #check if image already exists
    if not os.path.exists(file_path):
        if not os.path.exists(path):
            os.makedirs(path)
        files_endpt = "https://api.gdc.cancer.gov/files"

        context = {"op":"=",
                    "content":{
                        "field":"file_name",
                        "value":file_name
                        }
                    }
        params = {
        "filters":json.dumps(context),
        "fields" : "file_id",
        "format" :"json",
        "pretty" :"true",
        "size"   : "5",
        }

        response = requests.get(files_endpt, params = params)
        file_data = json.loads(response.content.decode("utf-8"))
        file_id = file_data['data']['hits'][0]['file_id']

        data_endpt = "https://api.gdc.cancer.gov/data/{}".format(file_id)
        print("downloading image {} to path {}".format(file_name,file_path))
        response = requests.get(data_endpt, headers = {"Content-Type": "application/json"})

        with open(file_path, "wb") as output_file:
            output_file.write(response.content)
    else:
        print("{} already exists, not downloading anything".format(file_path))

def get_image_labels(case_set, data):
    set_images = list(map( lambda x : data['case to images'][x], case_set))
    #flatten
    set_images = [img for case in set_images for img in case]
    
    #get sample ids to match with labels
    sample_ids = list(map( lambda x : data['image to sample'][x], set_images))
    set_labels = pd.DataFrame()
    for sample in sample_ids:
        set_labels = set_labels.append(data['labels'][data['labels']['sample.barcode']==sample])

    return set_labels, set_images