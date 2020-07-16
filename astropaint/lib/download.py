# @purpose: Get the necessary DataSets

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"

import requests
from tqdm.auto import tqdm
import os



#########################################################
#           Get the necessary datasets
#########################################################


def _get_dataset(ds_name, ds_url):
    """
        @purpose: Connect and Download the necessary files
        @input:
            ds_name: str: DataSet Names
            ds_url: str: link to the datasets

    """
    response = requests.get(ds_url, stream=True)

    with tqdm.wrapattr(open(f'{ds_name}.csv', "w"), "write", miniters=1,
                       total=int(response.headers.get('content-length', 0)),
                       desc=ds_name) as fout:
        for chunk in response.iter_content(chunk_size=4096):
            fout.write(chunk)


def download_dataset(dataset_name:str='All'):
    """
        @purpose: Download datasets mannualy if not exists into data folder.
        @Input: Input a string:
                "All" : Downloads all the relevant datasets.
                "WebSky_lite" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/WebSky_lite.csv', 
                "WebSky_2x2" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/WebSky_2x2111.csv',
                "Sehgal" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/Sehgal.csv',
                "MICE" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/MICE.csv'

    
    """
    dict_datasets = {
        "WebSky_lite" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/WebSky_lite.csv', 
        "WebSky_2x2" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/WebSky_2x2111.csv',
        "Sehgal" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/Sehgal.csv',
        "MICE" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/MICE.csv'
    }
    
    if dataset_name == 'All':
        for ds_name, ds_url in dict_datasets.items():
            file_location = os.path.join(os.path.dirname(__file__), os.path.pardir,'data',f'{ds_name}.csv')
            if ~ os.path.exists(file_location):
                _get_dataset(ds_name,ds_url)
            else:
                print(f'{ds_name} already exists!')
    else:
        file_location = os.path.join(os.path.dirname(__file__), os.path.pardir,'data',f'{dataset_name}.csv')
        if ~ os.path.exists(file_location):
            _get_dataset(ds_name,ds_url, file_location)
        else:
            print(f'{dataset_name} already exists!')
    print('END!')

