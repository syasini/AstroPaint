# @purpose: Get the necessary DataSets

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"

import requests
from tqdm import tqdm
import os



#########################################################
#           Get the necessary datasets
#########################################################


def _get_dataset(file_location, ds_url, dataset_name=None):
    """
        @purpose: Connect and Download the necessary files
        @input:
            file_location: str: file location to download dataset
            ds_url: str: link to the datasets
            dataset_name: str: name of the dataset

    """
    response = requests.get(ds_url, stream=True)

    if dataset_name is None:
        dataset_name = file_location

    with tqdm.wrapattr(open(file_location, "wb"), "write", miniters=1,
                       total=int(response.headers.get('content-length', 0)),
                       desc=dataset_name) as fout:
        for chunk in response.iter_content(chunk_size=4096):
            fout.write(chunk)


def download(dataset_name:str='All'):
    """
        @purpose: Download datasets mannualy if not exists into data folder.
        @Input: Input a string:
                "All" : Downloads all the relevant datasets.
                "WebSky_lite" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/WebSky_lite.csv', 
                "WebSky_2x2" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/WebSky_2x2.csv',
                "Sehgal" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/Sehgal.csv',
                "MICE" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/MICE.csv'

    
    """
    datasets_url_dict = {
        "WebSky_lite" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/WebSky_lite.csv', 
        "WebSky_2x2" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/WebSky_2x2.csv',
        "Sehgal" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/Sehgal.csv',
        "MICE" : 'https://media.githubusercontent.com/media/syasini/AstroPaint/master/astropaint/data/MICE.csv',
    }

    # make a list of all the datasets to be downloaded
    download_list = []

    if dataset_name == 'All':
        download_list = [ds_name for ds_name in datasets_url_dict.keys()]
    elif type(dataset_name) == list:
        download_list = dataset_name
    elif type(dataset_name) == str:
        download_list.append(dataset_name)
    # TODO: add an else statement here to raise an error 

    for dataset in download_list:
        file_location = os.path.join(os.path.dirname(__file__), os.path.pardir,'data',f'{dataset}.csv')
        _get_dataset(file_location, datasets_url_dict[dataset], dataset)
    #TODO: add overwrite=True arg to the function 

    print('Done!')

