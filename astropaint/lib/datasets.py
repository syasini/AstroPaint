# @purpose: Get the necessary DataSets

__author__ = "Siavash Yasini, Shobeir K. S. Mazinani"
__email__ = "yasini@usc.edu , shobeir.mazinani@gmail.com"
__credits__ = ['Siavash Yasini', 'Shobeir K. S. Mazinani']


import requests
from tqdm import tqdm
import os



#########################################################
#           Get the necessary datasets
#########################################################


def _get_dataset(file_location, ds_url, dataset_name=None):
    """
        @purpose: Connect and Download the necessary files and show progressbar
        @input:
            file_location: str: file location to download dataset
            ds_url: str: link to the datasets
            dataset_name: str: name of the dataset

    """
    
    response = requests.get(ds_url, stream=True)

    if dataset_name is None:
        dataset_name = file_location
    
    # Download and write the file in chuncks and show progressbar
    with tqdm.wrapattr(open(file_location, "wb"), "write", miniters=1,
                       total=int(response.headers.get('content-length', 0)),
                       desc=dataset_name) as fout:
        for chunk in response.iter_content(chunk_size=4096):
            fout.write(chunk)


def download(dataset_name:str='All', overwrite = True):
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
    
    # Generate a list of datasets that needs to be downloaded

    if dataset_name == 'All':
        download_list = [ds_name for ds_name in datasets_url_dict.keys()]
    elif type(dataset_name) == list:
        download_list = dataset_name
    elif type(dataset_name) == str:
        download_list.append(dataset_name)
    else:
        print("Please check the list of available options and try again")
        raise NameError(f'{dataset_name} is not valid!')

    for dataset in download_list:
        file_location = os.path.join(os.path.dirname(__file__), os.path.pardir,'data',f'{dataset}.csv')
        if overwrite:
            _get_dataset(file_location, datasets_url_dict[dataset], dataset)
        else:
            if os.path.getsize(file_location) >= 1000000:
                _get_dataset(file_location, datasets_url_dict[dataset], dataset)


    print('Done!')

