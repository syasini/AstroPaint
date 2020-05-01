from setuptools import setup
from astropaint import __version__

with open("requirements.txt", "r") as f:
    reqs = [line.rstrip("\n") for line in f if line != "\n"]

setup(
    name='AstroPaint',
    version=__version__,
    packages=['astropaint', 'astropaint.lib'],
    install_requires=reqs,
    url='https://github.com/syasini/AstroPaint',
    license='',
    author='Siavash Yasini',
    author_email='yasini@usc.edu',
    description='a python package for creating mock maps of astrophysical signals from halo catalog'
    )
