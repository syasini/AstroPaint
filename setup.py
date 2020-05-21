import os
from setuptools import setup, find_packages


with open("requirements.txt", "r") as f:
    reqs = [line.rstrip("\n") for line in f if line != "\n"]


def get_version():
    """get the package version from __init__.py
    modified from this source: https://packaging.python.org/guides/single-sourcing-package-version/
    """
    package_dir = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(package_dir, "astropaint", "__init__.py"), "r") as f:
        for line in f.readlines():
            if line.startswith('__version__'):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]
    raise RuntimeError("Package version not found!")


setup(
    name='AstroPaint',
    version=get_version(),
    packages=find_packages(),
    install_requires=reqs,
    url='https://github.com/syasini/AstroPaint',
    license='',
    author='Siavash Yasini',
    author_email='yasini@usc.edu',
    description='a python package for creating mock maps of astrophysical signals from halo catalog'
    )
