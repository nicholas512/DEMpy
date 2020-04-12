from setuptools import setup

setup(name='dempy',
      version='1.0',
      description='Functions for downloading DEM tiles',
      url='http://github.com/nicholas512/dempy',
      author='Nick Brown',
      packages=['dempy'], 
      install_requires=[
                        'numpy',
                        'pandas',
                        'requests'])