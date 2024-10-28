import sys
import setuptools
from decoil import __version__, _program

 

with open("README.md", "r", encoding="utf-8") as fh:
	long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as f:
	required = [x for x in f.read().splitlines() if not x.startswith("#")]

setuptools.setup(name='decoil',
				 version=__version__,
				 test_suite='pytest.collector',
				 tests_require=['pytest'],
				 author="Madalina Giurgiu",
				 author_email="giurgiumadalina25@gmail.com",
				 description="EcDNA reconstruction from long-read nanopore data",
				 long_description=long_description,
				 long_description_content_type="text/markdown",
				 license='MIT',
				 entry_points={
					 'console_scripts': [
						 'decoil-pipeline = decoil.command:entry_point',
						 'decoil = decoil.main:main'
					 ]
				 },
				 extra_require={
					 "dev": ["pytest>=7.0","wheel","setuptools","twine"],
					 "build": [
						 # Define environment variables for the build process
						 "IN_CONTAINER=False",
					 ],
				 },
				 package_data={'': ['Snakefile', '*.json', 'envs/*.yaml', 'rules/*.smk']},
				 install_requires=required,
				 conda_channels=['conda-forge', 'bioconda'],
				 conda_packages={
         				'survivor' : '1.0.7', 
             			'sniffles' : '1.0.12',
                		'deeptools': '3.5.5',
						'ngmlr' : '0.2.7',
						'samtools':'1.15.1',
						'python-dateutil':'2.8.0'},
				 python_requires=">=3.10.0",
				 packages=setuptools.find_packages(),
				 include_package_data=True,
				 keywords=[],
				 zip_safe=False)
