from setuptools import setup, find_packages

setup(
    name='bifrost_assemblatron',
    version='v2_2_16',
    description='Datahandling functions for bifrost (later to be API interface)',
    url='https://github.com/ssi-dk/bifrost_assemblatron',
    author="Kim Ng, Martin Basterrechea",
    author_email="kimn@ssi.dk",
    packages=find_packages(),
    install_requires=[
        'bifrostlib >= 2.1.9',
        'cyvcf2 >= 0.30.1'
    ],
    package_data={"bifrost_assemblatron": ['config.yaml', 'pipeline.smk']},
    include_package_data=True
)
