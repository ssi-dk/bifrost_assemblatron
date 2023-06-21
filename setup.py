from setuptools import setup, find_packages

setup(
    name='bifrost_assemblatron',
    version='2.2.20',
    description='Denovo assembly component for Bifrost',
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
