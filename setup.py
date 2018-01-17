
from setuptools import setup, find_packages

def readme():
    with open('README.md', 'r') as info:
        return info.read()

setup(name='datagen',
      version='0.1.0',
      description='A package to extract and analyse a dataset of all-beta protein structures from the CATH or SCOPe database',
      long_description=readme(),
      url='https://github.com/kls93/DataGen',
      author='Kathryn Shelley',
      author_email='kathryn.l.shelley@gmail.com',
      classifiers=['Programming Language :: Python :: 3.6'],
      packages=find_packages(),
      install_requires=['pandas', 'matplotlib', 'networkx'],
      package_data={'CATH_domains_desc': [''],
                    'SCOP_domains_desc': ['']},
      entry_points={'console_scripts': ['datagen.datagen:main']}
)
