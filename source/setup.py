from setuptools import setup, find_packages
import os 


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
    
setup(
    name='embedGen',
    version='0.0.1',    
    description='Reduced dimension embeddings for pathogen sequences',
    url='https://github.com/blab/cartography/',
    author='Sravani Nanduri <nandsra@cs.washington.edu> , John Huddleston <huddlej@gmail.com>', 
    author_email='nandsra@cs.washington.edu',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT License',
    #project_urls = {
    #    "Documentation": "" Do I need to publish my documentation somewhere in order to have this link work?
    #}
	packages=find_packages(exclude=['test']),
    #packages=['seaborn', 'scikit-learn', 'umap-learn', 'matplotlib', 'pandas', 'numpy', 'hdbscan'],
    install_requires=['numpy', 
                      'pandas'
                      ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',         
        'Programming Language :: Python :: 3.7',
    ],
)