from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='embedGen',
    version='0.0.1',    
    description='Reduced dimension embeddings for pathogen sequences',
    url='https://github.com/blab/cartography/',
    author='Sravani Nanduri, John Huddleston', 
    author_email='nandsra@cs.washington.edu, huddlej@gmail.com',
    license='MIT License',
    packages=['seaborn', 'scikit-learn', 'umap-learn', 'matplotlib', 'pandas', 'numpy', 'hdbscan'],
    install_requires=['augur>=12.0.0',
                      'numpy', 
                      'pandas' 
                      ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        'Intended Audience :: Science/Research',
        'License :: MIT License',         
        'Programming Language :: Python :: 3.7',
    ],
)