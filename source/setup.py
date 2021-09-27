from setuptools import setup, find_packages
import os 


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
    
setup(
    name='pathogen-embed',
    version='0.0.2',    
    description='Reduced dimension embeddings for pathogen sequences',
    url='https://github.com/blab/cartography/',
    author='Sravani Nanduri <nandsra@cs.washington.edu> , John Huddleston <huddlej@gmail.com>', 
    author_email='nandsra@cs.washington.edu',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT License',
    project_urls = {
        "Documentation": "https://blab.github.io/cartography/",
        "Bug Reports": "https://github.com/blab/cartography/issues",
        "Source Code": "https://github.com/blab/cartography/tree/master/source",
		"Change Log": "https://github.com/blab/cartography/tree/master/source/CHANGES.md",
    },
    package_dir={"": "src"},
	packages = find_packages(where="src", exclude=['test']),
    #packages=['seaborn', 'scikit-learn', 'umap-learn', 'matplotlib', 'pandas', 'numpy', 'hdbscan'],
    install_requires=['numpy', 
                      'pandas',
                      "biopython",
                      'seaborn', 
                      'scikit-learn', 
                      'umap-learn', 
                      'matplotlib', 
                      'hdbscan'
                      ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
    ],
    entry_points = {
        "console_scripts": [
            "embed = embed.__main__:main",
        ]
    }
)