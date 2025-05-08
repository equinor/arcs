
"""
arcs: Automated Reactions for CO2 Storage
"""

from os.path import abspath, dirname
from setuptools import find_packages, setup

setup(
    name='arcs',
    version='1.4.0',
    description='automated reactions for CO2 storage',
    url="https://github.com/badw/arcs",
    author="Benjamin A. D. Williamson",
    author_email="benjamin.williamson@ntnu.no",
    license='MIT',
    packages=['arcs','arcs/dash_app','app'],
    package_data={'arcs/dash_app':['./assets/style.css',
                                  './assets/images/logos.png'],
                  'app':['./data/*']
                  },
    install_requires=[
        'ase==3.23.0',
        'pymatgen==2024.4.13', 
        'scipy==1.14.0', 
        'numpy==1.26.4',
        'chempy==0.9.0',
        'tqdm',
        'networkx==3.2.1',
        'pathos==0.3.1',
        'dash==2.18.1',
        'monty==2024.2.26',
        'dash_bootstrap_templates==1.2.4',
        'dash_bootstrap_components==1.6.0',
        'dash-loading-spinners==1.0.3',
        'dash-core-components==2.0.0',
        'loguru==0.7.2',
        'pywebview',
        'setproctitle'
        ],
    python_requires=">=3.10",
    classifiers=[
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Development Status :: 3 - Beta",
            "Intended Audience :: Science/Research",
            "Intended Audience :: System Administrators",
            "Intended Audience :: Industry",
            "Intended Audience :: Information Technology",
            "Operating System :: OS Independent",
            "Topic :: Other/Nonlisted Topic",
            "Topic :: Scientific/Engineering",
        ],
    #entry_points = {
    #    "console_scripts":[
    #        "arcs-app = app.arcs_app:start"
    #        ]
    #    },
    #include_package_data=True
    )
