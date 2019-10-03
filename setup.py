import os
from setuptools import setup, find_packages

packages = [x for x in find_packages('.') if x.startswith('aisynphys')]

setup(
    name = "aisynphys",
    version = "0.0.1",
    author = "Luke Campagnola",
    author_email = "lukec@alleninstitute.org",
    description = ("Analysis and data management tools used for Allen Institute multipatch pipeline"),
    license = "Allen Institute Software License",
    url = "http://github.com/aiephys/aisynphys",
    packages=packages,
    classifiers=[
        "Development Status :: 3 - Alpha",
    ],
)


