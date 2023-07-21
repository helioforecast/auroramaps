# -*- coding: utf-8 -*-

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="auroramaps",
    version="0.3",
    author="Christian Moestl",
    author_email="chris.moestl@outlook.com",
    keywords=["geophysics", "heliophysics", "space weather"],
    description="Aurora model OVATION PRIME 2010 open source",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/helioforecast/auroramaps",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
)
