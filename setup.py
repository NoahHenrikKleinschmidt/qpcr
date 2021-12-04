import setuptools


long_description = "This project represents a python package that includes a number of functions useful for the analysis of qPCR data generated generated by Qiagen RotorGene® and is taylored to work with the Excel Spreadhseet exported from this device (or, more precisely, a csv copy of the same)."

setuptools.setup(
    name="qpcr", 
    version="0.1.7",
    author="Noah H. Kleinschmidt",
    author_email="noah.kleinschmidt@students.unibe.ch",
    description="A module to perform analysis of qPCR data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NoahHenrikKleinschmidt/qpcr.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License ::  GPL-3.0 License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)