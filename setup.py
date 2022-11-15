import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="qpcr", 
    version="4.1.0",
    author="Noah H. Kleinschmidt",
    author_email="noah.kleinschmidt@students.unibe.ch",
    description="A package to perform analysis of qPCR data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NoahHenrikKleinschmidt/qpcr.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.8',
)