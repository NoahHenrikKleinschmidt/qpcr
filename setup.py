import setuptools

with open("README.md", "r") as f:
    long_description = f.read()
    
setuptools.setup(
    name="qpcr", 
    version="0.2.1",
    author="Noah H. Kleinschmidt",
    author_email="noah.kleinschmidt@students.unibe.ch",
    description="A module to perform analysis of qPCR data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NoahHenrikKleinschmidt/qpcr.git",
    packages=setuptools.find_packages(),
    package_dir={"": "qpcr"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License ::  GPL-3.0 License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)