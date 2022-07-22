.. _index:

.. qpcr documentation master file, created by
   sphinx-quickstart on Wed Jul 20 17:13:20 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to qpcr's documentation!
================================

``qpcr`` is python package designed to facilitate the analysis of qPCR data through established Delta-Delta-Ct analysis.
It supports a number of ways to import data from a various of source files, pre-process them, analyze them, and finally visualize the results.
Hence, ``qpcr`` aims to provide a maximally user-friendly analysis pipeline for qPCR data, with a lot of costumizability.


Installation
============
This package can directly be installed via ``pip``.

.. code-block:: bash

   pip install qpcr


Example usage
=============


A very basic use case may be that you have an excel file that contains the qPCR data of some qPCR assays you want to analyse.
Assuming you are happy with default settings, the entire analysis could look like this:

.. code-block:: python

   import qpcr

   myfile = "my_datafile.xlsx"

   # read the datafile to get "assays-of-interest" and "normalisers"
   assays, normalisers = qpcr.read_multi_assay(  myfile, decorator = True )

   # now we can compute first Delta-Ct values on all our loaded assays
   assays = qpcr.delta_ct( assays )
   normalisers = qpcr.delta_ct( normalisers )

   # and now we can normalise our assays-of-interest using our normalisers
   results = qpcr.normalise( assays, normalisers )

   # and we can save our results and have a look at them
   results.save( "my_results.csv" )
   results.preview()


To facilitate data analysis, common workflows have been implemented in pre-defined *pipelines* that allow for quick data analysis with minimal user effort. 
You can learn more about the pipelines in `this tutorial <https://github.com/NoahHenrikKleinschmidt/qpcr/blob/main/Examples/2_pipeline_tutorial.ipynb>`_ 
or `this one <https://github.com/NoahHenrikKleinschmidt/qpcr/blob/main/Examples/5_customisable_piplines.ipynb>`_.
A more advanced example of such an analysis using the pre-defined ``BasicPlus`` pipeline is this:

.. code-block:: python

   from qpcr.Pipes import BasicPlus

   # get some datafiles
   normaliser_files = [
                     "./Examples/Example Data/28S.csv",
                     "./Examples/Example Data/actin.csv"
                  ]

   assay_files = [
                  "./Examples/Example Data/HNRNPL_nmd.csv",
                  "./Examples/Example Data/HNRNPL_prot.csv",
                  "./Examples/Example Data/SRSF11_nmd.csv",
                  "./Examples/Example Data/SRSF11_prot.csv",
               ]

   # define our experimental parameters
   reps = 6
   group_names = ["WT (-)", "WT (+)", "KO (-)", "KO (+)"] 

   # setting up the pipeline
   pipeline = BasicPlus()
   pipeline.save_to("./Example Results")

   pipeline.replicates(reps)
   pipeline.names(group_names)

   # feed in our data
   pipeline.add_assays(assay_files)
   pipeline.add_normalisers(normaliser_files)

   # run the pipeline
   pipeline.run()

   # and at this point the results are already saved...

   # show a preview of our results with some customization
   results = pipeline.results()

   results.drop_rel()
   fig = results.preview( mode = "interactive", title = "Normalized to 28S+Actin" )


.. raw:: html
   :file: resources/preview_index.html


If you find ``qpcr`` useful for your research, please cite it ^^
   Kleinschmidt, N. (2022). qpcr - a python module for easy and versatile qPCR data analysis for small-scale datasets and high-throughput (Version 3.1.5) [Computer software]. https://github.com/NoahHenrikKleinschmidt/qpcr.git


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   qpcr
   gettingstarted
   tutorials
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
