# New planned features of *some future* release

New planned features for any future release are:
### Readers

New support to read "regular" files that contain multiple columns. 
- Implement support of this for the BigTableReader where there is an id_col that specifies the samples but the assays are in side-by-side columns, the input here would be a list of column headers for the assay_col argument.


### 1. Integration of `Qlipper` functionality 
- 1.5 while we're at it, we could also implement a new pipeline that starts from absorption data instead of Ct values (we name it `QlipperPipe`, or something like that...)

- 1.6 support reading `rex` files directly by `Qlipper`. Should be easy enough locally, but might be difficult to do via streamlit... We'll need to figure out a way to get a FileIO from the UploadedFile ...



### 3. Statistical evaluation of results 
Currently, no options are available for preforming t-tests and such on the results. It would be nice to have that functionality at some point. However, the main tricky part here is to decide which groups of replicates to compare against each other. One option would be to simply compare all and draw a comparison heatmap for all groups. But this is invariably computing more than necessary or desired. A sample-index file of sorts would be required here probably...

> UPDATE: 
> Let's go first with a "baseline"-comparison multiple-T-Test for each Assay internally, where we compare against the anchor. 
> Then let's do the same as a pairwise comparison between a baseline assay and a second assay.
> Then, let's also do an ANOVA implementation where we merge the assays and check if there is significant variance between the groups. 

### 4. Double-normalisation
Currently the workflow is only designed for a single normalisation. Comparing the levels of transcript A from condition X against the same of condition Y (i.e. get the fold-change) is not directly difficult but currently not implemented due to the same problem for the statistical evaluation - missing knowledge which groups of replicates to pair together...

> IDEA: Just an idea about the index file implementation. Given that we can manually set ids after SampleReader has read in a file, we might simply implement a re-id at this point using the Basic pipeline workflow. Like this we would be able to pair groups of replicates together pretty straightforwardly (unless we loose that specificity again after the first normalisation step, but this should be no problem due to the `drop_rel` method of the qpcr.Results, so we should be able to re-store the original ids after the first normalisation pretty easily...)
> SPONTANEOUS UPDATE: such an index file should be a vertical big table 
> with:
> normaliser, assay
> normaliser, assay
> ... , ... 
> where the entries are the assay ids...

> IDEA 2: 
> Let the assays split their id to id + label attributes. So "HNRNPL nmd" and "HNRNPL prot" can be split into id "HNRNPL" and label "nmd" + "prot", respectively. Then develop a class / function, whatever, that will pair up assays according to the same id but different labels. The splitting should be available through both simple `.split` and `regex`...


### 8. __str__ methods 
At some point we should add some `__str__` methdos to all classes ...



