
# New features of *this* release

New planned features for this next release are: 

### 1. Integration of `Qlipper` functionality 
to allow starting from absorption curves already instead of Ct values

- 1.1 write a proper `Qlipper` module (sub-module)

- 1.2 that needs an automated function to find the optimal threshold where all curves are most linear

- 1.3 at that point we will probably also re-write the `Qlipper` stand-alone web-app to work with the new `qpcr.Qlipper` integration...

- 1.4 a visualisation of the absorption curves alongside the threshold would be in order here as well...

| Point | Done | In progress | Stuck |
| ----- | ---- | ----------- | ----- |
|    1.1   |      |             |       |
|    1.2   |      |             |       |
|    1.3   |      |             |       |
|    1.4   |      |             |       |


# New features of *some future* release

New planned features for any future release are:


### 1. Integration of `Qlipper` functionality 
- 1.5 while we're at it, we could also implement a new pipeline that starts from absorption data instead of Ct values (we name it `QlipperPipe`, or something like that...)

### 2. Merging Filter before-after figures together
Currently before and after filtering figures are separate. It would be nice to have the before filtering figure (maybe with different color or with less opacity) be plotted and then to overlay the after-filtering on-top of that - all in the same figure...

### 3. Statistical evaluation of results 
Currently, no options are available for preforming t-tests and such on the results. It would be nice to have that functionality at some point. However, the main tricky part here is to decide which groups of replicates to compare against each other. One option would be to simply compare all and draw a comparison heatmap for all groups. But this is invariably computing more than necessary or desired. A sample-index file of sorts would be required here probably...

### 4. Double-normalisation
Currently the workflow is only designed for a single normalisation. Comparing the levels of transcript A from condition X against the same of condition Y (i.e. get the fold-change) is not directly difficult but currently not implemented due to the same problem for the statistical evaluation - missing knowledge which groups of replicates to pair together...

> IDEA: Just an idea about the index file implementation. Given that we can manually set ids after SampleReader has read in a file, we might simply implement a re-id at this point using the Basic pipeline workflow. Like this we would be able to pair groups of replicates together pretty straightforwardly (unless we loose that specificity again after the first normalisation step, but this should be no problem due to the `drop_rel` method of the qpcr.Results, so we should be able to re-store the original ids after the first normalisation pretty easily...)

### 5. New pipelines with more versatile customization
Two new pipelines are supposed to join the current Basic ones: The `Blueprint` pipeline and the `Framework` pipeline. They are just again both doing essentially the same as the Basic ones, except that Blueprint allows to link an externally set up SampleReader, Analyser, and Normaliser as desired, while using default settings for anything that was not explicitly set. Framework on the otherhand really requires the user to explicitly set all of these, it only provides the workflow framework.


