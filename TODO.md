
# New features of *this* release

New planned features for this next release are: 

### 1. Qupid 
New Qupid web-app for easy access to the analysis pipeline for non-experienced users. 

- 1.1 Make a new Qupid web-app
- 1.2 Make functional Cores within `qpcr` for Qupid

| Point | Done | In progress | Stuck |
| ----- | ---- | ----------- | ----- |
|    1.1   |  x    |            |       |
|    1.2   |   x   |             |       |

### 2. Anchor 
- 2.1 This new release has fixed the "grouped" anchor default settings.
- 2.2 Also it added a new feature that allows linking a custom `function` as anchor instead of an externally computed value or the default "first" or "grouped" arguments.

| Point | Done | In progress | Stuck |
| ----- | ---- | ----------- | ----- |
|    2.1   |  x    |            |       |
|    2.2   |   x   |             |       |


### 3. Blueprint pipeline
We added a new pipeline `Blueprint` which allows customization of Analyser, Normaliser, and SampleReader

- 3.1 New Blueprint Pipeline

| Point | Done | In progress | Stuck |
| ----- | ---- | ----------- | ----- |
|    3.1   |  x    |            |       |
### 4. Infer group names
We add a method to the SampleReader which will adopt Replicate Group names based on their given sample column.  We will start with a simple inference which will only adopt group names if all replicates share the same name.

- 4.1 identical group identifier

| Point | Done | In progress | Stuck |
| ----- | ---- | ----------- | ----- |
|    4.1   |    x  |            |       |


### 5. Support Excel Files
Now users can upload an excel file containing the replicates and their ct values.
The qpcr.Reader will parse the excel sheet to find the replciates and their ct values and generate a pandas dataframe from them

- 5.1 Support Excel File Reading
- 5.2 Enable Excel Files for Qupid

| Point | Done | In progress | Stuck |
| ----- | ---- | ----------- | ----- |
|    5.1   |    x  |            |       |
|    5.2   |    x (I think?)  |            |       |

# New features of *some future* release

New planned features for any future release are:


### 1. Integration of `Qlipper` functionality 
- 1.5 while we're at it, we could also implement a new pipeline that starts from absorption data instead of Ct values (we name it `QlipperPipe`, or something like that...)

- 1.6 support reading `rex` files directly by `Qlipper`. Should be easy enough locally, but might be difficult to do via streamlit... We'll need to figure out a way to get a FileIO from the UploadedFile ...

### 2. Merging Filter before-after figures together
Currently before and after filtering figures are separate. It would be nice to have the before filtering figure (maybe with different color or with less opacity) be plotted and then to overlay the after-filtering on-top of that - all in the same figure...

### 3. Statistical evaluation of results 
Currently, no options are available for preforming t-tests and such on the results. It would be nice to have that functionality at some point. However, the main tricky part here is to decide which groups of replicates to compare against each other. One option would be to simply compare all and draw a comparison heatmap for all groups. But this is invariably computing more than necessary or desired. A sample-index file of sorts would be required here probably...

### 4. Double-normalisation
Currently the workflow is only designed for a single normalisation. Comparing the levels of transcript A from condition X against the same of condition Y (i.e. get the fold-change) is not directly difficult but currently not implemented due to the same problem for the statistical evaluation - missing knowledge which groups of replicates to pair together...

> IDEA: Just an idea about the index file implementation. Given that we can manually set ids after SampleReader has read in a file, we might simply implement a re-id at this point using the Basic pipeline workflow. Like this we would be able to pair groups of replicates together pretty straightforwardly (unless we loose that specificity again after the first normalisation step, but this should be no problem due to the `drop_rel` method of the qpcr.Results, so we should be able to re-store the original ids after the first normalisation pretty easily...)

### 5. New pipelines with more versatile customization
Two new pipelines are supposed to join the current Basic ones: The `Blueprint` pipeline and the `Framework` pipeline. They are just again both doing essentially the same as the Basic ones, except that Blueprint allows to link an externally set up SampleReader, Analyser, and Normaliser as desired, while using default settings for anything that was not explicitly set. Framework on the otherhand really requires the user to explicitly set all of these, it only provides the workflow framework.


### 6. Estimation of qPCR efficiency
There are some nice papers describing how they estimated qPCR amiplification efficiency based on linreg around the linear window (which we use currently just to get the optimal threshold through R^2). The slope of the optimal window should also be the corresponding efficiency. Let's test this out and add a method to replace the default `efficiency = 2` with a properly computed one. Let's try to get the slopes of our window ranges as well during the optimal threshold search, and check if we can truly use these. In the paper they keep using a log-scale to do things, so, maybe we'll have to do the same. Let's see...