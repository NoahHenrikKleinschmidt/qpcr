
# New features of *this* release

New planned features for this next release are: 

### 1. Qupid 
New Qupid web-app for easy access to the analysis pipeline for non-experienced users. 

- 1.1 Make a new Qupid web-app
- 1.2 Make functional Cores within `qpcr` for Qupid
- 1.3 factor out _Qupid classes into _Qupid submodule
- 1.4 Make Qupid work with Parsers (currently they don't like the StringIO as self._src, so we need to maybe some _QupidParser classes that will use similar hacks as the _Qupid_SampleReader...)

| Point | Done | In progress | Stuck |
| ----- | ---- | ----------- | ----- |
|    1.1   |      |     x       |       |
|    1.2   |   x   |             |       |
|    1.3   |      |             |       |
|    1.4   |      |             |       |


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
- 4.2 infer number of replicates in case of identical group identifiers

| Point | Done | In progress | Stuck |
| ----- | ---- | ----------- | ----- |
|    4.1   |    x  |            |       |
|    4.2   |   x   |            |       |


### 5. Support Excel Files
Now users can upload an excel file containing the replicates and their ct values.
The qpcr.Reader will parse the excel sheet to find the replciates and their ct values and generate a pandas dataframe from them

- 5.1 Support Excel File Reading
- 5.2 Enable Excel Files for Qupid

| Point | Done | In progress | Stuck |
| ----- | ---- | ----------- | ----- |
|    5.1   |    x  |            |       |
|    5.2   |    x (I think?)  |            |       |


### 6. Formula for unequal group lengths
If the replicate groups are of unequal size (i.e. we have both triplicates and duplicates or whatever) then so far each group had to be manually specified in a tuple. Now users can input a string formula to automatically generate the corresponding tuple.
For example, if we have four triplicates, four unicates, three duplicates, and one nonaplicates, then our tuple would have to look like this: `(3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 9)`. Users can either enter this manually or they can use a comprehensive formula now to describe: `"3:4,1:4,2:3,9"` which will be translated into the above tuple. 
The formula is always `n:m` where `n` is the number of replicates in the group and `m` is the number of times this kind of group size is repeated afterward. Another example with triplicates, interspersed with a duplicate, and then triplicates again: `"3:5,2,3:3` (you get the idea...)

- 6.1 a string formula method to supply unequal group sizes

| Point | Done | In progress | Stuck |
| ----- | ---- | ----------- | ----- |
|    6.1   |    x  |            |       |

### 8. Support irregular input files
Support irregular csv files where the assay must actually be extracted through parsing. The core of this is already done but still has to be implemented into `qpcr` proper.

- 8.1 implement `qpcr.Parsers` for irregular files
- 8.2 integrate Parsers with `qpcr.Reader`

| Point | Done | In progress | Stuck |
| ----- | ---- | ----------- | ----- |
|    8.1   |   x   |            |       |
|    8.2   |  x    |            |       |

### 7. Support MULTI-ASSAY FILES
We simply let people add a little decorator in the cell above the assay declaration
of the excel / csv file that will tell the program if it's a assay of interest or a normaliser. We'll use something like `@qpcr:assay` and `@qpcr:normaliser`.

- 7.1 make a decorator extraction method to add to the _CORE_Parser
- 7.2 add decorator support for multi-assay csv files
- 7.3 add decorator support for multi-assay excel files
- 7.4 make a new `qpcr.MultiReader` class to read multi-assay files
      this will probably be some kind of `SampleReader`-level class



| Point | Done | In progress | Stuck |
| ----- | ---- | ----------- | ----- |
|    7.1   |      |            |       |
|    7.2   |      |            |       |
|    7.3   |      |            |       |
|    7.4  |      |            |       |


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


### 7. Multi-Assay Excel support
We want to be able to read and split multi-assay containing excel files into individual csv files that adhere to the structure of input files for the `qpcr` module. We want to make 
a stand-alone web-app for this as well... 