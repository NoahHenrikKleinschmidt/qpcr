
# New implemented features of *this* release


# New planned features of *some future* release

New planned features for any future release are:
### Readers

New support to read "regular" files that contain multiple columns. 
- Implement support of this for the BigTableReader where there is an id_col that specifies the samples but the assays are in side-by-side columns, the input here would be a list of column headers for the assay_col argument.


### 1. Integration of `Qlipper` functionality 
- 1.5 while we're at it, we could also implement a new pipeline that starts from absorption data instead of Ct values (we name it `QlipperPipe`, or something like that...)

- 1.6 support reading `rex` files directly by `Qlipper`. Should be easy enough locally, but might be difficult to do via streamlit... We'll need to figure out a way to get a FileIO from the UploadedFile ...

### 2. Merging Filter before-after figures together
Currently before and after filtering figures are separate. It would be nice to have the before filtering figure (maybe with different color or with less opacity) be plotted and then to overlay the after-filtering on-top of that - all in the same figure...

> UPDATE Idea: so, we let the ReplicateBoxplot be what it is right now, but instead add a new FilterBoxPlot plotter that makes subplots n x 2 where n is the number of assays, which will be on rows above each other. Then left to right is before and after. This is then just one big figure, which should hopefully more nicely display the before-after relations within the assays, but also respond to the same legend so we should hopefully be able to modulate all before and after replicates boxes with one click... 

> Also, we should re-think if we maybe not only want to set the Ct values to NaN instead of truly dropping the indices entirely...



### 3. Statistical evaluation of results 
Currently, no options are available for preforming t-tests and such on the results. It would be nice to have that functionality at some point. However, the main tricky part here is to decide which groups of replicates to compare against each other. One option would be to simply compare all and draw a comparison heatmap for all groups. But this is invariably computing more than necessary or desired. A sample-index file of sorts would be required here probably...

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


### 5. New pipelines with more versatile customization
Two new pipelines are supposed to join the current Basic ones: The `Blueprint` (<< CHECK) pipeline and the `Framework` pipeline. They are just again both doing essentially the same as the Basic ones, except that Blueprint allows to link an externally set up SampleReader, Analyser, and Normaliser as desired, while using default settings for anything that was not explicitly set. Framework on the otherhand really requires the user to explicitly set all of these, it only provides the workflow framework.


### 6. Estimation of qPCR efficiency
There are some nice papers describing how they estimated qPCR amiplification efficiency based on linreg around the linear window (which we use currently just to get the optimal threshold through R^2). The slope of the optimal window should also be the corresponding efficiency. Let's test this out and add a method to replace the default `efficiency = 2` with a properly computed one. Let's try to get the slopes of our window ranges as well during the optimal threshold search, and check if we can truly use these. In the paper they keep using a log-scale to do things, so, maybe we'll have to do the same. Let's see...


### 7. Multi-Assay Excel support << CHECK ^^
We want to be able to read and split multi-assay containing excel files into individual csv files that adhere to the structure of input files for the `qpcr` module. We want to make 
a stand-alone web-app for this as well... 
### 1. Qupid 
New Qupid web-app for easy access to the analysis pipeline for non-experienced users. 

- 1.1 Make a new Qupid web-app
- 1.2 Make functional Cores within `qpcr` for Qupid
- 1.3 factor out _Qupid classes into _Qupid submodule
- 1.4 Make Qupid work with Parsers (currently they don't like the StringIO as self._src, so we need to maybe some _QupidParser classes that will use similar hacks as the _Qupid_SampleReader...)

| Point | Done | In progress | Stuck |
| ----- | ---- | ----------- | ----- |
|    1.1   |      |     x       |       |
|    1.2   |      |    x (currently broken, haven't looked at it yet)         |       |
|    1.3   |   dropped   |             |       |
|    1.4   |      |             |       |

### 8. __str__ methods 
At some point we should add some `__str__` methdos to all classes ...

### 11. Add vectorized support or pipe methods
Pipe methods and such should also be able to work directly with a list of input files ...



