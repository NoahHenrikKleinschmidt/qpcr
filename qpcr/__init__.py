"""
This module is designed to provide functions to analyse qPCR data. 
It is designed for maximal user-friendliness and streamlined data-visualisation.

## Terminology
---

Let's talk some terminology first. There are a number of important terms that you will meet when using the `qpcr` module that you may feel have nothing to do with qPCR. 
Read on to learn what these terms mean and why they are important. 

### `qpcr` vs qPCR
Throughout this documentation, we will refer to the python module as `qpcr` in all lowercase, while we will talk about the qPCR experiment itself 
by using `qPCR` with uppercase. 

### File (or "datafile")
Let's start simple. A `file` is simply one of your input datafiles (shocker). For the `qpcr` module this means either a `csv` or an `excel` file. 
Files are the at the very basis of our data pipeline. 
However, `qpcr` does not strictly require that these files correspond to anything in particular from "qPCR vocabulary" _per se_ 
(although the basic assumption underlying the default settings
is that a single "regular" datafile contains exactly one single qPCR assay). 

The main `qpcr` classes are designed to work with datafiles that contain Ct values and replicate identifiers (see below for "replicates").
In the very most basic case, a simple datafile contains extactly an `"id"` or `"Name"` or whatever column and a `"Ct"` column (can also be named differently),
that store values from exactly _one single qPCR assay_. 

However, if your experimental setup looks differently, there are ways to adapt `qpcr` to handle different setups.
But if your data follows the above mentioned standard arrangements, you can think of a file as identical to a "qPCR assay".

### "Regular" vs "irregular" datafiles
`Regular` datafiles only contain the `id` and `Ct` columns and nothing else, storing values of exactly one qPCR assay. 
`Irregular` datafiles also have to have these two columns, but they allow for more irregular stuff
around these columns. 
Irregular files may contain multiple datasets / assays (see below). Check out the documentation of the `qpcr.Parsers` and `qpcr.Readers` for 
details on how to work with irregular and multi-assay datafiles.

### An `qpcr.Assay` and a "dataset" / "assay"
Let's start by the term "dataset". Well, a "dataset" is simply a collection or replicate identifiers and corresponding Ct values belonging together. 
By default this usually corresponds to a single qPCR assay, hence you will often encounter the term "assay" used instead of "dataset" to be more 
intuitive. However, the two terms are really interchangeable. The `qpcr.Assay` class is used to store the Ct values and replicate identifiers 
from a single dataset extracted from a datafile. 

### Replicates
As far as the `qpcr` module is concerned, each row within your datasets corresponds to one replicate. 
Hence, a _replicate_ is just a single pair of some identifier and a corresponding Ct value. 
That means that `qpcr` does not terminologically distinguish between multiplets (i.e. "true replicates" as an experimentor would understand them) and unicates. 
If you feel more comfortable with a term like "measurement" then you can also think of the replicates like that. 

### Groups (of Replicates)
This is one of the most important terms. In the previous paragraphs we already started to talk about replicates "belonging together". 
Groups of replicates are, as their name already implies, well, a number of replicates that somehow do belong together.
For most cases, if your datafiles contain one assay each, then the groups of replicates are most likely your different qPCR samples / experimental conditions. 
We talk about _groups_ of replicates instead of samples or conditions, primarily, because there might be different data setups so that these terms might not be always appropriate.
If your data follows default arrangements, however, then a _group of replicates_ is just what you would think of as a "qPCR sample". 
Groups are assigned a numeric index starting from 0, which is how they are identified by the classes. 
However, they also come with a text label called the `group_name` (you can manually set and re-set the group names as you like). 
Many classes such as the `qpcr.DataReader` will actually just use the term `names` instead of the full `group_names`. 
Whenever you see anything "names"-related it is (super-duper most likely) a reference to the `group_names`.

### Specifying (Groups of) Replicates
Here's the best part: usually, we don't necessarily need to do anything because `qpcr.Assay` objects are able to infer the groups of replicates in your data 
automatically from the replicate identifiers (yeah!). However, you will be asked to manually provide replicate settings in case this fails. 
In case you want to / have to manually specify replicate settings, an `qpcr.Assay` accepts an input `replicates` which is where you can specify this information. 

`replicates` can be either an `integer`, a `tuple`, or a `string`. Why's that? Well, normally we perform experiments as "triplicates", or "duplicates", or whatever multiplets.
Hence, if we always have the same number of replicates in each group (say all triplicates) we can simply specify this number as `replicates = 3`. 
However, some samples might only be done in unicates (such as the diluent sample), while others are triplicates.

In these cases your dataset does not have uniformly sized groups of replicates and a single number will not do to describe the groups of replicates. 
For these cases you can specify the number of replicates in each group separately as a `tuple` such as `replicates = (3,3,3,3,1)` or as a `string` "formula"
which allows you to avoid repeating the same number of replicates many times like `replicates = "3:4,1"`, which will translate into the same tuple as we specified manually. 
Check out the documentation of `qpcr.Assay.replicates` for more information on this. 

### `Delta-Ct` vs `Delta-Delta-Ct` vs `normalisation`
The default analysis workflow in $\Delta \Delta Ct$ analysis is to first calculate a $\Delta Ct$ using an intra-assay reference and then calculate the $\Delta \Delta Ct$ using a normaliser assay. 
The first $\Delta Ct$ step is performed by a class called `qpcr.Analyser` using its native method `qpcr.Analyser.DeltaCt`. 
Why just `DeltaCt` and not `DeltaDeltaCt`? Well, we call this second `Delta`-step in `DeltaDeltaCt` differently. 
We name it `normalisation`, and it is handled by a class called `qpcr.Normaliser` using its native method `qpcr.Normaliser.normalise`. 
So, as far as the `qpcr` module is concerned there is only `qpcr.Analyser.DeltaCt` which performs the first $\Delta Ct$, and `qpcr.Normaliser.normalise` which later handles the second "delta"-step to get to $\Delta \Delta Ct$. 
Of course, this means that the `qpcr.Normaliser` will need to have knowledge about which `qpcr.Assay` objects contain actual assays-of-interest and which ones contain normaliser-assays (specifying that is easy, though, so don't worry about that).


### The `anchor` and the "reference group"
Next to the "groups of replicates", this is probably one of the most important terms. The `anchor` is simply the intra-dataset reference used by the `qpcr.Analyser` to perform its first $\Delta Ct$. 
If your datafiles contain one assay each, and your groups of replicates are your qPCR samples, then you will likely have some "wildtype", "untreated", or "control" sample, right? 
Well, in `qpcr` terms that would be your _reference group_.
Usually your `anchor` is part of or generated from the Ct values of your _reference group_ (like their `mean` for instance).
By default it is assumed that your reference group is the _very first_ group of replicates. However, it's not a big problem if this is not the case, as you can specify different anchors easily.
So, again, the `anchor` is the dataset-internal reference value used for the first $\Delta Ct$.

### "assays" vs "normalisers"
You will likely encounter methods and/or arguments that speak of "assays" and "normalisers", especially with the `qpcr.Normaliser`. 
For all intents and purposes, an "assay" is simply one of your datasets (we know this already).
However, in practice "assays" are the short notation for specifically "assays-of-interest" 
(or more formally "datasets-of-interest"), while "normalisers" refer to your normaliser-assays (from housekeeping genes like ActinB for instance). 
But again, if your datafiles do not conform to standard data arrangements, do not be distracted from the terminology here.

You will also find that the term "assays" is used within the final results dataframe (when using the summary-statistics mode). 
In this setting "assays" refers to the assay-of-interst whose data was analysed according to the provided normaliser-assays. 
In fact, this is a new "hybrid" assay identifier taht includes the names of all the normaliser-assays used during computation (check out what the final results look like and it'll be immediately clear).

### "Samples"
You may find that there is also a term "sample" within `qpcr`'s vocabulary. 
As far as the `qpcr` module is concerned, the term "sample" is not very important in itself and usually appears in the context of "sample assays".
In this setting it is used interchangeably with "assays-of-interest". 
Actually, we try to phase out the term "sample" and it currently mainly appears in hidden auxiliary functions which have retained the term from earlier development versions.


## Some more basics 
----

#### The values stored by `qpcr`
Please, note at this point that, as described in more detail in the documentation of the `qpcr.Analyser`, Delta-Ct values are directly computed as 
exponential values $ \Delta Ct' = 2^{ - \Delta Ct}$, while normalisation later performs $ \mathrm{norm. } \Delta\Delta Ct = \\frac{  \Delta Ct'_s  }{  \Delta Ct'_n  }$, where $s$ is an assay of interest's Delta-Ct ($\Delta Ct'$) value of some replicate, 
and $n$ is the corresponding value of the normaliser assay. This is based on the mathemathical equivalence of $n^{  a - b  } \equiv \\frac{  n^{ a } } {  n^{ b } }$. 
Hence, while the documentation will continuously use the terms $\Delta Ct$ and $\Delta\Delta Ct$, they are in fact the exponential deriviative of the conventional values.

#### Modes of `normalisation`
The `qpcr.Normaliser` supports custom functions for normalisation. However, it also comes with three built-in methods to normalise sample-assays against normalisers.
These are accessible via the `mode` argument of the `qpcr.Normaliser.normalise` method, which can be set to `"pair-wise"` (default), `"combinatoric"`, or `"permutative"`. 
The default option `"pair-wise"` is computationally the fastest and will rigidly normalise replicates against their corresponding partner from the normaliser. I.e. first against first, 
second against second, etc. This mode is appropriate for multiplex qPCR experiments. For qPCR reactions that were pipetted individually, there is not reason to strictly only pair first
with first, second with second etc. For these cases there are other two options `"combinatoric"` and `"permutative"`. `"combinatoric"` normalisation will calculate all possible group-wise 
combinations of a sample-assay replicate with all available normaliser replicates of the same group. I.e. first against first, and against second, etc. This will generate $n^2$ values where 
$n$ is the number of replicates within a group. This mode is appropriate for small-scale datasets but will substantially increase computation times for larger datasets.   
`"permutative"` on the other hand will reflect the equivalence of replicates within a group through random permutations wihtin the normaliser replicates. Hence, 
first may be normalised against first, or second, etc. This normalisation method may by used iteratively to increase the accuracy. Replacement during permutations are allowed (although disabled by default).
If replacement is desired by the user, the probability of each replicate to be chosen will be weighted based on a fitted normal distribution. This method is appropriate for larger datasets for which combinatoric normalisation
is not desired. 

### `pipeline`s 
A `pipeline` is essentially any workflow that starts from one or multiple input datafiles and ultimately pops out some results table you are happy with.
Pipelines can be manually created by assembling the main `qpcr` classes, usually starting with a Reader, passing to an Analyser, to an Normaliser, and you're good to go.
When manually assembling your workflow you can extract your data at any point and perform your own computations on it as you like. However, if you wish to "just do some good ol' Delta-Delta-Ct"
there are pre-defined pipelines that will handle writing the workflow and only require a very basic setup. You can find these in the `qpcr.Pipes` submodule.

### "Results" = my final results?
Yes and no. Anything that is computed through any of the `qpcr` classes is called a "result" of some kind.
In practice as soon as you pass your data through a `qpcr.Normaliser` you generate "results" which are actually stored in a separate class called `qpcr.Results` (that's not so important, though).
So, when using the term "result" we usually mean just anything that was explicitly computed by `qpcr`. 

### `get`ting your data
Too many classes and objects? Well, no worries, the underlying data is stored as `pandas DataFrames`. To get your data from the clutches of the `qpcr` classes you can always use the `get()` method. 
`get` is pretty universal in the `qpcr` module, so whenever you want to extract your data, there's a `get()` method to help you.

### `link` vs `add` vs `pipe`
Different classes have slightly different methods of adding data to them. Classes that only accept one single data input (such as a single `qpcr.Assay` object or a single filepath as `string`)
usually have a `link()` method that, well, links the data to them. After that the classes are ready to perform whatever actions they can perform (an `qpcr.Analyser` would perform `DeltaCt()` for instance). 
Some classes such as the `qpcr.Analyser` have a wrapper that will call both their `link()` as well as their actual core-functional method together in one go. This wrapper is called `pipe()`. 
So for the `qpcr.Analyser` you could either manually use `link()` and then `DeltaCt()`, or simply call `pipe()` which does both for you. It is noteworthy that `pipe` methods actually _return_ whatever
their output is, which is *not* normally the case otherwise (normally you'd use the `get()` method to extract your data, see above). Most `qpcr.Readers` and `qpcr.Parsers` are also equipped with `pipe` methods.
Alright, we know about `link` and `pipe` now, what about `add`?  Classes that accept multiple inputs have `add` methods, which tells the class where exactly to store the input data. 
`add`-methods are especially implemented within the pre-defined analysis pipelines of the `qpcr.Pipes` submodule. You will probably often use the methods `add_assays()` and `add_normalisers()` if you plan on using these predefined pipelines.
However, these classes usually still have a `link()` method somewhere that you can use as well. For instance, the `qpcr.Normaliser` uses a `link` method to add both assays-of-interest and normaliser-assays simultaneously.

## Getting started
----
You can find useful tutorials and applied examples as `jupyter notebooks` [on GitHub](https://github.com/NoahHenrikKleinschmidt/qpcr/tree/main/Examples).

"""

from qpcr.main import * 
from qpcr.defaults import init_log_format

import logging
logging.basicConfig( format = init_log_format )