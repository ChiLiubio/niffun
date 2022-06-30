# niffun
An R package for functional prediction of diazotrophic community based on nifH sequences.

![](https://img.shields.io/badge/Test-0.0.1-red.svg)

## Install

Install the latest development version from github.

```r
# If devtools package is not installed, first install it
install.packages("devtools")
devtools::install_github("ChiLiubio/niffun")
```

Please download blast tools from "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+". 
For windows system, please download v2.5.0 as some errors can come from the latest versions because of memory issue. 
Then, place the blast into the computer system env path or provide the directory (such as 'ncbi-blast-2.5.0+/bin') to blast_tool_path parameter of cal_blast function.

## Example

```r
library(niffun)
# use the example data dataset_nifH
data("dataset_nifH")
n1 <- trans_niffun$new(dataset = dataset_nifH)
# create test directory for the temporary files
dir.create("test")
# blast
n1$cal_blast(path_to_temp_folder = "test", evalue = 1e-05)
n1$cal_pathway()
# further analysis
library(microeco)
data(Tax4Fun2_KEGG)
# create a microtable object for pathways
func1 <- microtable$new(otu_table = n1$res_KEGG_pathway, tax_table = Tax4Fun2_KEGG$ptw_desc, sample_table = dataset_nifH$sample_table)
func1$tidy_dataset()
func1$cal_abund()
```


## Contributing

We welcome any contribution, including but not limited to code and idea.
Please report errors and questions on github [Issues](https://github.com/ChiLiubio/niffun/issues).
Any contribution via [Pull requests](https://github.com/ChiLiubio/niffun/pulls) or Email(liuchi0426@126.com) will be appreciated.
By participating in this project you agree to abide by the terms outlined in the [Contributor Code of Conduct](CODE_OF_CONDUCT.md).

