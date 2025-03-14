# DeeDeeExperiment

`DeeDeeExperiment` is an S4 class extending the `SummarizedExperiment` framework to
facilitate the integration and management of transcriptomic analysis results.
It introduces two dedicated slots to store Differential Expression (DE) analysis
results and Functional Enrichment analysis outcomes, providing a structured approach
for downstream analysis.

## Installation

You can install the development version of `DeeDeeExperiment` from GitHub with

``` r
library("remotes")
remotes::install_github("imbeimainz/DeeDeeExperiment",
                        dependencies = TRUE,
                        build_vignettes = TRUE)
```

## Strucutre and Usage

The `DeeDeeExperiment` extends `SummarizedExperiment` and contains additional attributes:

* `dea` : A slot for storing DE-related information (currently supported formats:
results from `DESeq2`, `edgeR`, `limma`).

* `fea` : A slot for storing Functional Enrichment related information

TODO: later attach a schematic representation of the class

## Example

TODO

## Development

If you encounter a bug, have usage questions, or want to share ideas and functionality
to make this package better, feel free to file an
[issue](https://github.com/imbeimainz/DeeDeeExperiment/issues).

## License

MIT
