# `muscat_res`

A small example object generated with
[`muscat::pbDS()`](https://rdrr.io/pkg/muscat/man/pbDS.html) ona subset
Kang et al. dataset. The results represent DE between stimulated vs
control samples.

## Format

A `list` structured as returned by
[`muscat::pbDS()`](https://rdrr.io/pkg/muscat/man/pbDS.html)

## Value

A named `list` of DE results tables compatible with
`muscat_list_for_dde`

## Details

The original data were obtained from Kang et al. dataset in
`ExperimentHub` and processed using `muscat` workflows. Only a small
subset of genes and clusters was retained to reduce object size.
