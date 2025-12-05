# Display available FEA formats

Display available FEA formats

## Usage

``` r
supported_fea_formats()
```

## Value

a data.frame of possible FEA input formats

## Examples

``` r
supported_fea_formats()
#>         Format         Package
#> 1   data.frame           topGO
#> 2 enrichResult clusterProfiler
#> 3   gseaResult clusterProfiler
#> 4  fgseaResult           fgsea
#> 5   data.frame      gprofiler2
#> 6   data.frame         enrichR
#> 7   data.frame           DAVID
#> 8   data.frame       GeneTonic
```
