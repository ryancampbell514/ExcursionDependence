# ExcursionDependence

This repository is for reproduction of the results and figures in the manuscript "Analysis of spatial extreme events via excursion set geometry" by Ryan Campbell, Elena Di Bernardino, and Thomas Opitz.


## Dependencies

After cloning this repository, install the code needed for calculating the lower extremal range (LER) of excursion sets by installing the ExcursionSets R package. 
This can be done with
```r
install.packages("ExcursionDependence/ExcursionSets_0.0.4.tar.gz")
```
The ExcursionSets_0.0.4.tar.gz file was initially made available in the [Local-Excursion-Set-Coefficients](https://github.com/pony-helping-treat/Local-Excursion-Set-Coefficients) repository.
Furthermore, obtain the code needed for unbiased excursion set perimeter estimation by cloning the [excursion-sets](https://github.com/RyanCotsakis/excursion-sets) repository. 
Lastly, if not already installed, the following should be executed
```r
install.packages(c("this.path",
                   "evd",
                   "truncnorm",
                   "rmutil",
                   "imager",
                   "terra"))
```


## Rainfall data

The raster files associated with the rainfall data application are too large to include in this repository. To obtain them, visit [this website](https://www.data.gouv.fr/datasets/reanalyses-comephore) and download the precipitation observations from the desired months.
