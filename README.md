## RFASTfMRI
## Fast Adaptive Smoothing and Thresholding

Functional Magnetic Resonance Imaging is a noninvasive tool used to study brain function. Detecting activation is challenged by many factors, and even more so in low-signal scenarios that arise in the performance of high-level cognitive tasks. We provide a fully automated and fast adaptive smoothing and thresholding (FAST) algorithm that uses smoothing and extreme value theory on correlated statistical parametric maps for thresholding. 

## Authors: Israel A Almodovar-Rivera and Ranjan Maitra

## Installation

RFASTfMRI requires

    R version 3.0.0 or higher.
    R package fftw.

The package can be installed via the devtools package:
```R
library(devtools)
install_github("ialmodovar/RFASTfMRI")
```
