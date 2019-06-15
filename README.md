# MSqRobSum
Robust differential protein expression analysis for label-free quantitative proteomics and robust peptide expression summarization

Label-Free Quantitative mass spectrometry based workflows for differential expression (DE) analysis of proteins is often challenging due to peptide-specific effects and context-sensitive missingness of peptide intensities. Peptide-based workflows, like MSqRob, test for DE directly from peptide intensities and outperform summarisation methods which first aggregate MS1 peptide intensities to protein intensities before DE analysis. However, they are computationally expensive, often hard to understand for the non-specialised end-user, and they do not provide protein summaries, which are important for visualisation or downstream processing. We propose a novel summarisation strategy, MSqRobSum, which estimates MSqRob’s model parameters in a two-stage procedure circumventing the drawbacks of peptide-based workflows. MSqRobSum maintains MSqRob’s superior performance, while providing useful protein expression summaries for plotting and downstream analysis. Summarising peptide to protein intensities considerably reduces the computational complexity, the memory footprint and the model complexity. Moreover, MSqRobSum renders the analysis framework to become modular, providing users the flexibility to develop workflows tailored towards specific applications.

This R package allows for both MSqRob and MSqRobSum analysis.
See the vignette for more information.

## Installation

 You first need to install the [devtools](https://cran.r-project.org/package=devtools) package.

```r
install.packages("devtools")
```

Then install saas using the `install_github` function in the
[devtools](https://cran.r-project.org/package=devtools) package.
```r
library(devtools)
install_github("statOmics/MSqRobSum")
```


## If you use this package, please cite:

Robust summarization and inference in proteome-wide label-free quantification  
Adriaan Sticker, Ludger Goeminne, Lennart Martens, Lieven Clement  
bioRxiv 668863; doi: https://doi.org/10.1101/668863 
