# MSqRobSum

Robust differential protein expression analysis for label-free quantitative proteomics and robust peptide expression summarization

Label-Free Quantitative mass spectrometry based workflows for differential expression (DE) analysis of proteins is often challenging due to peptide-specific effects and context-sensitive missingness of peptide intensities. Peptide-based workflows, like MSqRob, test for DE directly from peptide intensities and outperform summarisation methods which first aggregate MS1 peptide intensities to protein intensities before DE analysis. However, they are computationally expensive, often hard to understand for the non-specialised end-user, and they do not provide protein summaries, which are important for visualisation or downstream processing. We propose a novel summarisation strategy, MSqRobSum, which estimates MSqRob’s model parameters in a two-stage procedure circumventing the drawbacks of peptide-based workflows. MSqRobSum maintains MSqRob’s superior performance, while providing useful protein expression summaries for plotting and downstream analysis. Summarising peptide to protein intensities considerably reduces the computational complexity, the memory footprint and the model complexity. Moreover, MSqRobSum renders the analysis framework to become modular, providing users the flexibility to develop workflows tailored towards specific applications.


# IMPORTANT

- NOTE THAT THE PACKAGE IS NO LONGER MAINTAINED. 
- All functionalities are ported to our novel bioconductor tool [msqrob2](https://www.bioconductor.org/packages/release/bioc/html/msqrob2.html).
