# ggr-project
GGR Analysis

# Python requirements
- Anaconda - v5+
- GPy - 1.8.4 (for DP_GP)

# Bioinformatics requirements
- bedtools/2.25.0+
- homer
- DP_GP (Englehardt) - commit: eec12e74219f916aa86e253783905f7b5e30f6f4
- deeptools

# R requirements - v3.4
- seqLogo # check if this is used
- deseq2
- rGREAT
- RDAVIDWebService
- biomaRt
- qvalue
- gplots
- ggplot2
- reshape
- RColorBrewer
- fastcluster
- dendextend # check if this is used
- gridGraphics

NOTE: packages as of 2017-11-08

# Other notes

```
pip install weblogo #3.6.0
pip install deeptools
```

Install homer according to homer instructions on UCSD website

Note that conda curl has issues for deeptools - do conda remove curl to use normal curl, also remove curl-config