# BenFound (a.k.a. "Let's see your digits: Anomalous-state detection using Benford's Law")

Thanks for visiting! This repository contains the Knitr source for the paper "Let's see your digits: Anomalous-state detection using Benford's Law" by Samuel Maurus and Claudia Plant (currently under review).

Please familiarize yourself with the concepts of Knitr before you start here. Importantly, the idea is that all code and experiments required to reproduce the results in the paper is included in the source of the paper itself. This is done by knitr, which is something like a precompiler for LaTex which includes the ability to embed R code.

Reproducing our results is thus "theoretically" only a matter of building our document using knitr. This will execute all experiments on your local machine, generate the relevant figures and a *.tex file for you to build using e.g. pdflatex.

The knitr document that you need to build is: bilb.Rnw.

We say "theoretically" because some results, especially those to do with running the comparison techniques in other languages (e.g. C++ or Java), can't be embedded using knitr (some would also be too intractable due to computational effort). In these cases, we have included detailed information in bilb.Rnw directly in the relevant sections of the paper. Open bilb.Rnw in a text editor and you will find all the embedded R code, and detailed instructions for reproducing the results that were not able to be done in R.

The data that is used by the paper is in the data directory. The R code that can be used independently of the paper for your own experiments is in func.R. Be sure to look at load.R, and also load it using source("load.R"). It contains a bunch of R packages that are dependencies of this work.

In func.R, all functions are commented. The most important for this work would be the function extractKSSignal, which extracts the "Benfordness" signal for a given time series vector.
