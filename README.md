What is NPEB?
======
`NPEB` is a Non-Parametric Empirical Bayes estimation framework for compound estimation in the discrete linear exponential family, which includes a wide class of discrete distributions frequently arising from modern big data applications. The proposed framework directly estimates the Bayes shrinkage factor in the generalized Robbins' formula via solving a scalable convex
program, which is carefully developed based on a RKHS representation of the Stein's discrepancy
measure. The new NEB estimation framework is exible for incorporating various structural constraints into the data driven rule, and provides a unified approach to compound estimation with both regular and scaled squared error losses. 

How to use this repository?
-----
This repository holds the scripts that reproduce the analysis in the paper [1]. Send me an email if anything does not work as expected. If you are looking for the associated R-package `npeb`, please [visit this page](https://github.com/trambakbanerjee/npeb#npeb).

References
=======
[1.] A General Framework for Empirical Bayes Estimation in Discrete Linear Exponential Family _(under review)_     
Banerjee, T., Liu, Q., Mukherjee, G. and Sun, W.

