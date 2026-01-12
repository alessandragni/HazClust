# HazClust: Clustering of Hazards with Frailties

The implementation builds on theoretical developments from

> **Penalized Likelihood Optimization for Adaptive Neighborhood Clustering in Time-to-Event Data with Group-Level Heterogeneity**  
> [arXiv preprint](https://arxiv.org/abs/ADDNUMBERS)
> by A. Ragni, L. Cavinato, F. Ieva  

The core survival and frailty estimation routines
rely on the `parfm` R package [1], which provides
maximum likelihood estimation for parametric frailty models. Our contribution
extends `parfm` by embedding it within an iterative graph-based clustering framework [2,3].


## `R` folder
This folder contains the main functions to run the proposed method, with main function being `HazClust()`.


## `simulations` folder
This folder contains the simulation study described in Section 3 of the paper.


## References
[1] Munda, M., Rotolo, F., & Legrand, C. (2012). parfm: Parametric frailty models in R. Journal of Statistical Software, 51, 1-20.

[2] Nie, F., X. Wang, and H. Huang (2014). Clustering and projected clustering with adaptive neighbors.
In Proceedings of the 20th ACM SIGKDD international conference on Knowledge discovery and data
mining, pp. 977–986.

[3] Liu, C., W. Cao, S. Wu, W. Shen, D. Jiang, Z. Yu, and H.-S. Wong (2020). Supervised graph clustering
for cancer subtyping based on survival analysis and integration of multi-omic tumor data. IEEE/ACM
Transactions on Computational Biology and Bioinformatics 19 (2), 1193–1202
