# HazClust: Clustering of Hazards with Frailties

The implementation builds on theoretical developments from

> **Penalized Likelihood Optimization for Adaptive Neighborhood Clustering in Time-to-Event Data with Group-Level Heterogeneity**  
> [arXiv preprint](https://arxiv.org/abs/2601.07446)
> by A. Ragni, L. Cavinato, F. Ieva  

The core survival and frailty estimation routines
rely on the `parfm` R package [1], which provides
maximum likelihood estimation for parametric frailty models. Our contribution
extends `parfm` by embedding it within an iterative graph-based clustering framework [2,3].


## `R` folder
This folder contains the main functions to run the proposed method, with main function being `HazClust()`.


## `simulations` folder
This folder contains the simulation study described in Section 3 of the paper.
 
- `Tab1_Silhouette.ipynb` reproduces the results reported in Table 1 (Section 3.1.1). 
   The outputs are stored in the `output/tab1` folder.  

   To generate the results, run:
   ```bash
   chmod +x submit_Tab1_simulations.sh
   ./submit_Tab1_simulations.sh
   ```
   This script relies on enhanced_simulation.job, enhanced_simulation_main.R, and enhanced_simulate_data.R.
  
- `Tab2_AnalysisAccuracyARI.ipynb` reproduces the results reported in Table 2 (referenced in Section 3.1.2). 
   The outputs are stored in the `output/tab2` folder.
  
   To generate the results, run:
   ```bash
   chmod +x submit_Tab2_simulations.sh
   ./submit_Tab2_simulations.sh
   ```
   This script relies on enhanced_simulation.job, enhanced_simulation_main.R, and enhanced_simulate_data.R.
  
   
- `TabS61_AccuracyARISilhouette.ipynb` reproduces the results reported in Table S6.1 (referenced in Section S6.1). 
   The outputs are stored in the `output/tabS61` folder.
  
   To generate the results, run:
   ```bash
   chmod +x submit_TabS61_simulations.sh
   ./submit_TabS61_simulations.sh
   ```
   This script relies on enhanced_simulation.job, enhanced_simulation_main.R, and enhanced_simulate_data.R.
 
- `TabS62_Figs_Variabandestimation.ipynb` reproduces 
   - the descriptive results in Section 3.1.1
   - Table S6.2 (referenced in Web Supplementary Materials), 
   - Figures 1 and S6.1 in Web Supplementary (both referenced in Section 3.1.3).
   
   This notebook relies on the outputs stored in the `output/tab2` folder.  

## `casestudy` folder
This folder contains the main code used for the case study described in Section 4 of the paper, along with the corresponding output and results presented in the manuscript.

## References
[1] Munda, M., Rotolo, F., & Legrand, C. (2012). parfm: Parametric frailty models in R. Journal of Statistical Software, 51, 1-20.

[2] Nie, F., X. Wang, and H. Huang (2014). Clustering and projected clustering with adaptive neighbors.
In Proceedings of the 20th ACM SIGKDD international conference on Knowledge discovery and data
mining, pp. 977–986.

[3] Liu, C., W. Cao, S. Wu, W. Shen, D. Jiang, Z. Yu, and H.-S. Wong (2020). Supervised graph clustering
for cancer subtyping based on survival analysis and integration of multi-omic tumor data. IEEE/ACM
Transactions on Computational Biology and Bioinformatics 19 (2), 1193–1202
