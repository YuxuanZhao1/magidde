# MAGIDDE: Inference for Delay Differential Equations Using Manifold-Constrained Gaussian Processes

## Install
Use `install.packages()` to install the `magidde` R package:

`install.packages("magidde_0.1.tar.gz", repos = NULL, type="source")`.

## Simulation Study

### Numerical Validation of Approximation Schemes

The R scripts are located in `Simulation_study/Rscript_numerical_validation`.

We investigate three schemes to handle the historical outputs, including:
1. Linear interpolation: folder `Linear_interpolation/`, run `Linear_interpolation.R`
2. Conditional expectation: folder `Conditional_expectation/`, run `conditional.R`
3. Fully Bayesian construction: folder `Fully/`, run `fully_specified.R`

### Method Comparison

The R scripts are located in `Simulation_study/Rscript_method_comparison`.

The following four methods are implemented, using Hutchinson's equation as a benchmark system:
1. Our proposed approach (MAGIDDE): `MAGI.R`
2. debinfer [1]: `debinfer.R`
3. SMCDE [2]: folder `smcde/`. See `smcde/README`. 
4. NLS: `nls_Hutchinson.R`.

Furthermore, run `MAGI_discretization_level_change.R` to investigate the effect of the denseness of discretization points $\boldsymbol|I|$ using MAGIDDE.

### $Lac$ Operon System
The R scripts are located in `Simulation_study/Rscript_Lac_Operon`

We compare the following methods:
1. NLS: `nls_lac.R`
2. Our proposed approach (MAGIDDE): `lacoperon_dde_MAGI.R`


## Application

The R scripts and processed data for the SIRD model are in the folder `Application/`.

+ The processed Ontario COVID data for Jan. 24 2022 to Feb. 23 2022 are provided in `sir_data.csv`.
+ To conduct parameter inference using the 30 days of observations referenced in the main text (Section 5), run the first portion of `ird_model.R`.
+ To fit the model using the observations from the first period only ($t = \{0,1,\cdots, 15\}$) and generate predictions for the second period ($t = \{16,17,\cdots, 30\}$), run the second portion of `ird_model.R`.  
+ To demonstrate how to generate predictions with MAGIDDE (using simulated data under a correctly specified model), run `simulate_data_prediction.R`.


## References
<b id="1">[1]</b> 
P. H. Boersch-Supan, S. J. Ryan, and L. R. Johnson. deBInfer:
Bayesian inference for dynamical models of biological systems in R.
Methods in Ecology and Evolution, 8(4):511–518, 2017.

<b id="2">[2]</b> 
S. Wang, S. Ge, R. Doig, and L. Wang. Adaptive semiparametric
Bayesian differential equations via sequential Monte Carlo. Journal
of Computational and Graphical Statistics, 31(2):600–613, 2022

