# Rules of Life
Note: Input and output directories are excluded from this repo as they are *very large.*

## Boilerplate

#### include

Libraries. All/most files in this repo include from `include/R/general.R`. This file contains most data parsing and filtering functions.

Other files in directory `R`:
* `data_transform.R`: largely unused older functions for filtering counts and wrapper functions for transforming between log ratio representations
* `GP.R`: code for fitting Gaussian process model via stray::basset
* `metagenomics.R`: functions for Piphillin data and exploratory analyses
* `VC.R`: code for (exploratory) variance component analysis
* `visualization.R`: most plotting functions; also includes `calc_autocorrelation` function

Files in directory `cpp`:
* `dens_optim.cpp`: calculation of (marginal) matrix-t density used in optimization of variance components
* `fast_corr.cpp`: just a fast C++ implementation of Pearson correlation for huge vectors
* `MatDist.h`: stolen (temporarily) from stray package; fast samping functions on which `uncollapse_simple.cpp` depends; included for testing, likely to be removed in future
* `Riemann_dist.cpp`: fast, paralleled calculation of Riemannian distance over pair of PD matrices or set of matrices joined with `cbind`
* `uncollapse_simple.cpp`: posterior sampling for Gaussian process; included for testing, likely to be removed in future

## Exporatory data analysis

#### exploratory_analyses

All code associated with exploratory data analyses.

Of interest:
* `apply_variance_components.R`: Fits variance component model to estimate (rough) contribution of time, season, group, etc. to variation across samples
* `explore_replicates.R`: visualization of between- and across-replicate sample variation
* `find_dense_windows.R`: match "windows" in time (2001-2013) with counts of densely sampled (> some threshold) individuals
* `visualize_autocorrelation.R`: wrapper for autocorrelation visualizations
* `visualize_dense_windows.R`: visualizes the output of `find_dense_windows.R`

## Fitting the model and exploring the posterior fits

#### preprocessing

Data preprocessing wrappers.

Of interest:
* `agglomerate.R`: performs agglomeration at desired taxonomic level using phyloseq
* `apply_filter.R`: filters low abundance taxa into "Other" cohort using specified thresholds
* `formalize_parameters.R`: makes empirical estimates of variance parameters for Gaussian process kernels (currently periodic seasonal kernel and squared exponential daily time kernel)

#### fit_model

Single file `run.R` is a wrapper for `stray::basset()` call to fit Gaussian process.

#### postprocessing

Assortment of posterior visualizations.

Of interest:
* `calculate_posterior_distances.R`: wrapper to calculate Riemannian distance over posterior samples assumed to exist in model\_dir directory specified in `include/R/general.R`
* `embed_posteriors.R`: embeds posterior samples using pre-calculated distances and multidimensional scaling; UMAP version available now
* `find_universal_microbes.R`: plot hosts x interactions heatmaps to visually identify "universal" interactions across hosts and use thresholding (on MAP-estimate interactions) to identify "universal" microbes across hosts (hote: the reality is few/no bacterial genera or families appear universal by this method)
* `visualize_embedding.R`: plots embedding and labels according to specified annotation (e.g. host, group)

## Other analyses

#### metagenomics\_PAM

Code associated with metagenomics sample selection.

Of interest:
* `apply_PAM.R`: filters on well-annotated individuals and calculates K near optimal exemplar samples using partitioning around medoids (PAM) algorithm
* `apply_PAM_replicates.R`: applies PAM algorithm to replicates for sanity checking
* `pull_PAM_samples.R`: pulls PAM-selected sample IDs from STDOUT output of `apply_PAM.R`
* `visualize_PAM_samples.R`: ordinates PAM-selected samples as a sanity check that they span the space of variation

#### fit_joint_model

Code associated with finding "modes" of microbial dynamics using cluster and the elbow method.

Also includes code to fit joint model of the 10 best-sampled hosts via stray::basset. This model estimates a single set of microbial covariances over the included hosts. The purpose of this test was to evaluate the predictive error of the jointly fit version of the model; assuming it wasn't drastically worse that the individually fit models, this gives us an option to evaluate an effective number of individuals by systematic exclusion of individuals in the joint model.

Of interest:
* `cluster_elbow_plot.R`: iterate through K means fits and calculate SS error associated with each K; output to STDOUT
* `fit_individual_models.R`: fits stray::basset model over 10 best sampled individuals with parameter cocktail used by joint model
* `fit_joint_model.R`: fits joint version of stray::basset model over 10 best sampled individuals
* `plot_elbow_scores.R`: parses and plots error associated with a sweep over clusterings

#### testing_simulation

Code associated with model simulations.

Of interest:
* `downsample_real_data.R`: downsamples DUI's data and repeatedly fits stray::basset model so we can evaluate how average estimates change with decreasing sample number
* `downsample_simulated_data.R`: simulate data (from the stray::basset model), fits, and evaluates how average estimates change with decreasing sample number
* `simulate_Riemannian_ATA.R`: simulates the change in Riemannian distance induced by adding a constant amount of noise to a base (random) PD matrix
* `simulate_basset.R`: simulates K hosts over three conditions and visually evaluates accuracy of model fit (1) same baseline composition x different dynamics (2) different baseline compositions x same dynamics (3) same baseline composition x same dynamics (control; sampling/technical noise only model)
