# Roulette
[A mutation rate model at the basepair resolution identifies the mutagenic effect of Polymerase III transcription](https://doi.org/10.1101/2022.08.20.504670)

Rate estimates in VCF format can be downloaded here: http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/

## Project components
### adding_mutation_rates
Scripts are provided to facilitate annotating variant data with estimates from Roulette and the two models analyzed as comparison.
### analysis_of_mutational_effects

### mutation_rate_model

### population_genetics
Notebooks used to analyze the distribution of allele frequencies (SFS) under recurrent mutation using Roulette estimates.
- Use a demographic model of faster-than-exponential growth to predict the SFS in each Roulette bin.
- Estimate the distribution mutation rates among observed SNVs in RNU, tRNA, and IGK genes.
- Perform simulations to estimate the probability that rare minor alleles are ancestral.

### validation
Notebooks used to estimate the residual variance in mutation rates and compare Roulette estimates to other models.
- Calculate pseudo-R^2 for different models and use bootstrap samples to estimate confidence intervals.
- Compare the *de novo* mutation rates at SNV and non-SNV sites to estimate the residual variance in each mutation rate bin for each model.
