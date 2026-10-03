# quantum_sequential_sampler_public

## Description
Public Repository for the Quantum Sequential Sampler Model.

## Repository structure 

```
├── models: matlab codes for all the models
├── model_BIC_fits: BIC_fits
├── model_generalization_test: generalization_tests
├── data: the 1162-subject dataset for probablity judgment operated on a rating scale + Zhu et al. (2020)'s dataset
├── figure: figures in the publication (to be updated soon)
└── analysis_code: analysis codes for empirical result, and model comparisons
└── BS dataset: models and analysis codes for the dataset in Zhu et al. (2020)
```
## Note
The model fits uploaded before 2026 came from an older version of the model code that contained a small error. The error has been fixed in the updated code, and the fits have been rerun with it. The difference in the results is subtle (1/300 of the fitted likelihood).

## Reference
Huang, J., Busemeyer, J., Ebelt, Z., & Pothos, E. (2023). Quantum Sequential Sampler: a dynamical model for human probability reasoning and judgments. In Proceedings of the Annual Meeting of the Cognitive Science Society (Vol. 45, No. 45).
<br />
Zhu, J. Q., Sanborn, A. N., & Chater, N. (2020). The Bayesian sampler: Generic Bayesian inference causes incoherence in human probability judgments. Psychological review, 127(5), 719.
<br />
Huang, J., Busemeyer, J. R., Ebelt, Z., & Pothos, E. M. (2025). Bridging the Gap Between Subjective Probability and Probability Judgments: The Quantum Sequential Sampler. Psychological Review, 132(4), 916-955.

