# Markov State Modelling with Optimal Transport

The project cosists of two parts. One regarding the Dixit and Dill approach to building Markov models with optimal transport [1]. The other part is about Optimal Transport Transfer Operators (OTTO), based on the work of Koltai et al. [2]

## Abstract

Markov state models (MSM) provide a framework for analysing molecular dynamics data and model the kinetics of the studied system. Sampling long-time kinetics with molecular simulations is often expensive. A novel approach to MSM inference incorporates information from many short simulations by combining path entropy maximisation (Max Cal) with Optimal Transport Theory (OT). [1] In this Master thesis, we further examine the novel model estimator. We compare the method to the classical Maximum Likelihood (ML) estimator by approximating a given MSM based on limited data. We found the OT-based estimator, in its current form, does not improve ML estimation. However, the connection between MSMs and OT allows for an OT-based approach to analyse molecular simulation data directly. In this approach, an optimal transport plan yields an approximation of the transfer operator, enabling us to circumvent state space decomposition. A recent pre-print introduced the method for coherent set detection [2] but did not consider timescale estimation. We found the approach obtains timescale estimates comparable to those of MSMs.

[1] Dixit, P. D., & Dill, K. A. (2019). Building Markov state models using optimal transport theory. Journal of Chemical Physics, 150(5). https://doi.org/10.1063/1.5086681

[2] Koltai, P., von Lindheim, J., Neumayer, S., & Steidl, G. (2020). Transfer Operators from Optimal Transport Plans for Coherent Set Detection. arXiv preprint arXiv:2006.16085. http://arxiv.org/abs/2006.16085
