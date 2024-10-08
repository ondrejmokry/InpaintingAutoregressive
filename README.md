# InpaintingAutoregressive

This is the accompanying repository for the article *On the Use of Autoregressive Methods for Audio Inpainting* authored by Ondřej Mokrý and Pavel Rajmic, submitted to ICASSP 2025.

> The paper presents an evaluation of popular audio inpainting methods based on autoregressive modeling, namely the extrapolation-based and Janssen methods. A novel variant of the Janssen method suitable for inpainting of gaps is also proposed. The main differences between the particular popular approaches are pointed out. In the experimental part, the importance of the choice of the AR model estimator is confirmed by objective metrics. Then, a mid-scale computational experiment is presented, and its results are confirmed by a listening test. All the experiments demonstrate the superiority of the new gap-wise Janssen method.

The preprint is available at [arXiv](http://arxiv.org/abs/2403.04433).

## Contents

The repository includes the MATLAB source codes needed to reproduce the research:

- **plotting** – Matlab scripts that load the results and replicate the figures used in the paper + some more.
- **references** – Implementation of reference methods, namely SPAIN and SPAIN-MOD, with the help of [InpaintingRevisited](https://github.com/ondrejmokry/InpaintingRevisited) and [Dictionary learning for sparse audio inpainting](https://www.oeaw.ac.at/isf/forschung/fachbereiche-teams/mathematik/dictionary-learning-for-sparse-audio-inpainting).
- **results** – `.mat` files with all the SDR and ODG values from testing the methods. Recovered signals are *not* included due to file size limits. However, these can be reproduced using the `maintest.m` and `maintest_spain.m` scripts.
- **utils** – All the functions needed to run the main files, except for the codes for the [Psychoacoustically motivated evaluation](#psychoacoustically-motivated-evaluation).
- `gaps_table.mat` – Source signals and masks, taken from [TestSignals repository](https://github.com/ondrejmokry/TestSignals).
- `maintest.m` – Main code running the test of the AR-based methods.
- `maintest_spain.m` – Main code running the test of the SPAIN variants.

For supplementary material (graphs, audio), see the [accompanying website](https://ondrejmokry.github.io/InpaintingAutoregressive/).

## Psychoacoustically motivated evaluation

Note that the codes to compute the psychoacoustically motivated metrics (PEMO-Q, PEAQ) are not provided.

The PEAQ package can be acquired from [TSP Lab of McGill University](http://www-mmsp.ece.mcgill.ca/Documents/Software/), the PEMO-Q software is available through [University of Oldenburg](https://uol.de/en/mediphysics/downloads/pemo-q). Because the processing is very time- and space-demanding, the provided .mat files with the results include all the data precomputed.

## Further dependencies

The experiments were run in Matlab R2023a using Signal Processing Toolbox (Version 9.2), Parallel Computing Toolbox (Version 7.8) and [Large Time-Frequency Analysis Toolbox](http://ltfat.org/) (Version 2.4.0).

For running SPAIN with dictionary learning as a reference method, the [CVX toolbox](http://cvxr.com/cvx/) is necessary.
