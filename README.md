# InpaintingAutoregressive

This is the accompanying repository for the article *On the Use of Autoregressive Methods for Audio Inpainting* authored by Ondřej Mokrý and Pavel Rajmic, submitted to EUSIPCO 2024.

> The paper presents an evaluation of popular audio inpainting methods based on autoregressive modelling, namely, extrapolation-based and Janssen methods. A novel variant of the Janssen method suitable for gap inpainting is also proposed. The main differences between the particular popular approaches are pointed out, and a mid-scale computational experiment is presented. The results demonstrate the importance of the choice of the AR model estimator and the suitability of the new gap-wise Janssen method.

The preprint is available at [arXiv](http://arxiv.org/abs/2403.04433).

## Contents

The repository includes the MATLAB source codes needed to reproduce the research

- **references** – Implementation of reference methods, namely SPAIN and SPAIN-MOD, with the help of [InpaintingRevisited](https://github.com/ondrejmokry/InpaintingRevisited) and [Dictionary learning for sparse audio inpainting](https://www.oeaw.ac.at/isf/forschung/fachbereiche-teams/mathematik/dictionary-learning-for-sparse-audio-inpainting).
- **utils** – All the functions needed to run the main files, except for the codes for the [Psychoacoustically motivated evaluationk](#psychoacoustically-motivated-evaluation).
- `gaps_table.mat` - Source signals and masks, taken from [TestSignals repository](https://github.com/ondrejmokry/TestSignals).
- `maintest.m` - Main code running the test of the AR-based methods.
- `maintest_spain.m` - Main code running the test of the SPAIN variants.

## Psychoacoustically motivated evaluation

Note that the codes to compute the psychoacoustically motivated metrics (PEMO-Q, PEAQ) are not provided.

The PEAQ package ican be acquired from [TSP Lab of McGill University](http://www-mmsp.ece.mcgill.ca/Documents/Software/). The PEMO-Q software is no longer publicly available. For that reason, and because the processing is very time- and spacedemanding, the provided .mat files with the results include all the data precomputed.
