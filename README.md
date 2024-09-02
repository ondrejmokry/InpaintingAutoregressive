# InpaintingAutoregressive

This is the accompanying repository for the article *On the Use of Autoregressive Methods for Audio Inpainting* authored by Ondřej Mokrý and Pavel Rajmic, submitted to ICASSP 2025.

> The paper presents an evaluation of popular audio inpainting methods based on autoregressive modelling, namely, extrapolation-based and Janssen methods. A novel variant of the Janssen method suitable for gap inpainting is also proposed. The main differences between the particular popular approaches are pointed out, and a mid-scale computational experiment is presented. The results demonstrate the importance of the choice of the AR model estimator and the suitability of the new gap-wise Janssen method.

The preprint is available at [arXiv](http://arxiv.org/abs/2403.04433).

## Contents

The repository includes the MATLAB source codes needed to reproduce the research

- **plotting** - Matlab scripts that load the results and replicate the figures used in the paper + some more.
- **references** – Implementation of reference methods, namely SPAIN and SPAIN-MOD, with the help of [InpaintingRevisited](https://github.com/ondrejmokry/InpaintingRevisited) and [Dictionary learning for sparse audio inpainting](https://www.oeaw.ac.at/isf/forschung/fachbereiche-teams/mathematik/dictionary-learning-for-sparse-audio-inpainting).
- **results** - `.mat` files with all the SDR and ODG values from testing the methods. Recovered signals are *not* included due to file size limits. However, these can be reproduced using the `maintest.m` and `maintest_spain.m` scripts.
- **utils** – All the functions needed to run the main files, except for the codes for the [Psychoacoustically motivated evaluation](#psychoacoustically-motivated-evaluation).
- `gaps_table.mat` - Source signals and masks, taken from [TestSignals repository](https://github.com/ondrejmokry/TestSignals).
- `maintest.m` - Main code running the test of the AR-based methods.
- `maintest_spain.m` - Main code running the test of the SPAIN variants.

## Psychoacoustically motivated evaluation

Note that the codes to compute the psychoacoustically motivated metrics (PEMO-Q, PEAQ) are not provided.

The PEAQ package ican be acquired from [TSP Lab of McGill University](http://www-mmsp.ece.mcgill.ca/Documents/Software/), the PEMO-Q software is available through [University of Oldenburg](https://uol.de/en/mediphysics/downloads/pemo-q). Because the processing is very time- and spacedemanding, the provided .mat files with the results include all the data precomputed.

## Further dependencies

The experiments were run in Matlab R2023a using Signal Processing Toolbox (Version 9.2), Parallel Computing Toolbox (Version 7.8) and [Large Time-Frequency Analysis Toolbox](http://ltfat.org/) (Version 2.4.0).

If you want to try also the variant of SPAIN with dictionary learning as a reference method, the [CVX toolbox](http://cvxr.com/cvx/) is necessary.

## Statistical testing of LPC versus Burg algorithm

For the assessment of statistical significance of the effect of the AR model estimator (LPC versus Burg algorithm), the [Wilcoxon signed rank test](https://www.mathworks.com/help/stats/signrank.html) was employed using the following hypotheses:
- H0: Burg algorithm and LPC lead to results with the same median,
- HA: Burg algorithm leads to results with higher median.
This was performed separately for each inpainting method, AR model order, and evaluation metric (SDR, PEMO-Q ODG).

The p-value displayed in the tables below indicate the rejection of the null hypothesis, i.e., p-value < 0.05 implies that the data feature enough evidence to reject the equality of the medians in favor of the alternative hypothesis HA (at the significance level of 5%). On the other hand, p-value > 0.05 means the test is inconclusive.

| evaluation by SDR     | 256      | 512      | 1024     | 2048     | 3072     |
|-----------------------|----------|----------|----------|----------|----------|
| extrapolation-based   | 1.06e-21 | 2.12e-22 | 2.12e-22 | 2.12e-22 | 3.58e-20 |
| Janssen, gap-wise     | 2.16e-09 | 3.89e-13 | 5.54e-10 | 0.03     | 0.97     |
| Janssen, Hann window  | 7.82e-05 | 6.19e-4  | 0.90     | 1.00     | 1.00     |
| Janssen, rect. window | 6.69e-11 | 4.67e-06 | 0.39     | 0.99     | 1.00     |

| evaluation by ODG     | 256      | 512      | 1024     | 2048     | 3072     |
|-----------------------|----------|----------|----------|----------|----------|
| extrapolation-based   | 6.35e-22 | 2.12e-22 | 2.12e-22 | 2.12e-22 | 3.57e-18 |
| Janssen, gap-wise     | 8.94e-13 | 1.25e-18 | 1.86e-20 | 3.42e-16 | 5.10e-06 |
| Janssen, Hann window  | 0.20     | 0.09     | 0.06     | 1.00     | 1.00     |
| Janssen, rect. window | 2.11e-22 | 1.84e-18 | 0.02     | 0.99     | 0.99     |
