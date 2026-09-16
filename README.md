# InpaintingAutoregressive

This repository contains the code, data, results, and supplementary material for two papers on autoregressive inpainting of gaps in audio signals.

## Papers

### Tweaking autoregressive methods for inpainting of gaps in audio signals

Ondřej Mokrý and Pavel Rajmic, presented at EUSIPCO 2025.

> A novel variant of the Janssen method for audio inpainting is presented and compared to other popular audio inpainting methods based on autoregressive (AR) modeling. Both conceptual differences and practical implications are discussed. The experiments demonstrate the importance of the choice of the AR model estimator, window/context length, and model order. The results show the superiority of the proposed gap-wise Janssen approach using objective metrics, which is confirmed by a listening test.

The preprint is available at [arXiv](http://arxiv.org/abs/2403.04433), the official version is published at [IEEE Xplore](https://ieeexplore.ieee.org/document/11226154).

- [Paper-specific code and results](papers/tweaking-autoregressive/README.md)
- [Supplementary material](https://ondrejmokry.github.io/InpaintingAutoregressive/)

### Reviving Etter method for autoregressive inpainting: Generalization, evaluation, implementation

Ondřej Mokrý, Matěj Hrdlička and Pavel Rajmic, submitted to ICUMT 2026.

> Audio inpainting aims to restore missing segments in an audio waveform, as encountered in dropouts and packet losses. This paper revisits the autoregressive (AR) interpolation method proposed by Etter, which combines forward and backward AR prediction through a structured linear system, yet lacks a widely used full implementation for audio signals. We provide an open implementation and propose two practical extensions: A formulation that allows high-order AR models even when the gap length is shorter than the model order, and a causal variant tailored to packet loss concealment. The method is evaluated on musical excerpts with gaps up to 80 ms, using SNR and perceptually motivated objective difference grades, complemented by a listening test. The results show that Etter inpainting matches the performance of most AR- and sparsity-based baseline methods, with the exception of the iterative, gap-wise Janssen approach.

- [Paper-specific code and results](papers/reviving-etter/README.md)
- [Supplementary material](https://ondrejmokry.github.io/InpaintingAutoregressive/)

The implementation of the Etter method is contained in `papers/reviving-etter/Etter.m`; its experiment is run with `maintest_Etter.m`. Both use the shared data and dependencies described below.

## Contents

## Shared repository contents

The repository includes the MATLAB source codes and data needed to reproduce the research:

- **utils** – Shared functions used by the experiments.
- **references** – Implementation of reference methods, namely SPAIN and SPAIN-MOD, with the help of [InpaintingRevisited](https://github.com/ondrejmokry/InpaintingRevisited) and [Dictionary learning for sparse audio inpainting](https://www.oeaw.ac.at/isf/forschung/fachbereiche-teams/mathematik/dictionary-learning-for-sparse-audio-inpainting).
- **irmas** – The IRMAS file list and preprocessing script used by the mid-scale experiment.
- `gaps_table.mat` – Source signals and masks, taken from [TestSignals repository](https://github.com/ondrejmokry/TestSignals).
- `papers/` – Paper-specific runners, plotting scripts, and result files.

Regarding the mid-scale experiment using the IRMAS dataset, the original files can be downloaded [here](https://www.upf.edu/web/mtg/irmas).

For supplementary material, including plots and audio examples, see the [accompanying website](https://ondrejmokry.github.io/InpaintingAutoregressive/).

## Psychoacoustically motivated evaluation

Note that the codes to compute the psychoacoustically motivated metrics (PEMO-Q, PEAQ) are not provided.

The PEAQ package can be acquired from [TSP Lab of McGill University](http://www-mmsp.ece.mcgill.ca/Documents/Software/), the PEMO-Q software is available through [University of Oldenburg](https://uol.de/en/mediphysics/downloads/pemo-q). Because the processing is very time- and space-demanding, the provided .mat files with the results include all the data precomputed.

## Further dependencies

The experiments for *Tweaking autoregressive methods for inpainting of gaps in audio signals* were run in Matlab R2023a using Signal Processing Toolbox (Version 9.2), Parallel Computing Toolbox (Version 7.8) and [Large Time-Frequency Analysis Toolbox](http://ltfat.org/) (Version 2.4.0).

For running SPAIN with dictionary learning as a reference method, the [CVX toolbox](http://cvxr.com/cvx/) is necessary.

The experiments for *Reviving Etter method for autoregressive inpainting: Generalization, evaluation, implementation* were run in Matlab R2026a using Signal Processing Toolbox (Version 26.1).