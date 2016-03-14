# streamgap-pepper

Exploring the effect of impacts from a CDM-like population of
dark-matter subhalos on tidal streams.

This repository contains the code associated with the paper Bovy,
Erkal, \& Sanders (2016, BES16). The ipython notebooks used to
generate the plots in this paper can be found in the top-level
directory of this repository. This code uses
[galpy](https://github.com/jobovy/galpy) and the
[streampepperdf.py](https://gist.github.com/jobovy/1be0be25b525e5f50ea3)
``galpy`` extension, which implements the fast calculation of the
perturbed stream structure.

There are various useful notebooks, in addition to others that were
used during code development that are not described further here.

## 1. [GD1Like-simulation-figures.ipynb](GD1Like-simulation-figures.ipynb)

(render this notebook on [nbviewer](http://nbviewer.ipython.org/github/jobovy/streamgap-pepper/blob/master/GD1Like-simulation-figures.ipynb), where you can toggle the code)

This notebook contains the figures in Section 2 with the properties of
the GD-1-like stream used throughout the paper.

## 2. [meanOperpAndApproxImpacts.ipynb](meanOperpAndApproxImpacts.ipynb)

This notebook contains figures testing the approximations used in the
fast line-of-parallel-angle approach to computing the perturbed stream
structure:

1. The sampling figure that shows that perturbations  in the perpendicular
frequency direction are much smaller than those in the parallel direction.

2. The test of the fast line-of-parallel-angle approach for a single
impact compared to a mock sampling and a direct numerical evaluation.

3. The test of the fast line-of-parallel-angle approach for four
impactscompared to a mock sampling and a direct numerical evaluation.

4. The example perturbed density and frequency tracks for different
mass ranges and full mass range.

5. The scaling of the computational time vs. the number of impacts at
different times.

## 3. [StreamPepperAnalysisGD1Like.ipynb](StreamPepperAnalysisGD1Like.ipynb)

(render this notebook on [nbviewer](http://nbviewer.ipython.org/github/jobovy/streamgap-pepper/blob/master/StreamPepperAnalysisGD1Like.ipynb), where you can toggle the code)

This notebook computes the power spectrum of the density and parallel
frequency for many different simulations with different
parameters. Also contains the bispectrum in these spaces.

## 4. [StreamPepperAnalysisGD1LikeObserved.ipynb](StreamPepperAnalysisGD1LikeObserved.ipynb)

(render this notebook on [nbviewer](http://nbviewer.ipython.org/github/jobovy/streamgap-pepper/blob/master/StreamPepperAnalysisGD1LikeObserved.ipynb), where you can toggle the code)

Power spectra and bispectra in observed space (density and the stream
track's location in the sky, distance, and line-of-sight velocity).

## 5. Simulations

All simulations are run using the code in
[simulate_streampepper.py](simulate_streampepper.py). This script has
a help function that explains its use.

## 6. [CompareNbodySimulations.ipynb](CompareNbodySimulations.ipynb)

(render this notebook on [nbviewer](http://nbviewer.ipython.org/github/jobovy/streamgap-pepper/blob/master/CompareNbodySimulations.ipynb), where you can toggle the code)

This notebook contains comparisons between the simplified action-angle
modeling used in this paper and full N-body simulations with single
and multiple impacts.
