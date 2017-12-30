# DR-PID

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<b> Multichannel Double Recursive Frequentist-Bayesian Particle Identification (pion, kaon, proton, electron etc.) for High Energy Physics </b>

A fully working proof of concept code, to be written into a C++ class.
<br/>

The code utilizes preprocessed ionization dE/dx and TOF (time of flight) information in terms of independent Gaussian nsigmas -> independent Gaussian likelihoods, but is extendable for any kind of (multidimensional) likelihood based input. The Expectation Maximization (EM) algorithm is first used for individual particle abundancies and then in second phase for final state channel abundancies (pi+pi-, K+K-, ppbar), differentially in transverse momentum (or total 3-momentum). This is a frequentist fit approach, but also a Bayesian approach is possible via prior distributions. Finally, event by event channel probabilities can be used with different schemes such as Maximum (Posteriori) Probability classification, event by event weights or as cuts, with suitable frequentist ROC-working (false positive, true positive) points.

<br/>

mikael.mieskolainen@cern.ch, 2017
