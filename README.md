# DR-PID

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A fully working proof of concept code, to be written into a C++ class.
<br/>

Multichannel Double Recursive Frequentist-Bayesian Particle Identification (pion, kaon, proton, electron etc.) in High Energy Physics. The code utilizes preprocessed dE/dX and TOF (time of flight) information in terms of independent Gaussian nsigmas -> independent Gaussian likelihoods, but is extendable for any kind of (multidimensional) likelihood based information. The Expectation Maximization (EM) algorithm is first used for individual particle abundancies and then in second phase for final state channel abundancies (pi+pi-, K+K-, ppbar), differentially in transverse momentum. Event by event channel probabilities can be used with different schemes such as Maximum (Posteriori) Probability classification, event by event weights or as cuts.


mikael.mieskolainen@cern.ch, 2017
