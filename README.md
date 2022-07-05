PrebioDivSim
====

Simulation model for prebiotic replicator (like RNAs) dynamics and evolution.

This simulation model is constructed mimicking translation-coupled RNA relpication system (TCRR system; Ichihashi _et al_., _Nat_. _Commun_., 2013).

## Simulations

This repository (PrebioDivSim) includes two simulations, `PrebioDivPreserv` and `PrebioDivEvo`.

`PrebioDivPreserv` is a simulation for determining whether compartmentalized replication systems are maintained or not.
Firstly a replicator network and its activity matrix are given, then the population dynamics is simulated fixed them.

`PrebioDivEvo` is a simulation for determining how compartmentalized replication systems evolved.

In these simulations, the concentrations of replicators in each compatment are defined as a type `Lives`. 
The entity of this type is `std::vector<double>`.
And, the activities of replications are defined as a type `Act`.
The entity of this type is `std::vector<std::vector<double>>`, then a matrix.

The entity of simulation process is a class `World`. For simulation, environmental parameters are unified as a class `Env`.

## Dependency

To use these simulations, `openmp` and `boost` must be installed in environment.

And, their operations are confirmed in `g++11` and `c++17 and later`.

## Licence

Simplified BSD License (read file `LICENSE`)

## Reference
[1] Kamiura, R., Mizuuchi, R. & Ichihashi, N. Plausible pathway for a host–parasite molecular replication network to increase its complexity through Darwinian evolution. bioRxiv 2022.01.17.476531 (2022).

[2] Ichihashi, N. et al. Darwinian evolution in a translation-coupled RNA replication system within a cell-like compartment. Nat. Commun. 4, (2013).