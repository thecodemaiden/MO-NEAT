MO-NEAT
=======

(Multi-Objective NeuroEvolution of Augmenting Topologies)
---

MO-NEAT is my attempt to adapt the general framework of [NEAT] as described in "Evolving Neural Networks Through Augmenting Topologies", Stanley & Miikkulainen, 2002. 

The most valuable aspect borrowed from this paper is the concept of genetic history, and the method of forming alignable genomes. Recombination and inter-individual distances are made simple by comparing genomes. Niching is also maintained, however not necessarily as the explicit species set up in Stanley & Miikkulainen.

The current implementation contains a variant of the NEAT algorithm mixed with [NSGA-II] ( Deb, Pratap , Agarwal, Meyarivan, 2000), named MONEAT, and another based on [SPEA2] ( Zitzler, Laumanns, Thiele, 2001). There is also a base class, BaseNEAT, that contains the framework for additional NEAT-based algorithms, whether single- or multi-objective.

The algorithms are run on a simple implementation of a neural network, BasicNN. Any network can be output in Graphviz dot format for visualization. The XOR network can be correctly evolved, as a simple test. As a test of multiple objectives, there is also a goal to produce a network with two outputs, one that identifies multiples of 2 and the other that identifies multiples of 3, where the input is a 4 digit binary number.

My next task is to test the performance with time-series input, i.e., recurrent networks.

Examples and results will be available soon.

[NEAT]:http://www.cs.ucf.edu/~kstanley/neat.html
[NSGA-II]:http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.18.7210
[SPEA2]:http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.112.5073
