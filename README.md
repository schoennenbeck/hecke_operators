# Hecke operators for unitary groups over imaginary quadratic number fields
## Introduction
This repository contains the code that was used to compute the Hecke operators 
presented in the article [Automorphic Forms on Feit's Hermitian Lattices](https://arxiv.org/abs/1804.05884) (N. Dummigan and S. Schönnenbeck, 2019, preprint). The implementation depends on the
'forms'-package of Markus Kirschmer's which can be found on his [homepage](http://www.math.rwth-aachen.de/~Markus.Kirschmer/) (see 'Magma-packages'). Anybody is welcome to use these algorithms for their
own research as long as the relevant articles are cited.

## Usage
Overall the usability is subpar since most of the algorithms were written for 
one time use. Nevertheless, we believe them to be of general interest. Use at your
own risk. 

Start any computation by attaching the 'forms'-package:

        AttachSpec("path/to/lat.spec");


### Constructing Feit's unimodular lattices

To construct the 20 unimodular lattices from Feit's original article ("Some lattices 
over Q(sqrt(-3))", W. Feit, J. Alg., 1978) simply call

        load "path/to/Magma/FeitLattices.m";

This will construct a list 'Unimod12' containing the Hermitian lattices in question
in the exact order we consider in our article. In addition some intermediate results
are constructed, e.g. you will get lists

        Lat6, Lat8, Lat9, Lat10, and Lat11

containing the indecomposable unimodular lattices in dimension 6, 8, 9, 10, and 11, respecitvely.