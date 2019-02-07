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

### Computing Hecke operators

For computing Hecke operators call

        load "path/to/Magma/Hecke.m";

This provides you with two functions

        IntertwiningOperator(GenusOfL, P);
        HeckeOperator(L,P);

The first one takes as argument a list of representatives of a genus of lattices 
and an ideal P and computes the intertwining operator S_P (see Proposition (3.5)
of the article) as well as a list of representatives of the corresponding secondary genus 
of lattices. These can be used in a straightforward manner to compute the Hecke operator
of the genus at the prime P. However, we also provide a convenience function
'HeckeOperator' that takes only a lattice L and a prime P and computes (in this order):
* The neighbouring operator of the genus of L at P
* The intertwining operator S_P
* The corresponding operator S_P'
* Representatives for the genus of L
* Representatives for the corresponding secondary genus

It is the user's responsibility to check that the conditions of Proposition (3.5) are 
fulfilled when using these functions as they will also run in other cases in which 
the results cannot be trusted.

The function 'HeckeOperator(L,P)' can be made quite a bit faster if we already know
as system of representatives of the genus of L (or the second genus) and the corresponding
automorphism groups. It should be relatively easy to adjust the code accordingly.

## Bugs
It is highly likely that these algorithms contain bugs. In case you find one please 
feel free to let me know or - if you manage to fix it yourself - to submit a pull request.