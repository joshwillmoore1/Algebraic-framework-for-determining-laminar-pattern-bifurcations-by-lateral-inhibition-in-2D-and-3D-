# Algebraic framework for determining laminar pattern bifurcations by lateral inhibition in 2D and 3D

This repository contains the code for both the static and dynamic simulations shown in the paper entitled: "Algebraic framework for determining laminar pattern bifurcations by lateral inhibition in 2D and 3D".

Static simulations:
The static simulations were conducted using MATLAB 2019b and are performing parameter sweeps over the signal strength parameters \\w_{1}\\ and \\w_{2}\\. The static simulation folders contains the parameter sweeps for the 2D and 3D lattice geometries considered.

Dynamic simulations:
The dynamic simulations were conducted using Chaste (Cancer, Heart and Soft Tissue Environment) v2019.1. The dynamic simulations folder contains both src and test files required to run the Collier et al. 1996 NDM with the fixed and adaptive signal strength mechanisms in both 2D and 3D. These files should be included within a project initiated by the user.
