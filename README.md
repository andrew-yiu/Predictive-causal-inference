# Predictive causal inference

This repository contains R code for running the methods in the paper "Causal predictive inference and target trial emulation" by Andrew Yiu, Edwin Fong, Stephen Walker, and Chris Holmes. The preprint can be found here.

The procedures use the **bartMachine** R package, which requires Java. Instructions for installing **bartMachine** can be found [here](https://github.com/kapelner/bartMachine). Java necessitates the specification of memory allocation before running the code, which is achieved through the command **options(java.parameters = )**. A default value has been provided, but the user should adjust this parameter accordingly if the memory is insufficient.


