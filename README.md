# Basilisk Starter Pack

This repository contains a set of examples to get you started with fluid simulations in Basilisk. These examples contain visualization and data output choices that we often use in the FluidLab, for example, outputting VTK files and visualizing in Paraview.

## Pre-requisites

You will need an installation of Basilisk in your machine. For the post-processing examples, you will also need an installation of Python. 

## Example 1: Shear flow of a Saramito fluid

This example simulates the flow of a Saramito fluid under simple shear. We output the polymeric stress over time and compare to a semi-analytical solution. VTK files are also created and these can be visualized in Paraview.

### Running a simulation
To run this example, the first thing to do is compile your code. You can do this using the Basilisk compiler via the following command:  
> qcc shear_evp.c -lm -o shear_evp
 
