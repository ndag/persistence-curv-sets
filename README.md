# Curvature sets over persistence diagrams
Implementation of some of the ideas in our paper https://arxiv.org/pdf/2103.04470.pdf

The software was developed by Mario Gomez and Facundo Memoli. We also use functions written by Misha Belkin (L2_distance.m), Vin de Silva (px_fps.m), and Uli Bauer (Ripser: https://github.com/Ripser/ripser). We have also benefited from Chris Tralie's Matlab wrapper for Ripser.
![Alt text](stack-D-to1.png?raw=true "Stacking of persistence diagrams")


## Installation Instructions:
1. From a terminal go to the Ripser folder and type 'make ripser-coeff'
2. Check your Matlab distribution. You should install the Parallel Computing Toolbox https://www.mathworks.com/products/parallel-computing.html
3. Edit the file Persistence_Parsave.m. Pay attention to the different flags and variables such as **n** and **d** and **ncores**. Note that **ncores** determines the number of workers that Matlab will recruit in order to run  *parfor*. This will depend on your system. 
4. Once you've decided on the above paramters, run Persistence_Parsave.m from the Matlab command window. This will produce a .mat file containing the results.
