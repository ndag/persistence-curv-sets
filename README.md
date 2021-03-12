# Curvature sets over persistence diagrams
This repository contains the implementation of some of the ideas in our paper https://arxiv.org/pdf/2103.04470.pdf. 

The implementation below is via Matlab's *parfor*.

The software was developed by Mario Gomez and Facundo Memoli. We also use functions written by Misha Belkin (L2_distance.m), Vin de Silva (px_fps.m), and Uli Bauer (Ripser: https://github.com/Ripser/ripser). We have also benefited from Chris Tralie's Matlab wrapper for Ripser (https://github.com/ctralie/Math412S2017).
![Alt text](stack-D-to1.png?raw=true "Stacking of persistence diagrams")


# Installation Instructions:
<!--- 1. From a terminal go to the ripser folder and run `make ripser-coeff`.
- Add support for a function that samples from a continuous metric space instead of a discrete one.
--->
This code was developed in Matlab version R2019a. It requires the Parallel Computing Toolbox (https://www.mathworks.com/products/parallel-computing.html).

## Example
We have included a file dX_example_circle.mat with 1000 points sampled uniformly at random from the unit circle. As long as you have the Parallel Computing Toolbox installed, the script should run out of the box. You can further personalize it as described in the **Usage** section.

## Usage
The input is the distance matrix $dX$ of a large metric space $X$ stored in a .mat file. The program will load this matrix and sample metric subspaces from it in order to compute $D_{n,k}^{\mathrm{VR}}(X)$. You need to set the following parameters:

1. Store a distance matrix in a .mat file. In order to guarantee meaningful results, the matrix should be square, symmetric, have diagonal 0, and satisfy the triangle inequality.
	- Keep in mind that a reasonable size for `dX` is at most 5000x5000, depending on your system.
2. Edit the file Persistence_parsave_dX.m to set the following *book-keeping* variables:
	- `filename_metric_space`: Name (and path) of the .mat file containing the distance matrix.
	- `results_file`: The name of the file where the results will be stored.
	- `save_to_file`: Boolean flag to decide if you want to save the output to results_file.mat or if you only want to keep it in Matlab's workspace.
	- `temp_file`: The ending of the file where the temporary results will be stored. The files are named checkpoints_temp_file_#.mat, where # runs from 1 to the number of cores utilized. This variable is usually not edited, unless you are running several copies of the script at the same time.
	- `nCores`: The number of cores that you want to use. The maximum is the number of physical cores (not logical cores) available in your system. You can comment out that line to use the default chosen by Matlab.
		- Matlab can only work with physical cores. See https://www.mathworks.com/matlabcentral/answers/80129-definitive-answer-for-hyperthreading-and-the-parallel-computing-toolbox-pct#answer_89845 for a more detailed explanation. In short, logical cores require hyperthreading, but this negatively impacts Matlab's performance more often than not.
3. Set the simulation parameters:
	- `dim`: The dimension of homology that you want to calculate. The program will automatically set `n=2*dim+2` to calculate the *principal* persistence set.
	- `nReps`: Number of $n$-point samples to take from $dX$.
4. Once you've decided on the above parameters, run Persistence_parsave_dX.m from the Matlab command window. This will produce a .mat file containing the results.

## Output
The script produces a graph of $D_{n,k}^\mathrm{VR}(X)$, as described in the paper and in the figure above. It also saves the following variables in results_file.mat:

* `confs`: Matrix of size `[n, nReps]`. Each column `I = confs(:,i)` is the set of indices a set of $n$ points sampled from $X$.
* `dms`: Array of size `[n, n, nReps]`. Each page `dm = dms(:,:,i)` is the distance matrix of a sample of $n$ points from $X$. It is obtained from `dX` and `confs` by setting `I=confs(:,i)` and `dm=dX(I,I)`.
* `bd_times`: An array of size `[nReps, 2]` with the persistence diagrams of the samples. The first column `bd_times(:,1)` contains the birth times and the second, `bd_times(:,2)`, the death times. A row `bd_times(i,:)` will be non-zero if, and only if, the configuration with distance matrix `dms(:,:,i)` produced persistence homology.

# TO DO
- Do we want to save the results to a .mat file directly from Persistence_parsave_dX.m, or should we gather the results distributed in the checkpoint_temp_file_#.mat files?
	- In any case, it would be useful to upload the code that gathers the results from the checkpoints.
- Write a Ripser version of Persistence_parfor.m. The use of parfor introduced bugs unseen in the sequential version, but this will allow us to compute more general persistence sets, where $n > 2k+2$.
- Set a random seed. This is not straightforward with parallel code because each worker has its own random seed.
