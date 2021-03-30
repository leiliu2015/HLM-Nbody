## HLM-Nbody

HLM-***Nbody*** is a polymer physics-based analysis method to predict chromatin *n*-body contact probability (*n>2*) from Hi-C data. As an updated version of [HLM](https://github.com/leiliu2015/HLM), we optimized the numeric algorithm by considering [PHi-C](https://github.com/soyashinkai/PHi-C) and the following two papers:
- G. Le Treut, F. Képès, and H. Orlands, [A polymer model for the quantitative reconstruction of chromosome architecture from HiC and GAM data](http://dx.doi.org/10.1016/j.bpj.2018.10.032), Biophys. J. 115, 2286-2294 (2018).
- P. Metha, M. Bukov, C.-H. Wang, A. G. R. Day, C. Richardson, C. K. Fisher, and D. J. Schwab, [A high-bias, low-variance introduction to Machine Learning for physicists](https://www.sciencedirect.com/science/article/pii/S0370157319300766?via%3Dihub), Phys. Rep. 810, 1-124 (2019).

Necessary codes written in Python and [Gnuplot](gnuplot.sourceforge.net) scripts to reproduce most results in our recent [work]() are archived in this repository, which have been tested on ubuntu 16.04/18.04 LTS. We recommend [Anaconda](https://www.anaconda.com/distribution/) to manage the Python environment (Python ***3.7***) and packages, such as H5py, NumPy, and Scipy. Unlike [HLM](https://github.com/leiliu2015/HLM), it does *not* reply on molecular dynamics simulations any more.

### File Description
- tino/
  - [tino_nan.py](tino/tino_nan.py) (A Python script which reads a Hi-C contact probability matrix to calculates {*k<sub>ij</sub>*})
  - [tino_K2P.py](tino/tino_K2P.py) (A Python script to calculate *n*-body contacts *p<sub>n</sub>*, for *n>=2*, based on {*k<sub>ij</sub>*})
  - [tino_exp.py](tino/tino_exp.py) (A Python script to calculate the expected *p<sub>ijk</sub>*, as defined by Eq. 2 in the [prepreint]())
  - [tino_ps.py](tino/tino_ps.py) (A Python script which merges {*p<sub>ij</sub>*} from Hi-C and from the model for visual comparison)
  - [averageLogs.py](tino/averageLogs.py) (A Python script to average the {*k<sub>ij</sub>*} training log files over independent replicas)
  - [tino_K2cfgs.py](tino/tino_K2cfgs.py) (A *auxiliary* Python script to generate three-dimensional chromatin structures based on {*k<sub>ij</sub>*})
  - [mdRandomSeeds](tino/mdRandomSeeds) (A TXT file including random number seeds)
  - xxx.pal (Color palettes used in Gnuplot scripts)
- tric/
  - [tric.sh](tric/tric.sh) (A BASH script to perform all the modeling and analyses in this directory)
  - ESC/ERY_alpha.txt (Contact frequency matrices, {*p<sub>ij</sub>*}, measured by [Oudelaar *et al.*](http://dx.doi.org/10.1038/s41467-020-16598-7) with Hi-C, based on which we inferred the values of model parameters {*k<sub>ij</sub>*})
  - tric-data/ (A subdirectory including three-body chromatin contact frequencies measured by [Oudelaar *et al.*](http://dx.doi.org/10.1038/s41588-018-0253-2) with Tri-C)
  - tric_xxx.py (Case specific analyzing scripts)
  - tric_x.gnu (Gnuplot scripts to visualize the modeling and prediction results)
- Other direcotries/
  - Similar to [tric/](tric/), in each directory, there is *only* one BASH script to carry out all tasks, and one or more Gnuplot scripts to plot the results.

### User Guide
Except the directory [tino/](tino/), each directory includes an application of HLM-Nbody which can be executed in a similar manner. Taking the directory [tric/](tric/) as an example, all you need to do is to run the BASH script in your terminal
```
$ cd ./tric
$ bash ./tric.sh
``` 
This takes less than 5 minutes to finish the modeling on our desktop with a Intel® Core™ i5-9500 CPU. Results can then be plotted with the Gnuplot scripts, e.g.,
```
$ gnuplot -persist tric_a.gnu
``` 
The approximate running time in each directory and their corresponding figure indices in our [preprint]() are listed in the following table for your reference. The running time is relatively long in the last two directories, as the same genomic region is modeled 32 times with different choices in the process of training {*k<sub>ij</sub>*}. Those 32 models are then compared to determine the best choice to train other examples.

| Directory | Running Time /min | Figure |
| --------- | ------------------ | -------|
|gaussianChain | 3 | Fig. 5 and Fig. S8
|chromatinTracing | 1 | Fig. S8|
|tric | 4 | Fig. 2|
|mc4c | 3 | Fig. 3|
|sprite | 2 | Fig. 4 and Fig. S7|
|numeric_mESC | 20 | Fig. S2|
|numeric_GM12878 | 80 | Fig. S3|

All the output files can be deleted by typing `$ bash ./clearAll.sh` at the repository root. For further questions and possible applications about HLM-Nbody, please contact Lei Liu (leiliu2015@163.com)
