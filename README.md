# Bayesian atmospheric tomography with application to data from the 2015 Ginninderra release experiment

This repository contains a curated dataset from the 2015 Ginninderra controlled-release experiment, and accompanying scripts required to reproduce results presented in "Bayesian atmospheric tomography for detection and quantification of methane emissions: Application to data from the 2015 Ginninderra release experiment", by Laura Cartwright, Andrew Zammit-Mangion, Sangeeta Bhatia, Ivan Schroder, Frances Phillips, Trevor Coates, Karita Negandhi, Travis Naylor, Martin Kennedy, Steve Zegelin, Nick Wokker, Nicholas M. Deutscher, and Andrew Feitz.

All scripts are written in [R](https://www.r-project.org/). The easiest way to run them (and the way described below) is via [RStudio](https://www.rstudio.com/). Note you must download and install R before RStudio.

## Downloading the repository

This can either be done by clicking the green "Clone or download" button, or by running the following in the command line/terminal:

git clone https://github.com/Lcartwright94/BayesianAT

## Data files

### CSV files used by the scripts

* GA-inversion-data-1: This is the data for the 5.8 g/min release period 
* GA-inversion-data-2: This is the data for the 5.0 g/min release period 

### Official curated data set

The official curated data set (including suitable meta data and dictionaries) can be accessed at [ http://dx.doi.org/10.26186/5cb7f14abd710](http://dx.doi.org/10.26186/5cb7f14abd710).

## To reproduce results

* Open either "Main_Ginninderra_Linux.R" or "Main_Ginninderra_Windows_Mac.R", depending on the operating system you run, in RStudio. Set the working directory to the folder containing the BayesianAT scripts, which you downloaded from Github (titled "Scripts" unless you have renamed it since downloading). This is done by going to Session -> Working Directory -> Choose Directory... in the top menu bar in RStudio. 

* Once the working directory is set, ensure the required packages are installed. The required packages for this script are dplyr, tidyr, lubridate, fdrtool, coda, Matrix, and if running a Linux operating system, parallel. To load an already installed package (or check if a package is already installed), type 

```diff
library(<package_name>)
```

into the console and hit the enter/return button. If the package is not yet installed, the Console will return a message saying that the package could not be found. To install a package, type 

```diff
install.packages("<package_name>")
```

into the Console, and hit the enter/return key. 

* After installing the necessary packages, press the "Source" button in the top, right corner of the script panel. This will execute all commands in the script from top to bottom, reproducing all of the results using the full model, with both upwind and downwind measurements, and when the methane-point-source is active.

## To reproduce plots in the paper

* Open "Plots.R" in RStudio, and set the working directory as described above.

* Once the working directory is set, ensure the required packages are installed via the instructions above. The required packages for this script are dplyr, lubridate, ggplot2, and ggpubr. 

* After installing the necessary packages, press the "Source" button in the top, right corner of the script panel. This will execute all commands in the script from top to bottom, reproducing all of the plots. If the results from the paper have not yet been reproduced, then the final plot will not work, and an error will be produced. This is because the final plot is of the results and needs results files. If the results have been reproduced, your final plot should look like this:

<a rel="results" href="https://github.com/Lcartwright94/BayesianAT/blob/master/IMG/Final_res.png"><img alt="Results plot" style="border-width:0" src="https://github.com/Lcartwright94/BayesianAT/blob/master/IMG/Final_res.png" /></a>



## Software license

Reproducible code for "Bayesian atmospheric tomography for detection and quantification of methane emissions: Application to data from the 2015 Ginninderra release experiment" by Cartwright et al.  
Copyright (c) 2019 Laura Cartwright  
Author: Laura Cartwright (lcartwri@uow.edu.au)

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.


## Data license

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />All data in this repository is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a> (Geoscience Australia).

