# Project Title

This repository contains the dataset from the 2015 Ginninderra controlled-release experiment, and accompanying scripts to reproduce results presented in ... 

All scripts are written in [R](https://www.r-project.org/).

## Downloading the repository



## Data files provided

### Main data file: 

To be updated.

### CSV files used by the scripts

* GA-inversion-data-1: This is the data during the 5.8 g/min release period 
* GA-inversion-data-2: This is the data during the 5.0 g/min release period 

The above csv files were created using the main data set (see above)

## To reproduce results

* Open either "Main_Ginninderra_Linux.R" or "Main_Ginninderra_Windows_Mac.R", depending on the operating system you run, in RStudio. Set the working directory to the folder containing the BayesianAT scripts, which you downloaded from Github (titled "Scripts" unless you have renamed it since downloading). This is done by going to Session -> Working Directory -> Choose Directory... in the top menu bar in RStudio. 

* Once the working directory is set, ensure the required packages are installed. The required packages for this script are dplyr, tidyr, lubridate, fdrtool, coda, Matrix, and if running a Linux operating system, parallel. To load an already installed package (or check if a package is already installed), type 

library(<package_name>)

into the console and hit the enter/return button. If the package is not yet installed, the Console will return a message saying that the package could not be found. To install a package, type 

install.packages("<package_name>")

into the Console, and hit the enter/return key. 

* After installing the necessary packages, press the "Source" button in the top, right corner of the script panel. This will execute all commands in the script from top to bottom, reproducing all of the results using the full model, with both upwind and downwind measurements, and when the methane-point-source is active.

## To reproduce plots in AMT paper

* Open "Plots.R" in RStudio, and set the working directory as described above.

* Once the working directory is set, ensure the required packages are installed via the instructions above. The required packages for this script are dplyr, lubridate, and ggpubr. 

* After installing the necessary packages, press the "Source" button in the top, right corner of the script panel. This will execute all commands in the script from top to bottom, reproducing all of the plots. If the results from the AMT paper have not yet been reproduced, then the final plot will not work, and an error will be produced. This is because the final plot is of the results and needs results files. 

## Contributing


## Versioning


## Authors

* **Laura Cartwright** 

## License



## Acknowledgments

* 

