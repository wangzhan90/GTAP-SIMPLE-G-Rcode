# GTAP-SIMPLEG: Auxiliary R code
By: Zhan Wang, Department of Agricultural Economics, Purdue University.

Contact: zhanwang@purdue.edu
# Introduction
This repository contains the R script to perform the parameter calibration and post-simulation analysis reported in the paper "GTAP-SIMPLEG: Integrating Gridded Land Use, Crop Production and Environment Impacts into Global General Equilibrium Model of Trade".

# Contents
## \Calibration
Please put the all contents within \Calibration folder in the model's folder with the mini database (requiring GTAP license), then run the \CalibrateETCP.R to calibrate parameter ETCP, run \CalibrateETMC.R to calibrate the parameter ETMC. Upon the completation of calibration, please run the \do_validate.bat first and \Validate.R later to visualize results after calibration, which will be stored under \visualization folder.

## \PostSim
Please put the entire folder under the model's folder with the full database (requiring GTAP license). After model simuluation, please run the \PostSim\AnalyzeGSG.R script. Simulation results will be stored under \PostSim\out folder.
