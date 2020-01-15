# Massie_etal_2018_PowerAnalysis
Framework for assessing the ability to detect a macroscale effect of fish growth

File descriptions:
Folder “EmpiricalAnalysis”: 
•	FHC_VonBgrowth.R : R script for quantifying the intraspecific spatial variability of flathead catfish (FHC) growth using a Bayesian hierarchical von Bertalanffy growth model
•	FHCLengthAtAgeData.csv : file containing observed flathead catfish length-at-ages

Folder “PowerAnalysis”:
•	PowerAnaylsis.R : R script for running the power analysis. All steps are commented out. The framework is encouraged to be adapted to investigate other effects [magnitude and direction], species, sampling scenarios, and estimating effects on other growth parameters
•	Age_structure.txt : text file of FHC proportion-at-ages calculated from the observed data 


Abbreviated steps for assessing the ability to detect a macroscale effect of fish growth
Step 1: Use FHC_VonBgrowth.R script to run the Bayesian hierarchical von Bertalanffy  growth model using observed FHC length-at-age data
Step 2: Save MCMC output from spatial growth model as RDS file
Step 3: Open PowerAnalysis.R script and read in the spatial growth model output from step 2 
Step 4: Edit input parameters as needed (e.g.; effect magnitude, number of fish sampled from each lake, and number of lakes sampled)
Step 5: Read in Age_structure.txt file which will be used to generate length-at-age data with the same proportion of ages as observed data 
Step 6:  Grab MCMC parameter estimates from the range of posterior parameter distributions estimated from the empirical analysis (step 1)
Step 7: Generate a population of length-at-age data for each lake using information in steps 5 and 6
Step 8: Randomly sample a specified number of fish from each lake population
Step 9:  Fit von Bertalanffy growth model (step 1) with added covariate to simulated datasets
Step 10: Steps 7-9 are repeated for every simulation run
Step 11: Calculate statistical power, and or, detection probability of the sampling design
*For a more detailed description see Massie et al. 2020
**Important to note that the power analysis process is computationally intensive and could take up to three weeks to complete (depending on the sampling design being investigated)
