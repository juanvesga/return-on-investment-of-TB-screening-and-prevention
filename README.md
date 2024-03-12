# Return-on-investment of TB-screening-and-prevention

This repository contains all code for the TB transmission models in support of the work published as **The Return on Investment of Scaling Tuberculosis Screening and Preventive Treatment: A Modelling Study in Brazil, Georgia, Kenya, and South Africa** 
The work considers analysis for four countries (i.e., Brazil, Georgia, Kenya and South Africa). Given that methods for all countries are the same we share here all files relating to South Africa. 
All models and procedures for mathematical modelling were produced using MATLAB&#174;

**How to use**
The master script ("master.m") sereves as a control platform to perform all procedures:
1) run a calibration: change "action" variable in master.m to "mcmc". Run the master script. This will call a number of MCMC functions, run and store a MCMC chain
2) plot calibrations: After running an MCMC, change the "action" variable to "plot_mcmc". This will display model fits to data.
3) Run interventions: a series of interventions can be run after model calibration. action variable can be changed to "interventions", "interventions_noTPT","interventions_noTST", "interventions_fast","interventions_slow"
