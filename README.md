# AntioxidantsODE
MATLAB and python files for ODE simulation and analysis of b-lapachone treatment in head and neck cancer

<b>scrnaseq_anal.ipynb</b>
jupyter notebook file that generates figure 2 and analyzes scRNAseq data from GEO (GEO Accession: GSE103322).

<b>Figure3.m</b>
MATLAB file that runs sensitivity analysis for figure 3

<b>Figure4.m</b>
MATLAB file that runs the single cell ODE simulations based on scRNAseq data

<b>ODE_result_anal.ipynb</b>
jupyter notebook file that analyzes the single cell ODE simulation results

<b>blap_ODE_system.m</b> 
MATLAB file that contains the function with the ODE system of antioxidants with the following parameters:<br>
endtime - end time in simulated minutes<br>
step - timestep in minutes<br>
celldens - number of cells to simulate<br>
protein_abund - vector of protein abundances used to adjust ODE parameters and initial values<br>
blap_conc - concentration of b-lapachone<br>
parameter_pctchange - matrix to adjust parameters for sensitivity analysis (set to ones for no effect)<br>
species_pctchange - matrix to adjust initial values of species (set to ones for no effect)<br>

