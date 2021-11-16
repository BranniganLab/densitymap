# densitymap
VMD and python tools for calculating and plotting 2D lipid density maps 

## Running: TCL
Start with the TCL folder. In the helpers directory you have three TCL scripts. BinTools has required functions for the overall script. Lipid_Saturation_Head sets up a series of macros to organize lipids by head group and saturation. asign_helices_2BG9_CG breaks down the protein's structure into the different chains and allows users to isolate specific alpha helicies. Not it will only work for _2BG9 structure_. This version is a modification from Grace Brannigan's. 

The file densitymap contains the TCL functions. It is suggested to build a run file to load coarse grained simulations and call the polarDensityBin function from there. To run:

> source polarDensity.tcl
> polarDensityBin <file name of your choosing> <lipid species by resname head group or saturation> <min radius from center> <max radius from center> <radial step size> <number of theta bins>.
  
Please run for one lipid at a time.
**Side Note**: This script uses qwrap by Jerome Henin (https://github.com/jhenin/qwrap). I have not been sucessful installing qwrap on a Mac.
  
## Running: Python

polarDensity_helper contains the  organization and plotting routines. This is the file you want to change paths in for your data files! Density_Analysis runs the analysis. Please note, in it's current stage it has been designed to work with a few very specific membranes.
  
To run:
> Density(<membrane type>, ddg=False, lipids=[list of lipid resnames, saturations, or head groups], enrich=<True/False>)
If enrichment is True, it will produce an a polar enrichemnt plot normalized to the bulk membrane. If False, it will produce a density plot. Please leave ddg=False, this was implemented later in Polar_Binning_DeltaG (https://github.com/BranniganLab/Polar_Binning_DeltaG).
