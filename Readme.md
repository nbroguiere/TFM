TFM3D - Readme.md

Copyright (C) 2018 Nicolas Broguiere

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# Overview
This Matlab series of scripts enables 3D traction force microscopy in biological applications. It is implemented with a simple tracking on a cubic mesh of matrix deformations by cross-correlation, followed by force reconstruction under the assumptions of linear elastic isotropic materials under small deformations.

A variety of graphic and numeric outputs are provided to visualize the experimental displacement field, reconstructed displacement and force field, and the intensity of the stress or of the normal stress at the surface of the cells. 

# Experimental data
For a typical usage, a cell is imaged with a confocal microscope, to obtain a 3D pixel-data stack or timelapse (xyczt). 

An additional identical stack is acquired after killing the cell to provide a reference with relaxed forces. 

The data should contain at least two color channels to be used for cell segmentation (e.g. cell tracker, GFP, or mask from manual segmentation) and matrix tracking (e.g. backscattered light or fluorescent particles) respectively.  

The live data and reference relaxed stack have to be combined (typically using Fiji) to form a single hyperstack (xyczt with one timepoint containing the reference data, typically the last), converted to 8 bit composite and exported as a .tif single file.  

# Software
- Matlab 2017a with:
	- Image Processing toolbox
	- Statistics and Machine Learning Toolbox
	- Control System Toolbox
- Bio-Formats Matlab toolbox:
	- http://www.openmicroscopy.org/bio-formats/downloads/

# Quick start
- Download the scripts and unzip them (in a folder denoted 'scripts' here) as well as the Bio-Formats Matlab toolbox (in a folder denoted 'LOCI' here). 
- Create a folder (denoted 'root' here).
- Copy the dataset to analyze inside 'root' or in a subfolder named 'data'.
- Copy 'master_template.m' inside 'root' and rename it. Denoted 'master' here.
- Open Matlab, navigate to 'root' and edit 'master'. 
- Fill in all the parameters as instructed, including the location of 'scripts' and of 'LOCI', that might be copy-pasted from windows explorer. 
- Run the script (press play or F5)
- See your outputs in 'root>output' (.tif stack files and summary matlab figures), 'root>output_3D' (.tif image sequence for 3D+t animations, import in fiji as an image sequence for visualizing or re-exporting as a compressed .avi movie),  and 'root>output_fig' (matlab figures of individual frames of 3D+t animations that can be analyzed with free mouse controlled rotation/translation/zoom).

# Notes
The script sub-sections are numbered from from 1 to 6, and can only be used in this order. 
Script sections with a number+letter, for example 3b, are only used to generate graphic outputs and are not essential for the execution of the following numbered functions. 
It is always possible to change some parameters and restart the computations half way, for example to try different parameters for force reconstruction without recomputing the experimental displacement. 

A script for saving the workspace in the 'root' folder is provided ('TFM_save_workspace.m'). To restart from a saved workspace, double click on the backup, and additionally execute again 'TFM_open_stack.m'. The work can then be continued from where it was left. 

All the scripts are abundantly commented and the helper functions stackshow and simple_quiver are provided to visualize intermediate data such as smoothed or segmented cells, tensor components, or vector fields such as normals. Checking these things can help to optimize parameters or understand the origin of eventual artifacts. 
