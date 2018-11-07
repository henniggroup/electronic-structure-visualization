# Full electronic structure visualization tool

This is an ongoing effort to develop an interactive tool for visualizing full electronic structures of materials.


## Installation:
1. First, you will need to have python installed. We suggest using [anaconda/miniconda](https://conda.io/docs/user-guide/install/index.html)
2. Make sure you have the following python packages installed: 
	* numpy
	* pymatgen
	* plotly
	* dash
	* mendeleev
3. Simply download the scripts in this repository and you are good to go!


## Running the python script:
1. Run `python dash_test.py` from the command line.
2. Open your local server. 


## Providing input files:
Currently:
* user provides the path to data stored locally

Future:
* user can provide mpid to query data from Materials Project (MP doesn't have projected band structure data?)
* incorporate into online databases e.g. Materials Project, MaterialsWeb


## General features of interactive figures:
* zoom
* pan
* rotate (3D figure)
* view data on hover
* save figure as png
* turn on and off individual traces by clicking on the legend


## Structure figure:
Note: this isn't meant to replace the structure figure which can be generated more beautifully by jmol or other visualization software. The purpose of this simple structure figure is to have a representation of the structure which can be easily connected to the band structure figure.

Currently:
* atoms colored using same color scheme as jmol
* clickable atoms (text below the figure indicates the selected atom)

Future:
* connect the structure figure to the band structure figure so that the user can select the atom to project onto by clicking on the atom in the structure figure


## Band structure figure:
Currently:
* user can select the element(s) and orbitals(s) they want to view projections onto
* new plot is generated only upon clicking "Generate Plot" button
* band structure element projections are represented using "fat bands" (there is an option in the script to plot the band structure element projections as colored bands, but the way it is currently being implemented is very slow...)
* band structure orbital projections are represented using "fat bands" 
* up-spin contributions plotted as filled circles, down-spin contributions (if computed) plotted as empty circles
* user can turn on and off individual traces by clicking on the legend on the right

Future:
* enable projections onto selected atoms and/or sub-orbitals
* suggestions for improving visualization of up-spin/down-spin, contributions from different elements, etc.?


## Density of states figure:
Currently:
* user can select the element(s) and orbitals(s) they want to view projections onto
* new plot is generated only upon clicking "Generate Plot" button
* up-spin contributions plotted as solid lines, down-spin contributions (if computed) plotted as dashed lines
* user can turn on and off individual traces by clicking on the legend on the right

Future:
* enable projections onto selected atoms and/or sub-orbitals


NEW as of 11/07/2018:
* user can choose to plot *only* the band structure *or* density of states by specifying only that one `vasprun.xml` file in the appropriate field. If `vasprun.xml` files are specified in both fields, the combined band structure and density of states plot is generated.
