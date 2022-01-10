# Full electronic structure visualization tool

This is a simple interactive tool for visualizing full electronic structures of materials. Unfortunately, it is no longer being actively maintained :frowning_face: Updates for compatibility with major package updates may come sporadically, but don't count on it...

An online version of this tool is also available on the [MaterialsWeb website](https://materialsweb.org/electronic_visualization).

![screenshot of visualization tool](https://github.com/henniggroup/electronic-structure-visualization/blob/master/screenshot2.png)

UPDATE as of 01/10/2022:
* Updated to be compatible with pymatgen v. 2020.*, dash v. 2.0.0, plotly v. 5.1.0

NEW as of 11/11/2019:
* Implemented sub-orbital projections for the density of states plot!

NEW as of 08/27/2019:
* Changed the method for selecting elements and orbitals for projection -- responsive dropdown lists enable users to select specific element + orbital combinations easily


## Installation:
1. First, you will need to have python installed. We suggest using [anaconda/miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
2. Install the following python packages via either conda or pip (the latest version of each package that this tool has been tested on is listed in parenthesis): 
	* numpy (v. 1.21.2)
	* pymatgen (v. 2022.0.14)
	* chart_studio (v. 1.0.0)
	* plotly (v. 5.1.0)
	* dash (v. 2.0.0)
	* mendeleev (v. 0.9.0)
	* sqlalchemy (v. 1.4.27)
3. Simply download the scripts in this repository and you are good to go!


## Running the python script:
1. Run `python dash_main.py` from the command line.
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
* clickable atoms

Future:
* connect the structure figure to the band structure figure so that the user can select the atom to project onto by clicking on the atom in the structure figure


## Band structure figure:
Currently:
* user can select the element(s) and orbitals(s) they want to view projections onto
* new plot is generated only upon clicking "Generate Plot" button
* band structure element projections are represented using "fat bands" (there is an option in the script to plot the band structure element projections as colored bands, but the way it is implemented is very slow so this option is currently disabled...)
* band structure orbital projections are represented using "fat bands" 
* up-spin contributions plotted as filled circles, down-spin contributions (if computed) plotted as empty circles
* user can turn on and off individual traces by clicking on the legend on the right

Future:
* enable projections onto selected atoms and/or sub-orbitals
* suggestions for improving visualization of up-spin/down-spin, contributions from different elements, etc.?


## Density of states figure:
Currently:
* user can select the element(s) orbitals(s) and/or sub-orbitals(s) they want to view projections onto
* new plot is generated only upon clicking "Generate Plot" button
* up-spin contributions plotted as solid lines, down-spin contributions (if computed) plotted as dashed lines
* user can turn on and off individual traces by clicking on the legend on the right

Future:
* enable projections onto selected atoms


## Authors:
Anne Marie Z. Tan

Richard G. Hennig


## How to cite:
BibTex entry for this Github repository::

```
   @misc{electronic-structure-visualization,
     title        = {Full electronic structure visualization tool},
     author       = {A. M. Z Tan and R. G. Hennig},
     year         = 2019,
     publisher    = {GitHub},
     journal      = {GitHub repository},
     howpublished = {\url{https://github.com/henniggroup/electronic-structure-visualization}},
     url          = {https://github.com/henniggroup/electronic-structure-visualization},
     doi          = {}
   }
```
