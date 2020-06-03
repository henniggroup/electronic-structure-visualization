# Full electronic structure visualization tool

This is an ongoing effort to develop an interactive tool for visualizing full electronic structures of materials. It is in its very early stages so there will probably be bugs, so please be patient! Please send any questions or feedback to Anne Marie at annemarietan@ufl.edu :smiley:.

![screenshot of visualization tool](https://github.com/henniggroup/electronic-structure-visualization/blob/master/screenshot2.png)

NEW as of 11/11/2019:
* Implemented sub-orbital projections for the density of states plot! Sub-orbital projections for the band structure plot coming soon...

NEW as of 08/27/2019:
* Changed the method for selecting elements and orbitals for projection -- responsive dropdown lists enable users to select specific element + orbital combinations easily
* compatible with Dash v1.0


## Installation:
1. First, you will need to have python installed. We suggest using [anaconda/miniconda](https://conda.io/docs/user-guide/install/index.html)
2. Have the following python packages installed (the latest version of each package that this tool has been tested on is listed in parenthesis): 
	* numpy (v. 1.18.4)
	* pymatgen (v. 2019.4.29)
	* chart_studio (v. 1.0.0)
	* plotly (v. 4.8.0)
	* dash (v. 1.12.0)
	* mendeleev (v. 0.6.0)
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
