import numpy as np

import dash, json
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State

from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.dos import CompleteDos, _get_orb_type
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine

import plotly.tools as tls
tls.set_credentials_file(username='annemarietan', api_key='373kEaPah9OkvR1HbBha')

from gen_structfig import StructFig
from gen_bandsfig import BandsFig


class MyEncoder(json.JSONEncoder):
    ## convert numpy types to regular types for conversion to json
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)

## query materials project
## but I can't seem to get the projected band info...
#mpr = MPRester("tnS76clmyw18JNre")
#dos = mpr.get_dos_by_material_id("mp-149")
#bs = mpr.get_bandstructure_by_material_id("mp-149") ## Si 

## use local test dataset
#folder = "../bandstructureplots_Vallverdu/Si_bands/"
#structure = Vasprun(folder+"Bandes/vasprun.xml").structures[-1]
#dos = Vasprun(folder+"DOS/vasprun.xml").complete_dos
#bands = Vasprun(folder+"Bandes/vasprun.xml", parse_projected_eigen = True)
#bs = bands.get_band_structure(folder+"Bandes/KPOINTS", line_mode=True, efermi=dos.efermi)


## now let's try putting it together with other dash components
app = dash.Dash()

## set fonts
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})
mathjax = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML'
app.scripts.append_script({ 'external_url' : mathjax })

## set colors
colors = {'background': '#FFFFFF',
          'text1': '#0021A5',
          'text2': '#FA4616'
          }


app.layout = html.Div([
    
    ## set header
    html.H3('This is a test page',
            style={'textAlign': 'center',
                   'color': colors['text1']
                   }
            ), 

    ## we'll define this line of text as its own section (division)
    html.Div('currently a work in progress by Anne Marie Tan :)',
             style={'textAlign': 'center',
                    'color': colors['text2'],
                    'marginBottom': '20'
                   }
            ),

    html.Div('From local data:',
             style={'marginLeft': '50'
                   }
            ),
    
    dcc.Input(id='vasprun_dos',
              type='text',
              value='',
              placeholder='input path to vasprun.xml file from DoS calculation',
              style={'display': 'block',
                     'width': '30%',
                     'height': '30px',
                     'marginLeft': '50',
                     'marginBottom': '5',
                     'borderWidth': '1px',
                     'textAlign': 'center'}
              ),
    
    dcc.Input(id='vasprun_bands',
              type='text',
              value='',
              placeholder='input path to vasprun.xml file from bands calculation',
              style={'display': 'block',
                     'width': '30%',
                     'height': '30px',
                     'marginLeft': '50',
                     'marginBottom': '5',
                     'borderWidth': '1px',
                     'textAlign': 'center'}
              ),
    
    dcc.Input(id='kpts_bands',
              type='text',
              value='',
              placeholder='input path to KPOINTS file from bands calculation',
              style={'display': 'block',
                     'width': '30%',
                     'height': '30px',
                     'marginLeft': '50',
                     'borderWidth': '1px',
                     'textAlign': 'center'}
              ),
    
#    html.Div('From Materials Project data:',
#             style={'marginLeft': '50'
#                   }
#            ),
#    
#    dcc.Input(id='mpid',
#              type='text',
#              value='',
#              placeholder='input Materials Project id: mp-###',
#              style={'display': 'block',
#                     'width': '30%',
#                     'height': '30px',
#                     'marginLeft': '50',
#                     'marginBottom': '5',
#                     'borderWidth': '1px',
#                     'textAlign': 'center'}
#              ),
    
    ## hidden divs used to store dos and bs objects
    html.Div(id='dos_object', 
             style={'display': 'none'}),
             
    html.Div(id='bs_object', 
             style={'display': 'none'}),

    ## our simple clickable structure figure!
    dcc.Graph(id='unitcell',
              figure={'data': []},
              style={'display': 'inline-block',
                     'width': '30%'
                     }
              ),

    ## our interactive graph!
    dcc.Graph(id='DOS_bands',
              figure={'data': []},
              style={'display': 'inline-block',
                     'width': '70%',
                     'height': '600'
                     }
              ),

    ## this div tells which atom was selected (temporary)
    html.Div(id='select_atom',
            style={'display': 'inline-block',
               'float': 'left',
               'width': '30%',
               'marginLeft': '100'
               }
        ),
#    html.Div(['Atom(s) to project onto:',
#        dcc.Checklist(
#                id='select_atom',
#                options=[{'label': i, 'value': i}
#                            for i in ['all combined', 'Si1', 'Si2']],
#                values=['all combined'],
#                labelStyle={'display': 'block'}
#            )
#        ],
#        style={'display': 'block',
#           'float': 'left',
#           'width': '30%'}
#        ),

    ## this div contains the checklist for selecting element projections
    html.Div(['Element(s) to project onto:',
        dcc.Checklist(
                id='select_elem',
                options=[],
                values=['all combined'],
                labelStyle={'display': 'block'}
            )
        ],
        style={'display': 'inline-block',
               'float': 'center',
               'width': '20%'
               }
        ),

    ## this div contains the checklist for selecting orbital projections
    html.Div(['Orbitals(s) to project onto:',
        dcc.Checklist(
                id='select_orb',
                options=[{'label': i, 'value': i}
                            for i in ['total', 's', 'p', 'd']],
                values=['total'],
                labelStyle={'display': 'block'}
            )
        ],
        style={'display': 'inline-block',
               'float': 'center',
               'width': '20%'
               }
        ),
        
    html.Button(id='submit_button',
                n_clicks=0,
                children='Generate plot',
                style={'display': 'inline-block',
                       'float': 'right',
                       'width': '20%'
                       }
        )
    ],
    style={'backgroundColor': colors['background']}
)


## The "inputs" and "outputs" of our application interface
## are described declaratively through the app.callback decorator
## The inputs and outputs are simply the properties of a particular component
## Whenever an input property changes,
## the function that the callback decorator wraps will get called. 
## Dash passes the new value of the input property to the function
## and updates the output property with whatever was returned by the function.
        
@app.callback(Output('dos_object', 'children'),
              [Input('vasprun_dos', 'value')])
def get_dos(vasprun_dos):       
    ## get CompleteDos object and "save" in hidden div in json format
    dos = Vasprun(vasprun_dos).complete_dos 
    return json.dumps(dos.as_dict())

#@app.callback(Output('dos_object', 'children'),
#              [Input('mpid', 'value')])
#def get_dos_mp(mpid):       
#    ## get CompleteDos object and "save" in hidden div in json format
#    mpr = MPRester("tnS76clmyw18JNre")
#    dos = mpr.get_dos_by_material_id(mpid)
#    return json.dumps(dos.as_dict())

@app.callback(Output('bs_object', 'children'),
              [Input('dos_object', 'children'),
               Input('vasprun_bands', 'value'),
               Input('kpts_bands', 'value')])
def get_bs(dos, vasprun_bands, kpts_bands):      
    ## get BandStructureSymmLine object and "save" in hidden div in json format
    dos = CompleteDos.from_dict(json.loads(dos))
    bands = Vasprun(vasprun_bands, parse_projected_eigen = True)
    bs = bands.get_band_structure(kpts_bands, line_mode=True, efermi=dos.efermi)     
    return json.dumps(bs.as_dict(), cls=MyEncoder)

#@app.callback(Output('bs_object', 'children'),
#              [Input('mpid', 'value')])
#def get_bs_mp(mpid):      
#    ## get BandStructureSymmLine object and "save" in hidden div in json format
#    mpr = MPRester("tnS76clmyw18JNre")
#    bs = mpr.get_bandstructure_by_material_id(mpid)   
#    return json.dumps(bs.as_dict(), cls=MyEncoder)

@app.callback(Output('DOS_bands', 'figure'),
              [Input('submit_button', 'n_clicks')],
              [State('dos_object', 'children'),
               State('bs_object', 'children'),
               State('select_elem', 'values'),
               State('select_orb', 'values')])
def update_bandsfig(n_clicks, dos, bs, elems, orbs):
    ## de-serialize dos and bs from json format to pymatgen objects 
    dos = CompleteDos.from_dict(json.loads(dos))
    bs = BandStructureSymmLine.from_dict(json.loads(bs))
    if n_clicks > 0: ## only update figure on click   
        ## update the band structure and dos figure
        dosbandfig = BandsFig().generate_fig(dos, bs, elems, orbs)     
    return dosbandfig

#@app.callback(Output('DOS_bands', 'figure'),
#              [Input('vasprun_dos', 'value'),
#               Input('vasprun_bands', 'value'),
#               Input('kpts_bands', 'value'),
#               Input('select_elem', 'values'),
#               Input('select_orb', 'values')])
#def update_bandsfig(vasprun_dos,vasprun_bands,kpts_bands,elemproj,orbproj):
#    ## This version calls pymatgen every time...will it be slower?
#    ## get dos and bs objects
#    dos = Vasprun(vasprun_dos).complete_dos
#    bands = Vasprun(vasprun_bands, parse_projected_eigen = True)
#    bs = bands.get_band_structure(kpts_bands, line_mode=True, efermi=dos.efermi) 
#    ## update the band structure and dos figure
#    dosbandfig = gen_bandsfig.generate_fig(dos, bs, orbproj)   
#    return dosbandfig

@app.callback(Output('unitcell', 'figure'),
              [Input('vasprun_dos', 'value')])  
def update_structfig(vasprun_dos):
    ## Generate our simple structure figure
    structure = Vasprun(vasprun_dos).structures[-1]
    structfig = StructFig().generate_fig(structure)
    return structfig

@app.callback(Output('select_elem', 'options'),
              [Input('vasprun_dos', 'value')])  
def update_elem_options(vasprun_dos):
    ## get list of elements in the system to display as options
    structure = Vasprun(vasprun_dos).structures[-1]  
    elems = []          
    for site in structure.sites:
        if site.specie not in elems:
            elems.append(str(site.specie))
    return [{'label': i, 'value': i} for i in elems]

@app.callback(Output('select_elem', 'values'),
              [Input('select_elem', 'options')])  
def update_elem_value(elems):
    ## set first element as default value for projection
    return [elems[0]['value']]

@app.callback(Output('select_atom', 'children'),
              [Input('unitcell', 'clickData'),
               Input('vasprun_dos', 'value')])
def update_select_atom(clickData,vasprun_dos):
    structure = Vasprun(vasprun_dos).structures[-1]
    if clickData == None:
        return 'You haven\'t selected an atom yet!'
    elif clickData['points'][0]['curveNumber'] < structure.num_sites:
        return 'You\'ve selected atom {}.'.format(clickData['points'][0]['curveNumber'])
    else:
        return 'You\'ve selected something that\'s not an atom!'


if __name__ == '__main__':
    
    app.run_server(debug=True)
    
    