import numpy as np
from mendeleev import element


#import plotly as pltly             ## plotting functions    
#import chart_studio.tools as tls   ## plotly tools
import plotly.graph_objs as go     ## plot and configuration tools : Scatter, Line, Layout


class StructFig:
        

    def get_atomTrace(self, site):
       
        [x,y,z],el = site        
        atomTrace = go.Scatter3d(
            x = [x],
            y = [y],
            z = [z],
            mode = 'markers',
            marker = dict(size=20,color=element(el).jmol_color),
            name = str(el),
            showlegend = False
        )
        
        return atomTrace
    
    
    def get_edgeTrace(self, edge):
        
        x,y,z = zip(*edge)
        edgeTrace = go.Scatter3d(
            x = x,
            y = y,
            z = z,
            mode = 'lines',
            line = dict(color="#666666"),
            showlegend = False,
            hoverinfo = 'none'
        )
        
        return edgeTrace
    
    
    def layout(self,structure):
        
        mat = structure.lattice.matrix
        ## trying to hack the aspect ratio...
        yscale = np.amax(mat,0)[1]/np.amax(mat,0)[0]
        zscale = np.amax(mat,0)[2]/np.amax(mat,0)[0]
      
        layout = go.Layout(
        scene=dict(
            camera=dict(eye=dict(x=1.5*zscale, y=1.5*zscale, z=0.0)),
            aspectratio=dict(x=1,
                             y=0.8*yscale,
                             z=0.8*zscale),
            xaxis=dict(
#                autorange=False,
                showspikes=False,
#                showgrid=False,
                zeroline=False,
#                showline=False,
#                title='',
                ticks='',
                showticklabels=False
                ),
            yaxis=dict(
#                autorange=False,
                showspikes=False,
#                showgrid=False,
                zeroline=False,
#                showline=False,
#                title='',
                ticks='',
                showticklabels=False
                ),
            zaxis=dict(
#                autorange=False,
                showspikes=False,
#                showgrid=False,
                zeroline=False,
#                showline=False,
#                title='',
                ticks='',
                showticklabels=False
                )
            )
        )
            
        return layout
    
    
    def generate_fig(self, structure):                
        
        data = []    
        ## draw atoms
        sites = [[site.coords,site.species_string] for site in structure.sites]
        for site in sites:
            data.append(self.get_atomTrace(site))
        
        ## draw cell
        v0 = [0,0,0]
        v1,v2,v3 = structure.lattice.matrix
        v4,v5,v6,v7 = v2+v3,v1+v3,v1+v2,v1+v2+v3
        edges = [[v0,v1],[v0,v2],[v0,v3],[v1,v5],[v1,v6],[v2,v4],
                 [v2,v6],[v3,v4],[v3,v5],[v4,v7],[v5,v7],[v6,v7]]    
        for edge in edges:
            data.append(self.get_edgeTrace(edge))
        
        structfig = go.Figure(data=data, layout=self.layout(structure))
        
        return structfig

    