from mendeleev import element
import numpy as np

import suborb_utils as sub

from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import Spin, Orbital, OrbitalType
from pymatgen.electronic_structure.plotter import BSPlotter

#import plotly as pltly             ## plotting functions    
import chart_studio.tools as tls   ## plotly tools
import plotly.graph_objs as go     ## plot and configuration tools : Scatter, Line, Layout
from plotly.subplots import make_subplots
import plotly.io as pio


class BandsFig:
    
    def __init__(self):
    
        ## define some styling choices??
        self.spins = {'up': Spin.up,'down': Spin.down}
        self.orbs = {'s': OrbitalType.s,
                     'p': OrbitalType.p,'px': Orbital.px,
                     'py': Orbital.py,'pz': Orbital.pz,
                     'd': OrbitalType.d,'dz2': Orbital.dz2,
                     'dx2-y2': Orbital.dx2,'dxy': Orbital.dxy,
                     'dxz': Orbital.dxz,'dyz': Orbital.dyz}  
        self.fillcolors = {'s':{'up':'rgba(255,153,0,1)','down':'rgba(255,153,0,0)'},
                           'p':{'up':'rgba(255,0,0,1)','down':'rgba(255,0,0,0)'},
                           'px':{'up':'rgba(153,0,0,1)','down':'rgba(153,0,0,0)'},
                           'py':{'up':'rgba(255,0,0,1)','down':'rgba(255,0,0,0)'},
                           'pz':{'up':'rgba(250,128,114,1)','down':'rgba(250,128,114,0)'},
                           'd':{'up':'rgba(0,0,255,1)','down':'rgba(0,0,255,0)'},
                           'dz2':{'up':'rgba(0,0,153,1)','down':'rgba(0,0,153,0)'},
                           'dx2-y2':{'up':'rgba(51,102,204,1)','down':'rgba(51,102,204,0)'},
                           'dxy':{'up':'rgba(102,255,204,1)','down':'rgba(102,255,204,0)'},
                           'dxz':{'up':'rgba(51,204,51,1)','down':'rgba(51,204,51,0)'},
                           'dyz':{'up':'rgba(0,102,0,1)','down':'rgba(0,102,0,0)'}}   
        self.edgecolors = {'s':'rgb(255,153,0)',
                           'p':'rgb(255,0,0)','px':'rgb(153,0,0)',
                           'py':'rgb(255,0,0)','pz':'rgb(250,128,114)',
                           'd':'rgb(0,0,255)','dz2':'rgb(0,0,153)',
                           'dx2-y2':'rgb(51,102,204)','dxy':'rgb(102,255,204)',
                           'dxz':'rgb(51,204,51)','dyz':'rgb(0,102,0)'}
        self.arrows = {'up':u'\u2191','down':u'\u2193'}
 

    def get_elemcolors(self, elems):
        
        ## get the element colors corresponding to jmol
        ## and convert from hex to rgb for mixing
        elemcolors = []
        for elem in elems:
            elhex = element(elem).jmol_color.lstrip('#')
            elemcolors.append([int(elhex[0:2],16),
                               int(elhex[2:4],16),
                               int(elhex[4:6],16)])
        self.elemcolors = np.array(elemcolors)

        return
    

    def get_erange_dos(self, dos):
        
        ## Determine energy ranges for the plots
        cbm, vbm = dos.get_cbm_vbm()
        if cbm-vbm < 0.01:
            self.erange_dos = [-10,10]
        else:
            self.erange_dos = [vbm - dos.efermi - 4,
                               cbm - dos.efermi + 4]
            
    
    def get_erange_bs(self, bs):
        
        ## Determine energy ranges for the plots    
        if bs.is_metal():
            self.erange_bs = [-10,10]
        else:
            self.erange_bs = [bs.get_vbm()["energy"] - bs.efermi - 4,
                              bs.get_cbm()["energy"] - bs.efermi + 4]
    

    def get_bandrange(self, bs):
        
        ## set the range of bands to plot
        ## based on which show up within the energy window of interest 
        if bs.is_metal():
            bmin = 0
            bmax = bs.nb_bands            
        else:
            [emin,emax] = self.erange_bs
            if bs.is_spin_polarized:
                bvbm = max(bs.get_vbm()["band_index"][Spin.up][0],
                           bs.get_vbm()["band_index"][Spin.down][0])
                bcbm = min(bs.get_cbm()["band_index"][Spin.up][0],
                           bs.get_cbm()["band_index"][Spin.down][0])
                bmin, bmax = bvbm, bcbm
                while (bmin > 0 and 
                       max(bs.bands[Spin.up][bmin]) - bs.efermi > emin or
                       max(bs.bands[Spin.down][bmin]) - bs.efermi > emin):
                    bmin -= 1
                while (bmax < bs.nb_bands and
                       min(bs.bands[Spin.up][bmax]) - bs.efermi < emax or
                       min(bs.bands[Spin.down][bmax]) - bs.efermi < emax):
                    bmax += 1
            else:
                bvbm = bs.get_vbm()["band_index"][Spin.up][0]
                bcbm = bs.get_cbm()["band_index"][Spin.up][0]                
                bmin, bmax = bvbm, bcbm
                while (bmin > 0 and 
                       max(bs.bands[Spin.up][bmin]) - bs.efermi > emin):
                    bmin -= 1
                while (bmax < bs.nb_bands and
                       min(bs.bands[Spin.up][bmax]) - bs.efermi < emax):
                    bmax += 1
                
        self.bandrange = [b for b in range(bmin+1,bmax)]


    def get_bandsxy(self, bs, bandrange):
        
        ## get coords of the band structure data points
        bsplot = BSPlotter(bs)
        data = bsplot.bs_plot_data()
        
        x = [k for kbranch in data["distances"] for k in kbranch]
        yu = [[e - bs.efermi for e in bs.bands[Spin.up][band]] for band in bandrange]
        if bs.is_spin_polarized:
            yd = [[e - bs.efermi for e in bs.bands[Spin.down][band]] for band in bandrange]
        else:
            yd = None
            
        return [x,yu,yd]


    def get_kpt_labels(self, bs):
        
        bsplot = BSPlotter(bs)
        
        ## get unique K-points
        labels = bsplot.get_ticks()["label"]
        labelspos = bsplot.get_ticks()["distance"]
        labels_uniq = [labels[0]]
        labelspos_uniq = [labelspos[0]]
        for i in range(1,len(labels)):
            if labels[i] != labels[i-1]:
                labels_uniq.append(labels[i])
                labelspos_uniq.append(labelspos[i])
            
        labels_uniq = [label.replace("$\mid$","|") for label in labels_uniq]
        ## hack for dash which can't display latex :(
        labels_uniq = [label.replace("$\Gamma$",u"\u0393") for label in labels_uniq]
            
        return labels_uniq, labelspos_uniq

   
#    def get_el_contrib(self,x,pbands,band,elems,spin):
#        
#        ## compute normalized contribution from each element        
#        contrib = np.zeros((len(x),len(elems)))
#        for k in range(len(x)):
#            temp = [pbands[self.spins[spin]][band][k][elem]**2 for elem in elems]
#            tot = np.sum(temp)
#            if tot != 0.0:
#                for e in range(len(temp)):
#                    contrib[k,e] = temp[e]/tot
#                        
#        return contrib
    

    def get_dosTraces_tot(self, dos):
        
        ## get line plots for the total DoS 
                
        dosTraces = list()
        ## plot up-spin contributions
        dosTraces.append(self.dosTrace_tot(dos,spin='up'))
        ## plot down-spin contributions if spin-polarized calculation
        if len(dos.densities) == 2:
            dosTraces.append(self.dosTrace_tot(dos,spin='down'))
        
        return dosTraces
    
    
    def dosTrace_tot(self,dos,spin):
        
        ## generates, formats, and returns a line trace for the total DoS
        
        linestyle = dict(color="#202020")
        if spin == 'down':
            linestyle['dash'] = 'dot'
        
        return go.Scatter(x = dos.densities[self.spins[spin]],
                          y = dos.energies - dos.efermi,
                          mode = "lines",
                          name = "DoS: total "+self.arrows[spin],
                          line = linestyle
#                          line = dict(color="#202020"),
#                          fill = "tozeroy"
                          )

    
    def get_dosTraces_el(self, dos, elem):
        
        ## get line plots for the element and orbital projections of DoS 

        el_dos = dos.get_element_dos()
                
        dosTraces = list()
        ## plot up-spin contributions
        dosTraces.append(self.dosTrace_el(dos,el_dos,elem,spin='up'))
        ## plot down-spin contributions if spin-polarized calculation
        if len(el_dos[Element(elem)].densities) == 2:
            dosTraces.append(self.dosTrace_el(dos,el_dos,elem,spin='down'))
        
        return dosTraces


    def dosTrace_el(self,dos,el_dos,elem,spin):
        
        ## generates, formats, and returns a line trace for the element projected DoS
        
        linestyle = dict(color=element(elem).jmol_color)
        if spin == 'down':
            linestyle['dash'] = 'dot'
        
        return go.Scatter(x = el_dos[Element(elem)].densities[self.spins[spin]],
                          y = dos.energies - dos.efermi,
                          mode = "lines",
                          name = "DoS: "+elem+" "+self.arrows[spin],
                          line = linestyle
                          )
        
    
    def get_dosTraces_spd(self, dos, elem, orb):
        
        ## get line plots for the element and orbital projections of DoS 
    
        el_spd_dos = dos.get_element_spd_dos(elem)
                
        dosTraces = list()
        ## plot up-spin contributions
        dosTraces.append(self.dosTrace_spd(dos,el_spd_dos,elem,orb,spin='up'))
        ## plot down-spin contributions if spin-polarized calculation
        if len(el_spd_dos[self.orbs[orb]].densities) == 2:
            dosTraces.append(self.dosTrace_spd(dos,el_spd_dos,elem,orb,spin='down'))
        
        return dosTraces


    def dosTrace_spd(self,dos,el_spd_dos,elem,orb,spin):
        
        ## generates, formats, and returns a line trace for the orbital projected DoS
        
        linestyle = dict(color=self.edgecolors[orb])
        if spin == 'down':
            linestyle['dash'] = 'dot'
        
        return go.Scatter(x = el_spd_dos[self.orbs[orb]].densities[self.spins[spin]],
                          y = dos.energies - dos.efermi,
                          mode = "lines",
                          name = "DoS: "+elem+" "+orb+" "+self.arrows[spin],
                          line = linestyle
                          )
 

    def get_dosTraces_suborb(self, dos, elem, suborb):
        
        ## get line plots for the element and sub-orbital projections of DoS 
    
#        el_spd_dos = dos.get_element_spd_dos(elem)
        el_suborb_dos = sub.get_element_pdos(dos,elem)
                
        dosTraces = list()
        ## plot up-spin contributions
        dosTraces.append(self.dosTrace_suborb(dos,el_suborb_dos,elem,suborb,spin='up'))
        ## plot down-spin contributions if spin-polarized calculation
        if len(el_suborb_dos[self.orbs[suborb]].densities) == 2:
            dosTraces.append(self.dosTrace_spd(dos,el_suborb_dos,elem,suborb,spin='down'))
        
        return dosTraces


    def dosTrace_suborb(self,dos,el_suborb_dos,elem,suborb,spin):
        
        ## generates, formats, and returns a line trace for the sub-orbital projected DoS
        
        linestyle = dict(color=self.edgecolors[suborb])
        if spin == 'down':
            linestyle['dash'] = 'dot'
        
        return go.Scatter(x = el_suborb_dos[self.orbs[suborb]].densities[self.spins[spin]],
                          y = dos.energies - dos.efermi,
                          mode = "lines",
                          name = "DoS: "+elem+" "+suborb+" "+self.arrows[spin],
                          line = linestyle
                          )
       

    def get_bandTraces(self, bs):
        
        ## get simple line plots for the total band structure

        [x,yu,yd] = self.get_bandsxy(bs, self.bandrange)
 
        ## Each band is plotted as a separate trace
        ## and appended to the list bandTraces    
        bandTraces = list()
        ## plot up-spin contributions
        for b,band in enumerate(self.bandrange):
            bandTraces.append(self.bandTrace(x,yu[b],band,spin='up'))
        bandTraces[-1].showlegend = True  ## show only 1 legend for all traces in the group        
        ## plot down-spin contributions if spin-polarized calculation
        if bs.is_spin_polarized:
            for b,band in enumerate(self.bandrange):
                bandTraces.append(self.bandTrace(x,yd[b],band,spin='down'))
            bandTraces[-1].showlegend = True
            
        return bandTraces
    
    
    def bandTrace(self,x,y,band,spin):
        
        ## generates, formats, and returns a simple line trace for the band structure
        
#        linestyle = dict(color=element(elem).jmol_color)
        linestyle = dict(color="#202020")
        if spin == 'down':
            linestyle['dash'] = 'dot'
        
        return go.Scatter(x = x,
                          y = y,
                          mode = "lines",
                          line = linestyle,
                          name = "BS: total "+self.arrows[spin],
                          text = "Band "+str(band),
                          showlegend = False,
                          legendgroup = "total "+spin
                          )
    
    
#    def get_colorbandTraces(self, bs, elems):
#        
#        ## get line plots for the total band structure
#        ## colored by relative contribution from each element
#
#        self.get_elemcolors(elems)
#        [x,yu,yd] = self.get_bandsxy(bs, self.bandrange)
#        pbands = bs.get_projection_on_elements()
# 
#        ## Each segment of each band is plotted as a separate trace
#        ## and appended to the list colorbandTraces    
#        colorbandTraces = list()
#        ## plot up-spin contributions
#        for b,band in enumerate(self.bandrange):
#            contrib = self.get_el_contrib(x,pbands,band,elems,spin='up')
#            colorbandTraces.append(self.colorbandTrace(x,yu[b],contrib,band,elems,spin='up'))
#        colorbandTraces[-1][-1].showlegend = True  ## show only 1 legend for all traces in the group        
#        ## plot down-spin contributions if spin-polarized calculation
#        if bs.is_spin_polarized:
#            for b,band in enumerate(self.bandrange):
#                contrib = self.get_el_contrib(x,pbands,band,elems,spin='down')
#                colorbandTraces.append(self.colorbandTrace(x,yd[b],contrib,band,elems,spin='down'))
#            colorbandTraces[-1][-1].showlegend = True
#            
#        return colorbandTraces
#
#
#    def colorbandTrace(self,x,y,contrib,band,elems,spin):
#        
#        ## generates, formats, and returns a list of traces
#        ## to plot a band with color gradient
#        ## for band # "band" of the projection onto element "elem"
#        
#        colorBand = list()   
#        for k in range(len(x)-1):
#            ## determine the color of each line segment
#            ## based on the average of the elemental contributions at the 2 endpoints
#            contrib_ave = [(contrib[k,i] + contrib[k+1,i])/2 for i in range(len(elems))]
#            rgb = np.dot(self.elemcolors.T,contrib_ave)
#            linestyle = dict(color="rgb({},{},{})".format(rgb[0],rgb[1],rgb[2]))
#            if spin is 'down':
#                linestyle['dash'] = 'dot'
#            colorBand.append(
#                go.Scatter(
#                    x = [x[k]+1E-04, x[k+1]],
#                    y = [y[k], y[k+1]],
#                    mode = "lines",
#                    name = "BS: el. proj. total "+self.arrows[spin],
#                    text = "Band "+str(band),
#                    line = linestyle,
#                    showlegend = False,
#                    legendgroup = "el. proj. total "+spin
#                )
#            )
#
#        return colorBand

    def get_fatbandTraces_el(self, bs, elem):
                
        ## get dot plots for the projected band structure
        ## dot sizes proportionl to the relative contribution from that element
        
        [x,yu,yd] = self.get_bandsxy(bs, self.bandrange) 
        pbands = bs.get_projection_on_elements()
        
        ## Each band is plotted as a separate scatter trace
        ## and appended to the list fatbandTraces 
        fatbandTraces = list()
        ## plot up-spin contributions
        for b,band in enumerate(self.bandrange):
            fatbandTraces.append(self.fatbandTrace_el(x,yu[b],pbands,band,elem,'up'))
        fatbandTraces[-1].showlegend = True  ## show only 1 legend for all traces in the group
        ## plot down-spin contributions if spin-polarized calculation
        if bs.is_spin_polarized:
            for b,band in enumerate(self.bandrange):
                fatbandTraces.append(self.fatbandTrace_el(x,yd[b],pbands,band,elem,'down'))
            fatbandTraces[-1].showlegend = True
        
        return fatbandTraces


    def fatbandTrace_el(self,x,y,pbands,band,elem,spin):
        
        ## generates, formats, and returns a fatband trace
        ## for band # "band" of the projection onto element "elem"
        
        markersize = [15*pbands[self.spins[spin]][band][k][elem] for k in range(len(x))]
        
        return go.Scatter(x = x,
                          y = y,
                          mode = "markers",
                          marker = dict(size=markersize,
                                        color=element(elem).jmol_color,
                                        opacity=1,
                                        line=dict(color=element(elem).jmol_color,width=1)),
                          name = "BS: "+str(elem)+" "+self.arrows[spin],
                          text = "Band "+str(band),
                          showlegend = False,
                          legendgroup = str(elem)+spin
                          )

    
    def get_fatbandTraces_spd(self, bs, elem, orb):
                
        ## get dot plots for the projected band structure
        ## dot sizes proportionl to the relative contribution from that orbital
        
        [x,yu,yd] = self.get_bandsxy(bs, self.bandrange) 
        pbands = bs.get_projections_on_elements_and_orbitals({elem: [orb]})
        
        ## Each band is plotted as a separate scatter trace
        ## and appended to the list fatbandTraces 
        fatbandTraces = list() 
        ## plot up-spin contributions
        for b,band in enumerate(self.bandrange):
            fatbandTraces.append(self.fatbandTrace_spd(x,yu[b],pbands,band,elem,orb,'up'))
        fatbandTraces[-1].showlegend = True  ## show only 1 legend for all traces in the group
        ## plot down-spin contributions if spin-polarized calculation
        if bs.is_spin_polarized:
            for b,band in enumerate(self.bandrange):
                fatbandTraces.append(self.fatbandTrace_spd(x,yd[b],pbands,band,elem,orb,'down'))
            fatbandTraces[-1].showlegend = True
        
        return fatbandTraces


    def fatbandTrace_spd(self,x,y,pbands,band,elem,orb,spin):
        
        ## generates, formats, and returns a fatband trace
        ## for band # "band" of the projection onto element "elem", orbital "orb"
        
        markersize = [15*pbands[self.spins[spin]][band][k][elem][orb] for k in range(len(x))]
        
        return go.Scatter(x = x,
                          y = y,
                          mode = "markers",
                          marker = dict(size=markersize,
                                        color=self.fillcolors[orb][spin],
                                        opacity=1,
                                        line=dict(color=self.edgecolors[orb],width=1)),
                          name = "BS: "+str(elem)+" "+orb+" "+self.arrows[spin],
                          text = "Band "+str(band),
                          showlegend = False,
                          legendgroup = str(elem)+orb+spin
                          )
    
    
    def layout_dos(self, dos):
        
        ## Customize axes and other aspects of the dos plot.
        
        dosxaxis = go.layout.XAxis(
            title = "Density of states",
            showgrid = True,
            showline = True,
            zeroline = False,
            range = [.01, np.max(dos.densities[Spin.up]) + 0.2],
            mirror = "ticks",
            ticks = "inside",
            linewidth = 2,
            tickwidth = 2
        )
        
        dosyaxis = go.layout.YAxis(
#            title = "$E - E_f \quad / \quad \\text{eV}$",
            title = "E - Ef / eV",
            showgrid = True,
            showline = True,
            mirror = "ticks",
            ticks = "inside",
            linewidth = 1,
            tickwidth = 2,
    #        zeroline = False,
            zerolinecolor = '#ababab',
            zerolinewidth = 2
        )
        
        doslayout = go.Layout(
            xaxis = dosxaxis,
            yaxis = dosyaxis
        )
        
        return doslayout
    
    
    def layout_bands(self, erange, labels_uniq, labelspos_uniq):
        
        ## Customize axes and other aspects of the bands plot.
        
        bandxaxis = go.layout.XAxis(
            title="Wavevector",
            range=[0, labelspos_uniq[-1]],
            showline = True,
            showgrid = True,
            gridcolor='#ababab',
            gridwidth=2,
            ticks = "",
            showticklabels = True,
            linewidth = 1,
            tickvals = labelspos_uniq,
            ticktext = labels_uniq,
            mirror = True
        )
                   
        bandyaxis = go.layout.YAxis(
    #        title = "$E - E_f \quad / \quad \\text{eV}$",
            title = "E - Ef / eV",
            range = erange,
            showline = True,
            showgrid = True,
            mirror = "ticks",
            ticks = "inside",
            linewidth = 1,
            tickwidth = 2,
    #        zeroline = False,
            zerolinecolor = '#ababab',
            zerolinewidth = 2
        )
                   
        bandlayout = go.Layout(
            xaxis=bandxaxis,
            yaxis=bandyaxis
        )
        
        return bandlayout
    

    def generate_fig(self, dos, bs, projlist):
                
        dosbandfig = make_subplots(rows=1, cols=2, shared_yaxes=True)
        
        if bs:
            self.get_erange_bs(bs)
            self.get_bandrange(bs)
            
            ## get total band structure
            bandTraces = self.get_bandTraces(bs)
            ## add the total band structure to subplot (1,1)
            for btrace in bandTraces:
                dosbandfig.append_trace(btrace, 1, 1)
            
            if projlist != None:
                for proj in projlist:
                    elem, orb = proj.split()[0], proj.split()[1]
                    if orb == "Total":
                        ## get fat band structure colored by element
                        fatbandTraces = self.get_fatbandTraces_el(bs, elem)
                        ## add the fat bands to subplot (1,1)
                        for fbtrace in fatbandTraces:
                            dosbandfig.append_trace(fbtrace, 1, 1)
                    elif orb[0] == 's' or orb[0] == 'p' or orb[0] == 'd':
                        ## get the fat bands colored by orbital
                        fatbandTraces = self.get_fatbandTraces_spd(bs, elem, orb[0])
                        ## add the fat bands to subplot (1,1)
                        for fbtrace in fatbandTraces:
                            dosbandfig.append_trace(fbtrace, 1, 1)                                         

            ## format axes
            labels_uniq, labelspos_uniq = self.get_kpt_labels(bs)
            bandLayout = self.layout_bands(self.erange_bs, labels_uniq, labelspos_uniq)            
            dosbandfig["layout"].update(
                go.Layout(
                    title="Band structure",
                    xaxis1=bandLayout["xaxis"],
                    yaxis1=bandLayout["yaxis"]
                )
            )
                        
        if dos:
            self.get_erange_dos(dos)
            
            ## add total density of states
            dosTraces = self.get_dosTraces_tot(dos)
            ## add the densities to subplot (1,2)
            for dosTrace in dosTraces:
                dosbandfig.append_trace(dosTrace, 1, 2)
            
            if projlist != None:
                for proj in projlist:
                    elem, orb = proj.split()[0], proj.split()[1]
                    if orb == "Total":
                        ## get element projected density of states
                        dosTraces = self.get_dosTraces_el(dos, elem)
                    elif orb == 's' or orb == 'p' or orb == 'd':            
                        ## get spd projected density of states
                        dosTraces = self.get_dosTraces_spd(dos, elem, orb)
                    else:
                        ## get full sub-orbital projected density of states
                        dosTraces = self.get_dosTraces_suborb(dos, elem, orb)
                    ## add the densities to subplot (1,2)
                    for dosTrace in dosTraces:
                        dosbandfig.append_trace(dosTrace, 1, 2)
            
            ## format axes
            dosLayout = self.layout_dos(dos)            
            dosbandfig["layout"].update(
                go.Layout(
                    title="Band structure and density of states",
                    xaxis2=dosLayout["xaxis"]
                )
            )            
            if not bs:
                dosLayout["yaxis"].range = self.erange_dos
                dosbandfig["layout"].update(
                    go.Layout(
                        title="Density of states",
                        yaxis1=dosLayout["yaxis"]
                    )
                )
    
        ## adjust size of subplots
        if bs and dos:
            dosbandfig["layout"]["xaxis1"]["domain"] = [0., 0.7]
            dosbandfig["layout"]["xaxis2"]["domain"] = [0.702, 1.]
        elif bs and not dos:
            dosbandfig["layout"]["xaxis1"]["domain"] = [0.2, 0.9]
        elif dos and not bs:
            dosbandfig["layout"]["xaxis2"]["domain"] = [0.4, 0.7]
        
        return dosbandfig


    