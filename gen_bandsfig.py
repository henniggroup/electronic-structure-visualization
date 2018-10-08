from mendeleev import element
import numpy as np

from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import Spin, OrbitalType
from pymatgen.electronic_structure.plotter import BSPlotter

import plotly.plotly as pltly      ## plotting functions
import plotly.tools as tls         ## plotly tools
import plotly.graph_objs as go     ## plot and configuration tools : Scatter, Line, Layout
from plotly import offline


class BandsFig:
    
    def __init__(self):
    
        ## define some styling choices??
        self.spins = {'up': Spin.up,'down': Spin.down}
        self.orbs = {'s': OrbitalType.s,'p': OrbitalType.p,'d': OrbitalType.d}  
        self.fillcolors = {'s':{'up':'rgba(255,0,0,1)','down':'rgba(255,0,0,0)'},
                           'p':{'up':'rgba(0,128,0,1)','down':'rgba(0,128,0,0)'},
                           'd':{'up':'rgba(0,0,255,1)','down':'rgba(0,0,255,0)'}}   
        self.edgecolors = {'s':'rgb(255,0,0)','p':'rgb(0,128,0)','d':'rgb(0,0,255)'}
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
    
    
    def get_erange(self, bs):
        
        ## Determine energy ranges for the plots    
        if bs.is_metal():
            return [-10,10]
        else:
            return [bs.get_vbm()["energy"] - bs.efermi - 4,
                    bs.get_cbm()["energy"] - bs.efermi + 4]
    
    
    def get_bandrange(self, bs, brange=6):
        
        ## set the range of bands to plot
        ## default is 6 bands above and below band gap         
        if bs.is_metal():
            bmin = 0
            bmax = bs.nb_bands        
        elif bs.is_spin_polarized:
            bmin = max(bs.get_vbm()["band_index"][Spin.up][0]-brange+1,
                       bs.get_vbm()["band_index"][Spin.down][0]-brange+1,0)
            bmax = min(bs.get_cbm()["band_index"][Spin.up][0]+brange,
                       bs.get_cbm()["band_index"][Spin.down][0]+brange,bs.nb_bands)
        else:
            bmin = max(bs.get_vbm()["band_index"][Spin.up][0]-brange+1,0)
            bmax = min(bs.get_cbm()["band_index"][Spin.up][0]+brange,bs.nb_bands)
        
        return [b for b in range(bmin,bmax)]


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

   
    def get_el_contrib(self,x,pbands,band,elems,spin):
        
        ## compute normalized contribution from each element        
        contrib = np.zeros((len(x),len(elems)))
        for k in range(len(x)):
            temp = [pbands[self.spins[spin]][band][k][elem]**2 for elem in elems]
            tot = np.sum(temp)
            if tot != 0.0:
                for e in range(len(temp)):
                    contrib[k,e] = temp[e]/tot
                        
        return contrib
    

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
        if spin is 'down':
            linestyle['dash'] = 'dot'
        
        return go.Scatter(x = dos.densities[self.spins[spin]],
                          y = dos.energies - dos.efermi,
                          mode = "lines",
                          name = "DoS: total "+self.arrows[spin],
                          line = linestyle
#                          line = dict(color="#202020"),
#                          fill = "tozeroy"
                          )

    
    def get_dosTraces_el(self, dos, elems):
        
        ## get line plots for the element and orbital projections of DoS 
    
        self.get_elemcolors(elems)
        el_dos = dos.get_element_dos()
                
        dosTraces = list()
        for elem in elems:
            ## plot up-spin contributions
            dosTraces.append(self.dosTrace_el(dos,el_dos,elem,spin='up'))
            ## plot down-spin contributions if spin-polarized calculation
            if len(el_dos[Element(elem)].densities) == 2:
                dosTraces.append(self.dosTrace_el(dos,el_dos,elem,spin='down'))
        
        return dosTraces


    def dosTrace_el(self,dos,el_dos,elem,spin):
        
        ## generates, formats, and returns a line trace for the element projected DoS
        
        linestyle = dict(color=element(elem).jmol_color)
        if spin is 'down':
            linestyle['dash'] = 'dot'
        
        return go.Scatter(x = el_dos[Element(elem)].densities[self.spins[spin]],
                          y = dos.energies - dos.efermi,
                          mode = "lines",
                          name = "DoS: "+elem+" "+self.arrows[spin],
                          line = linestyle
                          )
        
    
    def get_dosTraces_spd(self, dos, elem, orbs):
        
        ## get line plots for the element and orbital projections of DoS 
    
        el_spd_dos = dos.get_element_spd_dos(elem)
                
        dosTraces = list()
        for orb in orbs: 
            ## plot up-spin contributions
            dosTraces.append(self.dosTrace_spd(dos,el_spd_dos,elem,orb,spin='up'))
            ## plot down-spin contributions if spin-polarized calculation
            if len(el_spd_dos[self.orbs[orb]].densities) == 2:
                dosTraces.append(self.dosTrace_spd(dos,el_spd_dos,elem,orb,spin='down'))
        
        return dosTraces


    def dosTrace_spd(self,dos,el_spd_dos,elem,orb,spin):
        
        ## generates, formats, and returns a line trace for the orbital projected DoS
        
        linestyle = dict(color=self.edgecolors[orb])
        if spin is 'down':
            linestyle['dash'] = 'dot'
        
        return go.Scatter(x = el_spd_dos[self.orbs[orb]].densities[self.spins[spin]],
                          y = dos.energies - dos.efermi,
                          mode = "lines",
                          name = "DoS: "+elem+" "+orb+" "+self.arrows[spin],
                          line = linestyle
                          )
        

    def get_bandTraces(self, bs):
        
        ## get simple line plots for the total band structure

        bandrange = self.get_bandrange(bs)
        [x,yu,yd] = self.get_bandsxy(bs, bandrange)
 
        ## Each band is plotted as a separate trace
        ## and appended to the list bandTraces    
        bandTraces = list()
        ## plot up-spin contributions
        for b,band in enumerate(bandrange):
            bandTraces.append(self.bandTrace(x,yu[b],band,spin='up'))
        bandTraces[-1].showlegend = True  ## show only 1 legend for all traces in the group        
        ## plot down-spin contributions if spin-polarized calculation
        if bs.is_spin_polarized:
            for b,band in enumerate(bandrange):
                bandTraces.append(self.bandTrace(x,yd[b],band,spin='down'))
            bandTraces[-1].showlegend = True
            
        return bandTraces
    
    
    def bandTrace(self,x,y,band,spin):
        
        ## generates, formats, and returns a simple line trace for the band structure
        
#        linestyle = dict(color=element(elem).jmol_color)
        linestyle = dict(color="#202020")
        if spin is 'down':
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
#        bandrange = self.get_bandrange(bs)
#        [x,yu,yd] = self.get_bandsxy(bs, bandrange)
#        pbands = bs.get_projection_on_elements()
# 
#        ## Each segment of each band is plotted as a separate trace
#        ## and appended to the list colorbandTraces    
#        colorbandTraces = list()
#        ## plot up-spin contributions
#        for b,band in enumerate(bandrange):
#            contrib = self.get_el_contrib(x,pbands,band,elems,spin='up')
#            colorbandTraces.append(self.colorbandTrace(x,yu[b],contrib,band,elems,spin='up'))
#        colorbandTraces[-1][-1].showlegend = True  ## show only 1 legend for all traces in the group        
#        ## plot down-spin contributions if spin-polarized calculation
#        if bs.is_spin_polarized:
#            for b,band in enumerate(bandrange):
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
        
        bandrange = self.get_bandrange(bs)
        [x,yu,yd] = self.get_bandsxy(bs, bandrange) 
        pbands = bs.get_projection_on_elements()
        
        ## Each band is plotted as a separate scatter trace
        ## and appended to the list fatbandTraces 
        fatbandTraces = list()
        ## plot up-spin contributions
        for b,band in enumerate(bandrange):
            fatbandTraces.append(self.fatbandTrace_el(x,yu[b],pbands,band,elem,'up'))
        fatbandTraces[-1].showlegend = True  ## show only 1 legend for all traces in the group
        ## plot down-spin contributions if spin-polarized calculation
        if bs.is_spin_polarized:
            for b,band in enumerate(bandrange):
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

    
    def get_fatbandTraces_spd(self, bs, elem, orbs):
                
        ## get dot plots for the projected band structure
        ## dot sizes proportionl to the relative contribution from that orbital
        
        bandrange = self.get_bandrange(bs)
        [x,yu,yd] = self.get_bandsxy(bs, bandrange) 
        pbands = bs.get_projections_on_elements_and_orbitals({elem: orbs})
        
        ## Each band is plotted as a separate scatter trace
        ## and appended to the list fatbandTraces 
        fatbandTraces = list()
        for orb in orbs: 
            ## plot up-spin contributions
            for b,band in enumerate(bandrange):
                fatbandTraces.append(self.fatbandTrace_spd(x,yu[b],pbands,band,elem,orb,'up'))
            fatbandTraces[-1].showlegend = True  ## show only 1 legend for all traces in the group
            ## plot down-spin contributions if spin-polarized calculation
            if bs.is_spin_polarized:
                for b,band in enumerate(bandrange):
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
            title = "$E - E_f \quad / \quad \\text{eV}$",
            showgrid = True,
            showline = True,
            mirror = "ticks",
            ticks = "inside",
            linewidth = 2,
            tickwidth = 2,
    #        zeroline = False,
            zerolinewidth = 0.5
        )
        
        doslayout = go.Layout(
            title = "Density of states of Silicon",
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
            title="Band structure of Silicon",
            xaxis=bandxaxis,
            yaxis=bandyaxis
        )
        
        return bandlayout
    
    
    def generate_fig(self, dos, bs, elems, orbs):
        
        dosbandfig = tls.make_subplots(rows=1, cols=2, shared_yaxes=True)
        
        ## add total band structure
        bandTraces = self.get_bandTraces(bs)
        ## add the total band structure to subplot (1,1)
        for btrace in bandTraces:
            dosbandfig.append_trace(btrace, 1, 1)                
        ## add total density of states
        dosTraces = self.get_dosTraces_tot(dos)
        ## add the densities to subplot (1,2)
        for dosTrace in dosTraces:
            dosbandfig.append_trace(dosTrace, 1, 2)
        
        if 'total' in orbs:            
            ## add band structure
            for elem in elems:                
                ## add the fat bands
                fatbandTraces = self.get_fatbandTraces_el(bs, elem)
                ## add the fat bands to subplot (1,1)
                for fbtrace in fatbandTraces:
                    dosbandfig.append_trace(fbtrace, 1, 1)
            ## add colored band structure
#            if len(elems) > 1: ## if more than 1 element
#                colorbandTraces = self.get_colorbandTraces(bs, elems)
#                ## add the colored band structure to subplot (1,1)
#                for cbtrace in colorbandTraces:
#                    for cbt in cbtrace:
#                        dosbandfig.append_trace(cbt, 1, 1)                         
            ## add element projected density of states
            eldosTraces = self.get_dosTraces_el(dos, elems)
            ## add the densities to subplot (1,2)
            for dosTrace in eldosTraces:
                dosbandfig.append_trace(dosTrace, 1, 2)                
                
        ## if any orbital projections are selected
        orbs1 = [b for b in orbs if b != "total"] 
        if len(orbs1) > 0:
            for elem in elems:                
                ## add the fat bands
                fatbandTraces = self.get_fatbandTraces_spd(bs, elem, orbs1)
                ## add the fat bands to subplot (1,1)
                for fbtrace in fatbandTraces:
                    dosbandfig.append_trace(fbtrace, 1, 1)                    
                ## add projected density of states
                dosTraces = self.get_dosTraces_spd(dos, elem, orbs1)
                ## add the densities to subplot (1,2)
                for dosTrace in dosTraces:
                    dosbandfig.append_trace(dosTrace, 1, 2)
            
        ## get axes ranges, labels         
        erange = self.get_erange(bs)
        labels_uniq, labelspos_uniq = self.get_kpt_labels(bs)
        bandLayout = self.layout_bands(erange, labels_uniq, labelspos_uniq) 
        dosLayout = self.layout_dos(dos)
            
        ## Customize axes and other aspects of the plot
        ## using previously defined axis and layout options
        dosbandfig["layout"].update(
            go.Layout(
                title="Band structure and density of states",
                xaxis1=bandLayout["xaxis"],
                yaxis1=bandLayout["yaxis"],
                xaxis2=dosLayout["xaxis"]
            )
        )
    
        ## adjust size of subplots
        dosbandfig["layout"]["xaxis1"]["domain"] = [0., 0.7]
        dosbandfig["layout"]["xaxis2"]["domain"] = [0.702, 1.]   
        
        return dosbandfig   


if __name__ == '__main__':
    
    
    tls.set_credentials_file(username='annemarietan', api_key='373kEaPah9OkvR1HbBha')
    
    ## query materials project
    ## but I can't seem to get the projected band info...
#    mpr = MPRester("tnS76clmyw18JNre")
#    dos = mpr.get_dos_by_material_id("mp-2815")
#    bs = mpr.get_bandstructure_by_material_id("mp-2815") ## MoS2 
    
    ## use local test dataset
    folder = "testdata/MoS2/"
    dos = Vasprun(folder+"vasprun.xml").complete_dos
    bands = Vasprun(folder+"vasprun.xml", parse_projected_eigen = True)
    bs = bands.get_band_structure(folder+"KPOINTS", line_mode=True, efermi=dos.efermi)

    dosbandfig = BandsFig().generate_fig(dos, bs, ["Mo","S"], ["total"])
    
    offline.plot(dosbandfig, filename="DOS_bands_MoS2_color.html", auto_open=False)
 
    