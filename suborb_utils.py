import numpy as np

from pymatgen.electronic_structure.core import Orbital, OrbitalType, Spin
from pymatgen.electronic_structure.dos import Dos
from pymatgen.core import Site
from pymatgen.core.periodic_table import Element
    

def add_densities(densities):
    
    """
    Method to sum two or more densities.

    Parameters
    ----------
    densities: list of densities to sum.
    
    Returns
    -------
    Dict of {spin: density}.
        
    """

    return {spin: sum(np.array(dens[spin]) for dens in densities) 
            for spin in densities[0].keys()}
    

def get_element_pdos(dos,el):
    
    """
    Get element and sub-orbital projected Dos

    Parameters
    ----------
    dos : CompleteDos object
    el : (str) element to obtain projection on
    
    Returns
    -------
    Dict of {orbital: Dos object}.
    
    """
    
    el_dos = {}
    for site, atom_dos in dos.pdos.items(): 
        ## .items() return (key,value) pairs
        if site.specie == Element(el):
            for orb, pdos in atom_dos.items():
                if orb not in el_dos:
                    el_dos[orb] = pdos
                else:
                    el_dos[orb] = add_densities([el_dos[orb], pdos])

    return {orb: Dos(dos.efermi, dos.energies, densities)
                 for orb, densities in el_dos.items()}

    