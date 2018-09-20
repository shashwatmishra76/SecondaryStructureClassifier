from Bio.PDB import *
import numpy as np

parser = PDBParser()
structure = parser.get_structure("test" , "1a9n.pdb")
for model in structure:
    cha = model['A']
phi = []
psi = []
atom_list = []
for res in cha :
    atom_list.append(res['N'])
    atom_list.append(res['CA'])
    atom_list.append(res['C'])
i = 0
while(i < len(atom_list)-5):
    v1 = atom_list[i].get_vector()
    v2 = atom_list[i+1].get_vector()
    v3 = atom_list[i+2].get_vector()
    v4 = atom_list[i+3].get_vector()
    angle_psi = calc_dihedral(v1,v2,v3,v4)
    psi.append(angle_psi)
    i=i+3
i = 2
while(i < len(atom_list)-2):
    v1 = atom_list[i].get_vector()
    v2 = atom_list[i+1].get_vector()
    v3 = atom_list[i+2].get_vector()
    v4 = atom_list[i+3].get_vector()
    angle_phi = calc_dihedral(v1,v2,v3,v4)
    phi.append(angle_phi)
    i=i+3


"""Eliminating first element of psi list and last element of phi list as the first residue only has a psi dihedral defined and last residue has only a phi residue defined."""
del psi[0]
del phi[-1]

SecStructures = []  //List to hold the secondary structure info for residues

for angle_phi in phi:
    for angle_psi in psi:
        if((angle_phi>0 && angle_psi>0) || (angle_phi<0 && angle_psi<0)):
            SecStructures.append('H')
        elif(angle_phi<0 && angle_psi>0):
            SecStructures.append('B')
        else:
            SecStructures.append('O')
