#!/usr/bin/env python

import numpy as np
from Bio.PDB.PDBParser import PDBParser
import os, sys

p = PDBParser(PERMISSIVE=1)
filename = "cyc.pdb"
structure = p.get_structure("cyc", filename)
model = structure[0]

# resnam = []; chn = []; resn = []

for chain in model:
    for residue in chain:
       for atoms in residue:
          if atoms.get_name() == "NA":
             coNa = atoms.get_coord()
          elif atoms.get_name() == "C1A":
              co_a1 = atoms.getcoord()
          elif atoms.get_name() == "C2A":
              co_a2 = atoms.get_coord()
          elif atoms.get_name() == "C3A":
              co_a3 = atoms.get_coord()
          elif atoms.get_name() == "C4A":
              co_a4 = atoms.get_coord()
          elif atoms.get_name() == "NB":
              co_Nb = atoms.get_coord()
          elif atoms.get_name() == "C1B":
              co_b1 = atoms.get_coord()
          elif atoms.get_name() == "C2B":
              co_b2 = atoms.get_coord()
          elif atoms.get_name() == "C3B":
              co_b3 = atoms.get_coord()
          elif atoms.get_name() == "C4B":
              co_b4 = atoms.get_coord()
          elif atoms.get_name() == "NC":
              co_Nc = atoms.get_coord()
          elif atoms.get_name() == "C1C":
              co_c1 = atoms.get_coord()
          elif atoms.get_name() == "C2C":
              co_c2 = atoms.get_coord()
          elif atoms.get_name() == "C3C":
              co_c3 = atoms.get_coord()
          elif atoms.get_name() == "C4C":
              co_c4 = atoms.get_coord()
          elif atoms.get_name() == "ND":
              co_Nd = atoms.get_coord()
          elif atoms.get_name() == "C1D":
              co_d1 = atoms.get_coord()
          elif atoms.get_name() == "C2D":
              co_d2 = atoms.get_coord()
          elif atoms.get_name() == "C3D":
              co_d3 = atoms.get_coord()
          elif atoms.get_name() == "C4D":
              co_d4 = atoms.get_coord()



