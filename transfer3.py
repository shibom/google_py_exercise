# You will need numpy and biopython to run this script.
import numpy as np
from Bio.PDB.PDBParser import PDBParser
import os

def lin(z):
    x = (z - c_xz)/m_xz
    y = (z - c_yz)/m_yz
    return x,y

p = PDBParser(PERMISSIVE=1)
structure_id = "PE"
filename = "PEoriginal.pdb" #Input PDB file
structure = p.get_structure(structure_id, filename) # parsing PDB file
model=structure[0]

# Calculate number of chromophores
n=0
for chain in model:
 for residue in chain:
   if residue.get_resname() == "CYC" or "BLA" or "PUB":   
     for atoms in residue:
        if atoms.get_name() == "NA":
           n+=1

h=0
b = ".txt"
c = str(h+1) + b 

resnam = []
ch = []
resno = []
m=1
# Loop for calculating the coordinates of the atoms in the pyrrole of the chromophore

for chain in model: 
 for residue in chain:    
  if residue.get_resname() == "CYC":
     m+=1
     resnam.append("CYC")
     ch.append(chain)
     resno.append(m+1)
     c = str(h+1) + b 
     with open(c,'w') as xyz:
       i=0
       j=0
       k=0
       l=0
       for atoms in residue:
           if atoms.get_name() == "NA":
              coordinatesaN = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "C1A":
              coordinatesa1 = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "C2A":
              coordinatesa2 = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "C3A":
              coordinatesa3 = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "C4A":
              coordinatesa4 = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "NB":
              coordinatesbN = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "C1B":
              coordinatesb1 = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "C2B":
              coordinatesb2 = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "C3B":
              coordinatesb3 = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "C4B":
              coordinatesb4 = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "NC":
              coordinatescN = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "C1C":
              coordinatesc1 = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "C2C":
              coordinatesc2 = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "C3C":
              coordinatesc3 = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "C4C":
              coordinatesc4 = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "ND":
              coordinatesdN = atoms.get_coord()
              l+=1
           elif atoms.get_name() == "C1D":
              coordinatesd1 = atoms.get_coord()
              l+=1
           elif atoms.get_name() == "C2D":
              coordinatesd2 = atoms.get_coord()
              l+=1
           elif atoms.get_name() == "C3D":
              coordinatesd3 = atoms.get_coord()
              l+=1
           elif atoms.get_name() == "C4D":
              coordinatesd4 = atoms.get_coord()
              l+=1

           if i==5:
             ringA = (coordinatesaN + coordinatesa1 + coordinatesa2 + coordinatesa3 + coordinatesa4)/5
             xyz.write("%f\t%f\t%f\n" % (ringA[0],ringA[1],ringA[2]))
             i=0
           if j==5:
             ringB = (coordinatesbN + coordinatesb1 + coordinatesb2 + coordinatesb3 + coordinatesb4)/5
             xyz.write("%f\t%f\t%f\n" % (ringB[0],ringB[1],ringB[2]))
             j=0
           if k==5:
             ringC = (coordinatescN + coordinatesc1 + coordinatesc2 + coordinatesc3 + coordinatesc4)/5
             xyz.write("%f\t%f\t%f\n" % (ringC[0],ringC[1],ringC[2]))
             k=0
           if l==5:
             ringD = (coordinatesdN + coordinatesd1 + coordinatesd2 + coordinatesd3 + coordinatesd4)/5
             xyz.write("%f\t%f\t%f\n" % (ringD[0],ringD[1],ringD[2]))   
             l=0
     h+=1
     xyz.close()
  elif residue.get_resname() == "BLA":
     m+=1
     resnam.append("BLA")
     ch.append(chain)
     resno.append(m+1) 
     c = str(h+1) + b 
     with open(c,'w') as pqr:
       i=0
       j=0
       k=0
       l=0
       for atoms in residue:
           if atoms.get_name() == "NA":
              coordinatesaN = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "C1A":
              coordinatesa1 = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "C2A":
              coordinatesa2 = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "C3A":
              coordinatesa3 = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "C4A":
              coordinatesa4 = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "NB":
              coordinatesbN = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "C1B":
              coordinatesb1 = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "C2B":
              coordinatesb2 = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "C3B":
              coordinatesb3 = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "C4B":
              coordinatesb4 = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "NC":
              coordinatescN = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "C1C":
              coordinatesc1 = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "C2C":
              coordinatesc2 = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "C3C":
              coordinatesc3 = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "C4C":
              coordinatesc4 = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "ND":
              coordinatesdN = atoms.get_coord()
              l+=1
           elif atoms.get_name() == "C1D":
              coordinatesd1 = atoms.get_coord()
              l+=1
           elif atoms.get_name() == "C2D":
              coordinatesd2 = atoms.get_coord()
              l+=1
           elif atoms.get_name() == "C3D":
              coordinatesd3 = atoms.get_coord()
              l+=1
           elif atoms.get_name() == "C4D":
              coordinatesd4 = atoms.get_coord()
              l+=1

           if i==5:
             ringA = (coordinatesaN + coordinatesa1 + coordinatesa2 + coordinatesa3 + coordinatesa4)/5
             pqr.write("%f\t%f\t%f\n" % (ringA[0],ringA[1],ringA[2]))
             i=0
           if j==5:
             ringB = (coordinatesbN + coordinatesb1 + coordinatesb2 + coordinatesb3 + coordinatesb4)/5
             pqr.write("%f\t%f\t%f\n" % (ringB[0],ringB[1],ringB[2]))
             j=0
           if k==5:
             ringC = (coordinatescN + coordinatesc1 + coordinatesc2 + coordinatesc3 + coordinatesc4)/5
             pqr.write("%f\t%f\t%f\n" % (ringC[0],ringC[1],ringC[2]))
             k=0
           if l==5:
             ringD = (coordinatesdN + coordinatesd1 + coordinatesd2 + coordinatesd3 + coordinatesd4)/5
             pqr.write("%f\t%f\t%f\n" % (ringD[0],ringD[1],ringD[2]))   
             l=0
     h+=1
     pqr.close()

  elif residue.get_resname() == "PUB":
     m+=1
     resnam.append("PUB")
     ch.append(chain)
     resno.append(m+1)
     c = str(h+1) + b 
     with open(c,'w') as mno:
       i=0
       j=0
       k=0
       l=0
       for atoms in residue:
           if atoms.get_name() == "NA":
              coordinatesaN = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "C1A":
              coordinatesa1 = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "C2A":
              coordinatesa2 = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "C3A":
              coordinatesa3 = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "C4A":
              coordinatesa4 = atoms.get_coord()
              i+=1
           elif atoms.get_name() == "NB":
              coordinatesbN = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "C1B":
              coordinatesb1 = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "C2B":
              coordinatesb2 = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "C3B":
              coordinatesb3 = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "C4B":
              coordinatesb4 = atoms.get_coord()
              j+=1
           elif atoms.get_name() == "NC":
              coordinatescN = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "C1C":
              coordinatesc1 = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "C2C":
              coordinatesc2 = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "C3C":
              coordinatesc3 = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "C4C":
              coordinatesc4 = atoms.get_coord()
              k+=1
           elif atoms.get_name() == "ND":
              coordinatesdN = atoms.get_coord()
              l+=1
           elif atoms.get_name() == "C1D":
              coordinatesd1 = atoms.get_coord()
              l+=1
           elif atoms.get_name() == "C2D":
              coordinatesd2 = atoms.get_coord()
              l+=1
           elif atoms.get_name() == "C3D":
              coordinatesd3 = atoms.get_coord()
              l+=1
           elif atoms.get_name() == "C4D":
              coordinatesd4 = atoms.get_coord()
              l+=1

           if i==5:
             ringA = (coordinatesaN + coordinatesa1 + coordinatesa2 + coordinatesa3 + coordinatesa4)/5
             mno.write("%f\t%f\t%f\n" % (ringA[0],ringA[1],ringA[2]))
             i=0
           if j==5:
             ringB = (coordinatesbN + coordinatesb1 + coordinatesb2 + coordinatesb3 + coordinatesb4)/5
             mno.write("%f\t%f\t%f\n" % (ringB[0],ringB[1],ringB[2]))
             j=0
           if k==5:
             ringC = (coordinatescN + coordinatesc1 + coordinatesc2 + coordinatesc3 + coordinatesc4)/5
             mno.write("%f\t%f\t%f\n" % (ringC[0],ringC[1],ringC[2]))
             k=0
           if l==5:
             ringD = (coordinatesdN + coordinatesd1 + coordinatesd2 + coordinatesd3 + coordinatesd4)/5
             mno.write("%f\t%f\t%f\n" % (ringD[0],ringD[1],ringD[2]))   
             l=0
     h+=1
     mno.close()
  else:
     m+=1
     continue  
X = np.zeros(n)
Y = np.zeros(n)
Z = np.zeros(n)
avg = np.zeros((n,3))

# Loop for calculating the best fit line for all the pyrrole's in the chromophore
for I in range(n):
  a = str(I+1)
  abc = np.loadtxt(a+".txt")
  # this will find the slope and x-intercept of a plane
  # parallel to the y-axis that best fits the data
  A_xz = np.vstack((abc[0], np.ones(len(abc[0])))).T
  m_xz, c_xz = np.linalg.lstsq(A_xz, abc[2])[0]

  # again for a plane parallel to the x-axis
  A_yz = np.vstack((abc[1], np.ones(len(abc[1])))).T
  m_yz, c_yz = np.linalg.lstsq(A_yz, abc[2])[0]

  zz = np.linspace(0,1,2) 
  xx,yy = lin(zz)
  X[I] = xx[1]-xx[0]
  Y[I] = yy[1]-yy[0]
  Z[I] = zz[1]-zz[0]
  # Calculating mid-point of chromophore
  avg[I][0] = (abc[0][0] + abc[1][0] + abc[2][0] + abc[3][0])/4
  avg[I][1] = (abc[0][1] + abc[1][1] + abc[2][1] + abc[3][1])/4
  avg[I][2] = (abc[0][2] + abc[1][2] + abc[2][2] + abc[3][2])/4

t0 = 2.5 # Fluorescence lifetime of the PEB chromophore, in ns
R0 = 50  # in Angstroms
# Loop for calculating the energy transfer between two chromophores  
with open("Results.txt",'w') as pqr:
 pqr.write("%s %i %s" %("Number of chromophores is:",n,"\n"))
 pqr.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("Chain A","Name A","Number A","|","Chain B","Name B","Number B","Energy transfer (in ns-1)")) 
 for I in range(n):
   J=I+1
   while J<n:
     R1 = np.sqrt(X[I]*X[I]+Y[I]*Y[I]+Z[I]*Z[I])
     R2 = np.sqrt(X[J]*X[J]+Y[J]*Y[J]+Z[J]*Z[J])
     angle = np.arccos((X[I]*X[J]+Y[I]*Y[J]+Z[I]*Z[J])/(R1*R2))
     R = np.sqrt((avg[I][0]-avg[J][0])*(avg[I][0]-avg[J][0])+(avg[I][1]-avg[J][1])*(avg[I][1]-avg[J][1])+(avg[I][2]-avg[J][2])*(avg[I][2]-avg[J][2]))
     rm = R0/R
     ket = (angle*angle)/(t0*rm*rm*rm*rm*rm*rm) # Formula for calculating energy transfer
     pqr.write("%s\t%s\t%i\t%s\t%s\t%s\t%i\t%s\n" %(ch[I],resnam[I],resno[I],"|",ch[J],resnam[J],resno[J],str(ket)))
     J+=1
pqr.close()
# Removing temporary files of the chromophores
for i in range(n):
  c = str(i+1) + b
  os.remove(c)     
