# Import the libraries (python,numpy)
import numpy,sys,os,pdb

import numpy,sys,os
newMatrix = []
a = {}
nodes=[10,20,30,40,50,60,70,80,90]
matrix = [10,20,30,40,30,40,50,60,60,50,70,80,30,40,80,90]
for i in range(0,len(nodes)):
    a[nodes[i]]=i+1
print a
for i in range(0,len(matrix)):
  newMatrix.append(a[matrix[i]])
print newMatrix
