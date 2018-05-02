# Import the libraries (OpenCMISS,python,numpy,scipy)
import numpy,math,cmath,csv,time,sys,os,pdb

newsolidMatrix=[]

solidNodesList=[10,20,30,40,50,60,70,80,90]
solidMatrix=[[10,20,30,40],[30,40,50,60],[60,70,80,90]]
#newsolidMatrix=[[1,2,3,4],[3,4,5,6],[6,7,8,9]]

for j in range(3):
    row = []
    for k in range(4):
	a=solidNodesList.index(solidMatrix[j][k])+1
        row.append(a)
	#print a
	#print row
    newsolidMatrix.append(row)
print ('>>>>newsolidMatrix<<<'),newsolidMatrix
 





        



