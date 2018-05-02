# Import the libraries (OpenCMISS,python,numpy,scipy)
import numpy,math,cmath,csv,time,sys,os,pdb
from opencmiss.iron import iron
import numpy as np

CellMLList=[]
CellMLFaces=[10,30,40,110,100,60,40,90,200,50,120,110,70,170,30]		
fluidNodesList=[10,20,30,40,50,60,70,80,90,100,110,120,130,140,150]
solidNodesList=[10,20,50,150,160,170,180,190,110,60,200,120]
#CellMLList=[30,40,100,90,70]
for i in range(0,len(CellMLFaces)):
    if CellMLFaces[i] in fluidNodesList:
        if CellMLFaces[i] not in solidNodesList:
            if CellMLFaces[i] not in CellMLList:
                  CellMLList.append(CellMLFaces[i])
print CellMLList
