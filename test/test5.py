# Import the libraries (OpenCMISS,python,numpy,scipy)
import math,numpy,csv,time,sys,os,pdb
from opencmiss.iron import iron


xi0=[0,1,1]
xi1=[1,0,1]
xi2=[1,1,0]
xi3=[1,1,1]
s1=[]
solidXi=[]
interfaceMatrix=[[10,20,30],[10,30,40],[40,50,60]]
solidMatrix=[[10,30,70,20],[90,100,110,80],[30,10,40,200],[60,150,50,40]]
for i in range(0,3):
  for k in range(0,4):
    j=0
    if interfaceMatrix[i][j] in solidMatrix[k]:
       if interfaceMatrix[i][j+1] in solidMatrix[k]:
	  if interfaceMatrix[i][j+2] in solidMatrix[k]:
	    print solidMatrix[k]
	    s1.append(solidMatrix[k])
	    for j in range(0,3):
		a=solidMatrix[k].index(interfaceMatrix[i][j])
                #print a
	        if a==0:
	          solidXi=xi0
                elif a==1:
	          solidXi=xi1
                elif a==2:
	          solidXi=xi2
                elif a==3:
	          solidXi=xi3
	        print solidXi
print s1

