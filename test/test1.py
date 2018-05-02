# Import the libraries (OpenCMISS,python,numpy,scipy)
import numpy,math,cmath,csv,time,sys,os,pdb
from opencmiss.iron import iron

a=0
List=[]
List2=[]
SolidMatrix=[]
M1=[10,20,30,40,30,40,50,60,60,50,70,80,30,40,80,90]
#M2=[10,20,30,40,50,60]
print M1


#for j in range(0,6):
  #b=M2.index(M2[j])
  #print b

for i in range(0,16):
  if M1[i] not in List:
    a+=1
    List.append(M1[i])
    List2.append(a)
  else:
    b=(List.index(M1[i]))+1
    List2.append(b)
print List2



k=0
for Num in range(0,4):
    N1=[0]*4
    N1[0]=List2[k]
    N1[1]=List2[k+1]
    N1[2]=List2[k+2]
    N1[3]=List2[k+3]
    k+=4
    SolidMatrix.append(N1)
print SolidMatrix


    






