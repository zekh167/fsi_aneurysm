# Import the libraries (OpenCMISS,python,numpy,scipy)
import numpy,math,cmath,csv,time,sys,os,pdb
from opencmiss.iron import iron
import numpy as np
#================================================================================================================================
#  Mesh Reading
#================================================================================================================================

print ('Reading geometry from files ...') 

solidElementNumber =0
fluidElementNumber = 0
numberOfSolidElements = 0
numberOfFluidElements = 0
numberOfSolidNodes = 0
numberOfFluidNodes = 0
solidNodesList=[]
fluidNodesList=[]
nodesOfSolid=[]
nodesOfFluid=[]
solidNodes=[]
fluidNodes=[]
newsolidNodes=[]
newfluidNodes=[]
ListS1=[]
ListF1=[]
solidMatrix=[]
fluidMatrix=[]
newsolidMatrix=[]
newfluidMatrix=[]
a=0
solidNodePositions=[]
fluidNodePositions=[]
interfaceNodePositions=[]
solidElementNumberMatrix=[]
fluidElementNumberMatrix=[]
# Read the ele file
with open('hollowcylinder16.1.ele','r') as elefile:
    target = elefile.readlines()
    for lineNum,line in enumerate(target):
        target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(":"," ").split(' ')
        target[lineNum] = filter(None, target[lineNum]) 
        if (lineNum == 0):
	    # Read the total number of elements
	    totalNumberOfElements = int(target[lineNum][0])   
	    print ('...totalNumberOfElements = ') , totalNumberOfElements
        else:
	    # Read the number of nodes and elements
	    regionMarker=int(target[lineNum][5])
	    if regionMarker == -10:
	        elementNumber=int(target[lineNum][0])
	 	#solidElementNumberMatrix.append(elementNumber)
		solidElementNumber+=1
		numberOfSolidElements+=1
		solidElementNodes=[0]*4
                solidElementNodes[0]=int(target[lineNum][1])
                solidElementNodes[1]=int(target[lineNum][2])
                solidElementNodes[2]=int(target[lineNum][3])
                solidElementNodes[3]=int(target[lineNum][4])
		solidMatrix.append(solidElementNodes)
	        #print elementNumber,solidElementNumber,solidElementNodes,regionMarker 
		for Index in range(0,4):
                    if solidElementNodes[Index] in solidNodesList:
	                numberOfSolidNodes+=0
	            else:
                        solidNodesList.append(solidElementNodes[Index])
	                numberOfSolidNodes+=1
		#for i in range(1,5):
		    #nodesOfSolid=int(target[lineNum][i])
		    #solidNodes.append(nodesOfSolid)
	    elif regionMarker == -20:
	        elementNumber=int(target[lineNum][0])
		#fluidElementNumberMatrix.append(elementNumber)
		numberOfFluidElements+=1
		fluidElementNumber+=1
		fluidElementNodes = [0]*4
                fluidElementNodes[0]=int(target[lineNum][1])
                fluidElementNodes[1]=int(target[lineNum][2])
                fluidElementNodes[2]=int(target[lineNum][3])
                fluidElementNodes[3]=int(target[lineNum][4])
		fluidMatrix.append(fluidElementNodes)
                #print elementNumber,fluidElementNumber,fluidElementNodes,regionMarker
		for Index in range(0,4):
                    if fluidElementNodes[Index] in fluidNodesList:
	                numberOfFluidNodes+=0
	            else:
                        fluidNodesList.append(fluidElementNodes[Index])
	                numberOfFluidNodes+=1
		#for i in range(1,5):
		    #nodesOfFluid=int(target[lineNum][i])
		    #fluidNodes.append(nodesOfFluid)
#print solidElementNumberMatrix
#print len(solidElementNumberMatrix)
#print fluidElementNumberMatrix
#print len(fluidElementNumberMatrix)
#print solidMatrix
#print fluidMatrix
#print solidNodesList
#print fluidNodesList

for j in range(0,len(solidMatrix)):
    row = []
    for k in range(4):
        row.append(solidNodesList.index(solidMatrix[j][k])+1)
	#print row
    newsolidMatrix.append(row)
#print ('>>>>newsolidMatrix<<<'),newsolidMatrix



for j in range(0,len(fluidMatrix)):
    row = []
    for k in range(4):
        row.append(fluidNodesList.index(fluidMatrix[j][k])+1)
	#print row
    newfluidMatrix.append(row)
#print ('>>>>newfluidMatrix<<<'),newfluidMatrix

#for j in range(0,len(solidNodes)):
    #if solidNodes[j] not in ListS1:
        #a+=1
        #ListS1.append(solidNodes[j])
	#newsolidNodes.append(a)
    #else:
        #b=ListS1.index(solidNodes[j])+1
        #newsolidNodes.append(b)

#k=0
#for Num in range(0,numberOfSolidElements):
    #N1=[0]*4
    #N1[0]=newsolidNodes[k]
    #N1[1]=newsolidNodes[k+1]
    #N1[2]=newsolidNodes[k+2]
    #N1[3]=newsolidNodes[k+3]
    #k+=4
    #newsolidMatrix.append(N1)
#print ('>>>>newsolidMatrix<<<'),newsolidMatrix
#print ('max(newsolidNodes)'),max(newsolidNodes)
#print ('min(newsolidNodes)'),min(newsolidNodes)
#print ('len(newsolidNodes)'),len(newsolidNodes)

#a=0
#for j in range(0,len(fluidNodes)):
    #if fluidNodes[j] not in ListF1:
        #a+=1
      	#ListF1.append(fluidNodes[j])
	#newfluidNodes.append(a)
    #else:
        #b=ListF1.index(fluidNodes[j])+1
        #newfluidNodes.append(b)

#k=0
#for Num in range(0,numberOfFluidElements):
    #N2=[0]*4
    #N2[0]=newfluidNodes[k]
    #N2[1]=newfluidNodes[k+1]
    #N2[2]=newfluidNodes[k+2]
    #N2[3]=newfluidNodes[k+3]
    #k+=4
    #newfluidMatrix.append(N2)
#print ('>>>>newfluidMatrix<<<'),newfluidMatrix
#print ('max(newfluidNodes)'),max(newfluidNodes)
#print ('min(newfluidNodes)'),min(newfluidNodes)
#print ('len(newfluidNodes)'),len(newfluidNodes)

print ('...numberOfSolidNodes = ') , numberOfSolidNodes
print ('...numberOfSolidElements = '), solidElementNumber
print ('...numberOfFluidNodes = ') , numberOfFluidNodes
print ('...numberOfFluidElements = '), fluidElementNumber

#------------------

# Read the face file
interfaceElementNumber = 0
numberOfInterfaceElements = 0
numberOfInterfaceNodes = 0
interfaceNodesList=[]
nodesOfInterface=[]
interfaceNodes=[]
newinterfaceNodes=[]
newinterfaceMatrix=[]
interfaceMatrix=[]
fluidInletFaceMatrix=[]
inletFluidNodesist=[]
solidInletFaceMatrix=[]
inletSolidNodesist=[]
fluidOutletFaceMatrix=[]
outletFluidNodesist=[]
with open('hollowcylinder16.1.face','r') as facefile:
    target = facefile.readlines()
    for lineNum,line in enumerate(target):
        target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(":"," ").split(' ')
        target[lineNum] = filter(None, target[lineNum])
        if (lineNum == 0):
	    # Read the total number of faces
	    totalNumberOfFaces = int(target[lineNum][0])   
	    print ('...totalNumberOfFaces = ') , totalNumberOfFaces
        else:
	    # Read the interface's nodes
	    regionMarker=int(target[lineNum][4])
	    if regionMarker == -5:
		interfaceElementNodes = [0]*3
                faceNumber=int(target[lineNum][0])
		interfaceElementNumber+=1
                interfaceElementNodes[0]=int(target[lineNum][1])
                interfaceElementNodes[1]=int(target[lineNum][2])
                interfaceElementNodes[2]=int(target[lineNum][3])
		interfaceMatrix.append(interfaceElementNodes)
	        #print faceNumber,interfaceElementNumber,interfaceElementNodes,regionMarker
	        numberOfInterfaceElements+=1 
		for Index in range(0,3):
                    if interfaceElementNodes[Index] in interfaceNodesList:
	                numberOfInterfaceNodes+=0
	            else:
                        interfaceNodesList.append(interfaceElementNodes[Index])
	                numberOfInterfaceNodes+=1 
		#for i in range(1,4):
		    #nodesOfInterface=int(target[lineNum][i])
		    #interfaceNodes.append(nodesOfInterface)	 
	    elif regionMarker == -1:
		fluidInletFaceNodes=[0]*3
		fluidInletFaceNodes[0]=int(target[lineNum][1])
		fluidInletFaceNodes[1]=int(target[lineNum][2])
		fluidInletFaceNodes[2]=int(target[lineNum][3])
		fluidInletFaceMatrix.append(fluidInletFaceNodes)
	        for Index in range(0,3):
		    if fluidInletFaceNodes[Index] in fluidNodesList:
			if fluidInletFaceNodes[Index] not in solidNodesList:
                            if fluidInletFaceNodes[Index] not in inletFluidNodesist:
                                inletFluidNodesist.append(fluidInletFaceNodes[Index])
		solidInletFaceNodes=[0]*3
		solidInletFaceNodes[0]=int(target[lineNum][1])
		solidInletFaceNodes[1]=int(target[lineNum][2])
		solidInletFaceNodes[2]=int(target[lineNum][3])
		solidInletFaceMatrix.append(solidInletFaceNodes)
	        for Index in range(0,3):
		    if solidInletFaceNodes[Index] in solidNodesList:
			if solidInletFaceNodes[Index] not in fluidNodesList:
                            if solidInletFaceNodes[Index] not in inletSolidNodesList:
                                inletSolidNodesist.append(solidInletFaceNodes[Index])
	    elif regionMarker == -2:
		fluidOutletFaceNodes=[0]*3
		fluidOutletFaceNodes[0]=int(target[lineNum][1])
		fluidOutletFaceNodes[1]=int(target[lineNum][2])
		fluidOutletFaceNodes[2]=int(target[lineNum][3])
		fluidOutletFaceMatrix.append(fluidOutletFaceNodes)
	        for Index in range(0,3):
		    if fluidOutletFaceNodes[Index] in fluidNodesList:
			if fluidOutletFaceNodes[Index] not in solidNodesList:
                            if fluidOutletFaceNodes[Index] not in outletFluidNodesist:
                                outletFluidNodesist.append(fluidOutletFaceNodes[Index])
#print len(fluidInletFaceMatrix)
#print len(inletFluidNodesist)
#print inletFluidNodesist
print len(solidInletFaceMatrix)
print len(inletSolidNodesist)
print inletSolidNodesList
#print len(fluidOutletFaceMatrix)
#print len(outletFluidNodesist)
#print outletFluidNodesist
pdb.set_trace()
#print ('>>>interfaceMatrix<<<'), interfaceMatrix

for j in range(0,len(interfaceMatrix)):
    row = []
    for k in range(3):
        row.append(interfaceNodesList.index(interfaceMatrix[j][k])+1)
    newinterfaceMatrix.append(row)
#print ('>>>>newinterfaceMatrix<<<'),newinterfaceMatrix        
#a = {}
#for i in range(0,len(interfaceNodesList)):
    #a[interfaceNodesList[i]]=i+1
#print a
#for i in range(0,len(interfaceNodes)):
    #newinterfaceNodes.append(a[interfaceNodes[i]])
#print newinterfaceNodes

#k=0
#for Num in range(0,numberOfInterfaceElements):
    #N3=[0]*3
    #N3[0]=newinterfaceNodes[k]
    #N3[1]=newinterfaceNodes[k+1]
    #N3[2]=newinterfaceNodes[k+2]
    #k+=3
    #newinterfaceMatrix.append(N3)
#print ('>>>>newinterfaceMatrix<<<'),newinterfaceMatrix

#print ('max(newinterfaceNodes'),max(newinterfaceNodes)
#print ('min(newinterfaceNodes)'),min(newinterfaceNodes)
#print ('len(newinterfaceNodes)'),len(newinterfaceNodes)
print ('...numberOfInterfaceNodes = ') , numberOfInterfaceNodes
print ('...numberOfInterfaceElements = ') , numberOfInterfaceElements

#------------------

# Read the node file
solidNodeNumber = 0
fluidNodeNumber = 0
interfaceNodeNumber = 0
fluidInletNodes=[]
with open("hollowcylinder16.1.node","r") as nodefile:
    target = nodefile.readlines()
    for lineNum,line in enumerate(target):
        target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(":"," ").split(' ')
        target[lineNum] = filter(None, target[lineNum]) 
        if lineNum == 0:
	    # Read the total number of nodes
            totalNumberOfNodes = int(target[lineNum][0])
	    print ('...totalNumberOfNodes = '), totalNumberOfNodes
        else:
	    # Initialise the coordinates
	    nodeNumber = int(target[lineNum][0])
	    if nodeNumber in solidNodesList:
		solidNodeNumber+=1
		solidNodeCoordinates=[0]*4
		solidNodeCoordinates[0]= nodeNumber
                solidNodeCoordinates[1] = float(target[lineNum][1])
                solidNodeCoordinates[2] = float(target[lineNum][2])
                solidNodeCoordinates[3] = float(target[lineNum][3])
		solidNodePositions.append(solidNodeCoordinates)
	    if nodeNumber in fluidNodesList:
  		fluidNodeNumber+=1
		fluidNodeCoordinates=[0]*4
		fluidNodeCoordinates[0]= nodeNumber
                fluidNodeCoordinates[1] = float(target[lineNum][1])
                fluidNodeCoordinates[2] = float(target[lineNum][2])
                fluidNodeCoordinates[3] = float(target[lineNum][3])
		fluidNodePositions.append(fluidNodeCoordinates)
		#if nodeNumber not in interfaceNodesList:
		    #if fluidNodeCoordinates[3] == 0:
			#fluidInletNodes.append(nodeNumber)
	    if nodeNumber in solidNodesList: 
		if nodeNumber in fluidNodesList:
		    interfaceNodeNumber+=1 
		    interfaceNodeCoordinates=[0]*4
		    interfaceNodeCoordinates[0]= nodeNumber
                    interfaceNodeCoordinates[1] = float(target[lineNum][1])
                    interfaceNodeCoordinates[2] = float(target[lineNum][2])
                    interfaceNodeCoordinates[3] = float(target[lineNum][3])
		    interfaceNodePositions.append(interfaceNodeCoordinates)
			    
#print fluidInletNodes
#print len(fluidInletNodes)
#print solidNodePositions
#print fluidNodePositions
#print len(fluidNodePositions)
#print interfaceNodePositions

#------------------

# Creat the connectivity Matrixes

S2=[]
F2=[]
I2=[]
solidConnectElementNumber=0
fluidConnectElementNumber=0
solidElementsConnection=[]
fluidElementsConnection=[]


for i in range (0,numberOfSolidNodes):
    solidNodeIdxMatrix=[0]*2
    solidNodeIdxMatrix[0]=i+1
    solidNodeIdxMatrix[1]=solidNodesList[i]
    S2.append(solidNodeIdxMatrix)
#print S2

for i in range (0,numberOfFluidNodes):
    fluidNodeIdxMatrix=[0]*2
    fluidNodeIdxMatrix[0]=i+1
    fluidNodeIdxMatrix[1]=fluidNodesList[i]
    F2.append(fluidNodeIdxMatrix)
#print F2

for i in range (0,numberOfInterfaceNodes):
    interfaceNodeIdxMatrix=[0]*2
    interfaceNodeIdxMatrix[0]=i+1
    interfaceNodeIdxMatrix[1]=interfaceNodesList[i]
    I2.append(interfaceNodeIdxMatrix)
#print I2

for i in range(0,numberOfInterfaceElements):
    j=0
    for k in range(0,numberOfSolidElements):
        if interfaceMatrix[i][j] in solidMatrix[k]:
            if interfaceMatrix[i][j+1] in solidMatrix[k]:
	        if interfaceMatrix[i][j+2] in solidMatrix[k]:
	            SandICommonElement=[0]*2
            	    solidConnectElementNumber=solidMatrix.index(solidMatrix[k])+1
	    	    interfaceConnectElementNumber=interfaceMatrix.index(interfaceMatrix[i])+1
	   	    SandICommonElement[0]=interfaceConnectElementNumber
	   	    SandICommonElement[1]=solidConnectElementNumber
                    solidElementsConnection.append(SandICommonElement)
		    #print solidConnectElementNumber
		    #print SandICommonElement
	            #print solidMatrix[k]
#print solidElementsConnection
    for z in range(0,numberOfFluidElements):
        if interfaceMatrix[i][j] in fluidMatrix[z]:
            if interfaceMatrix[i][j+1] in fluidMatrix[z]:
	        if interfaceMatrix[i][j+2] in fluidMatrix[z]:
	            FandICommonElement=[0]*2
                    fluidConnectElementNumber=fluidMatrix.index(fluidMatrix[z])+1
	            interfaceConnectElementNumber=interfaceMatrix.index(interfaceMatrix[i])+1
	            FandICommonElement[0]=interfaceConnectElementNumber
	            FandICommonElement[1]=fluidConnectElementNumber
                    fluidElementsConnection.append(FandICommonElement)
	            #print fluidMatrix[z]
#print fluidElementsConnection

#print interfaceNodesList
newNodeConnectionMatrix=[]
ConnectNodesMatrix=[numberOfInterfaceNodes]*3
for i in range(0,numberOfInterfaceNodes):
    ConnectNodesMatrix=[0]*3
    ConnectNodesMatrix[0]=solidNodesList.index(interfaceNodesList[i])+1
    ConnectNodesMatrix[1]=interfaceNodesList.index(interfaceNodesList[i])+1
    ConnectNodesMatrix[2]=fluidNodesList.index(interfaceNodesList[i])+1
    newNodeConnectionMatrix.append(ConnectNodesMatrix)
#print newNodeConnectionMatrix


print ('Reading geometry from files ... Done ')   


# ================================================================================================================================
#  Mesh Connectivity
# ================================================================================================================================

#MAPPING xi
solidElementsOnInterface=[]
fluidElementsOnInterface=[]
for i in range(0,numberOfInterfaceElements):
    j=0
    for k in range(0,numberOfSolidElements):
        if interfaceMatrix[i][j] in solidMatrix[k]:
            if interfaceMatrix[i][j+1] in solidMatrix[k]:
	        if interfaceMatrix[i][j+2] in solidMatrix[k]:
		    solidElementsOnInterface.append(solidMatrix[k])
	            interfacelocalNodes=[0]*3
	            interfacelocalNodes[0]=interfaceNodesList.index(interfaceMatrix[i][j])+1
	            interfacelocalNodes[1]=interfaceNodesList.index(interfaceMatrix[i][j+1])+1
	            interfacelocalNodes[2]=interfaceNodesList.index(interfaceMatrix[i][j+2])+1
	            solidlocalNodes=[0]*3
	            solidlocalNodes[0]=solidNodesList.index(interfaceMatrix[i][j])+1
	            solidlocalNodes[1]=solidNodesList.index(interfaceMatrix[i][j+1])+1
	            solidlocalNodes[2]=solidNodesList.index(interfaceMatrix[i][j+2])+1
                    #print interfacelocalNodes
                    #print solidlocalNodes
    for z in range(0,numberOfFluidElements):
	if interfaceMatrix[i][j] in fluidMatrix[z]:
            if interfaceMatrix[i][j+1] in fluidMatrix[z]:
	        if interfaceMatrix[i][j+2] in fluidMatrix[z]:
		    fluidElementsOnInterface.append(fluidMatrix[z])
	            interfacelocalNodes=[0]*3
	            interfacelocalNodes[0]=interfaceNodesList.index(interfaceMatrix[i][j])+1
	            interfacelocalNodes[1]=interfaceNodesList.index(interfaceMatrix[i][j+1])+1
	            interfacelocalNodes[2]=interfaceNodesList.index(interfaceMatrix[i][j+2])+1
	            fluidlocalNodes=[0]*3
	            fluidlocalNodes[0]=fluidNodesList.index(interfaceMatrix[i][j])+1
	            fluidlocalNodes[1]=fluidNodesList.index(interfaceMatrix[i][j+1])+1
	            fluidlocalNodes[2]=fluidNodesList.index(interfaceMatrix[i][j+2])+1
                    #print interfacelocalNodes
                    #print solidlocalNodes,interfacelocalNodes,fluidlocalNodes
#print solidElementsOnInterface
#print ('len(solidElementsOnInterface)'),len(solidElementsOnInterface)
#print fluidElementsOnInterface
#print ('len(fluidElementsOnInterface)'),len(fluidElementsOnInterface)

#mapping xi

xi0=[0,1,1]
xi1=[1,0,1]
xi2=[1,1,0]
xi3=[1,1,1]
solidXi=[]
fluidXi=[]
for i in range(0,numberOfInterfaceElements):
    for j in range(0,3):
	a=solidElementsOnInterface[i].index(interfaceMatrix[i][j])
        #print a
	if a==0:
	     solidXi=xi0
        elif a==1:
	     solidXi=xi1
        elif a==2:
	     solidXi=xi2
        elif a==3:
	     solidXi=xi3
        #print solidXi
    #print solidElementsOnInterface[i]
	b=fluidElementsOnInterface[i].index(interfaceMatrix[i][j])
        #print b
	if b==0:
	     fluidXi=xi0
        elif b==1:
	     fluidXi=xi1
        elif b==2:
	     fluidXi=xi2
        elif b==3:
	     fluidXi=xi3
	#print fluidXi
    #print fluidElementsOnInterface[i]


# ================================================================================================================================
#  Geometric Field
# ================================================================================================================================


for i in range(0,numberOfSolidNodes):  
    nodeNumber=solidNodePositions[i][0]
    xPosition=solidNodePositions[i][1]
    yPosition=solidNodePositions[i][2]
    zPosition=solidNodePositions[i][3]
    #print nodeNumber,xPosition,yPosition,zPosition


for i in range(0,numberOfFluidNodes):  
    nodeNumber=fluidNodePositions[i][0]
    xPosition=fluidNodePositions[i][1]
    yPosition=fluidNodePositions[i][2]
    zPosition=fluidNodePositions[i][3]
    #print nodeNumber,xPosition,yPosition,zPosition


for i in range(0,numberOfInterfaceNodes):  
    nodeNumber=interfaceNodePositions[i][0]
    xPosition=interfaceNodePositions[i][1]
    yPosition=interfaceNodePositions[i][2]
    zPosition=interfaceNodePositions[i][3]
    print nodeNumber,xPosition,yPosition,zPosition
















