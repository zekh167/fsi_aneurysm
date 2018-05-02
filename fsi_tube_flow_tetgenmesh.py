# > This example program solves a FSI problem in a tube using OpenCMISS.
# >
# > By Chris Bradley
# >
# >

# ================================================================================================================================
#  Start Program
# ================================================================================================================================

# Import the libraries (OpenCMISS,python,numpy,scipy)
import math, numpy, csv, time, sys, os, pdb
from opencmiss.iron import iron


LINEAR = 1
QUADRATIC = 2
CONSTANT = 3

pipeRadius = 1.5
lengthSize = 12

uInterpolation = LINEAR
pInterpolation = CONSTANT

# Time stepping parameters
startTime = 0.0
stopTime = 15.0
timeStep = 0.1

## Inlet velocity parameters
A = 0.25
B = 0.25
C = 5.0

# Material parmeters
# Fluid
Re = 1000

fluidDensity = 1060.0
# fluidDynamicViscosity = fluidDensity*(A+B)*(2.0*pipeRadius)/Re
fluidDynamicViscosity = 0.0035

fluidPInit = 0.0

# Solid
solidDensity = 1000.0
mooneyRivlin1 = 0.1  # N / m^2
mooneyRivlin2 = 0.2

# solidPInit=-mooneyRivlin1
solidPInit = 0.0

# Moving mesh
movingMeshKParameter = 1.0  # default

RBS = False
# RBS = True

outputFrequency = 1  # Result output frequency

# Output flags
fluidEquationsSetOutputType = iron.EquationsSetOutputTypes.NONE
# fluidEquationsSetOutputType = iron.EquationsSetOutputTypes.PROGRESS
fluidEquationsOutputType = iron.EquationsOutputTypes.NONE
# fluidEquationsOutputType = iron.EquationsOutputTypes.TIMING
# fluidEquationsOutputType = iron.EquationsOutputTypes.MATRIX
# fluidEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
fluidDynamicSolverOutputType = iron.SolverOutputTypes.NONE
# fluidDynamicSolverOutputType = iron.SolverOutputTypes.PROGRESS
# fluidDynamicSolverOutputType = iron.SolverOutputTypes.MATRIX
fluidNonlinearSolverOutputType = iron.SolverOutputTypes.NONE
# fluidNonlinearSolverOutputType = iron.SolverOutputTypes.MONITOR
# fluidNonlinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
# fluidNonlinearSolverOutputType = iron.SolverOutputTypes.MATRIX
fluidLinearSolverOutputType = iron.SolverOutputTypes.NONE
# fluidLinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
# fluidLinearSolverOutputType = iron.SolverOutputTypes.MATRIX
solidEquationsSetOutputType = iron.EquationsSetOutputTypes.NONE
# solidEquationsSetOutputType = iron.EquationsSetOutputTypes.PROGRESS
solidEquationsOutputType = iron.EquationsOutputTypes.NONE
# solidEquationsOutputType = iron.EquationsOutputTypes.TIMING
# solidEquationsOutputType = iron.EquationsOutputTypes.MATRIX
# solidEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
movingMeshEquationsSetOutputType = iron.EquationsSetOutputTypes.NONE
# movingMeshEquationsSetOutputType = iron.EquationsSetOutputTypes.PROGRESS
movingMeshEquationsOutputType = iron.EquationsOutputTypes.NONE
# movingMeshEquationsOutputType = iron.EquationsOutputTypes.TIMING
# movingMeshEquationsOutputType = iron.EquationsOutputTypes.MATRIX
# movingMeshEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
interfaceConditionOutputType = iron.InterfaceConditionOutputTypes.NONE
# interfaceConditionOutputType = iron.InterfaceConditionOutputTypes.PROGRESS
interfaceEquationsOutputType = iron.EquationsOutputTypes.NONE
# interfaceEquationsOutputType = iron.EquationsOutputTypes.TIMING
# interfaceEquationsOutputType = iron.EquationsOutputTypes.PROGRESS
# interfaceEquationsOutputType = iron.EquationsOutputTypes.MATRIX
# interfaceEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
# (NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
movingMeshLinearSolverOutputType = iron.SolverOutputTypes.NONE
# movingMeshLinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
# movingMeshLinearSolverOutputType = iron.SolverOutputTypes.MATRIX
fsiDynamicSolverOutputType = iron.SolverOutputTypes.NONE
# fsiDynamicSolverOutputType = iron.SolverOutputTypes.PROGRESS
# fsiDynamicSolverOutputType = iron.SolverOutputTypes.MATRIX
# fsiNonlinearSolverOutputType = iron.SolverOutputTypes.NONE
fsiNonlinearSolverOutputType = iron.SolverOutputTypes.MONITOR
# fsiNonlinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
# fsiNonlinearSolverOutputType = iron.SolverOutputTypes.MATRIX
fsiLinearSolverOutputType = iron.SolverOutputTypes.NONE
# fsiLinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
# fsiLinearSolverOutputType = iron.SolverOutputTypes.MATRIX

# Set solver parameters
fsiDynamicSolverTheta = [0.5]
nonlinearMaximumIterations = 100000000  # default: 100000
nonlinearRelativeTolerance = 1.0E-9  # default: 1.0E-05
nonlinearSolutionTolerance = 1.0E-9  # default: 1.0E-05
nonlinearAbsoluteTolerance = 1.0E-8  # default: 1.0E-10
nonlinearMaxFunctionEvaluations = 10000
nonlinearLinesearchAlpha = 1.0
linearMaximumIterations = 100000000  # default: 100000
linearRelativeTolerance = 1.0E-6  # default: 1.0E-05
linearAbsoluteTolerance = 1.0E-6  # default: 1.0E-10
linearDivergenceTolerance = 1.0E5  # default: 1.0E5
linearRestartValue = 30  # default: 30

progressDiagnostics = True
debug = True

# ================================================================================================================================
#  Should not need to change anything below here.
# ================================================================================================================================

numberOfDimensions = 3

fluidCoordinateSystemUserNumber = 1
solidCoordinateSystemUserNumber = 2
interfaceCoordinateSystemUserNumber = 3

solidRegionUserNumber = 1
fluidRegionUserNumber = 2
interfaceUserNumber = 3

uBasisUserNumber = 1
pBasisUserNumber = 2
interfaceBasisUserNumber = 3

fluidMeshUserNumber = 1
solidMeshUserNumber = 2
interfaceMeshUserNumber = 3
movingMeshUserNumber = 4

fluidDecompositionUserNumber = 1
solidDecompositionUserNumber = 2
interfaceDecompositionUserNumber = 3

fluidGeometricFieldUserNumber = 11
fluidEquationsSetFieldUserNumber = 12
fluidDependentFieldUserNumber = 13
fluidMaterialsFieldUserNumber = 14
fluidIndependentFieldUserNumber = 15
bcCellMLModelsFieldUserNumber = 16
bcCellMLStateFieldUserNumber = 17
bcCellMLParametersFieldUserNumber = 18
bcCellMLIntermediateFieldUserNumber = 19

solidGeometricFieldUserNumber = 21
solidFibreFieldUserNumber = 22
solidEquationsSetFieldUserNumber = 23
solidDependentFieldUserNumber = 24
solidMaterialsFieldUserNumber = 25
solidSourceFieldUserNumber = 26

movingMeshEquationsSetFieldUserNumber = 31
movingMeshDependentFieldUserNumber = 32
movingMeshMaterialsFieldUserNumber = 33
movingMeshIndependentFieldUserNumber = 34

interfaceGeometricFieldUserNumber = 41
interfaceLagrangeFieldUserNumber = 42

fluidEquationsSetUserNumber = 1
solidEquationsSetUserNumber = 2
movingMeshEquationsSetUserNumber = 3

bcCellMLUserNumber = 1

interfaceConditionUserNumber = 1

fsiProblemUserNumber = 1

# ================================================================================================================================
#  Initialise OpenCMISS
# ================================================================================================================================

# Get the computational nodes info
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# iron.OutputSetOn("Testing")

# ================================================================================================================================
#  Coordinate Systems
# ================================================================================================================================

if (progressDiagnostics):
    print('Coordinate systems ...')

# Create a RC coordinate system for the fluid region
fluidCoordinateSystem = iron.CoordinateSystem()
fluidCoordinateSystem.CreateStart(fluidCoordinateSystemUserNumber)
fluidCoordinateSystem.DimensionSet(3)
fluidCoordinateSystem.CreateFinish()
# Create a RC coordinate system for the solid region
solidCoordinateSystem = iron.CoordinateSystem()
solidCoordinateSystem.CreateStart(solidCoordinateSystemUserNumber)
solidCoordinateSystem.DimensionSet(3)
solidCoordinateSystem.CreateFinish()
# Create a RC coordinate system for the interface region
interfaceCoordinateSystem = iron.CoordinateSystem()
interfaceCoordinateSystem.CreateStart(interfaceCoordinateSystemUserNumber)
interfaceCoordinateSystem.DimensionSet(3)
interfaceCoordinateSystem.CreateFinish()

if (progressDiagnostics):
    print('Coordinate systems ... Done')

# ================================================================================================================================
#  Regions
# ================================================================================================================================

if (progressDiagnostics):
    print('Regions ...')

# Create a fluid region
fluidRegion = iron.Region()
fluidRegion.CreateStart(fluidRegionUserNumber, iron.WorldRegion)
fluidRegion.label = 'FluidRegion'
fluidRegion.coordinateSystem = fluidCoordinateSystem
fluidRegion.CreateFinish()
# Create a solid region
solidRegion = iron.Region()
solidRegion.CreateStart(solidRegionUserNumber, iron.WorldRegion)
solidRegion.label = 'SolidRegion'
solidRegion.coordinateSystem = solidCoordinateSystem
solidRegion.CreateFinish()

if (progressDiagnostics):
    print('Regions ... Done')

# ================================================================================================================================
#  Mesh Reading
# ================================================================================================================================

if (progressDiagnostics):
    print ('Reading geometry from files ...') 

#------------------

# Read the ele file
solidElementNumber = 0
fluidElementNumber = 0
numberOfSolidElements = 0
numberOfFluidElements = 0
numberOfSolidNodes = 0
numberOfFluidNodes = 0
solidMatrix=[]
fluidMatrix=[]
solidNodesList=[]
fluidNodesList=[]

with open('hollowcylinder16.1.ele','r') as elefile:
    target = elefile.readlines()
    for lineNum, line in enumerate(target):
        target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(":"," ").split(' ')
        target[lineNum] = filter(None, target[lineNum]) 
        if (lineNum == 0):
	    # Read the total number of elements
	    totalNumberOfElements = int(target[lineNum][0])   
	    print ('...totalNumberOfElements = '), totalNumberOfElements
        else:
	    # Read the number of nodes and elements of solid and fluid
	    regionMarker=int(target[lineNum][5])
	    if regionMarker == -10:
	        elementNumber=int(target[lineNum][0])
		solidElementNumber+=1
		numberOfSolidElements+=1
		solidElementNodes=[0]*4
                solidElementNodes[0]=int(target[lineNum][1])
                solidElementNodes[1]=int(target[lineNum][2])
                solidElementNodes[2]=int(target[lineNum][3])
                solidElementNodes[3]=int(target[lineNum][4])
		solidMatrix.append(solidElementNodes)
		for Index in range(0,4):
                    if solidElementNodes[Index] in solidNodesList:
	                numberOfSolidNodes+=0
	            else:
                        solidNodesList.append(solidElementNodes[Index])
	                numberOfSolidNodes+=1
	    elif regionMarker == -20:
	        elementNumber=int(target[lineNum][0])
		numberOfFluidElements+=1
		fluidElementNumber+=1
		fluidElementNodes = [0]*4
                fluidElementNodes[0]=int(target[lineNum][1])
                fluidElementNodes[1]=int(target[lineNum][2])
                fluidElementNodes[2]=int(target[lineNum][3])
                fluidElementNodes[3]=int(target[lineNum][4])
		fluidMatrix.append(fluidElementNodes)
		for Index in range(0,4):
                    if fluidElementNodes[Index] in fluidNodesList:
	                numberOfFluidNodes+=0
	            else:
                        fluidNodesList.append(fluidElementNodes[Index])
	                numberOfFluidNodes+=1

print ('...numberOfSolidNodes = '), numberOfSolidNodes
print ('...numberOfSolidElements = '), numberOfSolidElements
print ('...numberOfFluidNodes = '), numberOfFluidNodes
print ('...numberOfFluidElements = '), numberOfFluidElements

#------------------

# Read the face file
interfaceElementNumber = 0
numberOfInterfaceElements = 0
numberOfInterfaceNodes = 0
interfaceMatrix=[]
interfaceNodesList=[]
solidInletNodes=[]
fluidInletNodes=[]
interfaceInletNodes=[]
fluidOutletNodes=[]
interfaceOutletNodes=[]

with open('hollowcylinder16.1.face', 'r') as facefile:
    target = facefile.readlines()
    for lineNum,line in enumerate(target):
        target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(":"," ").split(' ')
        target[lineNum] = filter(None, target[lineNum])
        if (lineNum == 0):
	    # Read the total number of faces
	    totalNumberOfFaces = int(target[lineNum][0])   
	    print ('...totalNumberOfFaces = ') , totalNumberOfFaces
        else:
	    # Read the number of nodes and elements of interface
	    regionMarker=int(target[lineNum][4])
	    if regionMarker == -5:
                faceNumber=int(target[lineNum][0])
	        numberOfInterfaceElements+=1
		interfaceElementNumber+=1
		interfaceElementNodes = [0]*3
                interfaceElementNodes[0]=int(target[lineNum][1])
                interfaceElementNodes[1]=int(target[lineNum][2])
                interfaceElementNodes[2]=int(target[lineNum][3])
		interfaceMatrix.append(interfaceElementNodes)
		for Index in range(0,3):
                    if interfaceElementNodes[Index] in interfaceNodesList:
	                numberOfInterfaceNodes+=0
	            else:
                        interfaceNodesList.append(interfaceElementNodes[Index])
	                numberOfInterfaceNodes+=1 
	    # Read the inlet nodes of solid, fluid and Interface	 
	    elif regionMarker == -1:
		inletFaceNodes=[0]*3
		inletFaceNodes[0]=int(target[lineNum][1])
		inletFaceNodes[1]=int(target[lineNum][2])
		inletFaceNodes[2]=int(target[lineNum][3])
	        for Index in range(0,3):
		    if inletFaceNodes[Index] in solidNodesList:
                        if inletFaceNodes[Index] not in solidInletNodes:
                            solidInletNodes.append(inletFaceNodes[Index])
			if inletFaceNodes[Index] in fluidNodesList:
			    if inletFaceNodes[Index] not in interfaceInletNodes:
                                interfaceInletNodes.append(inletFaceNodes[Index])
		    elif inletFaceNodes[Index] not in solidNodesList:
			if inletFaceNodes[Index] in fluidNodesList:
                            if inletFaceNodes[Index] not in fluidInletNodes:
                                fluidInletNodes.append(inletFaceNodes[Index])
	    # Read the outlet nodes of fluid and interface
	    elif regionMarker == -2:
		outletFaceNodes=[0]*3
		outletFaceNodes[0]=int(target[lineNum][1])
		outletFaceNodes[1]=int(target[lineNum][2])
		outletFaceNodes[2]=int(target[lineNum][3])
	        for Index in range(0,3):
		    if outletFaceNodes[Index] in fluidNodesList:
			if outletFaceNodes[Index] in solidNodesList:
			    if outletFaceNodes[Index] not in interfaceOutletNodes:
                                interfaceOutletNodes.append(outletFaceNodes[Index])
			elif outletFaceNodes[Index] not in solidNodesList:
                            if outletFaceNodes[Index] not in fluidOutletNodes:
                                fluidOutletNodes.append(outletFaceNodes[Index])	        

print ('...numberOfInterfaceNodes = ') , numberOfInterfaceNodes
print ('...numberOfInterfaceElements = ') , numberOfInterfaceElements

#------------------

# Read the node file
solidNodePositions=[]
fluidNodePositions=[]
interfaceNodePositions=[]

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
	    # read the coordinates of solid,fluid and interface nodes
	    nodeNumber = int(target[lineNum][0])
	    if nodeNumber in solidNodesList:
		solidNodeCoordinates=[0]*4
		solidNodeCoordinates[0]= nodeNumber
                solidNodeCoordinates[1] = float(target[lineNum][1])
                solidNodeCoordinates[2] = float(target[lineNum][2])
                solidNodeCoordinates[3] = float(target[lineNum][3])
		solidNodePositions.append(solidNodeCoordinates)
		if nodeNumber in fluidNodesList: 
		    interfaceNodeCoordinates=[0]*4
		    interfaceNodeCoordinates[0]= nodeNumber
                    interfaceNodeCoordinates[1] = float(target[lineNum][1])
                    interfaceNodeCoordinates[2] = float(target[lineNum][2])
                    interfaceNodeCoordinates[3] = float(target[lineNum][3])
		    interfaceNodePositions.append(interfaceNodeCoordinates)
	    if nodeNumber in fluidNodesList:
		fluidNodeCoordinates=[0]*4
		fluidNodeCoordinates[0]= nodeNumber
                fluidNodeCoordinates[1] = float(target[lineNum][1])
                fluidNodeCoordinates[2] = float(target[lineNum][2])
                fluidNodeCoordinates[3] = float(target[lineNum][3])
		fluidNodePositions.append(fluidNodeCoordinates)	

#------------------

# Creat the connectivity Matrixes
solidConnectElementNumber=0
fluidConnectElementNumber=0
solidElementsConnection=[]
fluidElementsConnection=[]
# solid element and interface element connection
for i in range(0,numberOfInterfaceElements):
    j=0
    for k in range(0,numberOfSolidElements):
        if interfaceMatrix[i][j] in solidMatrix[k]:
            if interfaceMatrix[i][j+1] in solidMatrix[k]:
	        if interfaceMatrix[i][j+2] in solidMatrix[k]:
	            solidConnectMatrix=[0]*2
                    solidConnectElementNumber=solidMatrix.index(solidMatrix[k])+1
	            interfaceConnectElementNumber=interfaceMatrix.index(interfaceMatrix[i])+1
	            solidConnectMatrix[0]=interfaceConnectElementNumber
	            solidConnectMatrix[1]=solidConnectElementNumber
                    solidElementsConnection.append(solidConnectMatrix)

# fluid element and interface element connection
for i in range(0,numberOfInterfaceElements):
    j=0
    for k in range(0,numberOfFluidElements):
        if interfaceMatrix[i][j] in fluidMatrix[k]:
            if interfaceMatrix[i][j+1] in fluidMatrix[k]:
	        if interfaceMatrix[i][j+2] in fluidMatrix[k]:
	            fluidConnectMatrix=[0]*2
                    fluidConnectElementNumber=fluidMatrix.index(fluidMatrix[k])+1
	            interfaceConnectElementNumber=interfaceMatrix.index(interfaceMatrix[i])+1
	            fluidConnectMatrix[0]=interfaceConnectElementNumber
	            fluidConnectMatrix[1]=fluidConnectElementNumber
                    fluidElementsConnection.append(fluidConnectMatrix)

if (progressDiagnostics):
    print ('Reading geometry from files ... Done ')
   
# ================================================================================================================================
#  Bases
# ================================================================================================================================

if (progressDiagnostics):
    print('Basis functions ...')

numberOfNodesXi = uInterpolation + 1
numberOfGaussXi = uInterpolation + 1

uBasis = iron.Basis()
uBasis.CreateStart(uBasisUserNumber)
uBasis.type = iron.BasisTypes.SIMPLEX
uBasis.numberOfXi = 3
if (uInterpolation == LINEAR):
    uBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX] * 3
elif (uInterpolation == QUADRATIC):
    uBasis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX] * 3
else:
    print('Invalid u interpolation')
    exit()
uBasis.quadratureOrderSet = 3
uBasis.CreateFinish()

pBasis = iron.Basis()
pBasis.CreateStart(pBasisUserNumber)
pBasis.type = iron.BasisTypes.SIMPLEX
pBasis.numberOfXi = 3
pBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX] * 3
pBasis.quadratureOrderSet = 3
pBasis.CreateFinish()

interfaceBasis = iron.Basis()
interfaceBasis.CreateStart(interfaceBasisUserNumber)
interfaceBasis.type = iron.BasisTypes.SIMPLEX
interfaceBasis.numberOfXi = 2
if (uInterpolation == LINEAR):
    interfaceBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX] * 2
elif (uInterpolation == QUADRATIC):
    interfaceBasis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX] * 2
else:
    print('Invalid u interpolation')
    exit()
interfaceBasis.quadratureOrderSet = 3
interfaceBasis.CreateFinish()

if (progressDiagnostics):
    print('Basis functions ... Done')

# ================================================================================================================================
#  Mesh
# ================================================================================================================================

if (progressDiagnostics):
    print('...Mesh Parameters=')
    print('...numberOfSolidNodes= %d' % (numberOfSolidNodes))
    print('...numberOfFluidNodes= %d' % (numberOfFluidNodes))
    print('...numberOfInterfaceNodes= %d' % (numberOfInterfaceNodes))
    print('...numberOfSolidElements= %d' % (numberOfSolidElements))
    print('...numberOfFluidElements= %d' % (numberOfFluidElements))
    print('...numberOfInterfaceElements= %d' % (numberOfInterfaceElements))

#the creation of fluid mesh
fluidNodes = iron.Nodes()
fluidNodes.CreateStart(fluidRegion, numberOfFluidNodes)
fluidNodes.UserNumbersAllSet(fluidNodesList)
fluidNodes.CreateFinish()

fluidMesh = iron.Mesh()
fluidMesh.CreateStart(fluidMeshUserNumber, fluidRegion, 3)
fluidMesh.NumberOfElementsSet(numberOfFluidElements)
fluidMesh.NumberOfComponentsSet(2)

if (debug):
    print('  Fluid Elements:')
fluidUElements = iron.MeshElements()
fluidUElements.CreateStart(fluidMesh, 1, uBasis)
fluidPElements = iron.MeshElements()
fluidPElements.CreateStart(fluidMesh, 2, pBasis)
for elementIdx in range(0, numberOfFluidElements):
    fluidUElements.NodesSet(elementIdx+1, fluidMatrix[elementIdx])
    fluidPElements.NodesSet(elementIdx+1, fluidMatrix[elementIdx])
fluidUElements.CreateFinish()
fluidPElements.CreateFinish()

fluidMesh.CreateFinish()

#the creation of solid mesh
solidNodes = iron.Nodes()
solidNodes.CreateStart(solidRegion, numberOfSolidNodes)
solidNodes.UserNumbersAllSet(solidNodesList)
solidNodes.CreateFinish()

solidMesh = iron.Mesh()
solidMesh.CreateStart(solidMeshUserNumber, solidRegion, 3)
solidMesh.NumberOfElementsSet(numberOfSolidElements)
solidMesh.NumberOfComponentsSet(2)

if (debug):
    print('  Solid Elements:')

solidUElements = iron.MeshElements()
solidUElements.CreateStart(solidMesh, 1, uBasis)
solidPElements = iron.MeshElements()
solidPElements.CreateStart(solidMesh, 2, pBasis)

for elementIdx in range(0,numberOfSolidElements):
    solidUElements.NodesSet(elementIdx+1, solidMatrix[elementIdx])
    solidPElements.NodesSet(elementIdx+1, solidMatrix[elementIdx])
solidUElements.CreateFinish()
solidPElements.CreateFinish()

solidMesh.CreateFinish()

if (progressDiagnostics):
    print('Meshes ... Done')

# ================================================================================================================================
#  Interface
# ================================================================================================================================

if (progressDiagnostics):
    print('Interface ...')

# Create an interface between the two meshes
interface = iron.Interface()
interface.CreateStart(interfaceUserNumber, iron.WorldRegion)
interface.LabelSet('Interface')
# Add in the two meshes
solidMeshIndex = interface.MeshAdd(solidMesh)
fluidMeshIndex = interface.MeshAdd(fluidMesh)
interface.CoordinateSystemSet(interfaceCoordinateSystem)
interface.CreateFinish()

if (progressDiagnostics):
    print('Interface ... Done')

# ================================================================================================================================
#  Interface Mesh
# ================================================================================================================================

if (progressDiagnostics):
    print('Interface Mesh ...')

# Create an interface mesh
interfaceNodes = iron.Nodes()
interfaceNodes.CreateStartInterface(interface, numberOfInterfaceNodes)
interfaceNodes.UserNumbersAllSet(interfaceNodesList)
interfaceNodes.CreateFinish()

interfaceMesh = iron.Mesh()
interfaceMesh.CreateStartInterface(interfaceMeshUserNumber, interface, 2)
interfaceMesh.NumberOfElementsSet(numberOfInterfaceElements)
interfaceMesh.NumberOfComponentsSet(1)

interfaceElements = iron.MeshElements()
interfaceElements.CreateStart(interfaceMesh, 1, interfaceBasis)
for elementIdx in range(0, numberOfInterfaceElements):
    interfaceElements.NodesSet(elementIdx+1, interfaceMatrix[elementIdx])
interfaceElements.CreateFinish()
interfaceMesh.CreateFinish()

if (progressDiagnostics):
    print('Interface Mesh ... Done')

# ================================================================================================================================
#  Mesh Connectivity
# ================================================================================================================================

if (progressDiagnostics):
    print('Interface Mesh Connectivity ...')

# Couple the interface meshes
interfaceMeshConnectivity = iron.InterfaceMeshConnectivity()
interfaceMeshConnectivity.CreateStart(interface, interfaceMesh)
interfaceMeshConnectivity.BasisSet(interfaceBasis)

numberOfLocalInterfaceNodes=3
interfaceNodes = [0]*(numberOfInterfaceNodes)
solidNodes = [0]*(numberOfInterfaceNodes)
fluidNodes = [0]*(numberOfInterfaceNodes)
interfacelocalNodes = [0] * numberOfLocalInterfaceNodes
solidlocalNodes = [0] * numberOfLocalInterfaceNodes
fluidlocalNodes = [0] * numberOfLocalInterfaceNodes
solidElementsOnInterface=[]
fluidElementsOnInterface=[]

for i in range(0, numberOfInterfaceElements):
    interfaceElementNumber = solidElementsConnection[i][0]
    solidElementNumber = solidElementsConnection[i][1]
    fluidElementNumber = fluidElementsConnection[i][1]
    print('  Interface Element %8d:' % (interfaceElementNumber))  

    # Map interface elements
    interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber, solidMeshIndex, solidElementNumber)
    interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber, fluidMeshIndex, fluidElementNumber)
    print('    Solid Element %8d;     Fluid Element %8d' % (solidElementNumber,fluidElementNumber))
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
    for z in range(0,numberOfFluidElements):
	if interfaceMatrix[i][j] in fluidMatrix[z]:
            if interfaceMatrix[i][j+1] in fluidMatrix[z]:
	        if interfaceMatrix[i][j+2] in fluidMatrix[z]:
		    fluidElementsOnInterface.append(fluidMatrix[z])
	            fluidlocalNodes=[0]*3
	            fluidlocalNodes[0]=fluidNodesList.index(interfaceMatrix[i][j])+1
	            fluidlocalNodes[1]=fluidNodesList.index(interfaceMatrix[i][j+1])+1
	            fluidlocalNodes[2]=fluidNodesList.index(interfaceMatrix[i][j+2])+1

# Map interface xi
    xi0=[0,1,1]
    xi1=[1,0,1]
    xi2=[1,1,0]
    xi3=[1,1,1]
    solidXi=[]
    fluidXi=[]
    for localNodeIdx in range(0,3):
        a=solidElementsOnInterface[i].index(interfaceMatrix[i][localNodeIdx])
        if a==0:
	    solidXi=xi0
        elif a==1:
	    solidXi=xi1
        elif a==2:
	    solidXi=xi2
        elif a==3:
	    solidXi=xi3
        b=fluidElementsOnInterface[i].index(interfaceMatrix[i][localNodeIdx])
        if b==0:
	    fluidXi=xi0
        elif b==1:
	    fluidXi=xi1
        elif b==2:
	    fluidXi=xi2
        elif b==3:
	    fluidXi=xi3

        interfaceNodes[interfacelocalNodes[localNodeIdx]-1] = interfacelocalNodes[localNodeIdx]
        solidNodes[interfacelocalNodes[localNodeIdx]-1] = solidlocalNodes[localNodeIdx]
        fluidNodes[interfacelocalNodes[localNodeIdx]-1] = fluidlocalNodes[localNodeIdx]
        interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,solidMeshIndex,solidElementNumber,
                                                           localNodeIdx+1,1,solidXi)
        interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber,
                                                           localNodeIdx+1,1,fluidXi)
        print('    Local node    %8d:' % (localNodeIdx ))
        print('      Interface node    %8d:' % (interfacelocalNodes[localNodeIdx]))
        print('      Solid node        %8d; Solid xi = [ %.2f, %.2f, %.2f ]' % \
                    (solidlocalNodes[localNodeIdx], solidXi[0], solidXi[1], solidXi[2]))
        print('      Fluid node        %8d; Fluid xi = [ %.2f, %.2f, %.2f ]' % \
                    (fluidlocalNodes[localNodeIdx], fluidXi[0], fluidXi[1], fluidXi[2]))


# Map interface nodes
interfaceMeshConnectivity.NodeNumberSet(interfaceNodes,solidMeshIndex,solidNodes,fluidMeshIndex,fluidNodes)

interfaceMeshConnectivity.CreateFinish()

if (progressDiagnostics):
    print('Interface Mesh Connectivity ... Done')

# ================================================================================================================================
#  Decomposition
# ================================================================================================================================

if (progressDiagnostics):
    print('Decomposition ...')

# Create a decomposition for the fluid mesh
fluidDecomposition = iron.Decomposition()
fluidDecomposition.CreateStart(fluidDecompositionUserNumber, fluidMesh)
fluidDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
fluidDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
fluidDecomposition.CalculateFacesSet(True)
fluidDecomposition.CreateFinish()

# Create a decomposition for the solid mesh
solidDecomposition = iron.Decomposition()
solidDecomposition.CreateStart(solidDecompositionUserNumber, solidMesh)
solidDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
solidDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
solidDecomposition.CalculateFacesSet(True)
solidDecomposition.CreateFinish()

# Create a decomposition for the interface mesh
interfaceDecomposition = iron.Decomposition()
interfaceDecomposition.CreateStart(interfaceDecompositionUserNumber, interfaceMesh)
interfaceDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
interfaceDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
interfaceDecomposition.CreateFinish()

if (progressDiagnostics):
    print('Decomposition ... Done')

# ================================================================================================================================
#  Geometric Field
# ================================================================================================================================

if (progressDiagnostics):
    print('Geometric Field ...')

# Start to create a default (geometric) field on the fluid region
fluidGeometricField = iron.Field()
fluidGeometricField.CreateStart(fluidGeometricFieldUserNumber, fluidRegion)
# Set the decomposition to use
fluidGeometricField.MeshDecompositionSet(fluidDecomposition)
# Set the scaling to use
fluidGeometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
fluidGeometricField.VariableLabelSet(iron.FieldVariableTypes.U, 'FluidGeometry')
# Set the domain to be used by the field components.
fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
# Finish creating the second field
fluidGeometricField.CreateFinish()

# Start to create a default (geometric) field on the solid region
solidGeometricField = iron.Field()
solidGeometricField.CreateStart(solidGeometricFieldUserNumber, solidRegion)
# Set the decomposition to use
solidGeometricField.MeshDecompositionSet(solidDecomposition)
solidGeometricField.meshDecomposition = solidDecomposition
# Set the scaling to use
solidGeometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
solidGeometricField.VariableLabelSet(iron.FieldVariableTypes.U, 'SolidGeometry')
# Set the domain to be used by the field components.
solidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
solidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
solidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
# Finish creating the first field
solidGeometricField.CreateFinish()

# Start to create a default (geometric) field on the Interface
interfaceGeometricField = iron.Field()
interfaceGeometricField.CreateStartInterface(interfaceGeometricFieldUserNumber, interface)
# Set the decomposition to use
interfaceGeometricField.MeshDecompositionSet(interfaceDecomposition)
# Set the scaling to use
interfaceGeometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
interfaceGeometricField.VariableLabelSet(iron.FieldVariableTypes.U, 'InterfaceGeometry')
# Set the domain to be used by the field components.
interfaceGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
interfaceGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
interfaceGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
# Finish creating the first field
interfaceGeometricField.CreateFinish()

if (progressDiagnostics):
    print('Geometric Field ... Done')

if (progressDiagnostics):
    print('Geometric Parameters ...')

if (debug):
    print('  Fluid Nodes:')
for nodeIdx in range(0, numberOfFluidNodes):
    nodeNumber=fluidNodePositions[nodeIdx][0]
    #nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber, 1)
    #if (nodeDomain == computationalNodeNumber):
    xPosition=fluidNodePositions[nodeIdx][1]
    yPosition=fluidNodePositions[nodeIdx][2]
    zPosition=fluidNodePositions[nodeIdx][3]
    fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                 iron.FieldParameterSetTypes.VALUES,
                                                                 1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                                                 nodeNumber, 1, xPosition)
    fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                 iron.FieldParameterSetTypes.VALUES,
                                                                 1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                                                 nodeNumber, 2, yPosition)
    fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                 iron.FieldParameterSetTypes.VALUES,
                                                                 1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                                                 nodeNumber, 3, zPosition)
    if (debug):
            print('      Node        %d:' % (nodeNumber))
            print('         Position         = [ %.2f, %.2f, %.2f ]' % (xPosition, yPosition, zPosition))

# Update fields
fluidGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
fluidGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

if (debug):
    print('  Solid Nodes:')
for nodeIdx in range(0, numberOfSolidNodes):
    nodeNumber=solidNodePositions[nodeIdx][0]
    #nodeDomain = solidDecomposition.NodeDomainGet(nodeNumber, 1)
    #if (nodeDomain == computationalNodeNumber):
    xPosition=solidNodePositions[nodeIdx][1]
    yPosition=solidNodePositions[nodeIdx][2]
    zPosition=solidNodePositions[nodeIdx][3]
    solidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES,
                                                             1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                                             nodeNumber, 1, xPosition)
    solidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES,
                                                             1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                                             nodeNumber, 2, yPosition)
    solidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES,
                                                             1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                                             nodeNumber, 3, zPosition)
    if (debug):
            print('      Node        %d:' % (nodeNumber))
            print('         Position         = [ %.2f, %.2f, %.2f ]' % (xPosition, yPosition, zPosition))

# Update fields
solidGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
solidGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

if (debug):
    print('  Interface Nodes:')
for nodeIdx in range(0, numberOfInterfaceNodes):
    nodeNumber=interfaceNodePositions[nodeIdx][0]
    xPosition=interfaceNodePositions[nodeIdx][1]
    yPosition=interfaceNodePositions[nodeIdx][2]
    zPosition=interfaceNodePositions[nodeIdx][3]
    interfaceGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES,
                                                             1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                                             nodeNumber, 1, xPosition)
    interfaceGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES,
                                                             1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                                             nodeNumber, 2, yPosition)
    interfaceGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES,
                                                             1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                                             nodeNumber, 3, zPosition)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Position         = [ %.2f, %.2f, %.2f ]' % (xPosition, yPosition, zPosition))

#Update fields
interfaceGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
interfaceGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Geometric Parameters ... Done')

# Export results
fluidFields = iron.Fields()
fluidFields.CreateRegion(fluidRegion)
fluidFields.NodesExport("FSITubeFluid", "FORTRAN")
fluidFields.ElementsExport("FSITubeFluid", "FORTRAN")
fluidFields.Finalise()
solidFields = iron.Fields()
solidFields.CreateRegion(solidRegion)
solidFields.NodesExport("FSITubeSolid", "FORTRAN")
solidFields.ElementsExport("FSITubeSolid", "FORTRAN")
solidFields.Finalise()
interfaceFields = iron.Fields()
interfaceFields.CreateInterface(interface)
interfaceFields.NodesExport("FSITubeInterface", "FORTRAN")
interfaceFields.ElementsExport("FSITubeInterface", "FORTRAN")
interfaceFields.Finalise()

# ================================================================================================================================
#  Equations Set
# ================================================================================================================================

if (progressDiagnostics):
    print('Equations Sets ...')

# Create the equations set for the fluid region - Navier-Stokes
fluidEquationsSetField = iron.Field()
fluidEquationsSet = iron.EquationsSet()
if RBS:
    fluidEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                                      iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                      iron.EquationsSetSubtypes.ALE_RBS_NAVIER_STOKES]
else:
    fluidEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                                      iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                      iron.EquationsSetSubtypes.ALE_NAVIER_STOKES]
fluidEquationsSet.CreateStart(fluidEquationsSetUserNumber, fluidRegion, fluidGeometricField,
                              fluidEquationsSetSpecification, fluidEquationsSetFieldUserNumber,
                              fluidEquationsSetField)
fluidEquationsSet.OutputTypeSet(fluidEquationsSetOutputType)
fluidEquationsSet.CreateFinish()

if RBS:
    # Set max CFL number (default 1.0)
    fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                       iron.FieldParameterSetTypes.VALUES, 2, 1.0E20)
    # Set time increment (default 0.0)
    fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                       iron.FieldParameterSetTypes.VALUES, 3, timeStep)
    # Set stabilisation type (default 1.0 = RBS)
    fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                       iron.FieldParameterSetTypes.VALUES, 4, 1.0)

# Create the equations set for the solid region
solidEquationsSetField = iron.Field()
solidEquationsSet = iron.EquationsSet()
solidEquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                  iron.EquationsSetTypes.FINITE_ELASTICITY,
                                  iron.EquationsSetSubtypes.MOONEY_RIVLIN]
solidEquationsSet.CreateStart(solidEquationsSetUserNumber, solidRegion, solidGeometricField,
                              solidEquationsSetSpecification, solidEquationsSetFieldUserNumber,
                              solidEquationsSetField)
solidEquationsSet.OutputTypeSet(solidEquationsSetOutputType)
solidEquationsSet.CreateFinish()

# Create the equations set for the moving mesh
movingMeshEquationsSetField = iron.Field()
movingMeshEquationsSet = iron.EquationsSet()
movingMeshEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.LAPLACE_EQUATION,
                                       iron.EquationsSetSubtypes.MOVING_MESH_LAPLACE]
movingMeshEquationsSet.CreateStart(movingMeshEquationsSetUserNumber, fluidRegion, fluidGeometricField,
                                   movingMeshEquationsSetSpecification, movingMeshEquationsSetFieldUserNumber,
                                   movingMeshEquationsSetField)
movingMeshEquationsSet.OutputTypeSet(movingMeshEquationsSetOutputType)
movingMeshEquationsSet.CreateFinish()

if (progressDiagnostics):
    print('Equations Sets ... Done')

# ================================================================================================================================
#  Dependent Field
# ================================================================================================================================

if (progressDiagnostics):
    print('Dependent Fields ...')

# Create the equations set dependent field variables for dynamic Navier-Stokes
fluidDependentField = iron.Field()
fluidEquationsSet.DependentCreateStart(fluidDependentFieldUserNumber, fluidDependentField)
fluidDependentField.VariableLabelSet(iron.FieldVariableTypes.U, 'FluidDependent')
# Set the mesh component to be used by the field components.
for componentIdx in range(1, 4):
    fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, componentIdx, 1)
    fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, componentIdx, 1)
fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 4, 2)
fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 4, 2)
# fluidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.NODE_BASED)
# fluidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.NODE_BASED)
# Finish the equations set dependent field variables
fluidEquationsSet.DependentCreateFinish()

# Initialise the fluid dependent field
for componentIdx in range(1, 4):
    fluidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                                    componentIdx, 0.0)
# Initialise pressure component
fluidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 4,
                                                fluidPInit)

# Create the equations set dependent field variables for the solid equations set
solidDependentField = iron.Field()
solidEquationsSet.DependentCreateStart(solidDependentFieldUserNumber, solidDependentField)
solidDependentField.VariableLabelSet(iron.FieldVariableTypes.U, 'SolidDependent')
solidDependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN, 'SolidTraction')
for componentIdx in range(1, 4):
    solidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, componentIdx, 1)
    solidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, componentIdx, 1)
solidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 4, 2)
solidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 4, 2)
solidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 4, iron.FieldInterpolationTypes.NODE_BASED)
solidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN, 4,
                                              iron.FieldInterpolationTypes.NODE_BASED)
solidDependentField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
solidEquationsSet.DependentCreateFinish()

# Initialise the solid dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
for componentIdx in range(1, 4):
    solidGeometricField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,
                                                                 iron.FieldParameterSetTypes.VALUES, \
                                                                 componentIdx, solidDependentField,
                                                                 iron.FieldVariableTypes.U,
                                                                 iron.FieldParameterSetTypes.VALUES, componentIdx)
solidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 4,
                                                solidPInit)

# Create the equations set dependent field variables for moving mesh
movingMeshDependentField = iron.Field()
movingMeshEquationsSet.DependentCreateStart(movingMeshDependentFieldUserNumber, movingMeshDependentField)
movingMeshDependentField.VariableLabelSet(iron.FieldVariableTypes.U, 'MovingMeshDependent')
# Set the mesh component to be used by the field components.
for componentIdx in range(1, numberOfDimensions + 1):
    movingMeshDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, componentIdx, 1)
    movingMeshDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, componentIdx, 1)
# Finish the equations set dependent field variables
movingMeshEquationsSet.DependentCreateFinish()

# Initialise dependent field moving mesh
for componentIdx in range(1, numberOfDimensions + 1):
    movingMeshDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, \
                                                         componentIdx, 0.0)

fluidDependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
solidDependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
movingMeshDependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
fluidDependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
solidDependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
movingMeshDependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Dependent Fields ... Done')

# ================================================================================================================================
#  Materials Field
# ================================================================================================================================

if (progressDiagnostics):
    print('Materials Fields ...')

# Create the equations set materials field variables for dynamic Navier-Stokes
fluidMaterialsField = iron.Field()
fluidEquationsSet.MaterialsCreateStart(fluidMaterialsFieldUserNumber, fluidMaterialsField)
# Finish the equations set materials field variables
fluidEquationsSet.MaterialsCreateFinish()
fluidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                                fluidDynamicViscosity)
fluidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2,
                                                fluidDensity)

# Create the solid materials field
solidMaterialsField = iron.Field()
solidEquationsSet.MaterialsCreateStart(solidMaterialsFieldUserNumber, solidMaterialsField)
solidMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U, 'SolidMaterials')
solidMaterialsField.VariableLabelSet(iron.FieldVariableTypes.V, 'SolidDensity')
solidEquationsSet.MaterialsCreateFinish()
# Set Mooney-Rivlin constants c10 and c01 respectively
solidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                                mooneyRivlin1)
solidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2,
                                                mooneyRivlin2)
solidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V, iron.FieldParameterSetTypes.VALUES, 1,
                                                solidDensity)

# Create the equations set materials field variables for moving mesh
movingMeshMaterialsField = iron.Field()
movingMeshEquationsSet.MaterialsCreateStart(movingMeshMaterialsFieldUserNumber, movingMeshMaterialsField)
# Finish the equations set materials field variables
movingMeshEquationsSet.MaterialsCreateFinish()

movingMeshMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, \
                                                     movingMeshKParameter)

if (progressDiagnostics):
    print('Materials Fields ... Done')

# ================================================================================================================================
# Independent Field
# ================================================================================================================================

if (progressDiagnostics):
    print('Independent Fields ...')

# Create fluid mesh velocity independent field
fluidIndependentField = iron.Field()
fluidEquationsSet.IndependentCreateStart(fluidIndependentFieldUserNumber, fluidIndependentField)
fluidIndependentField.VariableLabelSet(iron.FieldVariableTypes.U, 'FluidIndependent')
# Set the mesh component to be used by the field components.
for componentIdx in range(1, numberOfDimensions + 1):
    fluidIndependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, componentIdx, 1)
# Finish the equations set independent field variables
fluidEquationsSet.IndependentCreateFinish()

# Create the moving mesh independent field
movingMeshIndependentField = iron.Field()
movingMeshEquationsSet.IndependentCreateStart(movingMeshIndependentFieldUserNumber, movingMeshIndependentField)
movingMeshIndependentField.VariableLabelSet(iron.FieldVariableTypes.U, 'MovingMeshIndependent')
# Set the mesh component to be used by the field components.
for componentIdx in range(1, numberOfDimensions + 1):
    movingMeshIndependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, componentIdx, 1)
# Finish the equations set independent field variables
movingMeshEquationsSet.IndependentCreateFinish()

# Initialise independent field moving mesh
movingMeshIndependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
                                                       movingMeshKParameter)

if (progressDiagnostics):
    print('Independent Fields ... Done')

# ================================================================================================================================
#  Equations
# ================================================================================================================================

if (progressDiagnostics):
    print('Equations ...')

# Fluid equations
fluidEquations = iron.Equations()
fluidEquationsSet.EquationsCreateStart(fluidEquations)
fluidEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
fluidEquations.outputType = fluidEquationsOutputType
fluidEquationsSet.EquationsCreateFinish()

# Solid equations
solidEquations = iron.Equations()
solidEquationsSet.EquationsCreateStart(solidEquations)
solidEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
solidEquations.outputType = solidEquationsOutputType
solidEquationsSet.EquationsCreateFinish()

# Moving mesh equations
movingMeshEquations = iron.Equations()
movingMeshEquationsSet.EquationsCreateStart(movingMeshEquations)
movingMeshEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
movingMeshEquations.outputType = movingMeshEquationsOutputType
movingMeshEquationsSet.EquationsCreateFinish()

if (progressDiagnostics):
    print('Equations ... Done')

# ================================================================================================================================
#  CellML
# ================================================================================================================================

if (progressDiagnostics):
    print('CellML ...')

# Create CellML equations for the temporal boundary conditions
bcCellML = iron.CellML()
bcCellML.CreateStart(bcCellMLUserNumber, fluidRegion)
bcCellMLIdx = bcCellML.ModelImport("input/poiseuilleinlet.cellml")
bcCellML.VariableSetAsKnown(bcCellMLIdx, "main/pipeRadius")
bcCellML.VariableSetAsKnown(bcCellMLIdx, "main/length")
bcCellML.VariableSetAsKnown(bcCellMLIdx, "main/dynamicViscosity")
bcCellML.VariableSetAsKnown(bcCellMLIdx, "main/A")
bcCellML.VariableSetAsKnown(bcCellMLIdx, "main/B")
bcCellML.VariableSetAsKnown(bcCellMLIdx, "main/C")
bcCellML.VariableSetAsKnown(bcCellMLIdx, "main/x")
bcCellML.VariableSetAsKnown(bcCellMLIdx, "main/y")
bcCellML.VariableSetAsWanted(bcCellMLIdx, "main/inletx")
bcCellML.VariableSetAsWanted(bcCellMLIdx, "main/inlety")
bcCellML.VariableSetAsWanted(bcCellMLIdx, "main/inletz")
bcCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps
bcCellML.FieldMapsCreateStart()
# Map geometric field to x and y
bcCellML.CreateFieldToCellMLMap(fluidGeometricField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                                bcCellMLIdx, "main/x", iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateFieldToCellMLMap(fluidGeometricField, iron.FieldVariableTypes.U, 2, iron.FieldParameterSetTypes.VALUES,
                                bcCellMLIdx, "main/y", iron.FieldParameterSetTypes.VALUES)
# Map fluid velocity to ensure dependent field isn't cleared when the velocities are copied back
bcCellML.CreateFieldToCellMLMap(fluidDependentField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES,
                                bcCellMLIdx, "main/inletx", iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateFieldToCellMLMap(fluidDependentField, iron.FieldVariableTypes.U, 2, iron.FieldParameterSetTypes.VALUES,
                                bcCellMLIdx, "main/inlety", iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateFieldToCellMLMap(fluidDependentField, iron.FieldVariableTypes.U, 3, iron.FieldParameterSetTypes.VALUES,
                                bcCellMLIdx, "main/inletz", iron.FieldParameterSetTypes.VALUES)
# Map inletx, inlety and inletz to dependent field
bcCellML.CreateCellMLToFieldMap(bcCellMLIdx, "main/inletx", iron.FieldParameterSetTypes.VALUES,
                                fluidDependentField, iron.FieldVariableTypes.U, 1, iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateCellMLToFieldMap(bcCellMLIdx, "main/inlety", iron.FieldParameterSetTypes.VALUES,
                                fluidDependentField, iron.FieldVariableTypes.U, 2, iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateCellMLToFieldMap(bcCellMLIdx, "main/inletz", iron.FieldParameterSetTypes.VALUES,
                                fluidDependentField, iron.FieldVariableTypes.U, 3, iron.FieldParameterSetTypes.VALUES)
bcCellML.FieldMapsCreateFinish()

# Create the CellML models field
bcCellMLModelsField = iron.Field()
bcCellML.ModelsFieldCreateStart(bcCellMLModelsFieldUserNumber, bcCellMLModelsField)
bcCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U, "BCModelMap")
bcCellML.ModelsFieldCreateFinish()

# Only evaluate BC on inlet nodes
bcCellMLModelsField.ComponentValuesInitialiseIntg(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 0)
if (debug):
    print('  CellML Boundary Conditions:')
    print('    Inlet Model Set:')
for nodeIdx in range(0,len(fluidInletNodes)):
    nodeNumber = fluidInletNodes[nodeIdx]
    #nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber, 1)
    #if (nodeDomain == computationalNodeNumber):
    bcCellMLModelsField.ParameterSetUpdateNodeIntg(iron.FieldVariableTypes.U,
                                                               iron.FieldParameterSetTypes.VALUES,
                                                               1, iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                                               nodeNumber, 1, 1)
    if (debug):
        print('      Node        %d:' % (nodeNumber))

bcCellMLModelsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
bcCellMLModelsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

# Create the CellML state field
bcCellMLStateField = iron.Field()
bcCellML.StateFieldCreateStart(bcCellMLStateFieldUserNumber, bcCellMLStateField)
bcCellMLStateField.VariableLabelSet(iron.FieldVariableTypes.U, "BCState")
bcCellML.StateFieldCreateFinish()

# Create the CellML parameters field
bcCellMLParametersField = iron.Field()
bcCellML.ParametersFieldCreateStart(bcCellMLParametersFieldUserNumber, bcCellMLParametersField)
bcCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U, "BCParameters")
bcCellML.ParametersFieldCreateFinish()

# Get the component numbers
pipeRadiusComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx, iron.CellMLFieldTypes.PARAMETERS, "main/pipeRadius")
lengthComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx, iron.CellMLFieldTypes.PARAMETERS, "main/length")
dynamicViscosityComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx, iron.CellMLFieldTypes.PARAMETERS,
                                                             "main/dynamicViscosity")
AComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx, iron.CellMLFieldTypes.PARAMETERS, "main/A")
BComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx, iron.CellMLFieldTypes.PARAMETERS, "main/B")
CComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx, iron.CellMLFieldTypes.PARAMETERS, "main/C")
# Set up the parameters field
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, \
                                                    pipeRadiusComponentNumber, pipeRadius)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, \
                                                    lengthComponentNumber, lengthSize)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, \
                                                    dynamicViscosityComponentNumber, fluidDynamicViscosity)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, \
                                                    AComponentNumber, A)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, \
                                                    BComponentNumber, B)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, \
                                                    CComponentNumber, C)

# Create the CELL intermediate field
bcCellMLIntermediateField = iron.Field()
bcCellML.IntermediateFieldCreateStart(bcCellMLIntermediateFieldUserNumber, bcCellMLIntermediateField)
bcCellMLIntermediateField.VariableLabelSet(iron.FieldVariableTypes.U, "BCIntermediate")
bcCellML.IntermediateFieldCreateFinish()

if (progressDiagnostics):
    print('CellML ... Done')

# ================================================================================================================================
#  Interface Condition
# ================================================================================================================================

if (progressDiagnostics):
    print('Interface Conditions ...')

# Create an interface condition between the two meshes
interfaceCondition = iron.InterfaceCondition()
interfaceCondition.CreateStart(interfaceConditionUserNumber, interface, interfaceGeometricField)
# Specify the method for the interface condition
interfaceCondition.MethodSet(iron.InterfaceConditionMethods.LAGRANGE_MULTIPLIERS)
# Specify the type of interface condition operator
interfaceCondition.OperatorSet(iron.InterfaceConditionOperators.SOLID_FLUID)
# Add in the dependent variables from the equations sets
interfaceCondition.DependentVariableAdd(solidMeshIndex, solidEquationsSet, iron.FieldVariableTypes.U)
interfaceCondition.DependentVariableAdd(fluidMeshIndex, fluidEquationsSet, iron.FieldVariableTypes.U)
# Set the label
interfaceCondition.LabelSet("FSI Interface Condition")
# Set the output type
interfaceCondition.OutputTypeSet(interfaceConditionOutputType)
# Finish creating the interface condition
interfaceCondition.CreateFinish()

if (progressDiagnostics):
    print('Interface Conditions ... Done')

if (progressDiagnostics):
    print('Interface Lagrange Field ...')

# Create the Lagrange multipliers field
interfaceLagrangeField = iron.Field()
interfaceCondition.LagrangeFieldCreateStart(interfaceLagrangeFieldUserNumber, interfaceLagrangeField)
interfaceLagrangeField.VariableLabelSet(iron.FieldVariableTypes.U, 'InterfaceLagrange')
# Finish the Lagrange multipliers field
interfaceCondition.LagrangeFieldCreateFinish()

for componentIdx in range(1, numberOfDimensions + 1):
    interfaceLagrangeField.ComponentValuesInitialise(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                                                     componentIdx, 0.0)

interfaceLagrangeField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
interfaceLagrangeField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Interface Lagrange Field ... Done')

if (progressDiagnostics):
    print('Interface Equations ...')

# Create the interface condition equations
interfaceEquations = iron.InterfaceEquations()
interfaceCondition.EquationsCreateStart(interfaceEquations)
# Set the interface equations sparsity
interfaceEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
# Set the interface equations output
interfaceEquations.outputType = interfaceEquationsOutputType
# Finish creating the interface equations
interfaceCondition.EquationsCreateFinish()

if (progressDiagnostics):
    print('Interface Equations ... Done')

# ================================================================================================================================
#  Problem
# ================================================================================================================================

if (progressDiagnostics):
    print('Problems ...')

# Create a FSI problem
fsiProblem = iron.Problem()
if RBS:
    fsiProblemSpecification = [iron.ProblemClasses.MULTI_PHYSICS,
                               iron.ProblemTypes.FINITE_ELASTICITY_NAVIER_STOKES,
                               iron.ProblemSubtypes.FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE]
else:
    fsiProblemSpecification = [iron.ProblemClasses.MULTI_PHYSICS,
                               iron.ProblemTypes.FINITE_ELASTICITY_NAVIER_STOKES,
                               iron.ProblemSubtypes.FINITE_ELASTICITY_NAVIER_STOKES_ALE]
fsiProblem.CreateStart(fsiProblemUserNumber, fsiProblemSpecification)
fsiProblem.CreateFinish()

if (progressDiagnostics):
    print('Problems ... Done')

# ================================================================================================================================
#  Control Loop
# ================================================================================================================================

if (progressDiagnostics):
    print('Control Loops ...')

# Create the fsi problem control loop
fsiControlLoop = iron.ControlLoop()
fsiProblem.ControlLoopCreateStart()
fsiProblem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], fsiControlLoop)
fsiControlLoop.LabelSet('TimeLoop')
fsiControlLoop.TimesSet(startTime, stopTime, timeStep)
fsiControlLoop.TimeOutputSet(outputFrequency)
fsiProblem.ControlLoopCreateFinish()

if (progressDiagnostics):
    print('Control Loops ... Done')

# ================================================================================================================================
#  Solvers
# ================================================================================================================================

if (progressDiagnostics):
    print('Solvers ...')

# Create the problem solver
bcCellMLEvaluationSolver = iron.Solver()
fsiDynamicSolver = iron.Solver()
fsiNonlinearSolver = iron.Solver()
fsiLinearSolver = iron.Solver()
movingMeshLinearSolver = iron.Solver()

fsiProblem.SolversCreateStart()
# Solvers for coupled FiniteElasticity NavierStokes problem
# Get the BC CellML solver
fsiProblem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, bcCellMLEvaluationSolver)
bcCellMLEvaluationSolver.outputType = iron.SolverOutputTypes.PROGRESS
# Get the dynamic ALE solver
fsiProblem.SolverGet([iron.ControlLoopIdentifiers.NODE], 2, fsiDynamicSolver)
fsiDynamicSolver.OutputTypeSet(fsiDynamicSolverOutputType)
fsiDynamicSolver.DynamicThetaSet(fsiDynamicSolverTheta)
# Get the dynamic nonlinear solver
fsiDynamicSolver.DynamicNonlinearSolverGet(fsiNonlinearSolver)
fsiNonlinearSolver.NewtonLineSearchTypeSet(iron.NewtonLineSearchTypes.LINEAR)
# fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS) #(.FD/EQUATIONS)
fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)  # (.FD/EQUATIONS)
fsiNonlinearSolver.NewtonMaximumFunctionEvaluationsSet(nonlinearMaxFunctionEvaluations)
fsiNonlinearSolver.OutputTypeSet(fsiNonlinearSolverOutputType)
fsiNonlinearSolver.NewtonAbsoluteToleranceSet(nonlinearAbsoluteTolerance)
fsiNonlinearSolver.NewtonMaximumIterationsSet(nonlinearMaximumIterations)
fsiNonlinearSolver.NewtonRelativeToleranceSet(nonlinearRelativeTolerance)
fsiNonlinearSolver.NewtonLineSearchAlphaSet(nonlinearLinesearchAlpha)
# Get the dynamic nonlinear linear solver
fsiNonlinearSolver.NewtonLinearSolverGet(fsiLinearSolver)
# fsiLinearSolver.LinearTypeSet(iron.LinearSolverTypes.ITERATIVE)
# fsiLinearSolver.LinearIterativeMaximumIterationsSet(linearMaximumIterations)
# fsiLinearSolver.LinearIterativeDivergenceToleranceSet(linearDivergenceTolerance)
# fsiLinearSolver.LinearIterativeRelativeToleranceSet(linearRelativeTolerance)
# fsiLinearSolver.LinearIterativeAbsoluteToleranceSet(linearAbsoluteTolerance)
fsiLinearSolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
fsiLinearSolver.OutputTypeSet(fsiLinearSolverOutputType)
# Linear solver for moving mesh
fsiProblem.SolverGet([iron.ControlLoopIdentifiers.NODE], 3, movingMeshLinearSolver)
movingMeshLinearSolver.OutputTypeSet(movingMeshLinearSolverOutputType)
# Finish the creation of the problem solver
fsiProblem.SolversCreateFinish()

if (progressDiagnostics):
    print('Solvers ... Done')

# ================================================================================================================================
#  CellML Equations
# ================================================================================================================================

if (progressDiagnostics):
    print('CellML Equations ...')

# Create CellML equations and add BC equations to the solver
bcEquations = iron.CellMLEquations()
fsiProblem.CellMLEquationsCreateStart()
bcCellMLEvaluationSolver.CellMLEquationsGet(bcEquations)
bcEquationsIndex = bcEquations.CellMLAdd(bcCellML)
fsiProblem.CellMLEquationsCreateFinish()

if (progressDiagnostics):
    print('CellML Equations ... Done')

# ================================================================================================================================
#  Solver Equations
# ================================================================================================================================

if (progressDiagnostics):
    print('Solver Equations ...')

# Start the creation of the fsi problem solver equations
fsiProblem.SolverEquationsCreateStart()
# Get the fsi dynamic solver equations
fsiSolverEquations = iron.SolverEquations()
fsiDynamicSolver.SolverEquationsGet(fsiSolverEquations)
fsiSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
# fsiSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.FULL
fsiSolidEquationsSetIndex = fsiSolverEquations.EquationsSetAdd(solidEquationsSet)
fsiFluidEquationsSetIndex = fsiSolverEquations.EquationsSetAdd(fluidEquationsSet)
fsiInterfaceConditionIndex = fsiSolverEquations.InterfaceConditionAdd(interfaceCondition)
# Set the time dependence of the interface matrix to determine the interface matrix coefficient in the solver matrix
# (basicly position in big coupled matrix system)
interfaceEquations.MatrixTimeDependenceTypeSet(fsiSolidEquationsSetIndex, True, \
                                               [iron.InterfaceMatricesTimeDependenceTypes.STATIC, \
                                                iron.InterfaceMatricesTimeDependenceTypes.FIRST_ORDER_DYNAMIC])
interfaceEquations.MatrixTimeDependenceTypeSet(fsiFluidEquationsSetIndex, True, \
                                               [iron.InterfaceMatricesTimeDependenceTypes.STATIC, \
                                                iron.InterfaceMatricesTimeDependenceTypes.STATIC])

# Create the moving mesh solver equations
movingMeshSolverEquations = iron.SolverEquations()
# Get the linear moving mesh solver equations
movingMeshLinearSolver.SolverEquationsGet(movingMeshSolverEquations)
movingMeshSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
# Add in the equations set
movingMeshEquationsSetIndex = movingMeshSolverEquations.EquationsSetAdd(movingMeshEquationsSet)

# Finish the creation of the fsi problem solver equations
fsiProblem.SolverEquationsCreateFinish()

if (progressDiagnostics):
    print('Solver Equations ...')

# ================================================================================================================================
#  Boundary Conditions
# ================================================================================================================================

if (progressDiagnostics):
    print('Boundary Conditions ...')

# Start the creation of the fluid boundary conditions
fsiBoundaryConditions = iron.BoundaryConditions()
fsiSolverEquations.BoundaryConditionsCreateStart(fsiBoundaryConditions)

if (debug):
    print('  Fluid Boundary Conditions:')
    print('    Inlet Boundary conditions:')
# Set inlet boundary conditions on the left hand edge
for nodeIdx in range(0,len(fluidInletNodes)):
    nodeNumber = fluidInletNodes[nodeIdx]
    #nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber, 1)
    #if (nodeDomain == computationalNodeNumber):
    fsiBoundaryConditions.SetNode(fluidDependentField, iron.FieldVariableTypes.U, 1, \
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                              nodeNumber, 1, iron.BoundaryConditionsTypes.FIXED_INLET, 0.0)
    fsiBoundaryConditions.SetNode(fluidDependentField, iron.FieldVariableTypes.U, 1, \
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                              nodeNumber, 2, iron.BoundaryConditionsTypes.FIXED_INLET, 0.0)
    fsiBoundaryConditions.SetNode(fluidDependentField, iron.FieldVariableTypes.U, 1, \
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                              nodeNumber, 3, iron.BoundaryConditionsTypes.FIXED_INLET, 0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Velocity         = [ %.2f, %.2f, %.2f ]' % (0.0, 0.0, 0.0))

if (debug):
    print('    Wall Boundary conditions:')
# Set no slip boundary conditions on the wall
for nodeIdx in range(0, len(interfaceNodesList)):
    nodeNumber = interfaceNodesList[nodeIdx]
    #nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber, 1)
    #if (nodeDomain == computationalNodeNumber):
    fsiBoundaryConditions.SetNode(fluidDependentField, iron.FieldVariableTypes.U, 1, \
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                              nodeNumber, 3, iron.BoundaryConditionsTypes.FIXED, 0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Velocity         = [ ????, ????, %.2f ]' % (0.0))

if (debug):
    print('    No Pressure Boundary conditions:')

nodeNumber = interfaceOutletNodes[3]
fsiBoundaryConditions.SetNode(fluidDependentField, iron.FieldVariableTypes.U, 1, \
                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                              nodeNumber, 4, iron.BoundaryConditionsTypes.PRESSURE, 0.0)

# Set no pressure boundary conditions on the outlet
for nodeIdx in range(0, len(fluidOutletNodes)):
    nodeNumber = fluidOutletNodes[nodeIdx] 
    #nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber, 2)
    #if (nodeDomain == computationalNodeNumber):
    fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                        iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                        nodeNumber,4,iron.BoundaryConditionsTypes.FIXED_OUTLET,0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Pressure         = %.2f' % (0.0))

if (debug):
    print('  Solid Boundary conditions:')
# Set solid boundary conditions
 #Set nodes on the axis to only slide on the axis.
#Top y-axis node - fix in the x-direction
for nodeIdx in range (0,2):
    M1=[5,37]
    nodeNumber = M1[nodeIdx]
    #nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
    #if (nodeDomain == computationalNodeNumber):
    fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                      nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Displacement     = [ %.2f, ????, ???? ]' % (0.0))  
                   
#Right x-axis node - fix in the y-direction
for nodeIdx in range (0,2):
    M2=[1,33]
    nodeNumber = M2[nodeIdx]
    #nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
    #if (nodeDomain == computationalNodeNumber):
    fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                      nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Displacement     = [ ????, %.2f, ???? ]' % (0.0))   
                       
#Bottom y-axis node - fix in the x-direction
for nodeIdx in range (0,2):
    M3=[13,45]
    nodeNumber = M3[nodeIdx]
    #nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
    #if (nodeDomain == computationalNodeNumber):
    fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                      nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Displacement     = [ %.2f, ????, ???? ]' % (0.0)) 
                       
#Left x-axis node - fix in the y-direction
for nodeIdx in range (0,2):
    M4=[9,41]
    nodeNumber = M4[nodeIdx]
    #nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
    #if (nodeDomain == computationalNodeNumber):
    fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                      nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Displacement     = [ ????, %.2f, ???? ]' % (0.0))                         

# Fix all solid nodes at the beginning of the tube in the z-direction
for nodeIdx in range(0, len(solidInletNodes)):
    nodeNumber = solidInletNodes[nodeIdx]
    #nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber, 1)
    #if (nodeDomain == computationalNodeNumber):
    fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                      nodeNumber,3,iron.BoundaryConditionsTypes.FIXED, 0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Displacement     = [ ????, ????, %.2f ]' % (0.0))

if (debug):
    print('  Lagrange Boundary conditions:')
# Remove Lagrange multipliers where solid displacement and fluid velocity is zero
nodeNumber = 37
#nodeDomain = computationalNodeNumber
#if (nodeDomain == computationalNodeNumber):
fsiBoundaryConditions.SetNode(interfaceLagrangeField, iron.FieldVariableTypes.U, 1, \
                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                  nodeNumber, 1, iron.BoundaryConditionsTypes.FIXED, 0.0)
if (debug):
    print('      Node        %d:' % (nodeNumber))
    print('         Lagrange         = [ %.2f, ????, ???? ]' % (0.0))
nodeNumber = 41
#nodeDomain = computationalNodeNumber
#if (nodeDomain == computationalNodeNumber):
fsiBoundaryConditions.SetNode(interfaceLagrangeField, iron.FieldVariableTypes.U, 1, \
                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                  nodeNumber, 2, iron.BoundaryConditionsTypes.FIXED, 0.0)
if (debug):
    print('      Node        %d:' % (nodeNumber))
    print('         Lagrange         = [ ????, %.2f, ???? ]' % (0.0))
nodeNumber = 45
#nodeDomain = computationalNodeNumber
#if (nodeDomain == computationalNodeNumber):
fsiBoundaryConditions.SetNode(interfaceLagrangeField, iron.FieldVariableTypes.U, 1, \
                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                  nodeNumber, 1, iron.BoundaryConditionsTypes.FIXED, 0.0)
if (debug):
    print('      Node        %d:' % (nodeNumber))
    print('         Lagrange         = [ %.2f, ????, ???? ]' % (0.0))
nodeNumber = 33
#nodeDomain = computationalNodeNumber
#if (nodeDomain == computationalNodeNumber):
fsiBoundaryConditions.SetNode(interfaceLagrangeField, iron.FieldVariableTypes.U, 1, \
                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                  nodeNumber, 2, iron.BoundaryConditionsTypes.FIXED, 0.0)
if (debug):
    print('      Node        %d:' % (nodeNumber))
    print('         Lagrange         = [ ????, %.2f, ???? ]' % (0.0))

for nodeIdx in range(0, len(interfaceInletNodes)):
    nodeNumber = interfaceInletNodes[nodeIdx]
    #nodeDomain = computationalNodeNumber
    #if (nodeDomain == computationalNodeNumber):
    fsiBoundaryConditions.SetNode(interfaceLagrangeField, iron.FieldVariableTypes.U, 1, \
                                      iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                      nodeNumber, 3, iron.BoundaryConditionsTypes.FIXED, 0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Lagrange         = [ ????, ????, %.2f ]' % (0.0))

    # Finish fsi boundary conditions
fsiSolverEquations.BoundaryConditionsCreateFinish()


pdb.set_trace()
# Start the creation of the moving mesh boundary conditions
movingMeshBoundaryConditions = iron.BoundaryConditions()
movingMeshSolverEquations.BoundaryConditionsCreateStart(movingMeshBoundaryConditions)
if (debug):
    print('  Moving Mesh Boundary Conditions:')
    print('    Fixed Wall Boundary conditions:')
# Inlet nodes
for nodeIdx in range(0,len(fluidInletNodes)):
    nodeNumber = fluidInletNodes[nodeIdx]
    #nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber, 1)
    #if (nodeDomain == computationalNodeNumber):
    movingMeshBoundaryConditions.AddNode(movingMeshDependentField, iron.FieldVariableTypes.U, 1,
                                                     iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, nodeNumber, 3,
                                                     iron.BoundaryConditionsTypes.FIXED_WALL, 0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))

if (debug):
    print('    Moving Wall Boundary conditions:')
for nodeIdx in range(0, len(interfaceNodesList)):
    nodeNumber = interfaceNodesList[nodeIdx]
    #nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber, 1)
    #if (nodeDomain == computationalNodeNumber):
    movingMeshBoundaryConditions.AddNode(movingMeshDependentField, iron.FieldVariableTypes.U, 1, \
                                                     iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                     nodeNumber, 1, iron.BoundaryConditionsTypes.FIXED, 0.0)
    movingMeshBoundaryConditions.AddNode(movingMeshDependentField, iron.FieldVariableTypes.U, 1, \
                                                     iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                     nodeNumber, 2, iron.BoundaryConditionsTypes.FIXED, 0.0)
    movingMeshBoundaryConditions.AddNode(movingMeshDependentField, iron.FieldVariableTypes.U, 1, \
                                                     iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                     nodeNumber, 3, iron.BoundaryConditionsTypes.FIXED, 0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))

# Finish moving mesh boundary conditions
movingMeshSolverEquations.BoundaryConditionsCreateFinish()

if (progressDiagnostics):
    print('Boundary Conditions ... Done')

# ================================================================================================================================
#  Run Solvers
# ================================================================================================================================

# quit()

# Export results
fields = iron.Fields()
fields.CreateRegion(fluidRegion)
fields.NodesExport("FSITube", "FORTRAN")
fields.ElementsExport("FSITube", "FORTRAN")
fields.Finalise()

fsiLinearSolver.MumpsSetIcntl(14, 20000)

# Create output directories
if not os.path.exists("output/Fluid"):
    os.makedirs("output/Fluid")
if not os.path.exists("output/Solid"):
    os.makedirs("output/Solid")
if not os.path.exists("output/Interface"):
    os.makedirs("output/Interface")

# Solve the problem
print('Solving problem...')
start = time.time()
fsiProblem.Solve()
end = time.time()
elapsed = end - start
print('Calculation Time = %3.4f' % elapsed)
print('Problem solved!')
print('#')
