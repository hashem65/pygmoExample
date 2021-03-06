'''
Created on 7/06/2018

@author: rjag008
'''
from __future__ import print_function
import numpy as np
from opencmiss.iron import iron
import exfile
from scipy.interpolate import griddata
import math

humphreyModel = False
    
class TubeGrower(object):
    '''
    classdocs
    '''
    if (humphreyModel):
        c1 = 3.2
        c2 = 0.39
        c3 = 11.0
        c4 = 0.49
    maxSolvebleRate = 0.1
    
    def __init__(self, circumferentialElements,axialElements,wallElements,discret=10,length=2.0,innerRadius=1.5,outerRadius=2.0,fixBottom=True,fixTop=False, stage=0):
        '''
        Constructor
        '''
        self.circumferentialElements = circumferentialElements
        self.axialElements = axialElements
        self.wallElements = wallElements
        self.length=length
        self.innerRadius=innerRadius
        self.outerRadius=outerRadius
        self.fixBottom=fixBottom
        self.fixTop=fixTop
        self.stage = stage
        xi = np.linspace(0,1.0,discret)
        xi1,xi2,xi3 = np.meshgrid(xi,xi,xi)
        self.xi = np.c_[xi1.flatten(),xi2.flatten(),xi3.flatten()].T.tolist()
        self.numPoints = discret**3
                
    def createBoundaryConditions(self,nonlinearEquations):
        dependentField = self.dependentField
        boundaryConditions = iron.BoundaryConditions()
        nonlinearEquations.BoundaryConditionsCreateStart(boundaryConditions)
        numberOfCircumferentialElementsPerQuarter = int(self.circumferentialElements/4)
        numberOfCircumferentialElements = self.circumferentialElements
        numberOfLengthNodes = self.axialElements+1
        numberOfCircumferentialNodes = numberOfCircumferentialElements
        
        
        #nodelistx = [91,99,107,115,123,131,139]
        #nodelisty = [73,77,81,85,139]
        #nodelistz = [65,66,67,68,69,70,71,72,137,138,139,140,141,142,143,144]
        nodelistx = [155,211] 
        nodelisty = [145,149,211]
        nodelistz = [211,215]
        for nodeNumber in nodelistz:        
        # Fix S3 and Z direction
            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)
            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0) 

        for nodeNumber in nodelisty:
        # Fix S2 and Y direction
            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
            #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
            #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0) 

        for nodeNumber in nodelistx:
        # Fix S1 and X direction
            boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
            #boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)        
        nonlinearEquations.BoundaryConditionsCreateFinish()

    def setupGrowthRates(self,growthElementRate):
        numberOfCircumferentialElements = self.circumferentialElements
        numberOfLengthElements = self.axialElements
        numberOfWallElements = self.wallElements
        numberOfGaussXi = 3 
        numberOfElements = numberOfCircumferentialElements*numberOfLengthElements*numberOfWallElements
        growthCellMLParametersField = self.growthCellMLParametersField
        #growthCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,fibreGrowthRate[0])
        #growthCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,sheetGrowthRate[0])
        #growthCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,normalGrowthRate[0])
        for elementNumber in range (1, numberOfElements+1):
            for xiIdx3 in range(1,numberOfGaussXi+1):
                for xiIdx2 in range(1,numberOfGaussXi+1):
                    for xiIdx1 in range(1,numberOfGaussXi+1):
                        gaussPointNumber = xiIdx1 + (xiIdx2-1)*numberOfGaussXi + (xiIdx3-1)*numberOfGaussXi*numberOfGaussXi
                        growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,1,growthElementRate[elementNumber-1, xiIdx1-1, xiIdx2-1,xiIdx3-1,0])
                        growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,2,growthElementRate[elementNumber-1, xiIdx1-1, xiIdx2-1,xiIdx3-1,1])
                        growthCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                        iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,3,growthElementRate[elementNumber-1, xiIdx1-1, xiIdx2-1,xiIdx3-1,2])
        growthCellMLParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        growthCellMLParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)    


    def setupGeometry(self,geometricField,stage):
        stageNumber = self.stage
        exnode = exfile.Exnode("mesh" + str(stageNumber) + "-8x8.part0.exnode")       # considered from 0-24   stageNumber+2 is the next stage ...
        derivativeValues = []
        #adding coordinates of the nodes 
        #adding dervatives of the nodes 
        for node_num in range(1, exnode.num_nodes + 1):
            version = 1
            derivative = 1        
            for component in range(1, 3 + 1):
                component_name = ["x", "y", "z"][component - 1]
                value = exnode.node_value("Coordinate", component_name, node_num, derivative)
                geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, version, derivative, node_num, component, value)        
            derivative = 2
            for component in range(1, 3 + 1):
                component_name = ["x", "y", "z"][component - 1]
                derivativeValues = exnode.node_values("Coordinate", component_name, node_num)
                geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, version, derivative, node_num, component, derivativeValues[1])
            derivative = 3
            for component in range(1, 3 + 1):
                component_name = ["x", "y", "z"][component - 1]
                derivativeValues = exnode.node_values("Coordinate", component_name, node_num)
                geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, version, derivative, node_num, component, derivativeValues[2])        
            derivative = 5
            for component in range(1, 3 + 1):
                component_name = ["x", "y", "z"][component - 1]
                derivativeValues = exnode.node_values("Coordinate", component_name, node_num)
                geometricField.ParameterSetUpdateNode(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, version, derivative, node_num, component, derivativeValues[4])
        geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

		
    def setupProblem(self,showProgress=False):
        # Hydrostatic pressure
        pInit = -8.0 # The initial hydrostatic pressure
        
        # Fibre angle
        #fibreAngle = math.pi/2.0 # The fibre angle wrt the for anisotropic fibres
        fibreAngle = 0.0 # The fibre angle wrt the for anisotropic fibres
                
        # Number of Gauss points used
        numberOfGaussXi = 3
        numberOfCircumferentialElements = self.circumferentialElements
        numberOfLengthElements = self.axialElements
        numberOfLengthNodes = self.axialElements+1
        numberOfCircumferentialNodes = numberOfCircumferentialElements
        numberOfWallNodes = self.wallElements+1
        numberOfWallElements = self.wallElements
        
        coordinateSystemUserNumber = 1
        regionUserNumber = 1
        tricubicHermiteBasisUserNumber = 1
        trilinearLagrangeBasisUserNumber = 2
        meshUserNumber = 1
        decompositionUserNumber = 1
        geometricFieldUserNumber = 1
        originalGeometricFieldUserNumber = 20
        fibreFieldUserNumber = 2
        dependentFieldUserNumber = 3
        equationsSetUserNumber = 1
        equationsSetFieldUserNumber = 5
        growthCellMLUserNumber = 1
        growthCellMLModelsFieldUserNumber = 6
        growthCellMLStateFieldUserNumber = 7
        growthCellMLParametersFieldUserNumber = 8
        constitutiveCellMLUserNumber = 2
        constitutiveCellMLModelsFieldUserNumber = 9
        constitutiveCellMLParametersFieldUserNumber = 10
        constitutiveCellMLIntermediateFieldUserNumber = 11
        problemUserNumber = 1
        
        # Get the number of computational nodes and this computational node number
        numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
        #computationalNodeNumber = iron.ComputationalNodeNumberGet()
        
        # Create a 3D rectangular cartesian coordinate system
        coordinateSystem = iron.CoordinateSystem()
        coordinateSystem.CreateStart(coordinateSystemUserNumber)
        # Set the number of dimensions to 3
        coordinateSystem.DimensionSet(3)
        # Finish the creation of the coordinate system
        coordinateSystem.CreateFinish()
        
        # Create a region and assign the coordinate system to the region
        region = iron.Region()
        region.CreateStart(regionUserNumber,iron.WorldRegion)
        region.LabelSet("HeartTubeRegion")
        # Set the regions coordinate system to the 3D RC coordinate system that we have created
        region.coordinateSystem = coordinateSystem
        # Finish the creation of the region
        region.CreateFinish()
        
        # Define basis
        # Start the creation of a tricubic Hermite basis function
        tricubicHermiteBasis = iron.Basis()
        tricubicHermiteBasis.CreateStart(tricubicHermiteBasisUserNumber)
        tricubicHermiteBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
        tricubicHermiteBasis.numberOfXi = 3
        tricubicHermiteBasis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*3
        tricubicHermiteBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
        tricubicHermiteBasis.CreateFinish()
        # Start the creation of a trilinear Hermite basis function
        trilinearLagrangeBasis = iron.Basis()
        trilinearLagrangeBasis.CreateStart(trilinearLagrangeBasisUserNumber)
        trilinearLagrangeBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
        trilinearLagrangeBasis.numberOfXi = 3
        trilinearLagrangeBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
        trilinearLagrangeBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
        trilinearLagrangeBasis.CreateFinish()
        
        # Start the creation of a manually generated mesh in the region
        numberOfNodes = numberOfCircumferentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
        numberOfElements = numberOfCircumferentialElements*numberOfLengthElements*numberOfWallElements
               
        # Define nodes for the mesh
        nodes = iron.Nodes()
        nodes.CreateStart(region,numberOfNodes)
        nodes.CreateFinish()
        mesh = iron.Mesh()
        
        # Create the mesh. The mesh will have two components - 1. tricubic Hermite elements; 2. trilinear Lagrange elements
        mesh.CreateStart(meshUserNumber,region,3)
        mesh.NumberOfComponentsSet(2)
        mesh.NumberOfElementsSet(numberOfElements)
        
        tricubicHermiteElements = iron.MeshElements()
        tricubicHermiteElements.CreateStart(mesh,1,tricubicHermiteBasis)
        trilinearLagrangeElements = iron.MeshElements()
        trilinearLagrangeElements.CreateStart(mesh,2,trilinearLagrangeBasis)
        
        elementNumber = 0
        for wallElementIdx in range(1,numberOfWallElements+1):
            for lengthElementIdx in range(1,numberOfLengthElements+1):
                for circumferentialElementIdx in range(1,numberOfCircumferentialElements+1):
                    elementNumber = elementNumber + 1
                    localNode1 = circumferentialElementIdx + (lengthElementIdx-1)*numberOfCircumferentialNodes + \
                        (wallElementIdx-1)*numberOfCircumferentialNodes*numberOfLengthNodes
                    if circumferentialElementIdx == numberOfCircumferentialElements:
                        localNode2 = 1 + (lengthElementIdx-1)*numberOfCircumferentialNodes + \
                            (wallElementIdx-1)*numberOfCircumferentialNodes*numberOfLengthNodes
                    else:
                        localNode2 = localNode1 + 1
                    localNode3 = localNode1 + numberOfCircumferentialNodes
                    localNode4 = localNode2 + numberOfCircumferentialNodes
                    localNode5 = localNode1 + numberOfCircumferentialNodes*numberOfLengthNodes
                    localNode6 = localNode2 + numberOfCircumferentialNodes*numberOfLengthNodes
                    localNode7 = localNode3 + numberOfCircumferentialNodes*numberOfLengthNodes
                    localNode8 = localNode4 + numberOfCircumferentialNodes*numberOfLengthNodes
                    localNodes = [localNode1,localNode2,localNode3,localNode4,localNode5,localNode6,localNode7,localNode8]
                    tricubicHermiteElements.NodesSet(elementNumber,localNodes)
                    trilinearLagrangeElements.NodesSet(elementNumber,localNodes)
        
        tricubicHermiteElements.CreateFinish()
        trilinearLagrangeElements.CreateFinish()
        
        # Finish the mesh creation
        mesh.CreateFinish() 
        
        # Create a decomposition for the mesh
        decomposition = iron.Decomposition()
        decomposition.CreateStart(decompositionUserNumber,mesh)
        # Set the decomposition to be a general decomposition with the specified number of domains
        decomposition.type = iron.DecompositionTypes.CALCULATED
        decomposition.numberOfDomains = numberOfComputationalNodes
        # Finish the decomposition
        decomposition.CreateFinish()
        
        # Create a field for the geometry
        geometricField = iron.Field()
        geometricField.CreateStart(geometricFieldUserNumber,region)
        # Set the decomposition to use
        geometricField.MeshDecompositionSet(decomposition)
        geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
        # Set the field label
        geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
        # Set the domain to be used by the field components to be tricubic Hermite
        geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
        geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
        geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
        # Set the scaling type
        geometricField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
        # Finish creating the field
        geometricField.CreateFinish()

        originalGeometricField= iron.Field()
        originalGeometricField.CreateStart(originalGeometricFieldUserNumber,region)
        # Set the decomposition to use
        originalGeometricField.MeshDecompositionSet(decomposition)
        originalGeometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
        # Set the field label
        originalGeometricField.VariableLabelSet(iron.FieldVariableTypes.U,"OriginalGeometry")
        # Set the domain to be used by the field components to be tricubic Hermite
        originalGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
        originalGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
        originalGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
        # Set the scaling type
        originalGeometricField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
        # Finish creating the field
        originalGeometricField.CreateFinish()
        
        self.setupGeometry(geometricField,self.stage)     
        # Update the geometric field
        geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        #Copy the geometric field
        iron.Field.ParametersToFieldParametersComponentCopy(
            geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
            originalGeometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
        iron.Field.ParametersToFieldParametersComponentCopy(
            geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
            originalGeometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
        iron.Field.ParametersToFieldParametersComponentCopy(
            geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
            originalGeometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
        originalGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        originalGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        
        
        # Create a fibre field and attach it to the geometric field
        fibreField = iron.Field()
        fibreField.CreateStart(fibreFieldUserNumber,region)
        fibreField.TypeSet(iron.FieldTypes.FIBRE)
        # Set the decomposition 
        fibreField.MeshDecompositionSet(decomposition)
        # Set the geometric field
        fibreField.GeometricFieldSet(geometricField)
        # Set the field variable label
        fibreField.VariableLabelSet(iron.FieldVariableTypes.U,"Fibre")
        # Set the fibre field to use trilinear-Lagrange elements
        fibreField.NumberOfComponentsSet(iron.FieldVariableTypes.U,3)
        fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,2)
        fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,2)
        fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,2)
        # Finish creating the field
        fibreField.CreateFinish()
        #Initialise the fibre field
        for wallNodeIdx in range(1,numberOfWallNodes+1):
            for lengthNodeIdx in range(1,numberOfLengthNodes+1):
                for circumferentialNodeIdx in range(1,numberOfCircumferentialNodes+1):
                    nodeNumber = circumferentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumferentialNodes + \
                        (wallNodeIdx-1)*numberOfCircumferentialNodes*numberOfLengthNodes
                    # Set the fibre angle
                    fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                            1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,fibreAngle)
                    fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                            1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,0.0)
                    fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                            1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,0.0)
        # Update the fibre field
        fibreField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        fibreField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
                
        # Create the dependent field
        dependentField = iron.Field()
        dependentField.CreateStart(dependentFieldUserNumber,region)
        dependentField.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)
        # Set the decomposition
        dependentField.MeshDecompositionSet(decomposition)
        # Set the geometric field
        dependentField.GeometricFieldSet(geometricField) 
        dependentField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
        # Set the field variables for displacement, traction, strain, stress and growth
        dependentField.NumberOfVariablesSet(5)
        dependentField.VariableTypesSet([iron.FieldVariableTypes.U,iron.FieldVariableTypes.DELUDELN,
                                         iron.FieldVariableTypes.U1,iron.FieldVariableTypes.U2,iron.FieldVariableTypes.U3])
        dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
        dependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"del U/del n")
        dependentField.VariableLabelSet(iron.FieldVariableTypes.U1,"Strain")
        dependentField.VariableLabelSet(iron.FieldVariableTypes.U2,"Stress")
        dependentField.VariableLabelSet(iron.FieldVariableTypes.U3,"Growth")
        dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,4)
        dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,4)
        dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U1,6)
        dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U2,6)
        dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U3,3)
        # Set the hydrostatic pressure to use tri-linear Lagrange elements
        dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,2)
        dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,2)
        # Set the strain, stress and growth variables to be Gauss point based.
        for comp in range(1,7):
            dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,comp,
                                                 iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
            dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,comp,
                                                 iron.FieldInterpolationTypes.GAUSS_POINT_BASED)

        for comp in range(1,4):
            dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U3,comp,
                                                 iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
        # Set the field scaling
        dependentField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
        # Finish creating the field
        dependentField.CreateFinish()
        
        # Initialise dependent field from undeformed geometry
        for comp in range(1,4):
            iron.Field.ParametersToFieldParametersComponentCopy(
                geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,
                dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp)
        # Initialise the hydrostatic pressure
        iron.Field.ComponentValuesInitialiseDP(dependentField,iron.FieldVariableTypes.U,
                                               iron.FieldParameterSetTypes.VALUES,4,pInit)
        
        # Create the equations_set
        equationsSetField = iron.Field()
        equationsSet = iron.EquationsSet()
        # Specify a finite elasticity equations set with the growth and constitutive law in CellML
        equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
            iron.EquationsSetTypes.FINITE_ELASTICITY,
            iron.EquationsSetSubtypes.CONSTIT_AND_GROWTH_LAW_IN_CELLML]
        equationsSet.CreateStart(equationsSetUserNumber,region,fibreField,
                                     equationsSetSpecification,equationsSetFieldUserNumber,
                                     equationsSetField)
        equationsSet.CreateFinish()
        
        # Set up the equation set dependent field
        equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
        equationsSet.DependentCreateFinish()
        
        # Create equations
        equations = iron.Equations()
        equationsSet.EquationsCreateStart(equations)
        # Use sparse equations
        equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
        # Do not output any equations information
        equations.outputType = iron.EquationsOutputTypes.NONE
        # Finish creating the equations
        equationsSet.EquationsCreateFinish()
       
        # Set up the growth CellML model
        growthCellML = iron.CellML()
        growthCellML.CreateStart(growthCellMLUserNumber,region)
        
        # Create the CellML environment for the simple growth law
        growthCellMLIdx = growthCellML.ModelImport("simplegrowth.cellml")
        # Flag the CellML variables that OpenCMISS will supply
        growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/fibrerate")
        growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/sheetrate")
        growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/normalrate")
        # Finish the growth CellML
        growthCellML.CreateFinish()
        
        # Create CellML <--> OpenCMISS field maps
        growthCellML.FieldMapsCreateStart()
        growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda1",iron.FieldParameterSetTypes.VALUES,
            dependentField,iron.FieldVariableTypes.U3,1,iron.FieldParameterSetTypes.VALUES)
        growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda2",iron.FieldParameterSetTypes.VALUES,
            dependentField,iron.FieldVariableTypes.U3,2,iron.FieldParameterSetTypes.VALUES)
        growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda3",iron.FieldParameterSetTypes.VALUES,
                dependentField,iron.FieldVariableTypes.U3,3,iron.FieldParameterSetTypes.VALUES)
        growthCellML.FieldMapsCreateFinish()
        
        # Create the CELL models field
        growthCellMLModelsField = iron.Field()
        growthCellML.ModelsFieldCreateStart(growthCellMLModelsFieldUserNumber,growthCellMLModelsField)
        growthCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthModelMap")
        growthCellML.ModelsFieldCreateFinish()
        
        # Create the CELL parameters field
        growthCellMLParametersField = iron.Field()
        growthCellML.ParametersFieldCreateStart(growthCellMLParametersFieldUserNumber,growthCellMLParametersField)
        growthCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthParameters")
        growthCellML.ParametersFieldCreateFinish()
        
        
        # Create the CELL state field
        growthCellMLStateField = iron.Field()
        growthCellML.StateFieldCreateStart(growthCellMLStateFieldUserNumber,growthCellMLStateField)
        growthCellMLStateField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthState")
        growthCellML.StateFieldCreateFinish()
        
        # Create the CellML environment for the consitutative law
        constitutiveCellML = iron.CellML()
        constitutiveCellML.CreateStart(constitutiveCellMLUserNumber,region)
        constitutiveCellMLIdx = constitutiveCellML.ModelImport("mooneyrivlin.cellml")
        # Flag the CellML variables that OpenCMISS will supply
        constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/C11")
        constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/C12")
        constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/C13")
        constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/C22")
        constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/C23")
        constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/C33")
        if (humphreyModel):
            constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/c1")
            constitutiveCellML.VariableSetAsKnown(constitutiveCellMLIdx,"equations/c2")
        # Flag the CellML variables that OpenCMISS will obtain
        constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev11")
        constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev12")
        constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev13")
        constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev22")
        constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev23")
        constitutiveCellML.VariableSetAsWanted(constitutiveCellMLIdx,"equations/Tdev33")
        constitutiveCellML.CreateFinish()
        
        # Create CellML <--> OpenCMISS field maps
        constitutiveCellML.FieldMapsCreateStart()
        constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,1,iron.FieldParameterSetTypes.VALUES,
            constitutiveCellMLIdx,"equations/C11",iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,2,iron.FieldParameterSetTypes.VALUES,
            constitutiveCellMLIdx,"equations/C12",iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,3,iron.FieldParameterSetTypes.VALUES,
            constitutiveCellMLIdx,"equations/C13",iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,4,iron.FieldParameterSetTypes.VALUES,
            constitutiveCellMLIdx,"equations/C22",iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,5,iron.FieldParameterSetTypes.VALUES,
            constitutiveCellMLIdx,"equations/C23",iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,6,iron.FieldParameterSetTypes.VALUES,
            constitutiveCellMLIdx,"equations/C33",iron.FieldParameterSetTypes.VALUES)
        
        constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev11",iron.FieldParameterSetTypes.VALUES,
            dependentField,iron.FieldVariableTypes.U2,1,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev12",iron.FieldParameterSetTypes.VALUES,
            dependentField,iron.FieldVariableTypes.U2,2,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev13",iron.FieldParameterSetTypes.VALUES,
            dependentField,iron.FieldVariableTypes.U2,3,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev22",iron.FieldParameterSetTypes.VALUES,
            dependentField,iron.FieldVariableTypes.U2,4,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev23",iron.FieldParameterSetTypes.VALUES,
            dependentField,iron.FieldVariableTypes.U2,5,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.CreateCellMLToFieldMap(constitutiveCellMLIdx,"equations/Tdev33",iron.FieldParameterSetTypes.VALUES,
            dependentField,iron.FieldVariableTypes.U2,6,iron.FieldParameterSetTypes.VALUES)
        constitutiveCellML.FieldMapsCreateFinish()
        
        # Create the CELL models field
        constitutiveCellMLModelsField = iron.Field()
        constitutiveCellML.ModelsFieldCreateStart(constitutiveCellMLModelsFieldUserNumber,
                                                   constitutiveCellMLModelsField)
        constitutiveCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstitutiveModelMap")
        constitutiveCellML.ModelsFieldCreateFinish()
        
        # Create the CELL parameters field
        constitutiveCellMLParametersField = iron.Field()
        constitutiveCellML.ParametersFieldCreateStart(constitutiveCellMLParametersFieldUserNumber,
                                                       constitutiveCellMLParametersField)
        constitutiveCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstitutiveParameters")
        constitutiveCellML.ParametersFieldCreateFinish()
        if (humphreyModel):
            # Set up the materials constants
            c1ComponentNumber = constitutiveCellML.FieldComponentGet(constitutiveCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"equations/c1")
            c2ComponentNumber = constitutiveCellML.FieldComponentGet(constitutiveCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"equations/c2")

            #c3ComponentNumber = constitutiveCellML.FieldComponentGet(constitutiveCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"equations/c3")
            #print "could succefully get the rates ... "

            for gaussPointNumber in range (1,27+1):
                for elementNumber in range (1,64+1):
                    constitutiveCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,c1ComponentNumber,self.c1)
                    constitutiveCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,c2ComponentNumber,self.c2)

            #constitutiveCellMLParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
            #constitutiveCellMLParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

            for gaussPointNumber in range (1,27+1):
                for elementNumber in range (64+1,128+1):
                    constitutiveCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,c1ComponentNumber,self.c3)
                    constitutiveCellMLParametersField.ParameterSetUpdateGaussPointDP(iron.FieldVariableTypes.U,
                                iron.FieldParameterSetTypes.VALUES,gaussPointNumber,elementNumber,c2ComponentNumber,self.c4)
            constitutiveCellMLParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
            constitutiveCellMLParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

        # Create the CELL intermediate field
        constitutiveCellMLIntermediateField = iron.Field()
        constitutiveCellML.IntermediateFieldCreateStart(constitutiveCellMLIntermediateFieldUserNumber,
                                                         constitutiveCellMLIntermediateField)
        constitutiveCellMLIntermediateField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstitutiveIntermediate")
        constitutiveCellML.IntermediateFieldCreateFinish()
        
        # Define the problem
        problem = iron.Problem()
        problemSpecification = [iron.ProblemClasses.ELASTICITY,
                iron.ProblemTypes.FINITE_ELASTICITY,
                iron.ProblemSubtypes.FINITE_ELASTICITY_WITH_GROWTH_CELLML]
        problem.CreateStart(problemUserNumber,problemSpecification)
        problem.CreateFinish()
        
        # Create control loops
        timeLoop = iron.ControlLoop()
        problem.ControlLoopCreateStart()
        problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],timeLoop)
        problem.ControlLoopCreateFinish()
        
        # Create problem solvers
        odeIntegrationSolver = iron.Solver()
        nonlinearSolver = iron.Solver()
        linearSolver = iron.Solver()
        cellMLEvaluationSolver = iron.Solver()
        problem.SolversCreateStart()
        problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,odeIntegrationSolver)
        problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,nonlinearSolver)
        nonlinearSolver.outputType = iron.SolverOutputTypes.NONE
        if showProgress:
            nonlinearSolver.outputType = iron.SolverOutputTypes.PROGRESS
        nonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
        nonlinearSolver.NewtonCellMLSolverGet(cellMLEvaluationSolver)
        nonlinearSolver.NewtonLinearSolverGet(linearSolver)
        linearSolver.linearType = iron.LinearSolverTypes.DIRECT
        problem.SolversCreateFinish()
        
        # Create nonlinear equations and add equations set to solver equations
        nonlinearEquations = iron.SolverEquations()
        problem.SolverEquationsCreateStart()
        nonlinearSolver.SolverEquationsGet(nonlinearEquations)
        nonlinearEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
        nonlinearEquationsSetIndex = nonlinearEquations.EquationsSetAdd(equationsSet)
        problem.SolverEquationsCreateFinish()
        
        # Create CellML equations and add growth and constitutive equations to the solvers
        growthEquations = iron.CellMLEquations()
        constitutiveEquations = iron.CellMLEquations()
        problem.CellMLEquationsCreateStart()
        odeIntegrationSolver.CellMLEquationsGet(growthEquations)
        growthEquationsIndex = growthEquations.CellMLAdd(growthCellML)
        cellMLEvaluationSolver.CellMLEquationsGet(constitutiveEquations)
        constitutiveEquationsIndex = constitutiveEquations.CellMLAdd(constitutiveCellML)
        problem.CellMLEquationsCreateFinish()
        
        # Prescribe boundary conditions (absolute nodal parameters)
        self.dependentField = dependentField
        self.elementNumber = elementNumber
        self.growthCellMLParametersField = growthCellMLParametersField
        self.constitutiveCellMLIntermediateField = constitutiveCellMLIntermediateField
        self.constitutiveCellMLParametersField = constitutiveCellMLParametersField
        self.geometricField = geometricField
        self.originalGeometricField = originalGeometricField
        self.timeLoop = timeLoop
        self.problem = problem
        self.growthCellMLStateField = growthCellMLStateField
        self.createBoundaryConditions(nonlinearEquations)
        self.region = region
        self.equationsSet = equationsSet
        
    def resetFieldsForNewSolution(self):
        #Reset to initial geometry 
        for comp in range(1,4):
            iron.Field.ParametersToFieldParametersComponentCopy(
                self.originalGeometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,
                self.geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp)
            iron.Field.ParametersToFieldParametersComponentCopy(
                self.originalGeometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,
                self.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp)
            
        self.geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        self.geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)        
        #Reset the dependent field
        for comp in range(1,7):
            iron.Field.ComponentValuesInitialiseDP(self.dependentField,iron.FieldVariableTypes.U2,
                                               iron.FieldParameterSetTypes.VALUES,comp,0.0)
        for comp in range(1,4):
            iron.Field.ComponentValuesInitialiseDP(self.dependentField,iron.FieldVariableTypes.U3,
                                               iron.FieldParameterSetTypes.VALUES,comp,0.0)
            iron.Field.ComponentValuesInitialiseDP(self.dependentField,iron.FieldVariableTypes.DELUDELN,
                                               iron.FieldParameterSetTypes.VALUES,comp,0.0)

        iron.Field.ComponentValuesInitialiseDP(self.dependentField,iron.FieldVariableTypes.DELUDELN,
                                               iron.FieldParameterSetTypes.VALUES,4,0.0)

        
        #Reset constitutive intermediate and parameters field
        '''
        for comp in range(1,7):
            self.constitutiveCellMLIntermediateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                           iron.FieldParameterSetTypes.VALUES,comp,0.0)
            self.constitutiveCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                           iron.FieldParameterSetTypes.VALUES,2,1.0)
        self.constitutiveCellMLIntermediateField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        self.constitutiveCellMLIntermediateField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        self.constitutiveCellMLParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        self.constitutiveCellMLParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        '''             
        #Reset growth field
        for comp in range(1,4):
            self.growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                           iron.FieldParameterSetTypes.VALUES,comp,1.0)
        self.growthCellMLStateField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        self.growthCellMLStateField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        
    def growthRatesInterpolation(self,growthRateVariables):
        layers = 2  #  can be 1 or 2              1 for only cardiac muscle  and  2  for cardiac muscle and cardiac jelly
        division = 4     # can be 2 or 4              2 for coarse meshgrid 18 vars and 4 for refine meshgrid with 36 vars 
        numberOfCircumferentialPoints = division + 1
        numberOfLongitudinalPoints = division + 1 
        numberOfCircumferentialElements = self.circumferentialElements
        numberOfLengthElements = self.axialElements
        numberOfLengthNodes = self.axialElements+1
        numberOfCircumferentialNodes = numberOfCircumferentialElements
        numberOfWallNodes = self.wallElements+1
        numberOfWallElements = self.wallElements
        numberOfNodes = numberOfCircumferentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
        numberOfElements = numberOfCircumferentialElements*numberOfLengthElements*numberOfWallElements
        grid_x, grid_y = np.mgrid[0:1:1000j, 0:1:1000j]
        # distribute the rates to the points 

        values = np.zeros(((division+1)*(division+1),3,layers)) 
        for k in range (layers):
            for m in range (division+1):
                for n in range (division):
                    values[m*n+n,:,k] =  growthRateVariables [(division+1)*division*k+m*n+n,:]		
        # there needs to be a symmetry as of unrolling cylinder 	
        #print (values)
        for j in range (division+1):
            values[(division+1)*j+division,:,:] = values[(division+1)*j,:,:] 
        values1 = np.zeros(((division+1)*(division+1),3))
        values1 [:,:] = values [:,:,0]
        if (layers == 2 ):
            values2 = np.zeros(((division+1)*(division+1),3))
            values2 [:,:] = values [:,:,1]
        #print (values1)
        #print (values2)

        # number of points
        numberOfPoints = (division + 1)*(division + 1)
        points = np.zeros((numberOfPoints,2))
        if (division == 4):
            points[:,0] = 0,0,0,0,0,0.25,0.25,0.25,0.25,0.25,0.5,0.5,0.5,0.5,0.5,0.75,0.75,0.75,0.75,0.75,1,1,1,1,1
            points[:,1] = 0,0.25,0.5,0.75,1,0,0.25,0.5,0.75,1,0,0.25,0.5,0.75,1,0,0.25,0.5,0.75,1,0,0.25,0.5,0.75,1
        elif (division == 2):
            points[:,0] = 0,0,0,0.5,0.5,0.5,1,1,1
            points[:,1] = 0,0.5,1,0,0.5,1,0,0.5,1

        grid_z00 = griddata(points, values1[:,0], (grid_x, grid_y), method='cubic')
        grid_z10 = griddata(points, values1[:,1], (grid_x, grid_y), method='cubic')
        grid_z20 = griddata(points, values1[:,2], (grid_x, grid_y), method='cubic')
        #grid_z30 = griddata(points, values1[:,3], (grid_x, grid_y), method='cubic')
        if (layers == 2):
            grid_z01 = griddata(points, values2[:,0], (grid_x, grid_y), method='cubic')
            grid_z11 = griddata(points, values2[:,1], (grid_x, grid_y), method='cubic')
            grid_z21 = griddata(points, values2[:,2], (grid_x, grid_y), method='cubic')
            #grid_z31 = griddata(points, values2[:,3], (grid_x, grid_y), method='cubic')
        # the location of gauss points  in [-1,1]  -0.577, 0 , 0.577
        # the location of gauss points  in [0,1]  (1-0.577)/2, 0.5 , (1+0.577)/2   ===>  0.211, 0.5, 0.789  or 0.21 , 0.5 , 0.79

        elementNumber = 0
        elementOrigin = []
        gaussPointsLocations = [0.21 , 0.5 , 0.79]
        interpolatedRates = np.zeros((numberOfElements,3,3,3,3))
        for numLayers in range (layers):
            for lengthElems in range (8):
                for circumElems in range (8):
                    elementNumber = numLayers*64 + lengthElems*8 + circumElems
                    for xi1 in range (3):
                        for xi2 in range (3):
                            for xi3 in range (3):
                                s = circumElems/8 + gaussPointsLocations[xi1]/8
                                t = lengthElems/8 + gaussPointsLocations[xi2]/8
                                p = gaussPointsLocations[xi3]
                                if (layers == 1 ):
                                    interpolatedRates[elementNumber, xi1, xi2,xi3,0] = grid_z00[int(1000*t) , int(1000*s)]
                                    interpolatedRates[elementNumber, xi1, xi2,xi3,1] = grid_z10[int(1000*t) , int(1000*s)] 
                                    interpolatedRates[elementNumber, xi1, xi2,xi3,2] = grid_z20[int(1000*t) , int(1000*s)]
                                elif (layers == 2 ):
                                    interpolatedRates[elementNumber, xi1, xi2,xi3,0] = ((1-p)*(grid_z00[int(1000*t),int(1000*s)])  +  p*(grid_z00[int(1000*t),
                                                int(1000*s)] + grid_z01[int(1000*t) , int(1000*s)])/2)*math.fabs(layers-1)  +  (p*(grid_z01[int(1000*t),
                                                int(1000*s)])  +  (1-p)*(grid_z00[int(1000*t),int(1000*s)] + grid_z01[int(1000*t) , int(1000*s)])/2)*math.fabs(layers)
                                    interpolatedRates[elementNumber, xi1, xi2,xi3,1] = ((1-p)*(grid_z10[int(1000*t),int(1000*s)])  +  p*(grid_z10[int(1000*t),
                                                int(1000*s)] + grid_z11[int(1000*t) , int(1000*s)])/2)*math.fabs(layers-1)  +  (p*(grid_z11[int(1000*t),
                                                int(1000*s)])  +  (1-p)*(grid_z10[int(1000*t),int(1000*s)] + grid_z11[int(1000*t) , int(1000*s)])/2)*math.fabs(layers)
                                    interpolatedRates[elementNumber, xi1, xi2,xi3,2] = ((1-p)*(grid_z20[int(1000*t),int(1000*s)])  +  p*(grid_z20[int(1000*t),
                                                int(1000*s)] + grid_z21[int(1000*t) , int(1000*s)])/2)*math.fabs(layers-1)  +  (p*(grid_z21[int(1000*t),
                                                int(1000*s)])  +  (1-p)*(grid_z20[int(1000*t),int(1000*s)] + grid_z21[int(1000*t) , int(1000*s)])/2)*math.fabs(layers)   


        # numberOfPoints = int(numberOfLengthNodes*numberOfCircumferentialNodes*numberOfWallElements)                                 
        # #numberOfPoints = 144      # 9*8*2
        # numberOfVariables = int((numberOfLengthElements/2 + 1)*(numberOfCircumferentialElements/2)*(numberOfWallElements)) 
        # #numberOfVariables = 40      # 5*4*2
        # # the purpose of this map from 40x3 points to 144 elements 
        # growthRatePoints = np.zeros((numberOfPoints,3))                                                                	# (8*9*2)*(3)
        # for i in range(1,numberOfWallElements+1):                                                                        # (1,2)
            # for j in range(1,int(numberOfLengthElements/2)+2):                                                            # (1,5)
                # for k in range(1, int(numberOfCircumferentialElements/2)+1):     	                            			# (1,9)
                    # m = (2*k-1)+(2*(j-1)*8)+(i-1)*8*9 - 1
                    # n = k+(j-1)*4+(i-1)*4*5 - 1
                    # growthRatePoints [m ,:] = growthRateVariables[n,:]
        # # ==================================================
        # # interpolate in the empty spots in 8x9x2x3         
        # # part 1 ... interpolate in half-filled lines ...                            
        # for j in range (int(numberOfPoints/numberOfCircumferentialElements)):                                                     # (0:9*2)
            # for i in range (int(numberOfCircumferentialElements/2)):                                                                   # (0:8)
                # if not (i == 3):
                    # growthRatePoints [j*8+2*i+1,:] = (growthRatePoints [j*8+2*i,:] + growthRatePoints [j*8+2*i+2,:])/2
                # else:
                    # growthRatePoints [j*8+2*i+1,:] = (growthRatePoints [j*8+2*i,:] + growthRatePoints [j*8,:])/2
        # # part 2 .... interpolate in empty lines 
        # for i in range (numberOfWallElements):
            # for k in range (int (numberOfLengthNodes/2)):
                # for j in range (numberOfCircumferentialNodes):
                    # pointNumber =  (2*k+1)*8 +  j  + i*72
                    # growthRatePoints[pointNumber,:] = (growthRatePoints[pointNumber-8,:] + growthRatePoints[pointNumber+8,:])/2
        # # ==================================================
        # growthElementRate = np.zeros((numberOfElements,3))                                                                                        # (8*9*2)*(3)
        # # map from points to elements 
        # elementNumber = 0
        # for wallElementIdx in range(1,numberOfWallElements+1):                                                                             # (1:2)
            # for lengthElementIdx in range(1,numberOfLengthElements+1):
                # for circumferentialElementIdx in range(1,numberOfCircumferentialElements+1):
                    # elementNumber = elementNumber + 1
                    # localNode1 = circumferentialElementIdx + (lengthElementIdx-1)*numberOfCircumferentialNodes + (wallElementIdx-1)*numberOfCircumferentialNodes*numberOfLengthNodes
                    # if circumferentialElementIdx == numberOfCircumferentialElements:
                        # localNode2 = 1 + (lengthElementIdx-1)*numberOfCircumferentialNodes + (wallElementIdx-1)*numberOfCircumferentialNodes*numberOfLengthNodes
                    # else:
                        # localNode2 = localNode1 + 1
                    # localNode3 = localNode1 + numberOfCircumferentialNodes
                    # localNode4 = localNode2 + numberOfCircumferentialNodes
                    # #print (elementNumber, localNode1,localNode2,localNode3,localNode4)
                    # growthElementRate[elementNumber-1,0] = (growthRatePoints[localNode1-1,0] + growthRatePoints[localNode2-1,0] + growthRatePoints[localNode3-1,0] + growthRatePoints[localNode4-1,0])/4
                    # growthElementRate[elementNumber-1,1] = (growthRatePoints[localNode1-1,1] + growthRatePoints[localNode2-1,1] + growthRatePoints[localNode3-1,1] + growthRatePoints[localNode4-1,1])/4
                    # growthElementRate[elementNumber-1,2] = (growthRatePoints[localNode1-1,2] + growthRatePoints[localNode2-1,2] + growthRatePoints[localNode3-1,2] + growthRatePoints[localNode4-1,2])/4
        # return growthElementRate
        return interpolatedRates
        
    def solveAndGetSurfaceDescriptors(self,growthRate):
        numberOfCircumferentialElements = self.circumferentialElements
        numberOfLengthElements = self.axialElements
        numberOfLengthNodes = self.axialElements+1
        numberOfCircumferentialNodes = numberOfCircumferentialElements
        numberOfWallNodes = self.wallElements+1
        numberOfWallElements = self.wallElements
        numberOfNodes = numberOfCircumferentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
        numberOfElements = numberOfCircumferentialElements*numberOfLengthElements*numberOfWallElements

        growthElementRate = np.zeros((numberOfElements,3))    
        growthElementRateDivided = np.zeros((numberOfElements,3))    
        growthElementRate = self.growthRatesInterpolation(growthRate)
        maxRate = growthElementRate.max()
        exponent = maxRate/self.maxSolvebleRate
        if exponent>1.0:
            divi = np.power(10,int(np.log10(exponent))+1)
            nsteps = divi
        else:
            nsteps = 1
            divi = 1.0
        #print("Using ",nsteps," timesteps(s) maxgrowth rate is ",self.maxSolvebleRate," new maxrate is ",maxRate/divi," rates ",growthRate," new rates ",growthRate/divi)
        growthElementRateDivided = growthElementRate/divi        
        # Initialise the parameters field
        self.setupGrowthRates(growthElementRateDivided)
        self.resetFieldsForNewSolution()
        time = 0.0
        timeIncrement = 1.0/nsteps
        for _ in range(nsteps):
            self.timeLoop.TimesSet(time,time+timeIncrement,timeIncrement)
            # Solve the problem
            self.problem.Solve()
            # Set geometric field to current deformed geometry
            for comp in range(1,4):
                iron.Field.ParametersToFieldParametersComponentCopy(
                    self.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp,
                    self.geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,comp)
            # Reset growth state to 1.0
                self.growthCellMLStateField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                               iron.FieldParameterSetTypes.VALUES,comp,1.0)
            self.geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
            self.geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
                
            self.growthCellMLStateField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
            self.growthCellMLStateField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

            
        return self.getSampleCoordinates()

    def getVolume(self):
        vol = 0.0
        for eno in range(1,self.elementNumber+1):
            vol += self.geometricField.GeometricParametersElementVolumeGet(eno)
        return vol
        
    def getSampleCoordinates(self):        
        #Get coordinate values at discretized points
        coords = None
        for eno in range(1,self.elementNumber+1):
            gcoord = self.geometricField.ParameterSetInterpolateMultipleXiDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,\
                                                                      1,eno,self.xi,(3,self.numPoints))
            if not coords is None:
                coords = np.r_[coords,np.array(gcoord).T]
            else:
                coords = np.array(gcoord).T
        rcoords = (coords*1e5).astype('int')
        #print("Volume ",self.getVolume())
        return rcoords
    
    def saveResults(self,filename):
        fields = iron.Fields()
        fields.CreateRegion(self.region)
        fields.NodesExport(filename,"FORTRAN")
        fields.ElementsExport(filename,"FORTRAN")            
        fields.Finalise()
