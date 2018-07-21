'''
Created on 09/06/2018

@author: rjag008
'''
from __future__ import print_function
import numpy as np
from opencmiss.iron import iron
import json,os,sys
import traceback
from io import StringIO
from subprocess import Popen,PIPE
from sqlitedict import SqliteDict  
size = 1
myrank = 0
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    size = comm.Get_size()
except:
    raise ImportError('mpi4py is required for parallelization')
print ('out here 1')
my_env = os.environ.copy()
'''
Even when launching iron through a subprocess,
MPI may find out that it was a part of a distributed process through environment and filesystem based flags
Fortunately, in the case of mpich, it has a few environmental flagsm which if removed the subprocess is loaded
with its own communicator!!

https://stackoverflow.com/questions/21090085/calling-mpi-binary-in-serial-as-subprocess-of-mpi-application

Here are the flags 
'''
try:
    del my_env['PMI_RANK']
except:
    pass
try:    
    del my_env['PMI_SIZE']
except:
    pass
try:    
    del my_env['PMI_FD']
except:
    pass
print ('out here 2')
from Growth import TubeGrower
growthParameters = np.zeros((128,3))                                                                                        # (8*9*2)*(3)
circumferentialElements=8
axialElements=8
wallElements=2
discret=10
length=2.0
innerRadius=1.5
outerRadius=2.0
fixBottom=True
fixTop = False
stage = 2
growthLogic = TubeGrower(circumferentialElements,axialElements,wallElements,discret,length,innerRadius,outerRadius,fixBottom,fixTop,stage)
growthLogic.setupProblem()
targetCoordinates = growthLogic.solveAndGetSurfaceDescriptors(growthParameters)
#result = {'nodePositions':targetCoordinates.tolist()}
#finaltargetCoordinates = np.array(result[u'nodePositions'])
print (targetCoordinates)
print ('out here 3')

def findDeformedCoordinates(rates):
    try:
        p = Popen(['python','parallelOptimization.py'],stdin=PIPE,stderr=PIPE,env=my_env)
        inputData = dict()
        inputData['optmin'] = rates.tolist() 
        subin = json.dumps(inputData)
        #print('##Initializing problem for parameters %s ' % (' ,'.join(map(str,rates.tolist()))))
        _,perr = p.communicate(str.encode(subin))
        p.wait()
        perr = perr.decode()
        #Ensure that we are starting at a json prefix
        ps = perr.index('POSResSTART') + 11
        pe = perr.index('POSResEND')
        result = json.loads(perr[ps:pe])
        if u'error' in result:
            raise ValueError(result[u'error'])
        return np.array(result[u'nodePositions'])
    except:
        traceback.print_exc(file=sys.stdout)

print ('out here 4')


dictionaryLoc = 'solutionsReal.sql'
solutionsDictionary = SqliteDict(dictionaryLoc, autocommit=True)

def getSolutionsDictionary():
    return solutionsDictionary

import random
import pygmo as pg
print ('out here 5')

class GrowthOptimization(object):
    print ('out here 6')

    precision = 1e4
    numberOfLengthElements = 8 
    numberOfCircumferentialElements = 8 
    numberOfWallElements = 2
    numberOfLengthNodes = numberOfLengthElements + 1
    numberOfCircumferentialNodes = numberOfCircumferentialElements
    numberOfWallNodes = numberOfWallElements + 1
    numberOfElements = numberOfLengthElements*numberOfCircumferentialElements*numberOfWallElements
    numberOfNodes = numberOfLengthNodes*numberOfCircumferentialNodes*numberOfWallNodes
    division = 2                                                                                                                       # division = 4    
    numberOfPoints = int((numberOfLengthNodes)*(numberOfCircumferentialNodes)*(numberOfWallElements))                                  # 9*8*2
    numberOfPoints = 144  
    numberOfVariables = int((numberOfLengthElements/division + 1)*(numberOfCircumferentialElements/division)*(numberOfWallElements))   # 5*4*2
    print ('out here 7')

    def __init__(self):
        growthShape = np.zeros((40,3))
        self.grshape = growthShape.shape
        #print ('self.grshape = ', self.grshape)
        #self.targetCoordinates = findDeformedCoordinates(growthRate)
    print ('out here 8')

    def growthRates(self,growthRateVariables):
        # generate some random numbers as the objective variables in a 4x5x2x3 net  
        #growthRateVariables = np.array((self.numberOfVariables,3))                                                                         # (8*9*2)*(3)
        # ==================================================
        # map from 4x5x2x3 points to 8x9x2x3 points 
        #growthRateVariables = np.reshape(random_floats,(numberOfVariables,3))   
        growthRatePoints = np.zeros((self.numberOfPoints,3))                                                                         		# (8*9*2)*(3)
        #print ('self.numberOfPoints', self.numberOfPoints)
        #print ('shape 1', growthRatePoints.shape )
        #print ('shape 2', growthRateVariables.shape )
        for i in range(1,self.numberOfWallElements+1):                                                                                      # (1,2)
            for j in range(1,int(self.numberOfLengthElements/self.division)+2):                                                             # (1,10)
                for k in range(1, int(self.numberOfCircumferentialElements/self.division)+1):     	                            			# (1,9)
                    m = (2*k-1)+(2*(j-1)*8)+(i-1)*8*9 - 1
                    n = k+(j-1)*4+(i-1)*4*5 - 1
                    growthRatePoints [m ,:] = growthRateVariables[n,:]

        # ==================================================
        # interpolate in the empty spots in 8x9x2x3         
        # part 1 ... interpolate in half-filled lines ...                            
        for j in range (int(self.numberOfPoints/self.numberOfCircumferentialElements)):                                                     # (0:9*2)
            for i in range (int(self.numberOfCircumferentialElements/2)):                                                                   # (0:8)
                if not (i == 3):
                    growthRatePoints [j*8+2*i+1,:] = (growthRatePoints [j*8+2*i,:] + growthRatePoints [j*8+2*i+2,:])/2
                else:
                    growthRatePoints [j*8+2*i+1,:] = (growthRatePoints [j*8+2*i,:] + growthRatePoints [j*8,:])/2

        # part 2 .... interpolate in empty lines 
        for i in range (self.numberOfWallElements):
            for k in range (int (self.numberOfLengthNodes/2)):
                for j in range (self.numberOfCircumferentialNodes):
                    pointNumber =  (2*k+1)*8 +  j  + i*72
                    growthRatePoints[pointNumber,:] = (growthRatePoints[pointNumber-8,:] + growthRatePoints[pointNumber+8,:])/2

        # ==================================================
        growthRate = np.zeros((self.numberOfElements,3))                                                                                        # (8*9*2)*(3)
        # map from points to elements 
        elementNumber = 0
        for wallElementIdx in range(1,self.numberOfWallElements+1):                                                                             # (1:2)
            for lengthElementIdx in range(1,self.numberOfLengthElements+1):
                for circumferentialElementIdx in range(1,self.numberOfCircumferentialElements+1):
                    elementNumber = elementNumber + 1
                    localNode1 = circumferentialElementIdx + (lengthElementIdx-1)*self.numberOfCircumferentialNodes + (wallElementIdx-1)*self.numberOfCircumferentialNodes*self.numberOfLengthNodes
                    if circumferentialElementIdx == self.numberOfCircumferentialElements:
                        localNode2 = 1 + (lengthElementIdx-1)*self.numberOfCircumferentialNodes + (wallElementIdx-1)*self.numberOfCircumferentialNodes*self.numberOfLengthNodes
                    else:
                        localNode2 = localNode1 + 1
                    localNode3 = localNode1 + self.numberOfCircumferentialNodes
                    localNode4 = localNode2 + self.numberOfCircumferentialNodes
                    print (elementNumber, localNode1,localNode2,localNode3,localNode4)
                    growthRate[elementNumber-1,0] = (growthRatePoints[localNode1-1,0] + growthRatePoints[localNode2-1,0] + growthRatePoints[localNode3-1,0] + growthRatePoints[localNode4-1,0])/4
                    growthRate[elementNumber-1,1] = (growthRatePoints[localNode1-1,1] + growthRatePoints[localNode2-1,1] + growthRatePoints[localNode3-1,1] + growthRatePoints[localNode4-1,1])/4
                    growthRate[elementNumber-1,2] = (growthRatePoints[localNode1-1,2] + growthRatePoints[localNode2-1,2] + growthRatePoints[localNode3-1,2] + growthRatePoints[localNode4-1,2])/4
        return growthRate
    print ('out here 9')

    def checkSolution(self,key):
        '''
        Checks for solutions within precision decimal places
        '''
        solutions = getSolutionsDictionary()
        mkey = (np.array(key)*self.precision).astype('int')
        with StringIO() as sfile:
            print(mkey, file=sfile)
            skey = sfile.getvalue()
            if skey in solutions:
                print("Found solution for ",key," ",solutions[skey])
                return solutions[skey]
        return None
    
    def addSolution(self,key,value):
        solutions = getSolutionsDictionary()
        mkey = (np.array(key)*self.precision).astype('int')
        with StringIO() as sfile:
            print(mkey, file=sfile)
            print(value, file=sfile)
            skey = sfile.getvalue()
            solutions[skey]=value
            
    def printSolutions(self):
        solutions = getSolutionsDictionary()
        for k,v in solutions.items():
            print(k," ",v)
    print ('out here 10')
 
    def fitness(self,x):
        f = self.checkSolution(x)
        
        if f is None:
            #exnode = exfile.Exnode("mesh" + str(stageNumber) + "-8x8.part0.exnode")       # considered from 0-24   stageNumber+2 is the next stage ... 
            print ('x shape',x.shape)
            xRates = x.reshape(self.grshape)			 # a = np.arange(6).reshape((3,2))    >>> a  array([[0, 1],[2, 3],[4, 5]])
            print ('xRate shape',xRates.shape)
            xElements = self.growthRates(xRates)       
            coordinates = findDeformedCoordinates(xElements)
            f = np.sum(np.linalg.norm(coordinates-targetCoordinates,axis=1))
            self.addSolution(x, f)
        print('Rate ',x,' objective ',f)
        return [f]
    print ('out here 11')

    def get_bounds(self):
        return ([-0.1]*(self.numberOfVariables*3),[0.15]*(self.numberOfVariables*3))
    
runParallel = False
import random
class MonteCarlo(pg.algorithm):
    """
     Monte-Carlo (random sampling) algorithm implemented purely in Python.
    """

    def __init__(self,iterations = 10):
        """
            Constructs a Monte-Carlo (random sampling) algorithm

             USAGE: algorithm.my_algorithm(iter = 10)

             NOTE: At the end of each iteration, the randomly generated
                     point substitutes the worst individual in the population if better

             * iter: number of random samples
        """
        #We start calling the base constructor
        super(MonteCarlo,self).__init__()
        #We then define the algorithm 'private' data members
        self.__iter = iterations

    #This is the 'juice' of the algorithm, the method where the actual optimzation is coded.
    def evolve(self,pop):
        #If the population is empty (i.e. no individuals) nothing happens
        if len(pop) == 0:
            return pop

        #Here we rename some variables, in particular the problem
        prob = pop.problem
        #Its dimensions (total and continuous)
        dim, cont_dim = prob.get_nx(), prob.get_ncx()
        #And the lower/upper bounds for the chromosome
        lb, ub = prob.get_bounds()
        
        #The algorithm now starts manipulating the population
        for _ in range(self.__iter):
            #We create a random vector within the bounds ... first the continuous part
            tmp_cont = [random.uniform(lb[i],ub[i]) for i in range(cont_dim)]
            #then the integer part
            tmp_int = [float(random.randint(lb[i],ub[i])) for i in range(cont_dim,dim)]
            #and we assemble them into one decision vector
            tmp_x = tmp_cont + tmp_int
            #which we push back in the population
            pop.push_back(tmp_x)
            #to then remove the worst individual
            #pop.erase(pop.get_worst_idx())
            #at the end of it all we return the 'evolved' population
            return pop

    def get_name(self):
        return "Monte Carlo (Python)"

    def human_readable_extra(self):
        return "n_iter=" + str(self.__n_iter)
print ('out here 12')

if __name__ == '__main__':
    gp = GrowthOptimization()
    prob = pg.problem(gp)
    algo = pg.algorithm(pg.pso(gen = 200))
    #algo = MonteCarlo(100)
    print ('out here 13')
    try:
        if not runParallel:
            pop = pg.population(prob,20)
            print(dir(pop))
            pop = algo.evolve(pop)
            print(pop.champion_f) 
        else:
            archi = pg.archipelago(n = 16, algo = algo, prob = prob, pop_size = 20, seed = 32)
            archi.evolve()
            print(archi)
            archi.wait()
            res = archi.get_champions_f()
            print(res) 
    except:
        traceback.print_exc(file=sys.stdout)
    