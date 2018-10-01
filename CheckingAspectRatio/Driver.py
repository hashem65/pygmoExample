from __future__ import print_function
import numpy as np
import json,os,sys
import traceback
import datetime 

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
my_env = os.environ.copy()
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

dictionaryLoc = 'aspectRatio.sql'
solutionsDictionary = SqliteDict(dictionaryLoc, autocommit=True)

def getSolutionsDictionary():
    return solutionsDictionary

import random
import pygmo as pg

class GrowthOptimization(object):
    precision = 1e3
    bestAnswer = 1e16
    def __init__(self):
        growthRate = np.zeros((8,6))
        self.grshape = growthRate.shape
        growthRate[0,:3] = 0.05,0.08,0.14
        growthRate[1,:3] = 0.08,0.06,0.08
        growthRate[2,:3] = 0.13,-0.14,0.07
        growthRate[3,:3] = -0.04,-0.13,0.07

        growthRate[4,:3] = 0.18,0.23,0.06
        growthRate[5,:3] = 0.19,0.14,0.05
        growthRate[6,:3] = -0.13,-0.03,0.04
        growthRate[7,:3] = -0.02,-0.12,0.06	

        growthRate[0:4,3:]=0.0
        growthRate[4:8,3:]=0.0
        self.targetCoordinates = findDeformedCoordinates(growthRate)
    
    def checkSolution(self,key):
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
            skey = sfile.getvalue()
            solutions[skey]=value
            
    def printSolutions(self):
        solutions = getSolutionsDictionary()
        for k,v in solutions.items():
            print(k," ",v)
                       
    def fitness(self,x):
        f = self.checkSolution(x)
        if f is None:
            try: 
                growthFull = np.zeros((8,6))
                #growthRatesMain = np.zeros((8,3))
                growthRatesCross = np.zeros((8,3))
                growthRatesMain = x.reshape(8,3)
                growthFull[:,0:3] = growthRatesMain[:,:] 
                coordinates = findDeformedCoordinates(growthFull)
                f1 = np.sum(np.linalg.norm(coordinates-self.targetCoordinates,axis=1)) 
                f2 = np.sum(np.linalg.norm(x.reshape(8,3),axis=1))*4e5
                f = f1 + f2
                self.addSolution(x, f)
                currentAnswer = f
                if (currentAnswer<0.95*self.bestAnswer):
                    with open("answers.txt","a") as fileAnswer:
                        fileAnswer.write ('  f1=    %f  ' %f1)
                        fileAnswer.write ('   f2=    %f  ' %f2)
                        fileAnswer.write ('  totalObjctive=    %f \r\n' %currentAnswer)
                        now = datetime.datetime.now()
                        datetime.time(now.hour, now.minute, now.second)
                        fileAnswer.write ('  hour %f  ' %now.hour)
                        fileAnswer.write ('minute %f  ' %now.minute)                        
                        np.savetxt(fileAnswer, x, newline ='\n')
                        print ('Rate   ',x,'    objective',f)
                        fileAnswer.write ('    *************************	\n \n ')
                        fileAnswer.close()
                        self.bestAnswer = currentAnswer
            except:
                f = 1e12	
        if (f < 1e5):
            print ('Rate ',x,' objective ',f)
            sys.exit()
        return [f]

    def get_bounds(self):
        return ([-0.2]*8*3,[0.25]*8*3)
    
runParallel = False
import random
class MonteCarlo(pg.algorithm):
    def __init__(self,iterations = 10):
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

if __name__ == '__main__':
    gp = GrowthOptimization()
    prob = pg.problem(gp)
    algo = pg.algorithm(pg.pso(gen = 300))
    #algo = MonteCarlo(100)
    try:
        if not runParallel:
            pop = pg.population(prob,40)
            print(dir(pop))
            pop = algo.evolve(pop)
            
            print(pop.champion_f) 
        else:
            archi = pg.archipelago(n = 16, algo = algo, prob = prob, pop_size = 40, seed = 32)
            archi.evolve()
            print(archi)
            archi.wait()
            res = archi.get_champions_f()
            print(res) 
    except:
        traceback.print_exc(file=sys.stdout)
    
