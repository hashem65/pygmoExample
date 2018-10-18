from __future__ import print_function
import numpy as np
import json,os,sys
import traceback
import datetime 
import math 

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

# def findTargetCoordinates(rates):
    # try:
        # p = Popen(['python','readingTargetMesh.py'],stdin=PIPE,stderr=PIPE,env=my_env)
        # inputData = dict()
        # inputData['optmin'] = rates.tolist() 
        # subin = json.dumps(inputData)
        # #print('##Initializing problem for parameters %s ' % (' ,'.join(map(str,rates.tolist()))))
        # _,perr = p.communicate(str.encode(subin))
        # p.wait()
        # perr = perr.decode()
        # #Ensure that we are starting at a json prefix
        # ps = perr.index('POSResSTART') + 11
        # pe = perr.index('POSResEND')
        # result = json.loads(perr[ps:pe])
        # if u'error' in result:
            # raise ValueError(result[u'error'])
        # return np.array(result[u'nodePositions'])
    # except:
        # traceback.print_exc(file=sys.stdout)
		
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

dictionaryLoc = 'solutions431.sql'
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
        growthRate[0,:3] = 0.18,0.0,0.0
        growthRate[1,:3] = 0.16,0.0,0.0
        growthRate[2,:3] = 0.17,0.0,0.0
        growthRate[3,:3] = 0.19,0.0,0.0	
		
        growthRate[4,:3] = 0.17,0.0,0.0
        growthRate[5,:3] = 0.15,0.0,0.0
        growthRate[6,:3] = 0.16,0.0,0.0
        growthRate[7,:3] = 0.18,0.0,0.0	
	

        growthRate[0:4,3:]=0.0
        growthRate[4:8,3:]=0.0
        growthRate[8:,3:]=0.0
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
        #print ('##############$$$$$$$$$$$$$$$$$$$$$$$$')
        if f is None:
            try: 
                growthFull = np.zeros((8,6))
                continuity = np.zeros((8,1))
                exactAnswer = np.zeros((8,1))
                growthRatesMain = np.zeros((8,1))
                growthRatesCross = np.zeros((8,5))
                growthRatesMain[:,:] = x.reshape(8,1)
                #growthRatesMain[4:8,:] = x.reshape(4,3) 
                #growthRatesMain[8:,:] = x.reshape(4,3)
                growthFull[:,0:1] = growthRatesMain[:,:] 
                coordinates = findDeformedCoordinates(growthFull)
                exactAnswer[:,0] = 0.18,0.16,0.17,0.19,0.17,0.15,0.16,0.18
                #growthFull[:,0] = 0.18,0.16,0.17,0.19,0.17,0.15,0.16,0.18
                #x[:] =  0.18,0.16,0.17,0.19,0.17,0.15,0.16,0.18
                #print (exactAnswer.shape)
                f1 = np.sum(np.linalg.norm(coordinates-self.targetCoordinates,axis=1)) 
                fit = np.sum(np.linalg.norm(x.reshape(8,1)-exactAnswer,axis=1))
                #print (f1,'%%%%%%%%%%%%%%%%%%%%%%%%')
                #print (growthFull)
                #print (fit,'%%%%%%%%%%%%%%%%%%%%%%%%')
                continuity[0,0]=(math.pow(growthFull[1,0]-growthFull[0,0],2)+math.pow(growthFull[3,0]-growthFull[0,0],2)+math.pow(growthFull[4,0]-growthFull[0,0],2)+math.pow(growthFull[5,0]-growthFull[0,0],2)/2+math.pow(growthFull[7,0]-growthFull[0,0],2)/2)/4
                continuity[1,0]=(math.pow(growthFull[0,0]-growthFull[1,0],2)+math.pow(growthFull[2,0]-growthFull[1,0],2)+math.pow(growthFull[5,0]-growthFull[1,0],2)+math.pow(growthFull[4,0]-growthFull[1,0],2)/2+math.pow(growthFull[6,0]-growthFull[1,0],2)/2)/4
                continuity[2,0]=(math.pow(growthFull[1,0]-growthFull[2,0],2)+math.pow(growthFull[3,0]-growthFull[2,0],2)+math.pow(growthFull[6,0]-growthFull[2,0],2)+math.pow(growthFull[5,0]-growthFull[2,0],2)/2+math.pow(growthFull[7,0]-growthFull[2,0],2)/2)/4
                continuity[3,0]=(math.pow(growthFull[0,0]-growthFull[3,0],2)+math.pow(growthFull[2,0]-growthFull[3,0],2)+math.pow(growthFull[7,0]-growthFull[3,0],2)+math.pow(growthFull[4,0]-growthFull[3,0],2)/2+math.pow(growthFull[6,0]-growthFull[3,0],2)/2)/4
                # The middle ones 
                continuity[4,0]=(math.pow(growthFull[0,0]-growthFull[4,0],2)+math.pow(growthFull[5,0]-growthFull[4,0],2)+math.pow(growthFull[7,0]-growthFull[4,0],2)+math.pow(growthFull[1,0]-growthFull[4,0],2)/2+math.pow(growthFull[3,0]-growthFull[4,0],2)/2)/4
                continuity[5,0]=(math.pow(growthFull[1,0]-growthFull[5,0],2)+math.pow(growthFull[4,0]-growthFull[5,0],2)+math.pow(growthFull[6,0]-growthFull[5,0],2)+math.pow(growthFull[0,0]-growthFull[5,0],2)/2+math.pow(growthFull[2,0]-growthFull[5,0],2)/2)/4
                continuity[6,0]=(math.pow(growthFull[2,0]-growthFull[6,0],2)+math.pow(growthFull[5,0]-growthFull[6,0],2)+math.pow(growthFull[7,0]-growthFull[6,0],2)+math.pow(growthFull[1,0]-growthFull[6,0],2)/2+math.pow(growthFull[3,0]-growthFull[6,0],2)/2)/4
                continuity[7,0]=(math.pow(growthFull[3,0]-growthFull[7,0],2)+math.pow(growthFull[6,0]-growthFull[7,0],2)+math.pow(growthFull[4,0]-growthFull[7,0],2)+math.pow(growthFull[0,0]-growthFull[7,0],2)/2+math.pow(growthFull[2,0]-growthFull[7,0],2)/2)/4				
                #print ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',continuity)
                # f2 = np.sum(np.linalg.norm(x.reshape(4,3),axis=1))*4e6
                f2 = np.sum(np.linalg.norm(continuity,axis=1))    # *5e8      # 0.68 minimum 
                #print ('norm',f2)
                f= fit + 5*(f2 - 0.0019)   #   +f1
                #print (f1,f2,f1,f2,f1,f2,f1,f2)
                self.addSolution(x, f)
                currentAnswer = f
                if (currentAnswer<0.95*self.bestAnswer):
                    with open("answers.txt","a") as fileAnswer:
                        fileAnswer.write ('  f1=    %f  ' %f1)
                        fileAnswer.write ('   f2=    %f  ' %f2)
                        fileAnswer.write ('   fit=    %f  ' %fit)                        
                        fileAnswer.write ('  totalObjctive=    %f \r\n' %currentAnswer)
                        now = datetime.datetime.now()
                        datetime.time(now.hour, now.minute, now.second)
                        fileAnswer.write ('  hour %f  ' %now.hour)
                        fileAnswer.write ('minute %f  ' %now.minute)
                        #fileAnswer.write ('  rates =    %d \r\n' np.c_[x])
                        np.savetxt(fileAnswer, x, newline ='\n')
                        print ('Rate   ',x,'    objective',f)
                        fileAnswer.write ('    *************************	\n \n ')
                        fileAnswer.close()
                    self.bestAnswer = currentAnswer
            except:
                f = 1e12	
            if (f < 1e-4):
                print ('Rate ',x,' objective ',f)
                sys.exit()
        return [f]

    def get_bounds(self):
        return ([-0.2]*8,[0.25]*8)
    
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
            archi = pg.archipelago(n = 16, algo = algo, prob = prob, pop_size = 20, seed = 32)
            archi.evolve()
            print(archi)
            archi.wait()
            res = archi.get_champions_f()
            print(res) 
    except:
        traceback.print_exc(file=sys.stdout)
    
