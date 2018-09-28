from __future__ import print_function
import numpy as np
from Growth import TubeGrower
    
import sys
import json

#from opencmiss.iron import iron
if __name__ == '__main__':
    data = json.load(sys.stdin)
    growthParameters = np.array(data['optmin'])
    circumferentialElements=4
    axialElements=2
    wallElements=1
    discret=10
    length=4
    innerRadius=1.0
    outerRadius=2.0
    fixBottom=True
    fixTop = False
    DMBC = True
    humphrey = False
    neoHookean=True
#    stage = 0
#    division = 2
#    layers = 1
    # growthParameters[0,:3] = 0.05,0.08,0.14
    # growthParameters[1,:3] = 0.08,0.06,0.08
    # growthParameters[2,:3] = 0.13,-0.14,0.07
    # growthParameters[3,:3] = -0.04,-0.13,0.07

    # growthParameters[4,:3] = 0.18,0.23,0.06
    # growthParameters[5,:3] = 0.19,0.14,0.05
    # growthParameters[6,:3] = -0.13,-0.03,0.04
    # growthParameters[7,:3] = -0.02,-0.12,0.06	

    # growthParameters[0:4,3:]=0.0
    # growthParameters[4:8,3:]=0.0

    growthLogic = TubeGrower(circumferentialElements,axialElements,wallElements,discret,length,innerRadius,outerRadius,fixBottom,fixTop, DMBC,humphrey, neoHookean) # ,stage,division, layers)
    growthLogic.setupProblem()
    #filename = "HeartTubeGrowth_" + "Undeformed"
    #growthLogic.saveResults(filename)
    try:
        targetCoordinates = growthLogic.solveAndGetSurfaceDescriptors(growthParameters)
        #filename = "HeartTubeGrowth_" + "target"
        filename = "HeartTubeGrowth_" + "deformed"
        growthLogic.saveResults(filename)
        result = {'nodePositions':targetCoordinates.tolist()}
    except Exception as e:
        result={'error':str(e)}
    print('POSResSTART\n%s\nPOSResEND'%json.dumps(result),file=sys.stderr)
