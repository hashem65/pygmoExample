from __future__ import print_function
import numpy as np
from Growth import TubeGrower  
import sys
import json
if __name__ == '__main__':
    data = json.load(sys.stdin)
    growthParameters = np.array(data['optmin'])
    #print (growthParameters, '##########################')
    circumferentialElements=4
    axialElements=2
    wallElements=1
    discret=10
    length=6
    innerRadius=0.5
    outerRadius=2.0
    fixBottom=True
    fixTop = False
    DMBC = False
    humphrey = False
    neoHookean=True
#    stage = 0
#    division = 2
#    layers = 1
    growthLogic = TubeGrower(circumferentialElements,axialElements,wallElements,discret,length,innerRadius,outerRadius,fixBottom,fixTop, DMBC,humphrey, neoHookean) # ,stage,division, layers)
    growthLogic.setupProblem()
    #filename = "HeartTubeGrowth_" + "Undeformed"
    #growthLogic.saveResults(filename)
    try:
        targetCoordinates = growthLogic.solveAndGetSurfaceDescriptors(growthParameters)
        #filename = "HeartTubeGrowth_" + "target"
        #filename = "HeartTubeGrowth_" + "deformed"
        #growthLogic.saveResults(filename)		
        result = {'nodePositions':targetCoordinates.tolist()}
    except Exception as e:
        result={'error':str(e)}
    print('POSResSTART\n%s\nPOSResEND'%json.dumps(result),file=sys.stderr)