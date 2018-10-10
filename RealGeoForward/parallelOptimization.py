from __future__ import print_function
import numpy as np
from Growth import TubeGrower  
import sys
import json
if __name__ == '__main__':
    data = json.load(sys.stdin)
    growthParameters = np.array(data['optmin'])
    #print (growthParameters, '@$@$@$@$@$@$@$@$@$@$@$@$@$@$')
    circumferentialElements=8
    axialElements=8
    wallElements=1
    discret=10
    length=6
    innerRadius=0.5
    outerRadius=2.0
    fixBottom=True
    fixTop = False
    DMBC = True
    humphrey = False
    neoHookean=True
    stage = 20
    division = 2
    layers = 1
    # growthParameters[0,:3] = 0.18,0.23,0.06
    # growthParameters[1,:3] = 0.19,0.14,0.05
    # growthParameters[2,:3] = -0.13,-0.03,0.04
    # growthParameters[3,:3] = -0.02,-0.12,0.06	
		
    # growthParameters[4,:3] = 0.18,0.23,0.06
    # growthParameters[5,:3] = 0.19,0.14,0.05
    # growthParameters[6,:3] = -0.13,-0.03,0.04
    # growthParameters[7,:3] = -0.02,-0.12,0.06	

    # growthParameters[8,:3] = 0.18,0.23,0.06
    # growthParameters[9,:3] = 0.19,0.14,0.05
    # growthParameters[10,:3] = -0.13,-0.03,0.04
    # growthParameters[11,:3] = -0.02,-0.12,0.06	

    # growthParameters[0:4,3:]=0.0
    # growthParameters[4:8,3:]=0.0
    # growthParameters[8:,3:]=0.0
    growthLogic = TubeGrower(circumferentialElements,axialElements,wallElements,discret,length,innerRadius,outerRadius,fixBottom,fixTop, DMBC,humphrey, neoHookean,stage,division, layers)
    growthLogic.setupProblem()
    filename = "HeartTubeGrowth_" + "Undeformed"
    growthLogic.saveResults(filename)
    try:
        targetCoordinates = growthLogic.solveAndGetSurfaceDescriptors(growthParameters)
        filename = "HeartTubeGrowth_" + "deformed"
        growthLogic.saveResults(filename)		
        result = {'nodePositions':targetCoordinates.tolist()}
    except Exception as e:
        result={'error':str(e)}
    print('POSResSTART\n%s\nPOSResEND'%json.dumps(result),file=sys.stderr)

# if __name__ == '__main__':
    # obj = TubeGrower(4,3,1,DMBC=True, humphrey = False, neoHookean=True)
    # obj.setupProblem(True)
    # growthParameters= np.zeros((12,6))
    # # growthParameters[4:8,:3]=0.0
    # growthParameters[0,:3] = 0.05,0.08,0.14
    # growthParameters[1,:3] = 0.08,0.06,0.08
    # growthParameters[2,:3] = 0.13,-0.14,0.07
    # growthParameters[3,:3] = -0.04,-0.13,0.07
    # growthParameters[4,:3] = 0.18,0.23,0.06
    # growthParameters[5,:3] = 0.19,0.14,0.05
    # growthParameters[6,:3] = -0.13,-0.03,0.04
    # growthParameters[7,:3] = -0.02,-0.12,0.06	
    # growthParameters[8,:3] = 0.09,0.01,-0.06
    # growthParameters[9,:3] = 0.16,0.03,-0.15
    # growthParameters[10,:3] = 0.13,-0.12,0.04	
    # growthParameters[11,:3] = -0.08,-0.16,0.06
    # growthParameters[0:4,3:]=0.2
    # growthParameters[4:8,3:]=0.1
    # growthParameters[8:,3:]=0.1	
    # obj.solveAndGetSurfaceDescriptors(growthParameters)
    # obj.saveResults('result')