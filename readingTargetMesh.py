'''
Created on 20/09/2017

@author: rjag008
'''
from __future__ import print_function
import numpy as np
from Growth import TubeGrower
    
import sys
import json

#from opencmiss.iron import iron

if __name__ == '__main__':
    data = json.load(sys.stdin)
    growthParameters = np.array(data['optmin'])
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
    
	# write out the files ... 
	
    try:
        targetCoordinates = growthLogic.solveAndGetSurfaceDescriptors(growthParameters)
        result = {'nodePositions':targetCoordinates.tolist()}
    except Exception as e:
        result={'error':str(e)}
    print('POSResSTART\n%s\nPOSResEND'%json.dumps(result),file=sys.stderr)