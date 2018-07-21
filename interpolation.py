'''
16/07/2018
hashem 
'''
import random 
import sys,os
import numpy as np
from opencmiss.iron import iron


numberOfLengthElements = 8 
numberOfCircumferentialElements = 8 
numberOfWallElements = 2
numberOfLengthNodes = numberOfLengthElements + 1
numberOfCircumferentialNodes = numberOfCircumferentialElements
numberOfWallNodes = numberOfWallElements + 1
numberOfElements = numberOfLengthElements*numberOfCircumferentialElements*numberOfWallElements
numberOfNodes = numberOfLengthNodes*numberOfCircumferentialNodes*numberOfWallNodes


### Based on the number of points used in the interpolation, the distribution has been structured ... 
### number of points for growth net = number of nodes * 2/3 
division = 2    
# division = 4    
numberOfPoints = int((numberOfLengthNodes)*(numberOfCircumferentialNodes)*(numberOfWallElements))                                  # 9*8*2
numberOfVariables = int((numberOfLengthElements/division + 1)*(numberOfCircumferentialElements/division)*(numberOfWallElements))   # 5*4*2
print ('numberOfPoints = ', numberOfPoints)
print ('numberOfVariables = ', numberOfVariables)

# generate some random numbers as the objective variables in a 4x5x2x3 net  
def random_floats(low,high):
    return [random.uniform(low, high)]
low = -0.15
high = 0.35
growthRateVariables = np.zeros((numberOfVariables,3))                                                                               # (8*9*2)*(3)
for i in range (numberOfVariables):
    for j in range (3):
        a = random_floats(low,high)
        growthRateVariables[i,j] = a[0]
print (growthRateVariables)
print ('\n'*3)


# ==================================================
# map from 4x5x2x3 points to 8x9x2x3 points 
#growthRateVariables = np.reshape(random_floats,(numberOfVariables,3))   
growthRatePoints = np.zeros((numberOfPoints,3))                                                                                    # (8*9*2)*(3)
for i in range(1,numberOfWallElements+1):                                                                                          # (1,2)
    for j in range(1,int(numberOfLengthElements/division)+2):                                                                      # (1,10)
        for k in range(1, int(numberOfCircumferentialElements/division)+1):                                                        # (1,9)
            growthRatePoints [(2*k-1)+(2*(j-1)*8)+(i-1)*8*9 - 1 ,:] = growthRateVariables[k+(j-1)*4+(i-1)*4*5 - 1,:]
            print ((2*k-1)+(2*(j-1)*8)+(i-1)*8*9 -1 , k+(j-1)*4+(i-1)*4*5 - 1 )
print ('\n'*6)
#print (growthRatePoints)


# ==================================================
# interpolate in the empty spots in 8x9x2x3         
# part 1 ... interpolate in half-filled lines ...                             
for j in range (int(numberOfPoints/numberOfCircumferentialElements)):                                                              # (0:9*2)
    for i in range (int(numberOfCircumferentialElements/2)):                                                                       # (0:8)
        if not (i == 3):
            growthRatePoints [j*8+2*i+1,:] = (growthRatePoints [j*8+2*i,:] + growthRatePoints [j*8+2*i+2,:])/2
        else:
            growthRatePoints [j*8+2*i+1,:] = (growthRatePoints [j*8+2*i,:] + growthRatePoints [j*8,:])/2
print ('\n'*9)
print (growthRatePoints)

# part 2 .... interpolate in empty lines 
for i in range (numberOfWallElements):
    for k in range (int (numberOfLengthNodes/2)):
        for j in range (numberOfCircumferentialNodes):
            pointNumber =  (2*k+1)*8 +  j  + i*72
            growthRatePoints[pointNumber,:] = (growthRatePoints[pointNumber-8,:] + growthRatePoints[pointNumber+8,:])/2
print ('\n'*12)
print (growthRatePoints)


# ==================================================
growthRate = np.zeros((numberOfElements,3))                                                                                        # (8*9*2)*(3)
# map from points to elements 
elementNumber = 0
for wallElementIdx in range(1,numberOfWallElements+1):                                                                             # (1:2)
    for lengthElementIdx in range(1,numberOfLengthElements+1):
        for circumferentialElementIdx in range(1,numberOfCircumferentialElements+1):
            elementNumber = elementNumber + 1
            localNode1 = circumferentialElementIdx + (lengthElementIdx-1)*numberOfCircumferentialNodes + (wallElementIdx-1)*numberOfCircumferentialNodes*numberOfLengthNodes
            if circumferentialElementIdx == numberOfCircumferentialElements:
                localNode2 = 1 + (lengthElementIdx-1)*numberOfCircumferentialNodes + (wallElementIdx-1)*numberOfCircumferentialNodes*numberOfLengthNodes
            else:
                localNode2 = localNode1 + 1
            localNode3 = localNode1 + numberOfCircumferentialNodes
            localNode4 = localNode2 + numberOfCircumferentialNodes
            growthRate[elementNumber-1,:] = (growthRatePoints[localNode1-1,:] + growthRatePoints[localNode1-1,:] + growthRatePoints[localNode1-1,:] + growthRatePoints[localNode1-1,:])/4
print ('\n'*12)
print ('length Matrix growth rates=',  len(growthRate))
print ('\n'*12)
print (growthRate)

# the end  - fully checked 