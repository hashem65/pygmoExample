from scipy.interpolate import griddata
#from matplotlib.mlab import griddata
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from numpy.random import uniform
import sys,os
import math
from numpy.random import randint

layers = 1  #  can be 1 or 2              1 for only cardiac muscle  and  2  for cardiac muscle and cardiac jelly
division = 4     # can be 2 or 4              2 for coarse meshgrid 18 vars and 4 for refine meshgrid with 36 vars 

numberOfCircumferentialPoints = division + 1
numberOfLongitudinalPoints = division + 1 

numberOfLongitudinalElements = 8
numberOfCircumferentialElements = 8
totalNumberOfElements = 64*layers


grid_x, grid_y = np.mgrid[0:1:1000j, 0:1:1000j]

# generate some random numbers for the field
valuesFibre1 = randint(-20, 25, (division+1)*(division+1))*0.01
valuesSheet1 = randint(-20, 25, (division+1)*(division+1))*0.01
valuesNormal1 = randint(-20, 25, (division+1)*(division+1))*0.01
if (layers == 2 ):
    valuesFibre2 = randint(-20, 25, (division+1)*(division+1))*0.01
    valuesSheet2 = randint(-20, 25, (division+1)*(division+1))*0.01
    valuesNormal2 = randint(-20, 25, (division+1)*(division+1))*0.01	
values = np.zeros(((division+1)*(division+1)*layers,4,layers))
values [:,0,0] = valuesFibre1[:]
values [:,1,0] = valuesSheet1[:]
values [:,2,0] = valuesNormal1[:]
if (layers == 2 ):
    values [:,0,1] = valuesFibre2[:]
    values [:,1,1] = valuesSheet2[:]
    values [:,2,1] = valuesNormal2[:]
# there needs to be a symmetry as of unrolling cylinder 	
for j in range (division+1):
    values[(division+1)*j+division,:,:] = values[(division+1)*j,:,:] 
values [:,3,:] = values[:,0,:]*values[:,0,:]+values[:,1,:]*values[:,1,:]+values[:,2,:]*values[:,2,:]
values1 = np.zeros(((division+1)*(division+1)*layers,4))
values1 [:,:] = values [:,:,0]
if (layers == 2 ):
    values2 = np.zeros(((division+1)*(division+1)*layers,4))
    values2 [:,:] = values [:,:,1]



        #if (layers == 1):
        #    values[(division+1)*j+division,:] = values[(division+1)*j,:] 
        #elif (layers == 2):
        #    values[(division+1)*j+division,:] = values[(division+1)*j,:] 
        #    values[(division+1)*(division+1)+(division+1)*j+division,:] = values[(division+1)*(division+1)+(division+1)*j,:]         

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
grid_z30 = griddata(points, values1[:,3], (grid_x, grid_y), method='cubic')
if (layers == 2 ):
    grid_z01 = griddata(points, values2[:,0], (grid_x, grid_y), method='cubic')
    grid_z11 = griddata(points, values2[:,1], (grid_x, grid_y), method='cubic')
    grid_z21 = griddata(points, values2[:,2], (grid_x, grid_y), method='cubic')
    grid_z31 = griddata(points, values2[:,3], (grid_x, grid_y), method='cubic')
# the location of gauss points  in [-1,1]  -0.577, 0 , 0.577
# the location of gauss points  in [0,1]  (1-0.577)/2, 0.5 , (1+0.577)/2   ===>  0.211, 0.5, 0.789  or 0.21 , 0.5 , 0.79
elementNumber = 0
elementOrigin = []
gaussPoints = [0.21 , 0.5 , 0.79]

#xiLocation = []
#allXiPoints = np.zeros ((totalNumberOfElements*numberOfXiPoints*numberOfXiPoints,2))
#for longitudinalElementIndex in range (numberOfLongitudinalElements):    # 0-7
#    for CircumferentialElementIndex in range (numberOfCircumferentialElements):  # 0-7
#        elementNumber = longitudinalElementIndex*numberOfCircumferentialElements + CircumferentialElementIndex
#        elementOrigin =	[(longitudinalElementIndex)/numberOfLongitudinalElements,(CircumferentialElementIndex)/numberOfCircumferentialElements ]
#        for k in range (numberOfXiPoints):
#            for m in range (numberOfXiPoints):
#                xiLocation = [elementOrigin[0]+(gaussPoints[k]) ,elementOrigin[1]+gaussPoints[m]]
#                allXiPoints[elementNumber*numberOfXiPoints*numberOfXiPoints + k*numberOfXiPoints + m,: ]= xiLocation

interpolatedRates = np.zeros((totalNumberOfElements,3,3,3,3))
for numLayers in range (layers):
    for lengthElems in range (8):
        for circumElems in range (8):
            elementNumber = numLayers*64 + lengthElems*8 + circumElems
            for xi1 in range (3):
                for xi2 in range (3):
                    for xi3 in range (3):
                        s = circumElems/8 + gaussPoints[xi1]/8
                        t = lengthElems/8 + gaussPoints[xi2]/8
                        p = gaussPoints[xi3]
                        if (layers == 1 ):
                            interpolatedRates[elementNumber, xi1, xi2,xi3,0] = grid_z00[int(1000*t) , int(1000*s)]
                            interpolatedRates[elementNumber, xi1, xi2,xi3,1] = grid_z10[int(1000*t) , int(1000*s)] 
                            interpolatedRates[elementNumber, xi1, xi2,xi3,2] = grid_z20[int(1000*t) , int(1000*s)]
                        elif (layers == 2 ):
                            interpolatedRates[elementNumber, xi1, xi2,xi3,0] = ((1-p)*(grid_z00[int(1000*t),int(1000*s)])  +  p*(grid_z00[int(1000*t),
                                                int(1000*s)] + grid_z01[int(1000*t) , int(1000*s)])/2)*math.abs(layers-1)  +  (p*(grid_z01[int(1000*t),
                                                int(1000*s)])  +  (1-p)*(grid_z00[int(1000*t),int(1000*s)] + grid_z01[int(1000*t) , int(1000*s)])/2)*math.abs(layers)
                            interpolatedRates[elementNumber, xi1, xi2,xi3,1] = ((1-p)*(grid_z10[int(1000*t),int(1000*s)])  +  p*(grid_z10[int(1000*t),
                                                int(1000*s)] + grid_z11[int(1000*t) , int(1000*s)])/2)*math.abs(layers-1)  +  (p*(grid_z11[int(1000*t),
                                                int(1000*s)])  +  (1-p)*(grid_z10[int(1000*t),int(1000*s)] + grid_z11[int(1000*t) , int(1000*s)])/2)*math.abs(layers)
                            interpolatedRates[elementNumber, xi1, xi2,xi3,2] = ((1-p)*(grid_z20[int(1000*t),int(1000*s)])  +  p*(grid_z20[int(1000*t),
                                                int(1000*s)] + grid_z21[int(1000*t) , int(1000*s)])/2)*math.abs(layers-1)  +  (p*(grid_z21[int(1000*t),
                                                int(1000*s)])  +  (1-p)*(grid_z20[int(1000*t),int(1000*s)] + grid_z21[int(1000*t) , int(1000*s)])/2)*math.abs(layers)   

# geting the gradient and calculating their magnitudes  and later on plotting them to if there is enough continuity in the gradient 
#gx00,gy00 = np.gradient(grid_z00)
#gx10,gy10 = np.gradient(grid_z10)
#gx20,gy20 = np.gradient(grid_z20)
#gx30,gy30 = np.gradient(grid_z30)
#gz00 = (gx00*gx00 + gy00*gy00)
#gz01 = (gx01*gx01 + gy01*gy01)
#gz02 = (gx02*gx02 + gy02*gy02)
#gz03 = (gx03*gx03 + gy03*gy03)
#if (layers == 2 ):
#     gx01,gy01 = np.gradient(grid_z01)
#     gx11,gy11 = np.gradient(grid_z11)
#     gx21,gy21 = np.gradient(grid_z21)
#     gx31,gy31 = np.gradient(grid_z31)
#     gz01 = (gx10*gx10 + gy10*gy10)
#     gz11 = (gx11*gx11 + gy11*gy11)
#     gz21 = (gx21*gx21 + gy21*gy21)
#     gz31 = (gx31*gx31 + gy31*gy31)



origin = 'lower' #origin = 'upper'
fig, ax = plt.subplots(nrows=3, ncols=3)
plt.subplot(331)
plt.imshow(grid_z00.T, extent=(0.0,1.0,0.0,1.0), origin=origin)
plt.contourf(grid_x,grid_y,grid_z00.T,30)
plt.title('FibreRates')

plt.subplot(332)
plt.imshow(grid_z10.T, extent=(0.0,1.0,0.0,1.0), origin=origin)
plt.contourf(grid_x,grid_y,grid_z10.T,30) # np.arange(0,1.0,0.01))
plt.title('SheetRates')

plt.subplot(333)
plt.imshow(grid_z20.T, extent=(0.0,1.0,0.0,1.0), origin=origin)
plt.contourf(grid_x,grid_y,grid_z20.T,30)
plt.title('NormalRates')

plt.subplot(334)
#plt.imshow(values, points, extent=(0,1,0,1), origin='lower')
plt.imshow(grid_z30.T, extent=(0.0,1.0,0.0,1.0), origin=origin)
plt.contourf(grid_x,grid_y,grid_z30.T,30)
plt.plot(points[:,0], points[:,1], 'k.', ms=1)
plt.title('OriginalPoints')

plt.subplot(335)
plt.imshow(gz)
plt.contourf(grid_x,grid_y,gz,15)
plt.title("gx")

plt.subplot(336)
plt.imshow(gz)
plt.contourf(grid_x,grid_y,gz,15)
plt.title("gy")

plt.subplot(337)
plt.scatter(plottingXiPoints[:,0], plottingXiPoints[:,1], c=plottingXiPointsValuesFibre[:] , marker=(2, 1))
plt.plot(points[:,0], points[:,1], 'k.', ms=1)
plt.title('GaussPoints')

plt.subplot(338)
plt.scatter(plottingXiPoints[:,0], plottingXiPoints[:,1], c=plottingXiPointsValuesFibre[:] , marker=(2, 1))
plt.plot(points[:,0], points[:,1], 'k.', ms=1)
plt.title('GaussPoints')

plt.subplot(339)
plt.scatter(plottingXiPoints[:,0], plottingXiPoints[:,1], c=plottingXiPointsValuesFibre[:] , marker=(2, 1))
plt.plot(points[:,0], points[:,1], 'k.', ms=1)
plt.title('GaussPoints')


plt.gcf().set_size_inches(6,6)
plt.show()