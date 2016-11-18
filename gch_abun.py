# -*- coding: utf-8 -*-

import sys
import pdb
import math
from math import *
import numpy as np

#COMPUTES THEORETICAL DATA FOR STELLAR MODEL INPUTS USING GCH MODEL
#INPUT: Metallicity [Z] 
#OUTPUT: .txt with isotope abundances computed at input [Z]
#
# Output txt files already computed for:
# [Z]=-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0,0.1,0.2
# Can be found in isotopes directory

#LOAD DATA FILES
ions=open('data/ions.dat').read().split() #ion names for output

arr = np.loadtxt('data/x_matrix_10_13_2013.dat') #solar abundance pattern decomposition
for i in range(len(arr)):
    for j in range(len(arr[i])):
       arr[i][j]=float(arr[i][j])

isomasses=[]
mass_point = open('data/isomasses_numbers.dat','r') #isotope masses (in amu) for conversion to mass fractions
masses = mass_point.readlines()
for j in range(len(masses)-5):
    isomasses.append(float(masses[j+5]))

#VARIABLES
logz = input('Enter [Z]: ')
logz=float(logz)
numiso=287
z=float(10.0**logz) #z stands for z/z_solar

#Compute Scaling Array
matrix=[]
matrix.append(z**1.227)
matrix.append(z**1.509)
matrix.append(1.0)
matrix.append(-2e-11*(1-(tanh(200.0*z-0.23)/tanh(200.0-0.23))))
matrix.append(z**1.230)
matrix.append(z**.938)
matrix.append(z**.938)
matrix.append(z**((1.509+.938)/2))
matrix.append(z*((tanh(5.024*z-2.722))+tanh(2.722))/(tanh(5.024-2.722)+tanh(2.722)))
matrix.append(logz+2.4983277)
matrix.append(1.0)
matrix.append(z**1.509)
matrix.append(z**.938)
matrix.append(1.0)
matrix.append(z**.938)
matrix.append(1.0)
matrix.append(z)

#CUT OFF STRONG S-PROCESS SCALING IF ABUN HITS ZERO
if matrix[3]<0.0:
    matrix[3]=0

for x in range(9,76):
    arr[x][10]=math.log10(arr[x][10])

abu1=np.dot([i[0:9] for i in arr],matrix[0:9])  
abu2=np.dot([i[11:17] for i in arr],matrix[11:17])

base=float(10)
sne=np.zeros(287,dtype=float)
for x in range(9,76):
    sne[x]=np.power(base,(np.dot(matrix[9],arr[x][9]) + np.dot(matrix[10],arr[x][10])))

abu=abu1+sne+abu2

#convert abu to mass fractions
for i in range(len(abu)):
    abu[i]=abu[i]*isomasses[i]

#OUTPUT ABUNDANCE FILE
if logz < 0:
    name="m"+str(abs(logz)).split(".")[0]+"p"+str(abs(logz)).split(".")[1]
else:
    name="p"+str(abs(logz)).split(".")[0]+"p"+str(abs(logz)).split(".")[1]
file=open("isotopes/gch_abun_"+name+".txt","w")
file.write("GCH Model Abundances (mass fractions)\n")
file.write("[Z]= ")
file.write(str(logz)+"\n") 
for i in range(len(abu)):
    file.write('{}'.format(ions[i]))
    file.write(" ")
    file.write('{}\n'.format(abu[i]))
file.close()
