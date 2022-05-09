### This file produced plots which depicts the stretched cores during impact, isolated, and calculates their length ###


import math, glob
import os, statistics
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
from scipy.misc import derivative 
from scipy import integrate
import numpy as np
from sigfig import round

class target:
    radii = []
    center = []
class impactor:
    radii = []
    center = []

#Converts a number in scientific notation to a float
def sci2Float (sciNot):
    A = sciNot.split('e')
    return float(A[0])*math.pow(10, int(A[1]))

#returns the index of stretched particles outside of the target
def outStretch (target, x, y, leng):
    indices = []
    for i in range(0, leng):
        if (1 < (((target.center[0]-x[i])/(target.radii[0]))**2 + ((target.center[1] - y[i])/(target.radii[1]))**2)):
            indices.append(i)
    return(indices)

#returns the index of stretched particles inside of the target at the end
def acc (target, x, y, leng):
    indices = []
    for i in range(0, leng):
        if (1 > (((target.center[0]-x[i])/(target.radii[0]))**2 + ((target.center[1] - y[i])/(target.radii[1]))**2)):
            indices.append(i)
    return(indices)

def func(x, p0, p1, p2):
    P = p0 + p1*x + p2*0.5*(3*x**2 - 1)
    return P

#   Global variable for function coefficients
C = None

def funcWrap(x):
    return func(x, *C)

def integrand(x):
    return math.sqrt(1 + derivative(funcWrap, x, dx=1e-06)**2)



directory = input('Enter directory name for data set: ')
while(True):
    files = sorted(glob.glob('./'+directory+'/*.dat'))
    if len(files) == 0: 
        directory = input('Non-existent or empty directory. Please try again: ')
    else:
        break    

sep = int(input('Separation index between target and impactor: '))



trgt = target()
impctr = impactor()
time = []
accPart = None

files = [file for file in files if not(str(file).endswith("disk.dat") or str(file).endswith("diss.dat"))]
files = files[:-1]
files.insert(1, files[-1])

for file in files:
    with open(file) as r:
        lines = r.readlines()[2:]
    print(str(file)[-9:-4])
    curr_time = sci2Float(lines[0])
    lines = lines[1:]
    x_t, x_i = [],[]
    y_t, y_i = [],[]
    z_t, z_i = [],[]
    iron_t, iron_i = [],[]
    
    for line in lines:
        coord = sci2Float(line.split()[1])
        x_t.append(coord)

    for line in lines:
        coord = sci2Float(line.split()[2])
        y_t.append(coord)
        
    for line in lines:
        coord = sci2Float(line.split()[3])
        z_t.append(coord)
        
    for line in lines:
        iron = int(line.split()[7])
        if iron==3:
            iron_t.append(True)
        else:
            iron_t.append(False)
    
    x_i = x_t[sep:]
    y_i = y_t[sep:]
    z_i = z_t[sep:]
    iron_i = iron_t[sep:]
    
    x_t = x_t[:sep-1]
    y_t = y_t[:sep-1]
    z_t = z_t[:sep-1]
    iron_t = iron_t[:sep-1]
    
    trgt.center = [sum(x_t)/len(x_t), sum(y_t)/len(y_t), sum(z_t)/len(z_t)]
    impctr.center = [sum(x_i)/len(x_i), sum(y_i)/len(y_i), sum(z_i)/len(z_i)]
    time.append(curr_time)
    
    if file == files[0]:
        trgt.radii = [abs(trgt.center[0]-max(x_t)), abs(trgt.center[1]-max(y_t)), abs(trgt.center[2]-max(z_t))]
        impctr.radii = [abs(trgt.center[0]-max(x_i)), abs(trgt.center[1]-max(y_i)), abs(trgt.center[2]-max(z_i))]
        continue
    
    xlim = 2*(trgt.radii[1] + impctr.radii[1])
    ylim = 2*(trgt.radii[1] + impctr.radii[1])
    cntrH = trgt.center + impctr.center
    cntr = [i/2 for i in cntrH]
    
    plt.xlim(cntr[0] - xlim, cntr[0] + xlim)
    plt.ylim(cntr[1] - ylim, cntr[1] + ylim)
    
    x_i_rn = []
    x_i_mntl = []
    y_i_rn = []
    y_i_mntl = []

    for i in range(len(iron_i)):
        if iron_i[i]:
            x_i_rn.append(x_i[i])
            y_i_rn.append(y_i[i])
        else:
            x_i_mntl.append(x_i[i])
            y_i_mntl.append(y_i[i])

#    plt.scatter(x_i_mntl, y_i_mntl, c='red')
            
    x_t_rn = []
    x_t_mntl = []
    y_t_rn = []
    y_t_mntl = []

    for i in range(len(iron_t)):
        if iron_t[i]:
            x_t_rn.append(x_t[i])
            y_t_rn.append(y_t[i])
        else:
            x_t_mntl.append(x_t[i])
            y_t_mntl.append(y_t[i])

#   Get accreted particles
    if files.index(file) == 1:
        accPart = acc(trgt, x_i_rn, y_i_rn, len(x_i_rn))
        continue
            
#   Get coordinates for all accreted iron particles outside of the target 
    outPart = outStretch(trgt, x_i_rn, y_i_rn, len(x_i_rn))
    outPart_x = []
    outPart_y = []
    for p in outPart:
        if p in accPart:
            outPart_x.append(x_i_rn[p])
            outPart_y.append(y_i_rn[p])
        
#   Perform COM-crawl to get the "useful" particles that are actually stretched    
    try:
        outPart_COM = [statistics.mean(outPart_x), statistics.mean(outPart_y)]
        dist_x = [x - outPart_COM[0] for x in outPart_x]
        dist_y = [y - outPart_COM[1] for y in outPart_y]
        stdDev = [statistics.stdev(dist_x), statistics.stdev(dist_y)]
    except statistics.StatisticsError:
        continue
    
    stretchedPart_x = []
    stretchedPart_y = []
    for i in range(0,len(dist_x)):
        if dist_x[i] <= 1.75*stdDev[0] and dist_y[i] <= 1.75*stdDev[1]:
            stretchedPart_x.append(outPart_x[i])
            stretchedPart_y.append(outPart_y[i])
    
    if len(stretchedPart_x) >= 20:
        popt, pcov = curve_fit(func, stretchedPart_x, stretchedPart_y)
        C = popt
    else:
        continue
    
    plt.scatter(stretchedPart_x, stretchedPart_y, c='blue', marker=".")
    xdata = np.linspace(min(stretchedPart_x), max(stretchedPart_x), 100)
    # Calculate stretching length
    rawLength = integrate.quad(integrand, min(stretchedPart_x), max(stretchedPart_x))
    length = (str(round(rawLength[0], sigfigs=4))+' +/- '+str(round(rawLength[1], sigfigs=3))+' m')
    plt.plot(xdata, func(xdata, *popt), label=length, c='red')
    plt.legend(loc='best')
    print(length)

    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
    plt.title(str(curr_time)+' hrs')
    
    if (not os.path.exists('./Time_Evolution_Figs/'+directory+'SL')):
        os.mkdir('./Time_Evolution_Figs/'+directory+'SL')

    plt.savefig('./Time_Evolution_Figs/'+directory+'SL/'+str(file)[-9:-4]+'.png')

    plt.close()
