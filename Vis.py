#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math, glob, os 
import matplotlib
import matplotlib.pyplot as plt

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

#Calculates the number of iron particles inside of the target
def numIron (target, x, y, z, iron, leng):
    ironCount = 0
    for i in range(0, leng):
        if (1 > ((target.center[0]-x[i])/(target.radii[0]))**2 + ((target.center[1] - y[i])/(target.radii[1]))**2 + 
            ((target.center[2] - z[i])/(target.radii[2]))**2 and iron[i]):
            ironCount+=1
    return(ironCount)


# In[2]:


directory = input('Enter directory name for data set: ')
while(True):
    files = sorted(glob.glob('./'+directory+'/*.dat'))
    if len(files) == 0: 
        directory = input('Non-existent or empty directory. Please try again: ')
    else:
        break    

sep = int(input('Separation index between target and impactor: '))

# In[3]:


trgt = target()
impctr = impactor()
time = []

files = [file for file in files if not(str(file).endswith("disk.dat") or str(file).endswith("diss.dat"))]
for file in files:
    # if str(file).endswith("disk.dat"):
    #     break
    with open(file) as r:
        lines = r.readlines()[2:]
    print(str(file))
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
            
    plt.scatter(x_t_mntl, y_t_mntl, c='orange', marker=".")
    plt.scatter(x_t_rn, y_t_rn, c='green', marker=".")
    
    plt.scatter(x_i_mntl, y_i_mntl, c='red', marker=".")
    plt.scatter(x_i_rn, y_i_rn, c='blue', marker=".")

    plt.xlabel('x position (m)')
    plt.ylabel('y position (m)')
    plt.title(str(curr_time)+' hrs')
    
    if (not os.path.exists('./Time_Evolution_Figs/'+directory)):
        os.mkdir('./Time_Evolution_Figs/'+directory)
    
    plt.savefig('./Time_Evolution_Figs/'+directory+'/'+str(file)[-9:-4]+'.png')

    plt.close()

    
