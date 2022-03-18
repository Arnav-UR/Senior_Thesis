#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math, glob 
import matplotlib
import matplotlib.pyplot as plt

class target:
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
time = []
ironAmt = []
ironRat = []

for file in files:
    if str(file).endswith("disk.dat"):
        break
    with open(file) as r:
        lines = r.readlines()[2:]
#         if str(file) == './sphg9__94\d0990.dat':
#             break
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
    time.append(curr_time)
    if str(file).endswith("1.dat"):
        trgt.radii = [abs(trgt.center[0]-max(x_t)), abs(trgt.center[1]-max(y_t)), abs(trgt.center[2]-max(z_t))]
    ironAmt.append(numIron(trgt, x_i, y_i, z_i, iron_i, len(x_i)))
    

log = open('./Accretion_Ratios/'+directory+'-accRatios.txt', 'w')

tIrn = 0
for n in iron_i:
    if n==True:
        tIrn += 1
for i in ironAmt:
    ironRat.append(i/tIrn)
plt.plot(time, ironRat)

log.write("Time(hr)\tMass Accretion Ratio\n")
for i in range(len(ironRat)):
    log.write(str(time[i])+'\t'+str(ironRat[i])+'\n')

plt.xlabel('Time (hrs)')
plt.ylabel('Mass Ratio of Accreted Iron to total Impactor Core Iron')
plt.title('Iron in Target vs. Time')
plt.savefig('./Mass_Accretion_Figs/'+directory+'.png')

plt.show()
plt.close()
log.close()
    
