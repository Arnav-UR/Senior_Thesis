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


directory_numbers = input('Enter numbers for data sets, separated by space: ')
sep = int(input('Separation index between target and impactor: '))
directories = []
for number in directory_numbers.split():
    if int(number) > 100:
        directories.append("sphg9_"+number)
    else:
        directories.append("sphg9__"+number)

# Make sure input is ok
l = True
while(l):
    for directory in directories:
        files = sorted(glob.glob('./'+directory+'/*.dat'))
        if len(files) == 0: 
            directories = input('Non-existent or empty directory '+directory)
            l = True
            break
        else:
            l = False
            
effective_radii = []
angles = [0, 30, 60, 90, 45]
impactor_mass = float(input("Enter the impactor's total mass in kg: "))

r_max = input('Enter stretching lengths in km, separated by space: ').split()
stretch_ratios = []

log = open('./Stretching_Ratios/'+directories[0]+'-stretchsRatios.txt', 'w')
    
for directory in directories:
    files = sorted(glob.glob('./'+directory+'/*.dat'))
    files = [file for file in files if not(str(file).endswith("disk.dat") or str(file).endswith("diss.dat"))]
    files = [files[0], files[-1]]
    
    trgt = target()
    
    for file in files:
        # if str(file).endswith("disk.dat") or files.index(file) != 0 or files.index(file) != (len(files)-1):
        #     continue
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

        if file == files[0]:
            trgt.radii = [abs(trgt.center[0]-max(x_t)), abs(trgt.center[1]-max(y_t)), abs(trgt.center[2]-max(z_t))]
            continue
        
        ironAccRat = numIron(trgt, x_i, y_i, z_i, iron_i, len(x_i))/len(iron_i)
        r_eff = ((3/4)*(0.0001)*(ironAccRat*impactor_mass)/math.pi)**(1/3)
        effective_radii.append(r_eff)

for i in range(len(r_max)):
    if effective_radii[i] != float(0):
        stretch_ratios.append(float(r_max[i])/effective_radii[i])
    else:
        stretch_ratios.append(float(0))
    
log.write("Impact Angle\tEffective Accretion Radius (km)\tMax Stretching (km)\tRatio rmax/r0\n")
for i in range(len(angles)):
    log.write(str(angles[i])+'\t'+str(effective_radii[i])+'\t'+r_max[i]+'\t'+str(stretch_ratios[i])+'\n')


plt.scatter(angles, stretch_ratios)
plt.xlabel('Impact Angle (degrees)')
plt.ylabel('rMax/r0')
plt.title('Stretching Ratio vs Impact Angle')
plt.savefig('./Stretching_Ratios/'+directories[0])
plt.show()

plt.close()
log.close()
    
