### This file calculates the stretching ratio for all simulations and all impacts ###
import numpy as np
import math, glob 
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class target:
    radii = []
    center = []

def func(inp, p0, p1, p2, c0, c1, c2):
    x,g = inp
    t = np.cos(x*math.pi/180) # cosine of impact angle
    P = p0 + p1*t + p2*0.5*(3*t**2 - 1)
    G = c0 + c1*g + c2*0.5*(3*g**2 - 1)
    return G*P


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

with open('./Stretching_Data.txt') as rS:
    dataLines = rS.readlines()[1:]

#NB: "2" denotes values associated with a secondary impact

all_angles = []
all_angles2 = []
all_ratios = []
all_ratios2 =[]
all_error = []
all_error2 = []
gamma = []
gamma2 =[]

for dataLine in dataLines:
    directory_numbers = dataLine.split('\t')[0]
    sep = int(dataLine.split('\t')[3])
    
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
    effective_radii2 = []
    
    angles = [0, 30, 60, 90, 45]
    if dataLine.split('\t')[4] != 'False':
        for angle in dataLine.split('\t')[4].split(' '):
            angles.insert(angles.index(int(angle))+1, int(angle))
    angles2 = []
    
    #  impactor's total mass in kg      
    impactor_mass = sci2Float(dataLine.split('\t')[1])

    # Maximum stretching length
    r_max = dataLine.split('\t')[2].split(' ')
    r_max2 = []
    stretch_ratios = []
    stretch_ratios2 = []
    sr_error = []
    sr_error2 = []
    
    sl_error = dataLine.split('\t')[-3].split(' ')
    sl_error2 = []
    
    reff_error = dataLine.split('\t')[-1].split(' ')
    reff_error2 = []
    
    log = open('./Stretching_Ratios/'+directories[0]+'-stretchsRatios.txt', 'w')
        
    for directory in directories:
        files = sorted(glob.glob('./'+directory+'/*.dat'))
        files = [file for file in files if not(str(file).endswith("disk.dat") or str(file).endswith("diss.dat"))]
        
        if angles.count(60) > 1 and directories.index(directory) == 2:
            impact_cutoff = dataLine.split('\t')[5].split(' ')[0].strip()
            files = [files[0], './'+directory+'/'+impact_cutoff+'.dat', files[-1]]
            angles2.append(angles[3])
            del angles[3]
            r_max2.append(r_max[3])
            del r_max[3]
            sl_error2.append(sl_error[3])
            del sl_error[3]
            reff_error2.append(sl_error[3])
            del reff_error[3]
            
        elif angles.count(90) > 1 and directories.index(directory) == 3:
            impact_cutoff = dataLine.split('\t')[5].split(' ')[1].strip()
            files = [files[0], './'+directory+'/'+impact_cutoff+'.dat', files[-1]]
            angles2.append(angles[4])
            del angles[4]
            r_max2.append(r_max[4])
            del r_max[4]
            sl_error2.append(sl_error[4])
            del sl_error[4]
            reff_error2.append(sl_error[4])
            del reff_error[4]
        else:
            files = [files[0], files[-1]]
        
        trgt = target()
        
        for file in files:
            with open(file) as r:
                lines = r.readlines()[2:]
        
            #print(str(file))
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
            
            # calculate effective accretion radius
            ironAccRat = numIron(trgt, x_i, y_i, z_i, iron_i, len(x_i))/sum(i for i in iron_i)
            
            if files.index(file) > 1:
                iron_mass = iron_i.count(True)*impactor_mass
                prev_ironAccRat = (math.pi/(0.000075*iron_mass))*(effective_radii[-1])**3
                ironAccRat = ironAccRat - prev_ironAccRat
                r_eff = ((3/4)*(0.0000001)*(ironAccRat*iron_mass)/math.pi)**(1/3)
                effective_radii2.append(r_eff)
                continue
            
            iron_mass = iron_i.count(True)*impactor_mass
            r_eff = ((3/4)*(0.0000001)*(ironAccRat*iron_mass)/math.pi)**(1/3)
            effective_radii.append(r_eff)
    
    for i in range(len(r_max)):
        try:
            rMax = float(r_max[i])
            stretch_ratio = rMax/effective_radii[i]
            error = stretch_ratio*math.sqrt(((float(sl_error[i].strip())/rMax)**2)+((float(reff_error[i].strip())/effective_radii[i])**2))
            stretch_ratios.append(stretch_ratio)
            sr_error.append(error)
        except ZeroDivisionError:
            stretch_ratios.append(float(0))
            sr_error.append(float(0))
    
    for i in range(len(r_max2)):
        try:
            rMax = float(r_max2[i])
            stretch_ratio = rMax/effective_radii[i]
            error = stretch_ratio*math.sqrt(((float(sl_error2[i].strip())/rMax)**2)+((float(reff_error2[i].strip())/effective_radii2[i])**2))
            stretch_ratios2.append(stretch_ratio)
            sr_error2.append(error)
        except ZeroDivisionError:
            stretch_ratios2.append(float(0))
            sr_error2.append(float(0))
        
    log.write("Impact Angle\tEffective Accretion Radius (m)\tMax Stretching (m)\tRatio rmax/r0\n")
    for i in range(len(angles)):
        log.write(str(angles[i])+'\t'+str(effective_radii[i])+'\t'+r_max[i]+'\t'+str(stretch_ratios[i])+'\n')
    for i in range(len(angles2)):
        log.write(str(angles2[i])+'\t'+str(effective_radii2[i])+'\t'+str(r_max2[i])+'\t'+str(stretch_ratios2[i])+'\n')
    
    print(angles)
    print(stretch_ratios)
    plt.scatter(angles, stretch_ratios)
    plt.xlabel('Impact Angle (degrees)')
    plt.ylabel('rMax/r0')
    plt.title('Stretching Ratio vs Impact Angle')
    plt.savefig('./Stretching_Ratios/'+directories[0])
    plt.show()
    plt.close()
    
    plt.scatter(angles, stretch_ratios)
    plt.xlabel('Impact Angle (degrees)')
    plt.ylabel('rMax/rEff')
    plt.title('Stretching Ratio vs Impact Angle (2nd Impact)')
    plt.savefig('./Stretching_Ratios/'+directories[0]+'2')
    plt.show()
    
    log.close()
    
    # Update total data
    for i in range(len(angles)):
        gamma.append(float(dataLine.split('\t')[-4]))
        
    for i in range(len(angles2)):
        gamma2.append(float(dataLine.split('\t')[-4]))
        
    all_angles.extend(angles)
    all_ratios.extend(stretch_ratios)
    all_error.extend(sr_error)
    
    all_angles2.extend(angles2)
    all_ratios2.extend(stretch_ratios2)
    all_error2.extend(sr_error2)
    
for i in range(len(all_angles)):
    print(str(all_angles[i])+' '+str(all_error[i]))

cmap=plt.cm.get_cmap('viridis')
plt.scatter(x=all_angles, y=all_ratios, c=gamma, cmap=cmap)

plt.errorbar(x=all_angles, y=all_ratios, yerr=all_error, linestyle="None", color='red')
plt.colorbar()

popt, pcov = curve_fit(func, [all_angles, gamma], all_ratios)
print(*popt)

angl_test = np.linspace(0, 90, 1000)
gamma_test1 = np.full(1000, 0.03)
gamma_test2 = np.full(1000, 0.1)
gamma_test3 = np.full(1000, 0.2)
gamma_test4 = np.full(1000, 0.3)
gamma_test5 = np.full(1000, 0.4)
gamma_test6 = np.full(1000, 0.5)

plt.plot(angl_test, func([angl_test, gamma_test1], *popt), label='γ=0.03', color=cmap(2*0.03))
plt.plot(angl_test, func([angl_test, gamma_test2], *popt), label='γ=0.1', color=cmap(2*0.1))
plt.plot(angl_test, func([angl_test, gamma_test4], *popt), label='γ=0.3', color=cmap(2*0.3))
plt.plot(angl_test, func([angl_test, gamma_test6], *popt), label='γ=0.5', color=cmap(2*0.5))

plt.legend()
plt.xlabel('Impact Angle (degrees)')
plt.ylabel('rMax/rEff')
plt.title('Stretching Ratio vs Impact Angle')
plt.savefig('./Stretching_Ratios/total.png')
plt.show()
plt.close()

cmap=plt.cm.get_cmap('viridis')
plt.scatter(x=all_angles2, y=all_ratios2, c=gamma2, cmap=cmap)

plt.errorbar(x=all_angles2, y=all_ratios2, yerr=all_error2, linestyle="None", color='red')
plt.colorbar()


plt.xlabel('Impact Angle (degrees)')
plt.ylabel('rMax/rEff')
plt.title('Stretching Ratio vs Impact Angle (2nd Impact)')
plt.savefig('./Stretching_Ratios/total2.png')
plt.show()
plt.close()
