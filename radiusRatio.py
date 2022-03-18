import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import math
import glob

# Functions whose parameters are to be minimized for the model 

# Legendre polynomial as a function of impact angle
def P(x, p0, p1, p2):
    t = np.cos(x*math.pi/180) # cosine of impact angle
    return p0 + p1*t + p2*0.5*(3*t**2 - 1)

# Polynomial as a function of mass ratio
def C(x, c0, c1, c2):
    c = c0 + c1*gamma + c2*(gamma**2)
    return c * P(x, *popt)

files = sorted(glob.glob('./radiusRat/*'))

for file in files:
    r = open(file)
    lines = r.readlines()
    
    # impactor to total mass ratio
    gamma = float(lines[0].split(":")[1].strip())
    
    #impact angles and plume to impactor radius ratios
    angle = []
    radRat = []
    
    for line in lines[2:]:
        angle.append(float(line.split('\t')[0]))
        radRat.append(float(line.split('\t')[1]))

#    angle = [0, 30, 60, 90]
#    rR = [1, 1, 0.87065987691, 0.2303460449]

    plt.scatter(angle, radRat)
    # curve fit wrt angle, then mass ratio
    popt, pcov = curve_fit(P, angle, radRat)
    popt2, pcov2 = curve_fit(C, angle, radRat)
    
#    print(popt)
#    print(popt2)

    w = open('./MR-PR.txt','a')
    w.write(str(gamma))
    for p in popt:
        w.write('\t'+str(p))
    w.write('\n')
    
    figname = str(file)[-9:-4]
    angl_test = np.linspace(0, 100, 1000)
    plt.plot(angl_test, C(angl_test, *popt2), label=str(popt)+'\n'+str(popt2))
    plt.legend()
    plt.xlabel('Impact Angle (degrees)')
    plt.ylabel('Radius Ratio')
    plt.title('Simulations '+figname+' (γ = '+str(gamma)+')')
    plt.savefig('./Fit_(angle_and_mass_ratio)/'+figname+'.png')
    plt.show()
    
    r.close()
    w.close()



