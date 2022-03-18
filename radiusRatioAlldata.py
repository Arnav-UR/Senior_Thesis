import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import math
import glob

# Functions whose parameters are to be minimized for the model 

def sigmoid(x, c0, c1, c2):
    return (c1/(1 + np.exp(x - c0))) + c2


def func(inp, p0, p1, p2, c0, c1, c2):
    x,g = inp
    t = np.cos(x*math.pi/180) # cosine of impact angle
    P = p0 + p1*t + p2*0.5*(3*t**2 - 1)
    return sigmoid(P, g*c0, c1, c2)

# Normalized to range from 0.0 to 1.0

def normFunc(inp, p0, p1, p2, c0, c1, c2):
    val = func(inp, p0, p1, p2, c0, c1, c2)
    val = [1.0 if v > 1.0 else v for v in val]
    val = [0.0 if v < 0.0 else v for v in val]
    return val

# impact angles, plume to impactor radius ratios, impactor to total mass ratio
angle = []
radRat = []
gamma = []

files = sorted(glob.glob('./radiusRat/*'))
    
for file in files:
    r = open(file)
    lines = r.readlines()
    
    # impactor to total mass ratio
    g = float(lines[0].split(":")[1].strip())
    

    for line in lines[2:]:
        angle.append(float(line.split('\t')[0]))
        radRat.append(float(line.split('\t')[1]))
        gamma.append(g)
        
    r.close()

print(gamma)
cmap=plt.cm.get_cmap('viridis')
plt.scatter(x=angle, y=radRat, c=gamma, cmap=cmap)
plt.colorbar()

popt, pcov = curve_fit(func, [angle, gamma], radRat)


angl_test = np.linspace(0, 90, 1000)
gamma_test1 = np.full(1000, 0.03)
gamma_test2 = np.full(1000, 0.1)
gamma_test3 = np.full(1000, 0.2)
gamma_test4 = np.full(1000, 0.3)
gamma_test5 = np.full(1000, 0.4)
gamma_test6 = np.full(1000, 0.5)

plt.plot(angl_test, normFunc([angl_test, gamma_test1], *popt), label='γ=0.03', color=cmap(2*0.03))
plt.plot(angl_test, normFunc([angl_test, gamma_test2], *popt), label='γ=0.1', color=cmap(2*0.1))
plt.plot(angl_test, normFunc([angl_test, gamma_test4], *popt), label='γ=0.3', color=cmap(2*0.3))
plt.plot(angl_test, normFunc([angl_test, gamma_test6], *popt), label='γ=0.5', color=cmap(2*0.5))

plt.legend()
plt.xlabel('Impact Angle (degrees)')
plt.ylabel('Radius Ratio')
plt.title('Radius Ratio vs Impact Angle at Various γ')
plt.savefig('./allDataFitTestsigGrad3.png')
plt.show()





