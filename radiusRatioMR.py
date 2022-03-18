import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np


def func(x, c0, c1, c2):
    return c0 + c1*x + c2*(x**2)


r = open('./MR-PR.txt')
lines = r.readlines()

gamma, p0, p1, p2 = [], [], [], []

for line in lines[1:]:
    values = line.split('\t')
    gamma.append(values[0])
    p0.append(values[1])
    p1.append(values[2])
    p2.append(values[3])

popt0, pcov0 = curve_fit(func, gamma, p0)
popt1, pcov1 = curve_fit(func, gamma, p1)
popt2, pcov2 = curve_fit(func, gamma, p2)

plt.scatter(gamma, p0)
gamma_test = np.linspace(0, 1)
plt.plot(gamma_test, func(gamma_test, *popt0), label=str(popt0))
plt.legend()
plt.xlabel('mass ratio')
plt.ylabel('p0')
plt.savefig('./p0MR.png')
plt.show()
plt.clf()

plt.scatter(gamma, p1)
gamma_test = np.linspace(0, 1)
plt.plot(gamma_test, func(gamma_test, *popt1), label=str(popt1))
plt.legend()
plt.xlabel('mass ratio')
plt.ylabel('p1')
plt.savefig('./p1MR.png')
plt.show()
plt.clf()

plt.scatter(gamma, p2)
gamma_test = np.linspace(0, 1)
plt.plot(gamma_test, func(gamma_test, *popt2), label=str(popt2))
plt.legend()
plt.xlabel('mass ratio')
plt.ylabel('p2')
plt.savefig('./p2MR.png')
plt.show()
plt.clf()
    

# files = sorted(glob.glob('./radiusRat/*'))

# for file in files:
#     r = open(file)
#     lines = r.readlines()
    
#     # impactor to total mass ratio
#     gamma = float(lines[0].split(":")[1].strip())
    
#     #impact angles and plume to impactor radius ratios
#     angle = []
#     radRat = []
    
#     for line in lines[2:]:
#         angle.append(float(line.split('\t')[0]))
#         radRat.append(float(line.split('\t')[1]))

# #    angle = [0, 30, 60, 90]
# #    rR = [1, 1, 0.87065987691, 0.2303460449]

#     plt.scatter(angle, radRat)
#     # curve fit wrt angle, then mass ratio
#     popt, pcov = curve_fit(P, angle, radRat)
    
# #    print(popt)
# #    print(popt2)
    
#     figname = str(file)[-9:-4]
#     angl_test = np.linspace(0, 100, 1000)
#     plt.plot(angl_test, P(angl_test, *popt), label=str(popt))
#     plt.legend()
#     plt.xlabel('Impact Angle (degrees)')
#     plt.ylabel('Radius Ratio')
#     plt.title('Simulations '+figname+' (γ = '+str(gamma)+')')
#     plt.savefig('./'+figname+'.png')
#     plt.show()



