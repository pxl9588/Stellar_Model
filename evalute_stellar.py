import numpy as np
import matplotlib.pyplot as plt

GAMMA = 5/3
stants = {"sigma": 5.67051e-5, "c": 2.99792458e+10, "a": 7.56591e-15, "G": 6.67259e-8, "k_B": 1.38065e-16,
                 "m_H": 1.673534e-24, "pi": np.pi, "gamma":  GAMMA, "gamrat":  GAMMA / (GAMMA - 1.0),
                 "kPad":  1.0, "g_ff": 1.0, "Rsun": 6.9599e+10, "Msun": 1.989e+33, "Lsun": 3.826e+33}
g = open("filenames.txt","r")
collection = []
for i in g.read().split("\n")[:-1]:
    f = open(i, 'r')

    temp_a = f.readline()[:-1].split(",")
    labels = temp_a
    temp_data = np.array([[0.0]*9])
    while len(temp_a) > 1:
        a = np.array([])
        temp_a = f.readline()[:-1].split(",")
        for i in temp_a:
            try:
                a = np.append(a, float(i))
            except:
                break
        if len(a) > 2:
            temp_data = np.insert(temp_data, 1, a, 0)
    collection.append(np.array(temp_data)[1:])
data0 = collection[0]
#data1 = collection[1]
data0_transpose = data0.transpose()
#data1_transpose = data1.transpose()
radius0 = data0_transpose[0][0:]
print("Radius: " + str(radius0[0]))
pres0 = data0_transpose[1][0:]
print("Pressure: " + str(pres0[0]))
temp0 = data0_transpose[2][0:]
print("Temperature: " + str(temp0[0]))
rho0 = data0_transpose[3][0:]
print("Rho: " + str(rho0[0]))
kappa0 = data0_transpose[4][0:]
print("Kappa: " + str(kappa0[0]))
eps0 = data0_transpose[5][0:]
qm0 = data0_transpose[6][0:]
lum0 = 4 * np.pi * radius0**2 * rho0 * eps0
print("Luminosity: " + str(lum0[0]))
dpt0 = data0_transpose[8][0:]
print("dlnP/lnT: " + str(dpt0[0]))
#x1 = data1_transpose[0][0:]
#y1 = data1_transpose[-1][0:]
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)
#fig, ax = plt.subplot()
ax1.plot(radius0, temp0)
ax1.set_ylabel("T(r)")
ax2.plot(radius0, pres0)
ax2.set_ylabel("P(r)")
ax3.plot(radius0, rho0)
ax3.set_ylabel("Rho(r)")
ax4.plot(radius0, kappa0)
ax4.set_ylabel("K(r)")
ax5.plot(radius0, lum0)
ax5.set_ylabel("L(r)")
ax6.plot(radius0, dpt0)
ax6.set_ylabel("dlnP/lnT(r)")
ax6.set_ylim(top=6, bottom=1)
ax1.set_xlabel("r/R0")
ax2.set_xlabel("r/R0")
ax3.set_xlabel("r/R0")
ax4.set_xlabel("r/R0")
ax5.set_xlabel("r/R0")
ax6.set_xlabel("r/R0")
#ax2.plot(x1/x1[0], y1)
#plt.ylim(top=6, bottom=0)
plt.xlabel('r/Ro')
plt.show()
