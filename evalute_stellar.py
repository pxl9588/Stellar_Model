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
data1 = collection[1]
data0_transpose = data0.transpose()
data1_transpose = data1.transpose()
x0 = data0_transpose[0][0:]
y0 = data0_transpose[-1][0:]
x1 = data1_transpose[0][0:]
y1 = data1_transpose[-1][0:]
fig,ax = plt.subplots()

line1, = ax.plot(x0/x0[0], y0, label='2M')
line2, = ax.plot(x1/x1[0], y1, label='1M')
ax.legend()
plt.ylim(top = 6, bottom = 0)
plt.show()
