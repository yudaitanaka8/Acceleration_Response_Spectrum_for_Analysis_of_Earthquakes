
# coding: utf-8

# In[ ]:


import numpy as np
import math
import matplotlib.pyplot as plt

T = 0.02
Max = 5.0
Num = int(Max/T)
dt = 0.01
NS = [] 
Spec = np.array([[[0.0 for i3 in range(Num)] for j in range(2)] for k in range(3)])
h_value_candidate = [0.01, 0.05, 0.1]

data = "" #Accelerogram of an earthquake which, for instance, can be downloaded from https://www.data.jma.go.jp/svd/eqev/data/kyoshin/jishin/
with open(data, "r", encoding="utf_8") as fileobj:
    for i, line in enumerate(fileobj):
        datalist = line.split(",")
        if i > 0:
            NS.append(float(datalist[0]))
NS = np.array(NS)
NS = NS*0.01

plt.figure(figsize=(8, 6,), dpi=200)

for i in range(3):
    h = h_value_candidate[i]
    for j in range(Num):
        Spec[i][0][j] = float(T*(j+1))
        w = 2*math.pi/(T*(j+1))
        a = 0
        v = 0
        x = 0
        k = 0
        for k in range(len(NS)-1):
            formula = -1*(NS[k+1] + 2*h*w*(v+a*dt/2) + w**2*(x+v*dt+a*dt**2/3))/(1 + h*w*dt + (w*dt)**2/6) #線形加速度法の式
            x = x + v*dt + (2*a+formula)*dt**2/6
            v = v + (a + formula) * dt / 2
            a = formula
            if abs(a + NS[k+1]) > Spec[i][1][j]:
                Spec[i][1][j] = a + NS[k+1]
            k += 1
    plt.plot(Spec[i, 0, :], Spec[i, 1, :], label =  "h=" + str(h))

plt.title("Acceleration Response Spectrum")
plt.xlabel("Time")
plt.ylabel("Acceleration")
plt.legend()
plt.savefig("filename.png")
plt.show()

