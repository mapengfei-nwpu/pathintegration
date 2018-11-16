import numpy as np
from scipy.interpolate import Rbf

f = open("161.01.txt")
x,y,z=[],[],[]
for line in f.readlines():
    data = line.split(" ")
    x.append(data[0])
    y.append(data[1])
    z.append(data[2])

f.close()
rbf = Rbf(x,y,z)
inteval = list(np.linspace(-5,5,51))
result = []
for xx in inteval:
    sum = 0.0
    for yy in inteval:
        sum += rbf([xx], [yy])[0]
    result.append(sum*0.2)

f = open('px.txt', "w") 
for i in range(len(inteval)):
        print("%.4f %.4f" % (inteval[i], result[i]), file = f)

f.close()