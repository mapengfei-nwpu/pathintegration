import numpy as np

f = open("149.23.txt")
xx,yy,errs=[],[],[]
for line in f.readlines():
    words=line.split(" ")
    x=float(words[0])
    y=float(words[1])
    appr = float(words[2])
    real = 0.0297*np.exp(-0.4*( -x*x + 0.1*x*x*x*x + + y*y ))
    err  = np.abs(appr-real)
    xx.append(x)
    yy.append(y)
    errs.append(err)

import matplotlib.pyplot as plt
plt.plot(xx,yy,errs)
plt.show()

file_name = 'errs.txt'
f = open(file_name, "w") 
for i in range(len(xx)):
    print("%.4f %.4f %.4f" % (xx[i],yy[i],errs[i]), file = f)
f.close()

err_norm1 = 0.0
for err in errs:
    if err>err_norm1:
        err_norm1=err

err_norm2 = 0.0
for err in errs:
    err_norm2 = err*err

err_norm2 = np.sqrt(err_norm2)
