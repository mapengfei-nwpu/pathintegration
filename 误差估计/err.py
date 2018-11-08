import numpy as np
f = open("149.23.txt")
errs=[]
for line in f.readlines():
    words=line.split(" ")
    x=float(words[0])
    y=float(words[1])
    appr = float(words[2])
    real = 0.0297*np.exp(-0.4*( -x*x + 0.1*x*x*x*x + + y*y ))
    err  = np.abs(appr-real)
    errs.append(err)

err_norm1 = 0.0
for err in errs:
    if err>err_norm1:
        err_norm1=err

err_norm2 = 0.0
for err in errs:
    err_norm2 = err*err

err_norm2 = np.sqrt(err_norm2)
