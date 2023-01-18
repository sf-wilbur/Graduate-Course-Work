
import numpy as np
x = [0, 0, 1, 2, 3, 4, 2, 0, -2, 0, 0]
y = []
for n,value in enumerate(x):
    print(n,value)
    y[n] = np.median(x[n-1], x[n], x[n+1])
    #
    # print(y)