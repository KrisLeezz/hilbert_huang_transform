import math
from matplotlib import pyplot as plt
import numpy as np
import csv
T=np.arange(0,10,0.001)
f1=5
f2=20
f3=50
S1=np.sin(2*f1*T*math.pi)#传入的是弧度
S2=np.sin(2*f2*T*math.pi)
S3=np.sin(2*f3*T*math.pi)
S4=S1+S2
S5=S2+S3
data=(np.vstack((S4,S5))).T
csvfile=open('sample.csv','w',newline ='')
csvwrite=csv.writer(csvfile)
csvwrite.writerows(data)
plt.figure(12)
plt.subplot(4,1,1)
plt.plot(T,S4)
plt.subplot(4,1,2)
plt.plot(T,S5)
plt.show()
csvfile.close()

