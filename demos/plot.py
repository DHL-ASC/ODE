import matplotlib.pyplot as plt
import numpy as np

#write a function to read and plot the file test_ode.txt with matplotlib

# read the file test_ode.txt
data = np.loadtxt('test_ode_ie.txt')
# plot the data
plt.plot(data[:,0],data[:,1],label='U_0(t)')
plt.plot(data[:,0],data[:,2],label='U_C(t)')
plt.xlabel('t')
plt.ylabel('y(t), y\'(t)')
plt.legend()
plt.show()
