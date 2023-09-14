#!/usr/bin/env python
# coding: utf-8

# In[4]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import numpy as np
from scipy import integrate 
import pylab as p
import matplotlib.pyplot as plt
import matplotlib.animation as animation

G=1
m1 = 1
m2 = 1
m3 = 1

Initial = [
1,0,0,0,-1,0.2,
0,0.5,0,-0.2,1,0,
-1,0,1,0,0,0] 
#x1,dx1,y1,dy1,z1,dz1,
#x2,dx2,y2,dy2,z2,dz2,
#x3,dx3,y3,dy3,z3,dz3

abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 30
numpoints = 3000
Time = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]


def ODE(F,t):
    x1,dx1,y1,dy1,z1,dz1,x2,dx2,y2,dy2,z2,dz2,x3,dx3,y3,dy3,z3,dz3 = F
    r12 = [x2-x1,y2-y1,z2-z1]
    r13 = [x3-x1,y3-y1,z3-z1]
    r23 = [x3-x2,y3-y2,z3-z2]
    magr12 = np.sqrt((r12[0])**2+(r12[1])**2+(r12[2])**2)
    magr13 = np.sqrt((r13[0])**2+(r13[1])**2+(r13[2])**2)
    magr23 = np.sqrt((r23[0])**2+(r23[1])**2+(r23[2])**2)
    return [dx1,G*(m2*r12[0]/(magr12**3)+m3*r13[0]/(magr13**3)),
            dy1,G*(m2*r12[1]/(magr12**3)+m3*r13[1]/(magr13**3)),
            dz1,G*(m2*r12[2]/(magr12**3)+m3*r13[2]/(magr13**3)),
            dx2,G*(-m1*r12[0]/(magr12**3)+m3*r23[0]/(magr23**3)),
            dy2,G*(-m1*r12[1]/(magr12**3)+m3*r23[1]/(magr23**3)),
            dz2,G*(-m1*r12[2]/(magr12**3)+m3*r23[2]/(magr23**3)),
            dx3,G*(-m1*r13[0]/(magr13**3)-m2*r23[0]/(magr23**3)),
            dy3,G*(-m1*r13[1]/(magr13**3)-m2*r23[1]/(magr23**3)),
            dz3,G*(-m1*r13[2]/(magr13**3)-m2*r23[2]/(magr23**3))]


#def ThreeDIM(Initial):
result = integrate.odeint(ODE,Initial,Time, atol=abserr, rtol=relerr)
x1 = result[:,0]
y1 = result[:,2]
z1 = result[:,4]                     
x2 = result[:,6]
y2 = result[:,8]
z2 = result[:,10] 
x3 = result[:,12]
y3 = result[:,14] 
z3 = result[:,16] 
          
dx1 = result[:,1]
dy1 = result[:,3]
dz1 = result[:,5]
dx2 = result[:,7]
dy2 = result[:,9]
dz2 = result[:,11]
dx3 = result[:,13]
dy3 = result[:,15]
dz3 = result[:,17]
        
v1 = np.sqrt(dx1**2+dy1**2+dz1**2)
v2 = np.sqrt(dx2**2+dy2**2+dz2**2)
v3 = np.sqrt(dx3**2+dy3**2+dz3**2)
                     
ax = plt.figure().add_subplot(projection='3d')
ax.plot(x1,y1,z1,color='red')
ax.plot(x2,y2,z2,color='green')
ax.plot(x3,y3,z3,color='blue')
plt.show()


# In[ ]:




