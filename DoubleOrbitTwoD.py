#!/usr/bin/env python
# coding: utf-8

# In[4]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import numpy as np
from scipy import integrate 
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def animate(i):
    Orbit1.set_xdata(x1[:i])
    Orbit1.set_ydata(y1[:i])
    Orbit2.set_xdata(x2[:i])
    Orbit2.set_ydata(y2[:i])
    return (Orbit1,Orbit2)

G=1
m1=1
m2=1
Initial = [-1,0,0,0,2,0,1,-0.5] #[x1,dx1,y1,dy1,x2,dx2,y2,dy2]

abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 40
numpoints = 1000
Time = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

def ODE(F,t):
    x1,dx1,y1,dy1,x2,dx2,y2,dy2 = F
    r12 = [x2-x1,y2-y1]
    magr = np.sqrt((x2-x1)**2+(y2-y1)**2)
    return [dx1,m2*G*r12[0]/(magr**3),
            dy1,m2*G*r12[1]/(magr**3),
            dx2,-m1*G*r12[0]/(magr**3),
            dy2,-m1*G*r12[1]/(magr**3)]


#def TwoDIM(Initial):
    
    
result = integrate.odeint(ODE,Initial,Time, atol=abserr, rtol=relerr)
x1 = result[:,0]
y1 = result[:,2]
x2 = result[:,4]
y2 = result[:,6]
dx1 = result[:,1]
dy1 = result[:,3]
dx2 = result[:,5]
dy2 = result[:,7]
    
v1 = np.sqrt(dx1**2+dy1**2)
v2 = np.sqrt(dx2**2+dy2**2)
KE1 = 0.5*m1*v1**2
KE2 = 0.5*m2*v2**2
GPE = -G*m1*m2/(np.sqrt((x1-x2)**2+(y1-y2)**2))
AngMom1 = m1*(x1*dy1 - y1*dx1)
AngMom2 = m2*(x2*dy2 - y2*dx2)
    
    
f0,ax = plt.subplots(figsize = (10,6))
Orbit1 = plt.plot(x1[0],y1[0],color='red',label = '1')[0]
Orbit2 = plt.plot(x2[0],y2[0],color='blue',label = '2')[0]
ani = animation.FuncAnimation(fig=f0, func=animate, frames=numpoints, interval=20)
ax.set(xlim=[-3, 3], ylim=[-3, 3], xlabel='x', ylabel='y')
plt.legend()
"""
f1 = p.subplots(figsize = (10,6))
plt.scatter(Time,v1,s=4,color='red',label = '1')
plt.scatter(Time,v2,s=4,color = 'blue',label = '2')             
p.xlabel('Time')
p.ylabel('Velocity')
plt.legend()
    
f2 = plt.subplots(figsize = (10,6))
plt.plot(Time,KE1,color='red',label = 'KE1')
plt.plot(Time,KE2,color = 'blue',label = 'KE2') 
plt.plot(Time,GPE,color = 'black',label = 'GPE') 
plt.plot(Time,KE1+KE2+GPE,color = 'purple',label = 'Total') 
plt.xlabel('Time')
plt.ylabel('Energy')
plt.legend()
    
f3 = plt.subplots(figsize = (10,6))
plt.plot(Time,AngMom1,color='red',label = 'AngMom1')
plt.plot(Time,AngMom2,color = 'blue',label = 'AngMom2')     
plt.plot(Time,AngMom1+AngMom2,color = 'Purple',label = 'Total')  
plt.xlabel('Time')
plt.ylabel('Angular Momentum')
plt.legend()
"""
plt.show()
    
#TwoDIM(TwoDInitial)


# In[ ]:




