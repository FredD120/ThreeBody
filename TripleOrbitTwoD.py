#!/usr/bin/env python
# coding: utf-8

# In[5]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import numpy as np
from scipy import integrate 
import pylab as p
import matplotlib.pyplot as plt
import matplotlib.animation as animation

G=1
m1 = 1
m2 = 1
m3 = 4

Initial = [-1,0,0,0.5,1,0,0,-0.5,0,0,0,0] #[x1,dx1,y1,dy1,x2,dx2,y2,dy2,x3,dx3,y3,dy3]

abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 15
numpoints = 3000
Time = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]


def animate(i):
    Orbit1.set_xdata(x1[:i])
    Orbit1.set_ydata(y1[:i])
    Orbit2.set_xdata(x2[:i])
    Orbit2.set_ydata(y2[:i])
    Orbit3.set_xdata(x3[:i])
    Orbit3.set_ydata(y3[:i])
    return (Orbit1,Orbit2,Orbit3)


def ODE(F,t):
    x1,dx1,y1,dy1,x2,dx2,y2,dy2,x3,dx3,y3,dy3 = F
    r12 = [x2-x1,y2-y1]
    r13 = [x3-x1,y3-y1]
    r23 = [x3-x2,y3-y2]
    magr12 = np.sqrt((r12[0])**2+(r12[1])**2)
    magr13 = np.sqrt((r13[0])**2+(r13[1])**2)
    magr23 = np.sqrt((r23[0])**2+(r23[1])**2)
    return [dx1,G*(m2*r12[0]/(magr12**3)+m3*r13[0]/(magr13**3)),
            dy1,G*(m2*r12[1]/(magr12**3)+m3*r13[1]/(magr13**3)),
            dx2,G*(-m1*r12[0]/(magr12**3)+m3*r23[0]/(magr23**3)),
            dy2,G*(-m1*r12[1]/(magr12**3)+m3*r23[1]/(magr23**3)),
            dx3,G*(-m1*r13[0]/(magr13**3)-m2*r23[0]/(magr23**3)),
            dy3,G*(-m1*r13[1]/(magr13**3)-m2*r23[1]/(magr23**3))]


#def ThreeDIM(Initial):
result = integrate.odeint(ODE,Initial,Time, atol=abserr, rtol=relerr)
x1 = result[:,0]
y1 = result[:,2]
x2 = result[:,4]
y2 = result[:,6]
x3 = result[:,8]
y3 = result[:,10]               
          
dx1 = result[:,1]
dy1 = result[:,3]
dx2 = result[:,5]
dy2 = result[:,7]
dx3 = result[:,9]
dy3 = result[:,11]
        
v1 = np.sqrt(dx1**2+dy1**2)
v2 = np.sqrt(dx2**2+dy2**2)
v3 = np.sqrt(dx3**2+dy3**2)
    
KE1 = 0.5*m1*v1**2
KE2 = 0.5*m2*v2**2
KE3 = 0.5*m3*v3**2
GPE12 = -G*m1*m2/(np.sqrt((x1-x2)**2+(y1-y2)**2)) 
GPE13 = -G*m1*m3/(np.sqrt((x1-x3)**2+(y1-y3)**2)) 
GPE32 = -G*m3*m2/(np.sqrt((x3-x2)**2+(y3-y2)**2))
AngMom1 = m1*(x1*dy1 - y1*dx1)
AngMom2 = m2*(x2*dy2 - y2*dx2)
AngMom3 = m3*(x3*dy3 - y3*dx3)
    
f0,ax = plt.subplots(figsize = (6,6))
Orbit1=plt.plot(x1[0],y1[0],color='red',label = '1')[0]
Orbit2=plt.plot(x2[0],y2[0],color ='green',label = '2')[0]
Orbit3=plt.plot(x3[0],y3[0],color ='blue',label = '3')[0]            
ani = animation.FuncAnimation(fig=f0, func=animate, frames=numpoints, interval=1)
ax.set(xlim=[-3, 3], ylim=[-3, 3], xlabel='x', ylabel='y')
plt.legend()
"""
f1 = p.subplots(figsize = (10,6))
plt.scatter(Time,v1,s=4,color='red',label = '1')
plt.scatter(Time,v2,s=4,color = 'green',label = '2')
plt.scatter(Time,v3,s=4,color = 'blue',label = '3')              
p.xlabel('Time')
p.ylabel('Velocity')
plt.legend()

f2 = p.subplots(figsize = (10,6))
plt.plot(Time,KE1,color='red',label = 'KE1')
plt.plot(Time,KE2,color = 'green',label = 'KE2') 
plt.plot(Time,KE3,color = 'blue',label = 'KE3') 
plt.plot(Time,GPE12,color = 'yellow',label = 'GPE12') 
plt.plot(Time,GPE13,color = 'magenta',label = 'GPE13') 
plt.plot(Time,GPE32,color = 'cyan',label = 'GPE32') 
#plt.plot(Time,KE1+KE2+KE3+GPE12+GPE13+GPE32,color = 'black',label = 'Total') 
p.xlabel('Time')
p.ylabel('Energy')
plt.legend()
    
f3 = p.subplots(figsize = (10,6))
plt.plot(Time,AngMom1,color='red',label = 'AngMom1')
plt.plot(Time,AngMom2,color = 'green',label = 'AngMom2')  
plt.plot(Time,AngMom3,color = 'blue',label = 'AngMom3') 
#plt.plot(Time,AngMom1+AngMom2+AngMom3,color = 'Purple',label = 'Total')  
p.xlabel('Time')
p.ylabel('Angular Momentum')
plt.legend()
"""
plt.show()
    
#ThreeDIM(ThreeDInitial)


# In[ ]:




