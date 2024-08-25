import matplotlib.pyplot as plt
from eccentric_cycloid_gear import EccentricCycloidGear
from math import pi,sin,cos
import numpy as np
from matplotlib.patches import Circle
from HideOnClick2 import hideOnClick

gear_pair = EccentricCycloidGear(
        z1=3,
        z2=6,
        a=100,
        rA_star=1.0,
        st_star=1.0,
        phi_j1=0.0,
        phi_As=60.0*pi/180,
        phi_Ae=170*pi/180,
        lambda_=0.97
    )


fig,ax=plt.subplots()

#A
t_range_circle=np.linspace(pi/(gear_pair.z2*4),3/4*pi/gear_pair.z2,num=1000)
t_range_contact=np.linspace(-0,-2*pi/gear_pair.z2,num=1000)

x,y=np.array((-np.sin(t_range_circle),np.cos(t_range_circle)))*gear_pair.ra2
ax.plot(x,y,label="Kreis",marker='',linestyle='-')

ra2=Circle((0,0),gear_pair.ra2,linestyle='-.',color='b',label=f"ra2",fill=False)
ax.add_patch(ra2)

x,y=[],[]

for a in t_range_contact:
    xc,yc=gear_pair.p_pOA(a)
    x.append(xc)
    y.append(yc)

A=gear_pair.A()
ax.plot(*A,label="A",marker='o')
ax.plot(x,y,label="Contact",linestyle='-.')


#E
t_range_circle=np.linspace(0,1*pi/gear_pair.z2,num=1000)
t_range_contact=np.linspace(-0,-pi/(gear_pair.z2*2),num=1000)
x,y=np.array((np.sin(t_range_circle),-np.cos(t_range_circle)))*gear_pair.ra1+np.array([[0],[gear_pair.a]])
ax.plot(x,y,label="Kreis ra1",marker='',linestyle='-')

x,y=[],[]

for a in t_range_contact:
    xc,yc=gear_pair.p_pOA(a)
    x.append(xc)
    y.append(yc)

E=gear_pair.E(steps=100)
ax.plot(*E,label="E",marker='o')
ax.plot(x,y,label="Contact2")



ax.grid()
ax.legend()
ax.set_aspect('equal')
hideOnClick(fig,ax)
plt.show()