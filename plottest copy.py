import matplotlib.pyplot as plt
from eccentric_cycloid_gear import EccentricCycloidGear
from math import pi,sin,cos,atan
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

def xi(t,lam=1):
    return atan(lam*sin(t)/(1-lam*cos(t)))
#A
x=np.linspace(-pi/(gear_pair.z2),pi/gear_pair.z2,num=1000)
y1=np.array([xi(a) for a in x])
ax.plot(x*180/pi,y1*180/pi,label="xi(1)",marker='',linestyle='-')
y1=np.array([xi(a,0.97) for a in x])
ax.plot(x*180/pi,y1*180/pi,label="xi(0.97)",marker='',linestyle='-')
y1=np.array([xi(a,0.9) for a in x])
ax.plot(x*180/pi,y1*180/pi,label="xi(0.9)",marker='',linestyle='-')


ax.grid()
ax.legend()
#ax.set_aspect('equal')
hideOnClick(fig,ax)
plt.show()