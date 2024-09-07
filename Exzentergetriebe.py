import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Circle, FancyArrow,Arc,Polygon
from matplotlib.transforms import Affine2D
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider
import numpy as np
from math import pi,sin,cos,asin,acos,atan,tan,atan2,hypot
from eccentric_cycloid_gear import EccentricCycloidGear
from HideOnClick2 import hideOnClick
from Utilities.Helper import calculate_arcPolygon
from matplotlib.animation import FuncAnimation

phi0=pi/2
z1=32
z2=30
alpha=20*pi/180
e=(z2-z1)/2


fig,ax=plt.subplots()
ax.set_aspect("equal")
ax.grid(True)
M=np.array((cos(phi0),sin(phi0)))*e

gear1=Circle(M,radius=z1/2,label="gear1", linestyle='-.', color='g', fill=False)
gear2=Circle((0,0),radius=z2/2,label="gear1",color="green",fill=False)
ax.add_patch(gear1)
ax.add_patch(gear2)
head=np.array((cos(phi0),sin(phi0)))*(e+z1/2)
angle=0
scale=5
dF=np.array((sin(phi0+angle+alpha),cos(phi0+angle+alpha)))*scale
kurbel=ax.arrow(x=head[0],y=head[1],dx=dF[0],dy=dF[1])

exzenter,=ax.plot((0,M[0]),(0,M[1]),color="black")

ax_slider = plt.axes([0.2, 0.05, 0.3, 0.02]) # [left, bottom, width, height]
#slider = Slider(ax_slider, 'Rotation', (0)*180/pi, (360+phi0)*180/pi, valinit=phi0,valstep=0.5)
slider = Slider(ax_slider, 'Rotation', (0), 2*pi, valinit=phi0,valstep=0.5)
def update_plot(val):
    angle=val
    M=np.array((cos(angle),sin(angle)))*e
    head=np.array((cos(phi0+angle),sin(phi0+angle)))*(e+z1/2)
    dF=np.array((sin(phi0+angle+alpha),cos(phi0+angle+alpha)))*scale
    kurbel.set_data(x=head[0],y=head[1],dx=dF[0],dy=dF[1])
    exzenter.set_data([0,M[0]],[0,M[1]])
    gear1.center=M

slider.on_changed(update_plot)
ax.set_xlim(-e-z1/1.5,e+z1/1.5)
ax.set_ylim(-e-z1/1.5,e+z1/1.5)
ax.legend()
plt.show()