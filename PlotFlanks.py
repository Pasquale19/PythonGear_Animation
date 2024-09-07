import matplotlib.pyplot as plt
from matplotlib.widgets import Slider as Slider
from matplotlib.gridspec import GridSpec
import matplotlib.patches as patches
from eccentric_cycloid_gear import EccentricCycloidGear
from Utilities.Converter import calculateFlankParameter,convert_to_rack
from Utilities.GearRack import GearRack2D,KreisbogenRack
from HideOnClick2 import hideOnClick
from math import pi,cos,sin,tan,atan,acos,asin,tau,atan2
import numpy as np


def plotFlanks(xg,yg,isRack:bool=False):
    # Create figure
    fig = plt.figure(figsize=(10, 8))
    arrowprops=dict(head_length=1.5, head_width=1.2, lw=1,length_includes_head=True)
    plt.rcParams['lines.markersize'] = 1
    # Create GridSpec
    gs = GridSpec(2, 2, width_ratios=[3, 1])

    # Create subplots
    ax1 = fig.add_subplot(gs[:, 0])  # Left plot, spans both rows
    ax2 = fig.add_subplot(gs[0, 1])  # Top right plot
    ax3 = fig.add_subplot(gs[1, 1])  # Bottom right plot

        
    if isRack:
        rack_x,rack_y=xg,yg
    else:
        rack_x,rack_y=convert_to_rack(xg,yg,33.333,center_x=0,center_y=100)
        c=patches.Circle(xy=(0,0),radius=33.333,label="r1",fill=False,linestyle=":")
        ax1.add_patch(c)
        ax1.scatter(xg, yg-100, c=range(len(xg)), cmap='viridis',label="gear")

    #num_points2=10000
    # xrack_interp=np.linspace(rack_x[0],rack_x[-1],num_points2,endpoint=False)
    # yrack_interp=np.interp(xrack_interp,xp=rack_x,fp=rack_y,right=rack_y[-1],left=rack_y[0])
    #arr=calculateFlankParameter(xrack_interp,yrack_interp,radius=gear.r1)
    arr=calculateFlankParameter(rack_x,rack_y,radius=gear.r1)
    scatter = ax1.scatter(rack_x,rack_y, c=range(len(rack_x)), cmap='viridis',label="rack")

    #scatter = ax1.scatter(xrack_interp, yrack_interp, c=range(len(xrack_interp)), cmap='viridis',label="rack")

    alpha=arr[1,:]
    alpha_deg=alpha*180/pi
    n=arr[2,:]
    x=arr[0,:]


    ax1.set_aspect("equal")
    ax1.set_xlim(min(rack_x)*1.1,max(rack_x)*1.1)
    ax1.set_ylim(min(rack_y*2),max(rack_y)*2)
    ax1.hlines(0,xmin=-30,xmax=30)
    ax1.grid()

    plt.colorbar(scatter,label='index')
    #ax1.plot(xr,yr,label=r"$rack$")

    #slider
    ax_slider = plt.axes([0.2, 0.05, 0.3, 0.02]) # [left, bottom, width, height]
    slider = Slider(ax_slider, 'Rotation', 0, len(x)-1, valinit=0,valstep=1)


    
    #ax1.plot(xg,yg-self.a,label=r"$gear$")

    ax3.scatter(x, n, c=range(len(n)), cmap='viridis',label="n")
    ax3.set_ylabel(r"normal $n(t)$")
    #ax2.plot(x,n,label="n")
    ax2.scatter(x, alpha*180/pi, c=range(len(x)), cmap='viridis',label=r"$\alpha$")
    ax2.set_ylabel(r"pressure angle $\alpha$ in °")
    #ax2.plot(x,alpha*180/pi,label=)


    #interactive elements
    i0=int(slider.val)
    i02=i0+1
    index_head=num_points//4,num_points//4+2
    index_A=index_head[1]+num_points//4

    description=ax1.text(0, 1, f'index of A={index_A}',
        horizontalalignment='left',
        verticalalignment='top',
        transform = ax1.transAxes)
    title_ax1=ax1.set_title(f"Gear Profile {i0}")
    title_ax2=ax2.set_title(f"pressure angle {alpha[i0]*180/pi}°")
    title_ax3=ax3.set_title(f"normal {n[i0]}")


    p_foot=np.array((x[i0],0))
    p_head=np.array((rack_x[i02],rack_y[i02]))
    d=p_head-p_foot
    n_arrow=ax1.arrow(x=p_foot[0],y=p_foot[1],dx=d[0],dy=d[1],**arrowprops,clip_on=True)
    n_text=ax1.text(x=p_foot[0],y=p_foot[1],s="n",ha="center",va="bottom")
    p2=np.array((rack_x[i02],rack_y[i02]))
    d2=p_head-p2
    tangent_arrow=ax1.arrow(x=p2[0],y=p2[1],dx=d2[0],dy=d2[1],**arrowprops,clip_on=True)
    m1,=ax1.plot(*p_foot,marker='o',markersize=4)
    m2,=ax2.plot(x[i0],alpha_deg[i0],marker='o',markersize=4)
    m3,=ax3.plot(x[i0],n[i0],marker='o',markersize=4)

    def update_plot(i):
        i=int(i)
        xi=x[i]
        n_i=n[i]
        alpha_i=alpha_deg[i]
        p_foot=np.array((xi,0))
        p_head=np.array((rack_x[i],rack_y[i]))
        d=p_head-p_foot
        p2=np.array((rack_x[i-1],rack_y[i-1]))
        d2=p_head-p2
        
        n_arrow.set_data(x=p_foot[0],y=p_foot[1],dx=d[0],dy=d[1])
        tangent_arrow.set_data(x=p2[0],y=p2[1],dx=d2[0],dy=d2[1])
        title_ax1.set_text(f"Gear Profile {i}")

        title_ax2.set_text(f"pressure angle {alpha_i:.1f}°")
        title_ax3.set_text(f"normal {n_i:.2f}")
        n_text.set_position((p_foot[0], p_foot[1]))
        m1.set_data(*p_foot)
        m2.set_data([xi],[alpha_i])
        m3.set_data([xi],[n_i])




    ax1.legend()
    ax2.legend()
    hideOnClick(fig,ax1)
    slider.on_changed(update_plot)
    plt.show()


if __name__ =="__main__":
    gear = EccentricCycloidGear(
        z1=3,
        z2=6,
        a=100,
        c_star=0.1,
        rA_star=1.0,
        st_star=1.0,
        phi_j1=0.0,
        phi_As=60.0*pi/180,
        phi_Ae=170*pi/180,
        lambda_=0.97
    )
    num_points=100
    #xr,yr=gear.get_rackProfile(num=num_points)
    x,y=gear.get_arc_profile3(rotation_angle=pi,num_points=num_points,full_gear=False)
    x,y=KreisbogenRack(ra=10,num_points=num_points)
    #x,y=GearRack2D(1)
    plotFlanks(x,y,isRack=True)
