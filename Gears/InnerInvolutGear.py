import matplotlib.pyplot as plt
import numpy as np
from math import pi,cos,sin,tan,atan2,atan,sqrt,acos,tau
import matplotlib as mpl
from matplotlib.patches import Circle
from BaseGear2 import InvoluteGear
from InvoluteGear2D import InvoluteGear2D,involute_coordinates,InvoluteGear2D2,InvoluteGear2D3
#from ..HideOnClick2 import hideOnClick

if __name__ =="__main__":
    hap=0.4
    hfp=0.5
    m=2.15
    z=-61
    alpha=35*pi/180
    M=np.array((23,12))*0
    x=0.0
    fig,ax=plt.subplots()
    center2=M
    r2=m*z/2
    rb2=m*z/2*cos(alpha)
    ra2=m*(z/2+x+hap)
    rf2=m*(z/2+x-hfp)

    def plotGear(zi=z,hapi=0.4,hfpi=0.6,xi=0):
        x2,y2=InvoluteGear2D3(z=zi,m=m,alpha=alpha,x=xi,rotation_angle=0,M=center2*1,hap=hapi,hfp=hfpi)
        coords= list(zip(x2, y2))
        ply= mpl.patches.Polygon(coords,facecolor='red', edgecolor='yellow',
                                                    linestyle='-', alpha=0.3, linewidth=2,label=f"gear $h_a$={hapi:.1f} $h_f$={hfpi:.1f} x={xi}")
        #ax.add_patch(ply)
    
    #plotGear(zi=z,hapi=0.4,hfpi=0.5,xi=0)
    plotGear(zi=abs(z),hapi=hap,hfpi=hfp,xi=x)

    x2,y2=InvoluteGear2D3(z=z,m=m,alpha=alpha,x=x,rotation_angle=0,M=center2*1,hap=hap,hfp=hfp,ax=ax)
    #ax.plot(x2,y2,marker=".", markersize=4,linestyle="",label=f"gear2 $x=${x} $hap=${hap}",color="orange")
    coords= list(zip(x2, y2))
    ply= mpl.patches.Polygon(coords,facecolor='lightblue', edgecolor='b',
                                                   linestyle='-', alpha=0.3, linewidth=2,label="gear2")

    ax.add_patch(ply)
 
    # alpha_f=acos(rb2/r2)
    # alpha_a=acos(rb2/ra2)
    # print(f"rb2={rb2:.2f} rf2={rf2:.2f} alpha_f={alpha_f*180/pi:.2f}")
    # delta_f=tan(alpha_f)-alpha_f
    # delta_a=tan(alpha_a)-alpha_a
    # t_max=delta_f+alpha
    # t_min=delta_a+alpha
    # print(f"t_max={t_max:.2f} delta_f={delta_f*180/pi:.2f} t_min={t_min:.2f} delta_a={delta_a*180/pi:.2f}")
    # t=np.linspace(t_min,t_max,num=100,endpoint=True)
    # delta0=tan(alpha)-alpha
    # offs=delta0+abs(tau/(z*4))
    # involute_x,involute_y=involute_coordinates(t,rb2,a=offs)
    # lbl="Invol $\delta_0$=" f"{delta0*180/pi:.1f}°" f"   $offs$={offs*180/pi:.0f}°"
    # ax.plot(involute_x,involute_y+center2[1],label=lbl,markersize=2,marker='o')
    #ax.plot(involute_x,-involute_y,label=lbl,markersize=2,marker='o')
    ax.grid(True)
    ax.axhline(y=center2[1])
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(0.2,0.7),loc="upper right",fontsize=6)
    #hideOnClick(fig=fig,ax=ax)
    plt.show()