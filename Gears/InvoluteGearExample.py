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
    z=51
    da=m*(z+hap*2)
    df=m*(z-hfp*2)
    alpha=35*pi/180
    db=z*m*cos(alpha)
    M=np.array((23,12))
    x2=-0.2
    gear_pair = InvoluteGear(
        z1=50,
        z2=-51,
        m=2.15,
        alpha=alpha,
        hap1=hap,hap2=hap,
        hfp1=0.5,hfp2=0.5,
        x1=0.2,x2=x2
    )
    gear_pair.print_Parameter()
    fig,ax=plt.subplots()

    
    def plot_reference_circles2():
        # Plot reference circles
        rw2_circle = Circle(gear_pair.center2, abs(gear_pair.rw2), linestyle='-.', color='b', fill=False,label=f"rw2={gear_pair.rw2:.2f}")  
        da2_circle = Circle(gear_pair.center2, abs(gear_pair.ra2), linestyle=':', color='g', fill=False,label=f"ra2={gear_pair.ra2:.2f}") 
        db2_circle = Circle(gear_pair.center2, abs(gear_pair.rb2), linestyle='--', color='g', fill=False,label=f"rb2={gear_pair.rb2:.2f} $a_w$={gear_pair.alpha_w*180/pi:.3f}° \t $alph$={gear_pair.alpha*180/pi:.3f}°") 
        df2_circle = Circle(gear_pair.center2, abs(gear_pair.rf2), linestyle='--', color='g', fill=False,label=f"rf2={gear_pair.rf2:.2f}")
        r2_circle = Circle(gear_pair.center2, abs(gear_pair.r2), linestyle='-.', color='r', fill=False,label=f"$r_2$={gear_pair.r2:.2f}")
       
        for circle in [rw2_circle, da2_circle,db2_circle,df2_circle,r2_circle]:
            ax.add_patch(circle)
    
    #x,y=gear_pair.gearGeometry2(rotation_angle=0)
    plot_reference_circles2()
    #ax.plot(x,y,marker=".", markersize=1,linestyle="",label="gear2")
    
    x2,y2=InvoluteGear2D3(z=gear_pair.z2,m=gear_pair.m,alpha=gear_pair.alpha,x=gear_pair.x2,rotation_angle=0,M=gear_pair.center2*1,hap=gear_pair.hap2,hfp=gear_pair.hfp2,ax=ax)
    ax.plot(x2,y2,marker=".", markersize=4,linestyle="",label=f"gear2 $x=${gear_pair.x2} $hap=${gear_pair.hap2}",color="orange")
    # coords= list(zip(x, y))
    # ply= mpl.patches.Polygon(coords,facecolor='lightblue', edgecolor='b',
    #                                                linestyle='-', alpha=0.5, linewidth=2,label="gear2")

    #ax.add_patch(ply)
    alpha_f=acos(gear_pair.rb2/gear_pair.rf2)
    alpha_a=acos(gear_pair.rb2/gear_pair.ra2)
    print(f"rb2={gear_pair.rb2:.2f} rf2={gear_pair.rf2:.2f} alpha_f={alpha_f*180/pi:.2f}")
    delta_f=tan(alpha_f)-alpha_f
    delta_a=tan(alpha_a)-alpha_a
    t_max=delta_f+alpha
    t_min=delta_a+alpha
    print(f"t_max={t_max:.2f} delta_f={delta_f*180/pi:.2f} t_min={t_min:.2f} delta_a={delta_a*180/pi:.2f}")
    t=np.linspace(t_min,t_max,num=100,endpoint=True)
    delta0=tan(gear_pair.alpha)-gear_pair.alpha
    offs=delta0+abs(tau/(gear_pair.z2*4))
    involute_x,involute_y=involute_coordinates(t,gear_pair.rb2,a=offs)
    lbl="Invol $\delta_0$=" f"{delta0*180/pi:.1f}°" f"   $offs$={offs*180/pi:.0f}°"
    ax.plot(involute_x,involute_y+gear_pair.center2[1],label=lbl,markersize=2,marker='o')
    #ax.plot(involute_x,-involute_y,label=lbl,markersize=2,marker='o')
    ax.grid(True)
    ax.axhline(y=gear_pair.center2[1])
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(0.2,1),loc="upper right",fontsize=6)
    #hideOnClick(fig=fig,ax=ax)
    plt.show()