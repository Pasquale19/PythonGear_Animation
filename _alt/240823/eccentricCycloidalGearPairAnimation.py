import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Circle,FancyArrow
from math import pi,sin,cos,asin,acos,atan,tan,atan2
from eccentric_cycloid_gear import EccentricCycloidGear
from HideOnClick2 import hideOnClick
import  numpy as np
import shapely
from shapely import Point,Polygon,MultiPolygon
from shapely.affinity import translate,rotate
from matplotlib.widgets import Slider


#Styles
mpl.rcParams['lines.linewidth'] = 2 #documentation https://matplotlib.org/stable/api/matplotlib_configuration_api.html



def plotReferenceCircle(ax,gear):
    r2_circle=Circle((0,0),gear.r2,linestyle='-.',color='r',label=f"r2={gear.r2:.1f}",fill=False)
    ax.add_patch(r2_circle)
    r1_circle=Circle((0,gear.a),gear.r1,linestyle='-.',color='b',label=f"r1={gear.r1:.1f}",fill=False)
    ax.add_patch(r1_circle)
    da1_circle=Circle((0,gear.a),gear.da1/2,linestyle=':',color='b',label=f"ra1={gear.ra1:.1f}",fill=False)
    ax.add_patch(da1_circle)
    da2_circle=Circle((0,0),gear.ra2,linestyle=':',color='r',label=f"ra2={gear.ra2:.1f}",fill=False)
    ax.add_patch(da2_circle)


def plotPoints(ax,gear,fontsize:int=12):
    A=gear.A(steps=100)
    if A:
        Ax,Ay=A
        ax.annotate('A',(Ax+0.1,Ay+0.1),color='g', fontsize=fontsize)
        line_A=ax.plot([Ax],[Ay],c='g',label="A",marker='o')
    ax.plot([0],[gear.r2],c='r',label="C",marker='o')
    ax.text(0+0.01,gear.r2-1.5,'C',c='r', fontsize=fontsize)
    Ex,Ey=gear.E(steps=100)
    ax.scatter(x=Ex,y=Ey,c='m',label="E")
    ax.text(Ex+1,Ey+0.1,'E',c='m', fontsize=fontsize)

    ax.scatter(x=0,y=gear.a,c='b',label="$O_1$")
    ax.text(0,gear.a+1,"$O_1$",c='b', fontsize=fontsize)
    ax.scatter(x=0,y=0,c='r',label="$O_2$")
    ax.text(0,0+1,"$O_2$",c='r', fontsize=fontsize)

def initForceVector(ax,gear,zeta,**kwargs):
    #allgemeines verzahnungsgesetz: Kraft muss immer durch den Wälzpunkt verlaufeb
    default_arrow_props = {
        'length_includes_head': True,
        'head_width':2,
        'head_length':3
    }
    M_length=kwargs.pop('M_length',1)
    mu=kwargs.pop('mu',0.10) 
    alpha=gear.alpha_t(zeta)
    xc,yc=gear.contact_Point(-zeta)
    F1_length=M_length/cos(alpha)*gear.r1
    dn=np.array((cos(alpha),-sin(alpha)))*F1_length
    Fn_arrow=ax.arrow(xc-dn[0],yc-dn[1],*dn,label="$F_n$",color='orangered',**default_arrow_props)
    dr=-np.array((sin(alpha),cos(alpha)))*F1_length*mu
    Fr_arrow=ax.arrow(xc-dr[0],yc-dr[1],*dr,label="$F_r$",color='r',**default_arrow_props)
    return Fn_arrow,Fr_arrow

def setForceVector(zeta,gear,Fn_arrow:FancyArrow=None,Fr_arrow=None,**kwargs):
    alpha=gear.alpha_t(zeta)
    xc,yc=gear.contact_Point(-zeta)
    #angle=atan2(yc,xc)
    M_length=kwargs.pop('M_length',1)
    mu=kwargs.pop('mu',0.10)
    F1_length=M_length/cos(alpha)*gear.r1
    if Fn_arrow:
        dn=np.array((cos(alpha),-sin(alpha)))*F1_length
        Fn_arrow.set_data(x=xc-dn[0],y=yc-dn[1],dx=dn[0],dy=dn[1])
    if Fr_arrow:
        dr=-np.array((sin(alpha),cos(alpha)))*F1_length*mu
        Fr_arrow.set_data(x=xc-dr[0],y=yc-dr[1],dx=dr[0],dy=dr[1])


def plot_gear_profiles(gear: EccentricCycloidGear):
    """
    Plots the profiles of both gears, their reference circles, and the path of contact.

    Args:
        gear (EccentricCycloidGear): An instance of the EccentricCycloidGear class.
    """
    length_Force=10
    initial_angle=-gear.phi_tooth/2-pi/2
        # Create the plot
    fig, ax = plt.subplots(figsize=(12, 12))
    plotPoints(ax,gear)

    initial_rotation=gear.zetaA
    initial_rotation=0
    
    # # Plot arc gear profile
    arc_x, arc_y = gear.get_arc_profile(rotation_angle=initial_rotation,num_points=1000,full_gear=True)
    line_arcGear,=ax.plot(arc_x, np.array(arc_y), 'b', label='Arc Gear',linestyle='-')
    polyPatch_Gear1=Polygon(list(zip(arc_x, arc_y)),edgecolor='b',linestyle='-',alpha=0.1,facecolor='b')
    ax.add_patch(polyPatch_Gear1)
    
    
    # # Get path of contact
    contact_x, contact_y = gear.calculate_path_of_contact(num_points=1000)
    ax.plot(contact_x, contact_y, 'g-', label='Path of Contact')
    contact_x, contact_y = gear.calculate_path_of_contact(num_points=1000,start=-2*pi/gear.z1,end=2*pi/gear.z1)#full path
    ax.plot(contact_x, contact_y, c='g',linestyle=':', label='Path of Contact full')
    
    # Plot cycloid gear profile
    cycloidal_gear2=gear.ShapelyCycloidal2()
    cycloid_x,cycloid_y=cycloidal_gear2.exterior.xy
    cycloidal_gear_line2,=ax.plot(np.array(cycloid_x) + 0, cycloid_y, 'm-', label='Cycloid Gear')


    poly=gear.ShapelyGear(initial_angle=initial_angle,center_x=0,center_y=gear.a)
    
    
    Fn_vec,Fr_vec=initForceVector(ax,gear,initial_rotation*1/gear.i)

    # Plot reference circles
    plotReferenceCircle(ax,gear)
    
    #add contact point
    zeta=0
    
    PC,=ax.plot(*gear.p_pOA(zeta),'yo',label="C2")


    # Set plot properties
    ax.set_aspect('equal', 'box')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Eccentric Cycloid Gear Pair')
    
    
    ax.legend()
    ax.grid(True)

    # Set axis limits to ensure both gears are visible
    ax.set_xlim(-gear.ra2 - 1, gear.ra2 + 1)
    ax.set_ylim(0 - 1, gear.a+gear.da1/2)

    # Text for displaying angles
    angle_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, verticalalignment='top')

    # Create slider
    
    ax_slider = plt.axes([0.2, 0.1, 0.6, 0.03])
    slider = Slider(ax_slider, 'Rotation', -1.5*pi/gear.z1, 1.5*pi/gear.z1, valinit=zeta)
    
    def update(kappa):
        zeta=-(kappa+initial_rotation)/gear.i
        #set gear 1
        arc_x, arc_y = gear.get_arc_profile2(rotation_angle=kappa+initial_rotation,num_points=100,full_gear=True)   
        line_arcGear.set_data(arc_x,arc_y)
        #Kontaktpunkt
        x,y=gear.p_pOA(-zeta)
        PC.set_data([x],[y])
        PC.set_label(f"C = {-zeta*180/pi:.1f}°")

        #line_gear1.set_data(x,y)
        #gearPoly1.set_data()
        angle_text.set_text(f"φ={kappa*180/pi:.1f}°")

        #set gear 2
        gear22=rotate(cycloidal_gear2,zeta,origin=(0,0),use_radians=True)
        cycloidal_gear_line2.set_data(*gear22.exterior.xy)
        
        #set arrows
        setForceVector(zeta,gear,Fn_arrow=Fn_vec,Fr_arrow=Fr_vec)
        ax.legend()
        
    slider.on_changed(update)

    hideOnClick(fig,ax)
    plt.show()

if __name__ == "__main__":
    # Create an instance of EccentricCycloidGear
    gear_pair = EccentricCycloidGear(
        z1=3,
        z2=6,
        a=100,
        rA_star=1.0,
        st_star=1.0,
        phi_j1=1.0*pi/180,
        phi_As=60*pi/180,
        phi_Ae=170*pi/180,
        lambda_=0.97
    )
    parameters_list = gear_pair.list_parameters()
    print(parameters_list)

    # Plot the gear profiles and path of contact
    plot_gear_profiles(gear_pair)