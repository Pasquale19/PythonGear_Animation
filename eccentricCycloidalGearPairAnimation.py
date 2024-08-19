import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from math import pi,sin,cos,asin,acos,atan
from eccentric_cycloid_gear import EccentricCycloidGear
from HideOnClick2 import hideOnClick
import  numpy as np
import shapely
from shapely import Point,Polygon,MultiPolygon
from shapely.affinity import translate,rotate
from matplotlib.widgets import Slider


def rotation_vector(gear: EccentricCycloidGear,rotation_angle:float=0):
    e=gear.e
    qF1=gear.q_F1
    z1=gear.z1
    angle2=-pi/z1+rotation_angle
    a=gear.a
    theta=rotation_angle-gear.phi_s1/2 #winkel des rihtungsvektors
    a2=-gear.phi_As+theta+pi
    xt,yt=np.array((sin(theta),-cos(theta)))*e+np.array((0,a))+np.array((sin(a2),-cos(a2)))*gear.rA
    a2=-gear.phi_As+theta
    xt,yt=np.array((sin(theta),-cos(theta)))*e+np.array((0,a))+np.array((-sin(a2),cos(a2)))*gear.rA
    x=[sin(rotation_angle)*e,0,sin(angle2)*qF1]
    y=[-cos(rotation_angle)*e+a,a,-cos(angle2)*qF1+a]
    x.insert(0, xt)
    y.insert(0, yt)
    return x,y


def plotReferenceCircle(ax,gear):
    r_start=Circle((0,gear.a),gear.r_start,linestyle='-.',color='b',label=f"r_start={gear.r_start:.2f}",fill=False)
    ax.add_patch(r_start)
    ecc_circle=Circle((0,gear.a),gear.e,linestyle='-.',color='b',label=f"e={gear.e:.2f}",fill=False)
    ax.add_patch(ecc_circle)
    r1_circle=Circle((0,gear.a),gear.r1,linestyle='--',color='b',label=f"r1={gear.r1:.1f}",fill=False)
    ax.add_patch(r1_circle)
    qF1_circle=Circle((0,gear.a),gear.q_F1,linestyle=':',color='b',label=f"qF1={gear.q_F1:.1f}",fill=False)
    ax.add_patch(qF1_circle)
    da1_circle=Circle((0,gear.a),gear.da1/2,linestyle=':',color='y',label=f"da1={gear.da1:.1f}",fill=False)
    ax.add_patch(da1_circle)
    df1_circle=Circle((0,gear.a),gear.df1/2,linestyle=':',color='y',label=f"df1={gear.df1:.1f}",fill=False)
    ax.add_patch(df1_circle)
    da2_circle=Circle((0,0),gear.ra2,linestyle=':',color='r',label=f"ra2={gear.ra2:.1f}",fill=False)
    ax.add_patch(da2_circle)
    df2_circle=Circle((0,0),gear.rf2,linestyle=':',color='r',label=f"rf2={gear.rf2:.1f}",fill=False)
    ax.add_patch(df2_circle)

def plot_gear_profiles(gear: EccentricCycloidGear):
    """
    Plots the profiles of both gears, their reference circles, and the path of contact.

    Args:
        gear (EccentricCycloidGear): An instance of the EccentricCycloidGear class.
    """
    # Get gear profiles
    initial_rotaion=-gear.phi_tooth/2
    initial_angle=-gear.phi_tooth/2-pi/2
        # Create the plot
    fig, ax = plt.subplots(figsize=(12, 12))
    arc_x, arc_y = gear.get_arc_profile2(rotation_angle=initial_rotaion,num_points=1000,full_gear=True)
    cycloid_x, cycloid_y = gear.get_cycloid_profile(num_points=1000)
    x_vec,y_vec=rotation_vector(gear,rotation_angle=0)
    line_vec,=ax.plot(x_vec,y_vec,label="rotation")
    
    # # Get path of contact
    contact_x, contact_y = gear.calculate_path_of_contact(num_points=1000)


    
    # # Plot arc gear profile
    line_arcGear,=ax.plot(arc_x, np.array(arc_y), 'b-', label='Arc Gear',linestyle='-')
    
    # # Plot cycloid gear profile
    ax.plot(np.array(cycloid_x) + 0, cycloid_y, 'r-', label='Cycloid Gear')
    # # Plot path of contact
    ax.plot(contact_x, contact_y, 'g-', label='Path of Contact')

    
    #initial_angle=-pi/2
    print(f"initial_angle={initial_angle:.1f}")
    poly=gear.ShapelyGear(initial_angle=initial_angle,center_x=0,center_y=gear.a)
    #line_gear1,=ax.plot(*poly.exterior.xy,color='b',label="gear1")



    # Plot reference circles
    plotReferenceCircle(ax,gear)
    
    #add contact point
    zeta=0
    ax.scatter(x=0,y=gear.r2,c='r',label="C")
    PC,=ax.plot(*gear.contact_Point(zeta),'yo',label="C2")
    Ax,Ay=gear.A()
    ax.annotate('A',(Ax+0.1,Ay+0.1),color='g', fontsize=14)
    line_A=ax.scatter(x=Ax,y=Ay,c='g',label="A")
    Ex,Ey=gear.E()
    #Ex,Ey=gear.contact_Point((gear.phi_Ae-gear.phi_As)/gear_pair.i)
    line_E=ax.scatter(x=Ex,y=Ey,c='m',label="E")
    ax.text(Ex+1,Ey+0.1,'E',c='m', fontsize=14)


    # Set plot properties
    ax.set_aspect('equal', 'box')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Eccentric Cycloid Gear Pair')
    
    
    ax.legend()
    ax.grid(True)
    params=gear.list_parameters()
    # Set axis limits to ensure both gears are visible
    #ax.set_xlim(-gear.r1 - 1, gear.a + gear.r2 + 1)
    #ax.set_ylim(0 - 1, gear.a+gear.da1/2)

    # Text for displaying angles
    angle_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, verticalalignment='top')

    # Create slider
    
    ax_slider = plt.axes([0.2, 0.1, 0.6, 0.03])
    slider = Slider(ax_slider, 'Rotation', -pi/gear.z1, pi/gear.z1, valinit=zeta)
    
    def update(val):
        arc_x, arc_y = gear.get_arc_profile2(rotation_angle=val+initial_rotaion,num_points=100,full_gear=True)
        line_arcGear.set_data(arc_x,arc_y)
        x,y=gear.contact_Point(val)
        PC.set_data([x],[y])
        poly2=rotate(poly,angle=val,origin=(0,gear.a),use_radians=True)
        x,y=poly2.exterior.xy
        #line_gear1.set_data(x,y)
        #gearPoly1.set_data()
        angle_text.set_text(f"φ={val*180/pi}°")

        x_vec,y_vec=rotation_vector(gear,val)
        line_vec.set_data(x_vec,y_vec)
        line_vec.set_label=f"rotation_vec{val*180/pi:.0f}°"
        
    slider.on_changed(update)

    hideOnClick(fig,ax)
    # Show the plot
    plt.show()

    # Print gear parameters
    print(params)

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