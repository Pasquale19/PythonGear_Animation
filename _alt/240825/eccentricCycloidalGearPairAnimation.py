import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Circle,FancyArrow,Polygon
from math import pi,sin,cos,asin,acos,atan,tan,atan2

from matplotlib.transforms import Affine2D
from matplotlib.gridspec import GridSpec
from eccentric_cycloid_gear import EccentricCycloidGear
from HideOnClick2 import hideOnClick
import  numpy as np
import shapely
from shapely.affinity import translate,rotate
from matplotlib.widgets import Slider


#Styles
mpl.rcParams['lines.linewidth'] = 2 #documentation https://matplotlib.org/stable/api/matplotlib_configuration_api.html
plt.rcParams['axes.spines.top']=False
plt.rcParams['axes.spines.right']=False

def plotReferenceCircle(ax,gear):
    r2_circle=Circle((0,0),gear.r2,linestyle='-.',color='g',fill=False)
    ax.add_patch(r2_circle)
    r1_circle=Circle((0,gear.a),gear.r1,linestyle='-.',color='b',fill=False)
    ax.add_patch(r1_circle)
    da1_circle=Circle((0,gear.a),gear.da1/2,linestyle=':',color='b',fill=False)
    ax.add_patch(da1_circle)
    da2_circle=Circle((0,0),gear.ra2,linestyle=':',color='g',fill=False)
    ax.add_patch(da2_circle)


def plotPoints(ax,gear,fontsize:int=12,color='black'):
    A=gear.A(steps=100)
    if A:
        Ax,Ay=A
        ax.annotate('A',(Ax+0.3,Ay+1),color=color, fontsize=fontsize)
        line_A=ax.plot([Ax],[Ay],c=color,label=f"A  $ζ={gear.zetaA*180/pi:.1f}$",marker='o')
    ax.plot([0],[gear.r2],c=color,label="C",marker='o')
    ax.text(0-0.01,gear.r2-2,'C',c=color, fontsize=fontsize,ha='center',va='top')
    Ex,Ey=gear.E(steps=100)
    ax.plot([Ex],[Ey],c=color,label=f"E  $ζ={gear.zetaE*180/pi:.1f}$",marker='o')
    ax.text(Ex-1,Ey-1.6,'E',c=color, fontsize=fontsize,ha='right',va='top')

    ax.scatter(x=0,y=gear.a,c='b',label="$O_1$")
    ax.text(0,gear.a+1,"$O_1$",c='b', fontsize=fontsize)
    ax.scatter(x=0,y=0,c='g',label="$O_2$")
    ax.text(0,0+1,"$O_2$",c='g', fontsize=fontsize)

def initForceVector(ax,gear,zeta,**kwargs):
    #allgemeines verzahnungsgesetz: Kraft muss immer durch den Wälzpunkt verlaufeb
    default_arrow_props = {
        'length_includes_head': True,
        'head_width':2,
        'head_length':3
    }
    M_length=kwargs.pop('M_length',0.5)
    mu=kwargs.pop('mu',0.20) 
    alpha=gear.alpha_t(zeta)
    xc,yc=gear.p_pOA(zeta)
    F1_length=M_length/cos(alpha)*gear.r1
    dn=np.array((cos(alpha),sin(alpha)))*F1_length
    Fn_arrow=ax.arrow(xc-dn[0],yc-dn[1],*dn,label="$F_n$",color='orangered',**default_arrow_props)
    dr=-np.array((sin(alpha),-cos(alpha)))*F1_length*mu
    Fr_arrow=ax.arrow(xc-dr[0],yc-dr[1],*dr,label="$F_r$",color='r',**default_arrow_props)
    return Fn_arrow,Fr_arrow

def setForceVector(zeta,gear,Fn_arrow:FancyArrow=None,Fr_arrow=None,**kwargs):
    if zeta==0:
        zeta=0.001
    alpha=gear.alpha_t(zeta)
    xc,yc=gear.p_pOA(zeta)
    #angle=atan2(yc,xc)
    M_length=kwargs.pop('M_length',0.5)
    mu=kwargs.pop('mu',0.20)
    F1_length=M_length/cos(alpha)*gear.r1
    for F in (Fn_arrow,Fr_arrow):
        F.set_visible(False)
    if Fn_arrow:
        dn=np.array((cos(alpha),sin(alpha)))*F1_length
        Fn_arrow.set_data(x=xc-dn[0],y=yc-dn[1],dx=dn[0],dy=dn[1])
    if Fr_arrow:
        dr=-np.array((sin(alpha),-cos(alpha)))*F1_length*mu
        Fr_arrow.set_data(x=xc-dr[0],y=yc-dr[1],dx=dr[0],dy=dr[1])

def initSpeedVector(ax,gear,zeta,**kwargs):
    #allgemeines verzahnungsgesetz: Kraft muss immer durch den Wälzpunkt verlaufeb
    default_arrow_props = {
        'length_includes_head': True,
        'head_width':2,
        'head_length':3
    }
    w_length=kwargs.pop('w_length',1)
    vt1,vt2,vg=gear.v_gear(zeta,w_length)
    xc,yc=gear.p_pOA(zeta)
    vt1_arrow=ax.arrow(xc-vt1[0],yc-vt1[1],*vt1,label="$v_{t1}$",color='gold',**default_arrow_props)
    vt2_arrow=ax.arrow(xc-vt2[0],yc-vt2[1],*vt2,label="$v_{t2}$",color='salmon',**default_arrow_props)
    vg_arrow=ax.arrow(xc-vg[0],yc-vg[1],*vg,label="$v_g$",color='goldenrod',**default_arrow_props)
    return vt1_arrow,vt2_arrow,vg_arrow

def setSpeedVector(zeta,gear,vt1_arrow:FancyArrow=None,vt2_arrow=None,vg_arrow=None,**kwargs):
    if zeta==0:
        zeta=0.001
    w_length=kwargs.pop('w_length',1)
    vt1,vt2,vg=gear.v_gear(zeta,w_length)
    alpha=gear.alpha_t(zeta)
    xc,yc=gear.p_pOA(zeta)
    for v in (vt1_arrow,vt2_arrow):
        v.set_visible(False)

    if vt1_arrow:
        vt1_arrow.set_data(x=xc-vt1[0],y=yc-vt1[1],dx=vt1[0],dy=vt1[1])
    if vt2_arrow:
        vt2_arrow.set_data(x=xc-vt2[0],y=yc-vt2[1],dx=vt2[0],dy=vt2[1])
    if vg_arrow:
        vg_arrow.set_data(x=xc-vg[0],y=yc-vg[1],dx=vg[0],dy=vg[1])

def setupPressureAnglePlot(ax,gear,zeta:float=0):
    # gear.A()
    # gear.E()
    zetaA=gear.zetaA
    zetaE=gear.zetaE
    angles=np.linspace(zetaA,zetaE,num=400)
    
    #xis=np.array([(gear.xi(a))*180/pi for a in angles])
    alphas=np.array([(pi/2-gear.xi(a))*180/pi for a in angles])
    
    x_coords,y_coords=gear.calculate_path_of_contact(start=zetaA,end=zetaE,num_points=len(angles))
    lengths = np.zeros(len(angles))
    for i in range(1, len(angles)):
        lengths[i] = lengths[i-1] + np.sqrt((x_coords[i] - x_coords[i-1])**2 + (y_coords[i] - y_coords[i-1])**2)
    ax.plot(lengths,alphas,label="α_t",color="black")
    ax.set_xlabel("Path of Contact [mm]")
    ax.set_title(r"Transverse pressure angle $\alpha_t$ [°]", loc='center',y=0.9)
    ax.text(0.5,0.9,r"Transverse pressure angle $\alpha_t$ [°]",
        horizontalalignment='center', 
        verticalalignment='top', 
        transform=ax.transAxes)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    alpha_zeta=gear.alpha_t(zeta)*180/pi

    yinterp = np.interp(-alpha_zeta, -alphas, lengths)
    print(f"alpha_point={alpha_zeta:.1f}° and length={yinterp:.1f}mm and zeta={zeta*180/pi:.1f}°")
    alpha_point,=ax.plot(lengths[0],alphas[0],marker='o',color="black",linestyle="",linewidth=4)
    return alpha_point

def setupSlidingVelocityPlot(ax,gear,zeta:float=0,**kwargs):
    # gear.A()
    # gear.E()
    zetaA=kwargs.get("start",gear.zetaA)
    zetaE=kwargs.get("end",gear.zetaE)
    angles=np.linspace(zetaA,zetaE,num=400)
    
    #xis=np.array([(gear.xi(a))*180/pi for a in angles])
    alphas=np.array([(pi/2-gear.xi(a))*180/pi for a in angles])
    
    x_coords,y_coords=gear.calculate_path_of_contact(start=zetaA,end=zetaE,num_points=len(angles))
    lengths = np.zeros(len(angles))
    for i in range(1, len(angles)):
        lengths[i] = lengths[i-1] + np.sqrt((x_coords[i] - x_coords[i-1])**2 + (y_coords[i] - y_coords[i-1])**2)
    ax.plot(lengths,alphas,label="v_g",color="black")
    ax.set_ylabel("Sliding Velocity [mm]")
    ax.set_xlabel("Path of Contact [mm]")
    ax.text(0.5,0.9,r"Transverse pressure angle $\alpha_t$ [°]",
        horizontalalignment='center', 
        verticalalignment='top', 
        transform=ax.transAxes)
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    alpha_zeta=gear.alpha_t(zeta)*180/pi

    yinterp = np.interp(-alpha_zeta, -alphas, lengths)
    print(f"alpha_point={alpha_zeta:.1f}° and length={yinterp:.1f}mm and zeta={zeta*180/pi:.1f}°")
    alpha_point,=ax.plot(lengths[0],alphas[0],marker='o',color="black",linestyle="",linewidth=4)
    return alpha_point


def addLines(ax,gear,zeta: float=0,**kwargs):
    alpha=gear.alpha_t(zeta)
    length=kwargs.get("length",90)
    props={"color":"black","linestyle":"--","linewidth":0.5,"alpha":0.8}
    C=gear.p_pOA(zeta)
    dt=np.array((-sin(alpha),cos(alpha)))*length*-1
    dn=np.array((cos(alpha),sin(alpha)))*length*-1
    P1_t=C-dt/2
    P1_n=C-dn/2
    tangent=ax.arrow(x=P1_t[0],y=P1_t[1],dx=dt[0],dy=dt[1],label="tangent",**props)
    normal=ax.arrow(x=P1_n[0],y=P1_n[1],dx=dn[0],dy=dn[1],label="line of action",**props)
    text_normal=ax.text(*P1_n,"flank normal",va="bottom",ha="left",fontsize=8)
    text_tangent=ax.text(*P1_t,"flank tangent",va="top",ha="left",fontsize=8,rotation=0)
    dct=np.array((cos(0),sin(0)))*length
    P1_ct=np.array((0,gear.r2))-2*dct/3
    ax.arrow(x=P1_ct[0],y=P1_ct[1],dx=dct[0],dy=dct[1],label="common tangent",**props)
    ax.text(*P1_ct,"common tangent",va="bottom",ha="left",fontsize=8)
    return tangent,normal,text_tangent,text_normal

def updateLines(ax,gear,zeta: float,tangent,normal,text_tangent,text_normal,**kwargs):
    alpha=gear.alpha_t(zeta)
    length=kwargs.get("length",90)
    C=gear.p_pOA(zeta)
    dt=np.array((-sin(alpha),cos(alpha)))*length*-1
    dn=np.array((cos(alpha),sin(alpha)))*length
    P1_t=C-dt/2
    P1_n=C-dn/2
    tangent.set_data(x=P1_t[0],y=P1_t[1],dx=dt[0],dy=dt[1])
    normal.set_data(x=P1_n[0],y=P1_n[1],dx=dn[0],dy=dn[1])
    text_normal.set_position(P1_n)
    text_tangent.set_position(P1_t)

def plot_gear_profiles(gear: EccentricCycloidGear):
    """
    Plots the profiles of both gears, their reference circles, and the path of contact.

    Args:
        gear (EccentricCycloidGear): An instance of the EccentricCycloidGear class.
    """
    # Create the plot
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)

    # Create GridSpec
    gs = GridSpec(2, 2, width_ratios=[3, 1])

    # Create subplots
    ax = fig.add_subplot(gs[:, 0])  # Left plot, spans both rows
    ax2 = fig.add_subplot(gs[0, 1])  # Top right plot
    ax3 = fig.add_subplot(gs[1, 1])  # Bottom right plot


    #fig, (ax,ax1) = plt.subplots(1, 2,figsize=(12, 12))
    plotPoints(ax,gear)
    


    initial_rotation=0
    #initial_rotation=14.7*180/pi*1/gear.i*0
    zeta=gear.zetaA
    kappa=gear.i*zeta+initial_rotation
    setupPressureAnglePlot(ax2,gear,zeta)
    setupSlidingVelocityPlot(ax3,gear,zeta)
    
    # Gear1
    arc_x, arc_y = gear.get_arc_profile(rotation_angle=kappa*0,num_points=1000,full_gear=True)
    polyPatch_Gear1=mpl.patches.Polygon(list(zip(arc_x, arc_y)),facecolor='lightblue',edgecolor='b',linestyle='-',alpha=0.5,linewidth=2)
    rotation = Affine2D().rotate_around(0, gear.a, kappa)
    polyPatch_Gear1.set_transform(rotation + ax.transData)
    ax.add_patch(polyPatch_Gear1)

    # Plot cycloid gear profile
    cycloidal_gear2=gear.ShapelyCycloidal2()
    cycloid_x,cycloid_y=cycloidal_gear2.exterior.xy
    polyPatch_Gear2=mpl.patches.Polygon(list(zip(cycloid_x, cycloid_y)),facecolor='lightgreen',edgecolor='g',linestyle='-',alpha=0.5,linewidth=2)
    ax.add_patch(polyPatch_Gear2)
    rotation2 = Affine2D().rotate_around(0, 0, -zeta)
    polyPatch_Gear2.set_transform(rotation2 + ax.transData)
    
    # Get path of contact
    contact_x, contact_y = gear.calculate_path_of_contact(num_points=1000)
    ax.plot(contact_x, contact_y, color='black',linestyle='-', label='Path of Contact')
    contact_x, contact_y = gear.calculate_path_of_contact(num_points=1000,start=-2*pi/gear.z1,end=2*pi/gear.z1)#full path
    ax.plot(contact_x, contact_y, c='black',linestyle=':', label='Path of Contact full')
    
    Fn_vec,Fr_vec=initForceVector(ax,gear,zeta)
    vt1_arrow,vt2_arrow,vg_arrow=initSpeedVector(ax,gear,zeta)
    tangent,normal,text_tangent,text_normal=addLines(ax,gear,zeta)
    # Plot reference circles
    plotReferenceCircle(ax,gear)
    
    #add contact point
    p_poa=gear.p_pOA(zeta)
    PC,=ax.plot(*p_poa,'yo',label="C2")
    p_poa_arrow=ax.arrow(0,0,*p_poa,label="$\overrightarrow{p}_{poa}$",color='black',linewidth=0.5)
    p_poa1,=ax.plot([0,p_poa[0]],[0,p_poa[1]],label="$\overrightarrow{p1}_{poa}$",color='black',linewidth=0.5)
    # Set plot properties
    ax.set_aspect('equal', 'box')
    ax.set_title('Eccentric Cycloid Gear Pair')
    
    ax.legend(loc='upper right',bbox_to_anchor=(0.2,1), fontsize=8)
    #ax.legend(loc='upper left',bbox_to_anchor=(1,1), prop={'size': 6})
    #ax.grid(True)

    # Set axis limits to ensure both gears are visible
    ax.set_xlim(-gear.ra2 - 1, gear.ra2 + 1)
    ax.set_ylim(0 - 5, gear.a+gear.da1/2)

    # Text for displaying angles
    angle_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, verticalalignment='top')

    # Create slider
    
    ax_slider = plt.axes([0.2, 0.1, 0.4, 0.03])
    slider = Slider(ax_slider, 'Rotation', gear.zetaA, gear.zetaE, valinit=zeta)
    
    def update(zeta):
        #zeta=(kappa+initial_rotation)/gear.i
        kappa=zeta*gear.i
        rotation = Affine2D().rotate_around(0, gear.a, kappa)
        polyPatch_Gear1.set_transform(rotation + ax.transData)
        #set gear 2
        rotation2 = Affine2D().rotate_around(0, 0, -zeta)
        polyPatch_Gear2.set_transform(rotation2 + ax.transData)
        #Kontaktpunkt
        x,y=gear.p_pOA(zeta)
        PC.set_data([x],[y])
        PC.set_label(f"C = {zeta*180/pi:.1f}°")

        angle_text.set_text(f"ζ={zeta*180/pi:.1f}°")
       
        updateLines(ax,gear,zeta,tangent,normal,text_tangent,text_normal)

        #set arrows
        setForceVector(zeta,gear,Fn_arrow=Fn_vec,Fr_arrow=Fr_vec)
        setSpeedVector(zeta,gear,vt1_arrow=vt1_arrow,vt2_arrow=vt2_arrow,vg_arrow=vg_arrow)
        p_poa=gear.p_pOA(zeta)
        p_poa_arrow.set_data(x=0,y=0,dx=p_poa[0],dy=p_poa[1])
        p_poa1.set_data([0,p_poa[0]],[gear.a,p_poa[1]])
        ax.legend()


    def invertVisibility(event):
        if event.key == 'v':
            vis_F=Fn_vec.get_visible()
            Fn_vec.set_visible(not vis_F)
            Fr_vec.set_visible(not vis_F)
            vis_tangent=tangent.get_visible()
            tangent.set_visible(not vis_tangent)
            normal.set_visible(not vis_tangent)
            text_normal.set_visible(not vis_tangent)
            text_tangent.set_visible(not vis_tangent)
            update(slider.val)

    
    slider.on_changed(update)
    fig.canvas.mpl_connect('key_press_event', invertVisibility)
    hideOnClick(fig,ax)
    fig.savefig("export/start.svg", format='svg', dpi=1200)
    ax.axis('off')
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

    # Plot the gear profiles and path of contact
    plot_gear_profiles(gear_pair)