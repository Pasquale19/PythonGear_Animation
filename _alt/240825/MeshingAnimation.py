import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Circle, FancyArrow,Arc
from matplotlib.transforms import Affine2D
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider
import numpy as np
from math import pi,sin,cos,asin,acos,atan,tan,atan2,hypot
from eccentric_cycloid_gear import EccentricCycloidGear
from HideOnClick2 import hideOnClick

#Styles
mpl.rcParams['lines.linewidth'] = 2 #documentation https://matplotlib.org/stable/api/matplotlib_configuration_api.html
plt.rcParams['axes.spines.top']=False
plt.rcParams['axes.spines.right']=False


def calculate_arc(center, radius, start_angle, end_angle, num_points=100):
    '''calculates the coordinates of an arc, the angles must be defined in deg'''
    theta = np.linspace(np.deg2rad(start_angle), np.deg2rad(end_angle), num_points)
    x = center[0] + radius * np.cos(theta)
    y = center[1] + radius * np.sin(theta)
    return x, y

# Define your class here
class GearMeshingPlot:
    def __init__(self, gear:EccentricCycloidGear):
        self.gear = gear
        self.fig = None
        self.ax = None
        self.ax_pressure_angle = None
        self.ax_sliding_velocity = None
        self.slider = None
        self.angle_text = None
        self.zeta0=gear.zetaA
        # Plot elements to be updated
        self.polyPatch_Gear1 = None
        self.polyPatch_Gear2 = None
        self.Fn_vec = None
        self.Fr_vec = None
        self.vt1_arrow = None
        self.vt2_arrow = None
        self.vg_arrow = None
        self.tangent = None
        self.normal = None
        self.text_tangent = None
        self.text_normal = None
        self.PC = None
        self.p_poa_arrow = None
        self.p_poa1 = None
        self.marker_2=None
        self.marker_pressureAngle=None
        self.arc_alpha=None

        self.arrowprops=dict(head_length=1.5, head_width=1.2, lw=1,length_includes_head=True)
        self.setup_plot()

    def setup_plot(self):
        # Create figure and subplots
        gear=self.gear
        self.fig = plt.figure(figsize=(19.2, 10.8), dpi=100,num="MeshingAnimation")
        gs = GridSpec(2, 2, width_ratios=[3, 1], left=0.00, right=0.9, top=0.99, bottom=0.01, wspace=0.1, hspace=0.6) #distance to the left side, right side -> left=0,right=1 is maximum

        self.ax = self.fig.add_subplot(gs[:, 0])
        self.ax_pressure_angle = self.fig.add_subplot(gs[0, 1])
        self.ax_sliding_velocity = self.fig.add_subplot(gs[1, 1])
        zeta0=self.zeta0
        # Initial plot configuration
        self.plotPoints(ax=self.ax,gear=self.gear)
        self.plot_initial_gear_profiles()
        self.setup_pressure_angle_plot(ax=self.ax_pressure_angle,gear=self.gear,zeta=zeta0)
        self.setup_sliding_velocity_plot(ax=self.ax_sliding_velocity,gear=self.gear,zeta=zeta0)
        self.setup_slider()
        self.init_angle(zeta0)

        # Connect event handlers
        self.fig.canvas.mpl_connect('key_press_event', self.invert_visibility)

        # Style the plot
        ax=self.ax
        ax.set_xlim(-gear.ra2 - 1, gear.ra2 + 1)
        ax.set_ylim(0 - 5, gear.a+gear.da1/2)
        self.ax.axis('off')
        ax.set_aspect('equal', 'box')
        ax.set_title('Eccentric Cycloid Gear Pair')
        ax.legend(loc='upper right',bbox_to_anchor=(0.2,1), fontsize=8)

        plt.show()

    def plotPoints(self,ax,gear,fontsize:int=12,color='black'):
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

    def plot_initial_gear_profiles(self):
        # Plot initial points and gear profiles
        self.plot_reference_circles()
        self.plot_gear_polygons()
        self.plot_contact_path()
        self.initialize_vectors()
        self.add_tangent_and_normal_lines()

        self.angle_text = self.ax.text(0.05, 0.95, '', transform=self.ax.transAxes, verticalalignment='top')

    def plot_reference_circles(self):
        # Plot reference circles
        r2_circle = Circle((0, 0), self.gear.r2, linestyle='-.', color='g', fill=False)
        r1_circle = Circle((0, self.gear.a), self.gear.r1, linestyle='-.', color='b', fill=False)
        da1_circle = Circle((0, self.gear.a), self.gear.da1 / 2, linestyle=':', color='b', fill=False)
        da2_circle = Circle((0, 0), self.gear.ra2, linestyle=':', color='g', fill=False)
        for circle in [r2_circle, r1_circle, da1_circle, da2_circle]:
            self.ax.add_patch(circle)

    def plot_gear_polygons(self):
        # Initial rotation
        initial_rotation = 0
        zeta = self.gear.zetaA
        kappa = self.gear.i * zeta + initial_rotation

        # Plot Gear1
        arc_x, arc_y = self.gear.get_arc_profile(rotation_angle=kappa * 0, num_points=1000, full_gear=True)
        self.polyPatch_Gear1 = mpl.patches.Polygon(list(zip(arc_x, arc_y)), facecolor='lightblue', edgecolor='b',
                                                   linestyle='-', alpha=0.5, linewidth=2)
        rotation = Affine2D().rotate_around(0, self.gear.a, kappa)
        self.polyPatch_Gear1.set_transform(rotation + self.ax.transData)
        self.ax.add_patch(self.polyPatch_Gear1)

        # Plot Gear2
        cycloidal_gear2 = self.gear.ShapelyCycloidal2()
        cycloid_x, cycloid_y = cycloidal_gear2.exterior.xy
        self.polyPatch_Gear2 = mpl.patches.Polygon(list(zip(cycloid_x, cycloid_y)), facecolor='lightgreen',
                                                   edgecolor='g', linestyle='-', alpha=0.5, linewidth=2)
        rotation2 = Affine2D().rotate_around(0, 0, -zeta)
        self.polyPatch_Gear2.set_transform(rotation2 + self.ax.transData)
        self.ax.add_patch(self.polyPatch_Gear2)

    def plot_contact_path(self,**kwargs):
        # Get path of contact
        contact_x, contact_y = self.gear.calculate_path_of_contact(num_points=1000)
        self.ax.plot(contact_x, contact_y, color='black', linestyle='-', label='Path of Contact')
        contact_x, contact_y = self.gear.calculate_path_of_contact(num_points=1000,
                                                                   start=-2 * np.pi / self.gear.z1,
                                                                   end=2 * np.pi / self.gear.z1)
        self.ax.plot(contact_x, contact_y, c='black', linestyle=':', label='maximum possible Path of Contact')
        zeta=kwargs.get("zeta",self.zeta0)
        p_poa=self.gear.p_pOA(zeta)
        self.PC,=self.ax.plot(*p_poa,'yo',label="C2")

    def initialize_vectors(self):
        zeta = self.gear.zetaA
        self.Fn_vec, self.Fr_vec = self.init_force_vector(zeta)
        self.vt1_arrow, self.vt2_arrow, self.vg_arrow = self.init_speed_vector(zeta)

    def add_tangent_and_normal_lines(self):
        zeta = self.zeta0
        self.tangent, self.normal, self.text_tangent, self.text_normal = self.add_lines(zeta)

    def setup_slider(self):
        zeta = self.zeta0
        ax_slider = plt.axes([0.2, 0.1, 0.4, 0.03])
        self.slider = Slider(ax_slider, 'Rotation', self.gear.zetaA, self.gear.zetaE, valinit=zeta)
        self.slider.on_changed(self.update_plot)

    def update_plot(self, zeta):
        kappa = zeta * self.gear.i
        rotation = Affine2D().rotate_around(0, self.gear.a, kappa)
        self.polyPatch_Gear1.set_transform(rotation + self.ax.transData)

        rotation2 = Affine2D().rotate_around(0, 0, -zeta)
        self.polyPatch_Gear2.set_transform(rotation2 + self.ax.transData)

        # Update contact point
        x, y = self.gear.p_pOA(zeta)
        self.PC.set_data([x], [y])
        self.PC.set_label(f"C = {zeta * 180 / np.pi:.1f}°")
        self.angle_text.set_text(f"ζ={zeta * 180 / np.pi:.1f}°")

        # Update lines and vectors
        self.update_lines(zeta)
        self.update_vectors(zeta)
        self.update_pressureAnglePlot(zeta)
        self.update_pressureAnglePlot(zeta)
        self.updateAngles(zeta)
        self.ax.legend()
    
    def update_pressureAnglePlot(self, zeta):
        m=self.marker_pressureAngle
        alpha=self.gear.alpha_t(zeta)*180/pi
        self.marker_pressureAngle.set_data([zeta],[alpha])

    def update_SlidingVelocityPlot(self, zeta):
        m=self.marker_2
        #m.set_data()

    def update_vectors(self,zeta:float,**kwargs):
        gear=self.gear
        vt1_arrow,vt2_arrow,vg_arrow=self.vt1_arrow,self.vt2_arrow,self.vg_arrow
        if zeta==0:
            zeta=0.001
        w_length=kwargs.pop('w_length',1)
        vt1,vt2,vg=gear.v_gear(zeta,w_length)
        alpha=gear.alpha_t(zeta)
        xc,yc=gear.p_pOA(zeta)

        if vt1_arrow:
            vt1_arrow.set_data(x=xc-vt1[0],y=yc-vt1[1],dx=vt1[0],dy=vt1[1])
        if vt2_arrow:
            vt2_arrow.set_data(x=xc-vt2[0],y=yc-vt2[1],dx=vt2[0],dy=vt2[1])
        if vg_arrow:
            vg_arrow.set_data(x=xc-vg[0],y=yc-vg[1],dx=vg[0],dy=vg[1])

        #Forces
        M_length=kwargs.pop('M_length',0.5)
        mu=kwargs.pop('mu',0.20)
        F1_length=M_length/cos(alpha)*gear.r1
        Fn_arrow,Fr_arrow=self.Fn_vec,self.Fr_vec

        if Fn_arrow:
            dn=np.array((cos(alpha),sin(alpha)))*F1_length
            Fn_arrow.set_data(x=xc-dn[0],y=yc-dn[1],dx=dn[0],dy=dn[1])
        if Fr_arrow:
            dr=-np.array((sin(alpha),-cos(alpha)))*F1_length*mu
            Fr_arrow.set_data(x=xc-dr[0],y=yc-dr[1],dx=dr[0],dy=dr[1])

    def updateAngles(self,zeta:float):
        angle_alpha=self.arc_alpha
        alpha=self.gear.xi(zeta)*180/pi
        x,y=calculate_arc(center=(0,self.gear.r2),radius=20,start_angle=0,end_angle=alpha)
        angle_alpha.set_data(x,y)

    def invert_visibility(self, event):
        if event.key == 'v':
            vis_F = self.Fn_vec.get_visible()
            self.Fn_vec.set_visible(not vis_F)
            self.Fr_vec.set_visible(not vis_F)
            vis_tangent = self.tangent.get_visible()
            self.tangent.set_visible(not vis_tangent)
            self.normal.set_visible(not vis_tangent)
            self.text_normal.set_visible(not vis_tangent)
            self.text_tangent.set_visible(not vis_tangent)
            self.update_plot(self.slider.val)
            self.fig.draw()

    def init_angle(self,zeta):
        ax=self.ax
        xt,yt=(0,self.gear.r2)
        width=25
        ax.text(xt+width*0.9,yt,r"$\alpha$",va="top",ha="right",fontsize=12)
        alpha=90+self.gear.xi(self.zeta0)*180/pi
        print(f"alpha0={alpha:.1f}")
        #self.arc_alpha=Arc((xt,yt),width=width*2,height=width*2,theta1=-alpha,theta2=0,color="black",label=rf"$\alpha$={alpha:.1f}")
        x,y=calculate_arc(center=(xt,yt),radius=20,start_angle=0,end_angle=alpha)
        self.arc_alpha,=ax.plot(x,y,color="black")
        #ax.add_patch(self.arc_alpha)

    def init_force_vector(self,zeta:float,**kwargs):
        #allgemeines verzahnungsgesetz: Kraft muss immer durch den Wälzpunkt verlaufeb
        default_arrow_props = self.arrowprops
        gear=self.gear
        ax=self.ax
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

    def init_speed_vector(self,zeta,**kwargs):
        #allgemeines verzahnungsgesetz: Kraft muss immer durch den Wälzpunkt verlaufeb
        default_arrow_props = self.arrowprops
        default_arrow_props2 = {
            'length_includes_head': True,
            'head_width':2,
            'head_length':3
        }
        gear=self.gear
        ax=self.ax
        w_length=kwargs.pop('w_length',1)
        vt1,vt2,vg=gear.v_gear(zeta,w_length)
        xc,yc=gear.p_pOA(zeta)
        # vt1_arrow=ax.arrow(xc-vt1[0],yc-vt1[1],*vt1,label="$v_{t1}$",color='gold',**default_arrow_props)
        vt1_arrow=ax.arrow(xc-vt1[0],yc-vt1[1],*vt1,label="$v_{t1}$",color='gold',head_width=1,head_length=2,length_includes_head=True)
        vt2_arrow=ax.arrow(xc-vt2[0],yc-vt2[1],*vt2,label="$v_{t2}$",color='salmon',**default_arrow_props)
        vg_arrow=ax.arrow(xc-vg[0],yc-vg[1],*vg,label="$v_g$",color='goldenrod',**default_arrow_props2)
        return vt1_arrow,vt2_arrow,vg_arrow

    def addLines(self,zeta: float=0,**kwargs):
        gear=self.gear
        ax=self.ax
        alpha=gear.alpha_t(zeta)
        length=kwargs.get("length",60)
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
    
    def update_lines(self,zeta: float,**kwargs):
        gear=self.gear
        alpha=gear.alpha_t(zeta)
        length=kwargs.get("length",90)
        C=gear.p_pOA(zeta)
        dt=np.array((-sin(alpha),cos(alpha)))*length*-1
        dn=np.array((cos(alpha),sin(alpha)))*length
        P1_t=C-dt/2
        P1_n=C-dn/2
        self.tangent.set_data(x=P1_t[0],y=P1_t[1],dx=dt[0],dy=dt[1])
        self.normal.set_data(x=P1_n[0],y=P1_n[1],dx=dn[0],dy=dn[1])
        self.text_normal.set_position(P1_n)
        self.text_tangent.set_position(P1_t)
    
    def setup_pressure_angle_plot(self,ax,gear,zeta:float=0):
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
        ax.plot(angles*180/pi,alphas,label="α_t",color="black")
        ax.set_xlabel("Rotation ζ °")
        ax.text(0.5,0.9,r"Transverse pressure angle $\alpha_t$ [°]",
            horizontalalignment='center', 
            verticalalignment='top', 
            transform=ax.transAxes)
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        alpha_zeta=gear.alpha_t(zeta)*180/pi

        yinterp = np.interp(-alpha_zeta, -alphas, lengths)
        print(f"marker_pressureAngle={alpha_zeta:.1f}° and length={yinterp:.1f}mm and zeta={zeta*180/pi:.1f}°")
        self.marker_pressureAngle,=ax.plot(angles[0]*180/pi,alphas[0],marker='o',color="black",linestyle="",linewidth=4)
        
    def setup_sliding_velocity_plot(self,ax,gear,zeta:float=0,**kwargs):
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
        print(f"marker_pressureAngle={alpha_zeta:.1f}° and length={yinterp:.1f}mm and zeta={zeta*180/pi:.1f}°")
        self.marker_2,=ax.plot(lengths[0],alphas[0],marker='o',color="black",linestyle="",linewidth=4)
    
    def add_lines(self,zeta: float=0,**kwargs):
        ax=self.ax
        gear=self.gear
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
    # Additional methods: `init_force_vector`, `init_speed_vector`, `add_lines`, `update_lines`, `update_vectors`

if __name__ == "__main__":
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
    gear_meshing_plot = GearMeshingPlot(gear_pair)
