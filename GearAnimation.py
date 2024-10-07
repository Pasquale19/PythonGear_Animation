import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Circle, FancyArrow,Arc,Polygon
from matplotlib.transforms import Affine2D
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider
import numpy as np
from Gears.BaseGear2 import InvoluteGear
from math import pi,sin,cos,asin,acos,atan,tan,atan2,hypot,tau
from HideOnClick2 import hideOnClick
from Utilities.Helper import calculate_arcPolygon
from matplotlib.animation import FuncAnimation
import time
#Styles
mpl.rcParams['lines.linewidth'] = 2 #documentation https://matplotlib.org/stable/api/matplotlib_configuration_api.html
plt.rcParams['axes.spines.top']=True
plt.rcParams['axes.spines.right']=True


class GearMeshingPlot:
    def __init__(self, gear:InvoluteGear):
        self.gear = gear
        self.fig = None
        self.ax = None
        self.ax_pressure_angle = None
        self.ax_sliding_velocity = None
        self.slider = None
        self.vis_slider=True
        self.angle_text = None
        self.phi_A=gear.phi_A


        # Plot elements to be updated
        self.polyPatch_Gear1 = None
        self.polyPatch_Gear2 = None
        self.Fn_vec = None
        self.Fr_vec = None
        self.arrow_Fres=None
        self.vt1_arrow = None
        self.vt2_arrow = None
        self.vg_arrow = None
        self.tangent = None
        self.normal = None
        self.common_tangent=None
        self.text_tangent = None
        self.text_normal = None
        self.text_common_tangent=None
        self.PC = None
        self.p_poa_arrow = None
        self.p_poa1 = None
        self.marker_2=None
        self.marker_pressureAngle=None
        self.text_pressureAngle=None
        self.arc_alpha=None
        self.text_alpha=None

        self.circle_F1=None
        self.marker_K=None

        #grouping elements
        self.group_lines=[self.tangent,self.normal,self.common_tangent,
            self.text_tangent,self.text_normal,self.text_common_tangent,self.arc_alpha,self.text_alpha]
        self.group_forces=[ self.Fn_vec,self.Fr_vec,self.arrow_Fres]
        self.group_speeds=[ self.vt1_arrow,self.vt2_arrow,self.vg_arrow]

        #style elements
        self.limits0=None
        self.arrowprops=dict(head_length=0.7, head_width=0.3, lw=1,length_includes_head=True)
        self.angle_radius=25


        #animation stuff
        self.frames_per_anim=100
        self.fps=24
        self.framenumber=0
        #self.setup_plot()

    def setup_plot(self):
        # Create figure and subplots
        gear=self.gear
        self.fig = plt.figure(figsize=(19.2, 10.8), dpi=100,num="MeshingAnimation")
        self.ax = self.fig.add_subplot(111)

        height_main_plot=gear.center2[1]+gear.ra2*0.6
        width=height_main_plot*19.2/10.8
        width_main_plot=max(gear.ra2,gear.ra1)*2+2

        #self.ax.set_xlim(-width_main_plot/2, -width_main_plot/2+width)
        lb_y=-5
        #self.ax.set_ylim(lb_y,lb_y+height_main_plot)
        margin=0.1
        width_subPlot=(1-width_main_plot/width-margin)
        height_subPlot=1/3
        xleft=width_main_plot/width+margin
        #self.ax_pressure_angle=self.ax.inset_axes([xleft,1-(margin+height_subPlot),width_subPlot,height_subPlot])
        #self.ax_sliding_velocity=self.ax.inset_axes([xleft,margin,width_subPlot,height_subPlot])
        # self.ax = self.fig.add_subplot(gs[:, 0])
        # self.ax_pressure_angle = self.fig.add_subplot(gs[0, 1])
        # self.ax_sliding_velocity = self.fig.add_subplot(gs[1, 1])


        zeta0=self.phi_A
        # Initial plot configuration
        self.plotPoints(ax=self.ax,gear=self.gear)
        self.setup_slider()
        self.plot_initial_gear_profiles()
        #self.setup_pressure_angle_plot(ax=self.ax_pressure_angle,gear=self.gear)
        # self.setup_sliding_velocity_plot(ax=self.ax_sliding_velocity,gear=self.gear,zeta=zeta0)
        self.initialize_vectors()
        self.init_angle(zeta0)
        #self.plotExtras()
        # Connect event handlers
        self.fig.canvas.mpl_connect('key_press_event', self.key_press)

        # Style the plot
        
        ax=self.ax
        ax.set_position([0.05, 0.05, 0.9, 0.9])
   
        self.limits0=(*ax.get_xlim(),*ax.get_ylim())
        #self.ax.axis('off') #hide axis
        ax.set_aspect('equal', 'box')
        self.title=ax.text((width_main_plot/2)/width, 1, 'Eccentric Cycloid Gear',transform=ax.transAxes,horizontalalignment='center',
            verticalalignment='top',fontsize=16,
            bbox=dict(facecolor='white', edgecolor='white', alpha=1))
        parameter_list:str=f"$z_1$={gear.z1}" "\t" f"$z_2$={gear.z2}" "\n" f"$x_1$={gear.x1:.2f}" "\t" f"$x_2$={gear.x2:.2f}" "\n" r"$\alpha$" f"={gear.alpha*180/pi:.1f}° ""\t\t" r"$\alpha_w$=" f"{gear.alpha_w*180/pi:.1f}°" "\n" r"$a$=" f"{gear.a:.1f}" "\n" r"$\varphi_A$" f"={gear.phi_A*180/pi:.1f}°" "\t" r"$\varphi_E$" f"={gear.phi_E*180/pi:.1f}°" "\n" r"$h_{ap1}$" f"={gear.hap1:.2f} \t" r"$h_{ap2}$" f"={gear.hap2:.2f}"  r"$\epsilon_{\alpha}=$" f"{gear.e_alpha:.2f}"
        self.parameter_list=ax.text(0.02, 0.05, parameter_list,transform=ax.transAxes,horizontalalignment='left',
            verticalalignment='bottom',fontsize=10, animated=False,
            bbox=dict(facecolor='yellow', edgecolor='red', alpha=0.8))
        ax.legend(loc='upper left',bbox_to_anchor=(0.0,1),framealpha=0, prop={'size': 10})
        ax.grid(True)
        hideOnClick(self.fig,self.ax)
        #self.createAnimation(frames_per_anim=10,fps=1)
        
        
        plt.show()
        #plt.close()
    


    def plotPoints(self,ax,gear,fontsize:int=12,color='black'):
        Ax,Ay=gear.A
        ax.annotate('A',(Ax,Ay),color=color, fontsize=fontsize)
        line_A=ax.plot([Ax],[Ay],c=color,label=f"A  start of meshing",marker='o')
        
        
        Ex,Ey=gear.E
        ax.plot([Ex],[Ey],c=color,label=f"E end of meshing",marker='o')
        ax.text(Ex,Ey,'E',c=color, fontsize=fontsize,ha='right',va='bottom')
        Cx,Cy=gear.C
        ax.plot([Cx],[Cy],c="magenta",label=f"C",marker='o')
        ax.text(0-0.01,Cy,'C',c=color, fontsize=fontsize,ha='center',va='top')
        pointLabels="$O_1$","$O_2$"
        ax.scatter(x=0,y=gear.a,c='g')
        ax.text(0,gear.a+1,pointLabels[1],c='g', fontsize=fontsize,clip_on = True)
        ax.scatter(x=0,y=0,c='b')
        ax.text(0,0+1,pointLabels[0],c='b', fontsize=fontsize,clip_on = True)
        E0x,E0y=gear.E0
        ax.plot([E0x],[E0y],c=color,label=f"$T_2$",marker='o')
        ax.annotate("$T_2$",(E0x,E0y),color=color, fontsize=fontsize)
        A0x,A0y=gear.A0
        ax.plot([A0x],[A0y],c=color,label=f"$T_1$",marker='o')
        ax.annotate("$T_1$",(A0x,A0y),color=color, fontsize=fontsize)

    def plot_initial_gear_profiles(self):
        # Plot initial points and gear profiles
        self.plot_reference_circles()
        self.plot_gear_polygons()
        self.plot_contact_path()

        self.angle_text = self.ax.text(0.05, 0.95, '', transform=self.ax.transAxes, verticalalignment='top')
        self.angle_text.set_visible(False)

    def plot_reference_circles(self):
        # Plot reference circles
        rw1_circle = Circle((0,0), self.gear.rw1, linestyle='-.', color='g', fill=False,label=f"rw1={self.gear.rw1:.2f}")
        rw2_circle = Circle(self.gear.center2, self.gear.rw2, linestyle='-.', color='b', fill=False,label=f"rw2={self.gear.rw2:.2f}")
        da1_circle = Circle((0, 0), self.gear.ra1 , linestyle=':', color='b', fill=False,label=f"ra1={self.gear.ra1:.2f}")
        da2_circle = Circle(self.gear.center2, self.gear.ra2, linestyle=':', color='g', fill=False,label=f"ra2={self.gear.ra2:.2f}")
        db1_circle = Circle((0, 0), self.gear.rb1 , linestyle='--', color='b', fill=False,label=f"rb1={self.gear.rb1:.2f}")
        db2_circle = Circle(self.gear.center2, self.gear.rb2, linestyle='--', color='g', fill=False,label=f"rb2={self.gear.rb2:.2f}")
        df1_circle = Circle((0, 0), self.gear.rf1 , linestyle='--', color='b', fill=False,label=f"rf1={self.gear.rf1:.2f}")
        df2_circle = Circle(self.gear.center2, self.gear.rf2, linestyle='--', color='g', fill=False,label=f"rf2={self.gear.rf2:.2f}")

        r2_circle = Circle(self.gear.center2, self.gear.r2, linestyle='-.', color='g', fill=False,label=f"$r_2$={self.gear.r2:.2f}")
        r1_circle = Circle((0, 0), self.gear.r1, linestyle='-.', color='b', fill=False,label=f"$r_1$={self.gear.r1:.2f}")
        for circle in [rw2_circle, rw1_circle, da1_circle, da2_circle,db1_circle,db2_circle,df1_circle,df2_circle,r1_circle,r2_circle]:
            self.ax.add_patch(circle)

    def plot_gear_polygons(self):
        # Initial rotation
        
        initial_rotation = pi/2-1/4*tau/self.gear.z1
        initial_rotation2=pi/2*np.sign(self.gear.i)+1/4*tau/self.gear.z2*np.sign(self.gear.i)
        #psi=-pi/2-1/4*tau/self.gear.z2
        # Plot Gear1
        gear1_x, gear1_y = self.gear.gearGeometry1( num=1000,rotation_angle=initial_rotation-self.gear.phi_fuss*1)
        self.polyPatch_Gear1 = mpl.patches.Polygon(list(zip(gear1_x, gear1_y)), facecolor='lightblue', edgecolor='b',
                                                   linestyle='-', alpha=0.5, linewidth=2,label="gear1")
        rotation = Affine2D().rotate_around(0, 0, theta=self.gear.phi_fuss*1)
        self.polyPatch_Gear1.set_transform(rotation + self.ax.transData)
        self.ax.add_patch(self.polyPatch_Gear1)

        # Plot Gear2
        gear2_x, gear2_y = self.gear.gearGeometry2( num=1000,rotation_angle= initial_rotation2)
        self.polyPatch_Gear2 = mpl.patches.Polygon(list(zip(gear2_x, gear2_y)), facecolor='lightgreen',
                                                   edgecolor='g', linestyle='-', alpha=0.5, linewidth=2)
        rotation2 = Affine2D().rotate_around(*self.gear.center2 ,0)
        self.polyPatch_Gear2.set_transform(rotation2 + self.ax.transData)
        self.ax.add_patch(self.polyPatch_Gear2)

    def plot_contact_path(self,**kwargs):
        # Get path of contact
        contact_x, contact_y = self.gear.calculate_path_of_contact()
        self.ax.plot(contact_x, contact_y, color='black', linestyle='-', label='Path of Contact')

    # def add_lines(self,zeta: float=0,**kwargs):
    #     ax=self.ax
    #     gear=self.gear
    #     alpha=gear.alpha_w
    #     length=kwargs.get("length",90)
    #     props={"color":"black","linestyle":"--","linewidth":0.5,"alpha":0.8}
    #     C=self.gear.C
    #     dt=np.array((-sin(alpha),cos(alpha)))*length*-1
    #     dn=np.array((cos(alpha),sin(alpha)))*length*-1
    #     P1_t=C-dt/2
    #     P1_n=C-dn/2
    #     tangent=ax.arrow(x=P1_t[0],y=P1_t[1],dx=dt[0],dy=dt[1],**props,clip_on=True)
    #     normal=ax.arrow(x=P1_n[0],y=P1_n[1],dx=dn[0],dy=dn[1],**props,clip_on=True)
    #     text_normal=ax.text(*P1_n,"flank normal",va="bottom",ha="left",fontsize=8,clip_on=True)
    #     text_tangent=ax.text(*P1_t,"flank tangent",va="top",ha="left",fontsize=8,rotation=0,clip_on=True)
    #     dct=np.array((cos(0),sin(0)))*length
    #     P1_ct=np.array((0,gear.r2))-2*dct/3
    #     self.common_tangent=ax.arrow(x=P1_ct[0],y=P1_ct[1],dx=dct[0],dy=dct[1],**props,clip_on=True)
    #     self.text_common_tangent=ax.text(*P1_ct,"common tangent",va="bottom",ha="left",fontsize=8,clip_on=True)
    #     return tangent,normal,text_tangent,text_normal

    def initialize_vectors(self,phi0:float=0):
        phi0 = self.gear.phi
        self.Fn_vec, self.Fr_vec = self.init_force_vector(phi0)
        self.vt1_arrow, self.vt2_arrow, self.vg_arrow = self.init_speed_vector(phi0)

    def init_angle(self,zeta:float,fontsize=10):
        zeta0=zeta
        ax=self.ax
        xt,yt=self.gear.C
        Ex,Ey=self.gear.E
        EC=self.gear.E-self.gear.C
        radius=hypot(EC[1],EC[0])*1.5
        self.text_alpha=ax.text(Ex,Ey,r"$\alpha$",va="top",ha="right",fontsize=fontsize,color="black")
        alpha=self.gear.alpha_w
        #print(f"alpha0={alpha:.1f}")
        #self.arc_alpha=Arc((xt,yt),width=width*2,height=width*2,theta1=-alpha,theta2=0,color="black",label=rf"$\alpha$={alpha:.1f}")
        xy=calculate_arcPolygon(center=(xt,yt),radius=radius,start_angle=pi-alpha,end_angle=pi,use_radians=True)
        self.arc_alpha=Polygon(xy, closed=True,fill="grey",facecolor="grey",edgecolor="black",alpha=0.3)
        ax.add_patch(self.arc_alpha)
        

        #rotations vectors
        r1=self.gear.rw1/2
        angle=60
        xy=calculate_arcPolygon(center=(0,0),radius=r1,start_angle=90-angle/2,end_angle=90+angle/2,use_radians=False,includeCenter=False)
        self.arc_w1=Polygon(xy, closed=False,fill=False,facecolor="none",edgecolor="black",alpha=1)
        ax.add_patch(self.arc_w1)
        ax.text(0,r1,r"$\omega_1$",va="top",ha="center",fontsize=fontsize,color="black",bbox=dict(facecolor='white', edgecolor='white', alpha=1))
        xy=calculate_arcPolygon(center=self.gear.center2,radius=r1,start_angle=270-angle/2,end_angle=270+angle/2,use_radians=False,includeCenter=False)
        self.arc_w2=Polygon(xy, closed=False,fill="black",facecolor="none",edgecolor="black",alpha=1)
        ax.add_patch(self.arc_w2)
        self.text_w2=ax.text(0,self.gear.a-r1,r"$\omega_2$",va="bottom",ha="center",fontsize=fontsize,color="black",bbox=dict(facecolor='white', edgecolor='white', alpha=1))
        #Pfeilspitzen
        Pfeil1=np.array((np.sin(angle/2*pi/180),np.cos(angle/2*pi/180)))*r1
        Pfeil2=-np.array((np.sin(angle/2*pi/180),np.cos(angle/2*pi/180)))*r1+self.gear.center2
        ax.plot(*Pfeil1,marker=(3, 0, -angle/2), markersize=6, linestyle='None',color="black")
        ax.plot(*Pfeil2,marker=(3, 0, -angle/2+180), markersize=6, linestyle='None',color="black")

    def init_force_vector(self,zeta:float,**kwargs):
        #allgemeines verzahnungsgesetz: Kraft muss immer durch den Wälzpunkt verlaufeb
        default_arrow_props = self.arrowprops
        gear=self.gear
        ax=self.ax
        M_length=kwargs.pop('M_length',0.5)
        mu=kwargs.pop('mu',0.20) 
        alpha=gear.alpha_w
        xc,yc=gear.P
        F1_length=M_length/cos(alpha)*gear.r1
        dn=np.array((cos(alpha),sin(alpha)))*F1_length
        Fn_arrow=ax.arrow(xc-dn[0],yc-dn[1],*dn,label="$F_n$",color='orangered',**default_arrow_props)
        dr=-np.array((sin(alpha),-cos(alpha)))*F1_length*mu
        #print("Fr not correct at init")
        Fr_arrow=ax.arrow(xc-dr[0],yc-dr[1],*dr,label="$F_r$",color='r',**default_arrow_props)
        dF2=-dn-dr
        self.arrow_Fres=ax.arrow(xc-dF2[0],yc-dF2[1],*dF2,label="$F_{res}$",color='magenta',**default_arrow_props)
        return Fn_arrow,Fr_arrow

    def init_speed_vector(self,zeta,**kwargs):
        #allgemeines verzahnungsgesetz: Kraft muss immer durch den Wälzpunkt verlaufeb
        default_arrow_props = self.arrowprops
        gear=self.gear
        ax=self.ax
        gear.calculate_state(self.slider.val/180*pi)
        w_length=kwargs.pop('w_length',1)
        vt1,vt2,vg=gear.vt1,gear.vt2,gear.vg
        xc,yc=gear.P
        # vt1_arrow=ax.arrow(xc-vt1[0],yc-vt1[1],*vt1,label="$v_{t1}$",color='gold',**default_arrow_props)
        vt1_arrow=ax.arrow(xc-vt1[0],yc-vt1[1],*vt1,color='gold',label="$v_{t1}$",**default_arrow_props)
        vt2_arrow=ax.arrow(xc-vt2[0],yc-vt2[1],*vt2,color='salmon',label="$v_{t2}$",**default_arrow_props)
        vt1_arrow.set_visible(True)
        vt2_arrow.set_visible(True)
        vg_arrow=ax.arrow(xc-vg[0],yc-vg[1],*vg,label="$v_g$",color='goldenrod',**default_arrow_props)
        return vt1_arrow,vt2_arrow,vg_arrow




    def setup_slider(self):
        phi_A = self.phi_A*180/pi
        ax_slider = plt.axes([0.2, 0.05, 0.3, 0.02]) # [left, bottom, width, height]
        #self.slider = Slider(ax_slider, 'Rotation', phi_A, self.gear.phi_E*180/pi, valinit=phi_A,valstep=0.5)
        self.slider = Slider(ax_slider, 'Rotation', -38, 38, valinit=phi_A*0,valstep=0.01)
        self.slider.on_changed(self.update_plot)
    



    #region update functions
    def update_plot(self, phi:float):
        artists = []
        state=self.gear.calculate_state(phi/180*pi)
        vt1=state["vt1"]
        vt2=state["vt2"]
        vg=state["vg"]
        P=state["P"]
        M_res=state["M_res"]
        F_res=state["F_res"]
        Fr2=state["Fr2"]
        Fn=state.pop("Fn",0)
        r_p=hypot(P[1],P[0])
        phi_P=atan2(P[1],P[0])
        #print(f"phi_P={phi_P*180/pi-90:.1f}°")
        delta_i=self.gear.delta_i(r_p)
        phi_gear=phi/180*pi+delta_i
        psi=phi_gear/self.gear.i
        psi=(phi_gear-self.gear.phi_fuss)/self.gear.i
        rotation = Affine2D().rotate_around(0, 0, phi_gear)
        self.polyPatch_Gear1.set_transform(rotation + self.ax.transData)
        artists.append(self.polyPatch_Gear1)
        rotation2 = Affine2D().rotate_around(0, self.gear.a, psi)
        self.polyPatch_Gear2.set_transform(rotation2 + self.ax.transData)
        artists.append(self.polyPatch_Gear2)

        
        self.title.set_text(r"$\varphi_{Gear}$=" f"{phi_gear*180/pi:.1f}\n" r"$\delta_{i}$=" r"$\varphi$=" f"{phi:.1f}\n" f"{delta_i*180/pi:.1f}\n $\psi_1$={psi*180/pi:.1f}")
        self.text_w2.set_text(r"$M_2$="f"{M_res:.2f}\n $M_1$={self.gear.M1:.1f}\nvs expected $M_2$={self.gear.M1*self.gear.i:.1f}")
        # Update contact point
        # p_pOA = self.gear.p_pOA(zeta)
        # self.PC.set_data([p_pOA[0]], [p_pOA[1]])
        # artists.append(self.PC)
        #self.PC.set_label(f"C = {zeta * 180 / np.pi:.1f}°")
        self.angle_text.set_text(r"$\varphi$" f"={phi_gear * 180 / np.pi:.1f}°")
        artists.append(self.angle_text)
        # Update lines and vectors
        # artists.extend(self.update_lines(zeta,alpha=alpha_rad,p_pOA=p_pOA))
        artists.extend(self.update_vectors(phi,alpha=self.gear.alpha_w,P=P,vt1=vt1,vt2=vt2,vg=vg,F_res=F_res,Fn=Fn,Fr2=Fr2))
        # artists.extend(self.update_pressureAnglePlot(zeta,alpha=alpha_deg))
        # artists.extend(self.update_SlidingVelocityPlot(zeta))
        # artist_alpha,=self.updateAlpha(alpha=float(alpha_deg))
        #self.updateExtras(zeta)
        #artists.append(artist_alpha)
        #self.fig.savefig(f"export/frame{self.framenumber:06d}.png", format='png', dpi=120)
        return artists
    
    def update_pressureAnglePlot(self,zeta:float, alpha:float): 
        self.marker_pressureAngle.set_data([zeta*180/pi],[alpha])
        self.text_pressureAngle.set_text(r"pressure angle  $\alpha$" f"={alpha:.1f}°")
        return self.marker_pressureAngle,

    def update_SlidingVelocityPlot(self, zeta):
        m=self.marker_2
        Kg=self.gear.Kg(zeta)
        m.set_data([zeta*180/pi],[Kg])
        return m,

    def update_vectors(self,zeta:float,alpha:float,P,**kwargs):
        gear=self.gear
        vt1_arrow,vt2_arrow,vg_arrow=self.vt1_arrow,self.vt2_arrow,self.vg_arrow
        if zeta==0:
            zeta=0.001
        vt1=kwargs.get("vt1",np.array((0,0)))
        vt2=kwargs.get("vt1",np.array((0,0)))
        vg=kwargs.get("vg",np.array((0,0)))
        
        #vt1,vt2,vg=gear.v_gear(zeta,w_length)
        xc,yc=P

        if vt1_arrow:
            vt1_arrow.set_data(x=xc-vt1[0],y=yc-vt1[1],dx=vt1[0],dy=vt1[1])
        if vt2_arrow:
            vt2_arrow.set_data(x=xc-vt2[0],y=yc-vt2[1],dx=vt2[0],dy=vt2[1])
        if vg_arrow:
            vg_arrow.set_data(x=xc-vg[0],y=yc-vg[1],dx=vg[0],dy=vg[1])

        #Forces
        M_length=kwargs.pop('M_length',0.5)
        mu=kwargs.pop('mu',0.20)
        F1_length=M_length/cos(-alpha)*gear.r1
        Fn_arrow,Fr_arrow=self.Fn_vec,self.Fr_vec
        Fn=kwargs.get("Fn",0)
        Fr2=kwargs.get("Fr2",0)
        F_res=kwargs.pop("F_res",0)
        if Fn_arrow:
            dn=np.array((-cos(alpha),sin(alpha)))*F1_length
            #print(f"dn={dn}\tFn{Fn}")
            dn=Fn
            Fn_arrow.set_data(x=xc-dn[0],y=yc-dn[1],dx=dn[0],dy=dn[1])
        if Fr_arrow:
            dr=np.array((dn[1],-dn[0]))*mu
            #print(f"dr={dr}\tFr2{Fr2}")
            dr=Fr2
            Fr_arrow.set_data(x=xc-dr[0],y=yc-dr[1],dx=dr[0],dy=dr[1])
        d_res=dn+dr
        self.arrow_Fres.set_data(x=xc-d_res[0],y=yc-d_res[1],dx=d_res[0],dy=d_res[1])
        artists=[self.arrow_Fres,self.vg_arrow,self.vt1_arrow,self.vt2_arrow,self.Fn_vec,self.Fr_vec]
        return artists

    def updateAlpha(self,alpha:float):
        angle_alpha=self.arc_alpha
        radius=self.angle_radius
        xy=calculate_arcPolygon(center=(0,self.gear.r2),radius=radius,start_angle=0,end_angle=-alpha)
        angle_alpha.set_xy(xy)
        return angle_alpha,

    # def update_lines(self,zeta: float,alpha:float,p_pOA,**kwargs):
    #     length=kwargs.get("length",90)
    #     C=p_pOA
    #     s_alpha,c_alpha=sin(alpha),cos(alpha)
    #     dt=np.array((-s_alpha,-c_alpha))*length*-1
    #     dn=np.array((-c_alpha,s_alpha))*length
    #     P1_t=C-dt/2
    #     P1_n=C-dn/2
        
    #     self.tangent.set_data(x=P1_t[0],y=P1_t[1],dx=dt[0],dy=dt[1])
    #     self.normal.set_data(x=P1_n[0],y=P1_n[1],dx=dn[0],dy=dn[1])
    #     self.text_normal.set_position(P1_n)
    #     self.text_tangent.set_position(P1_t)
    #     artists=[self.tangent,self.normal,self.text_normal,self.text_tangent]
    #     return artists
    #endregion
   
    def key_press(self, event,**kwargs):
        if event==None:
            key=key=kwargs.get("key",'2')
        else: key=event.key
        
        if key == 'v':
            vis_v = not self.vg_arrow.get_visible()
            self.showSpeeds(vis_v)
            if vis_v:
                self.title.set_text("sliding velocity $v_g$")   
        if key == 'n':
            F_vec=self.Fn_vec,self.Fr_vec
            vis_v = not F_vec[0].get_visible()
            self.showForces(vis_v)
            if vis_v:
                self.title.set_text("Forces")            
        if key == '2':
            lines=self.tangent,self.normal,self.text_tangent,self.text_normal,self.common_tangent,self.text_common_tangent,self.text_alpha,self.arc_alpha
            vis_v = not lines[0].get_visible()
            self.showLines(vis_v)
        if key=="5":
            self.ax_pressure_angle.set_visible(False)
            self.ax_sliding_velocity.set_visible(False)
        if key=="3":
            self.showSlider(not self.vis_slider)
        if key=="z":
            self.zoomIn(1)
        if key=="r":
            self.zoomIn(0)
        if key=="6":
            self.Animation1(True)
        if key=="7":
            self.Animation2(True)
        if key=="8":
            self.Animation3(True)
        if key=="9":
            self.Animation4(True)
        self.fig.canvas.draw()

    #region Setup side Plots
    def setup_pressure_angle_plot(self,ax,gear,zeta:float=0,**kwargs):
        default_arrow_props = self.arrowprops
        gear=self.gear
        M_length=kwargs.pop('M_length',gear.M1)
        mu=kwargs.pop('mu',gear.mu) 
        alpha=gear.alpha_w
        xc,yc=gear.P
        xc=0
        yc=0
        F1_length=M_length/cos(alpha)*gear.r1
        dn=np.array((-cos(alpha),sin(alpha)))*F1_length
        self.Fn_arrow_=ax.arrow(xc-dn[0],yc-dn[1],*dn,label="$F_n$",color='orangered',**default_arrow_props)
        dr=np.array((-sin(alpha),-cos(alpha)))*F1_length*mu
        #print("Fr not correct at init")
        self.Fr_arrow_=ax.arrow(xc-dr[0],yc-dr[1],*dr,label="$F_r$",color='r',**default_arrow_props)
        dF2=dn+dr
        self.arrow_Fres_=ax.arrow(xc-dF2[0],yc-dF2[1],*dF2,label="$F_res$",color='magenta',**default_arrow_props)

        self.text_pressureAngle=ax.text(0.5,0.9,r"pressure angle $\alpha$",
            horizontalalignment='center', 
            verticalalignment='top', 
            transform=ax.transAxes)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        
    def setup_sliding_velocity_plot(self,ax,gear,zeta:float=0,**kwargs):
        zetaA=kwargs.get("start",gear.zetaA)
        zetaE=kwargs.get("end",gear.zetaE)
        angles=np.linspace(zetaA,zetaE,num=400)
        
        #xis=np.array([(gear.xi(a))*180/pi for a in angles])
        Kgs=np.array([(gear.Kg(a)) for a in angles])
        
        #x_coords,y_coords=gear.calculate_path_of_contact(start=zetaA,end=zetaE,num_points=len(angles))
        # lengths = np.zeros(len(angles))
        # for i in range(1, len(angles)):
        #     lengths[i] = lengths[i-1] + np.sqrt((x_coords[i] - x_coords[i-1])**2 + (y_coords[i] - y_coords[i-1])**2)
        ax.plot(angles*180/pi,Kgs,label="sliding factor $K_g$",color="black")
        ax.set_ylabel("sliding factor $K_g$")
        ax.set_xlabel("Rotation ζ in °")
        ax.text(0.5,0.9,"sliding factor $K_g$",
            horizontalalignment='center', 
            verticalalignment='top', 
            transform=ax.transAxes)
        ax.grid(axis='y')
        self.marker_2,=ax.plot(angles[0]*180/pi,Kgs[0],marker='o',color="black",linestyle="",linewidth=4)
    
    #endregion
    
    #region Visibilty stuff
    def get_lines(self):
        self.group_lines=[self.tangent,self.normal,self.common_tangent,
            self.text_tangent,self.text_normal,self.text_common_tangent,self.arc_alpha,self.text_alpha]
        return self.group_lines
    
    def showAll(self,show:bool=False):
        groups=self.get_lines()
      
        group_speeds=[ self.vt1_arrow,self.vt2_arrow,self.vg_arrow]
        group_forces=[self.Fn_vec,self.Fr_vec,self.arrow_Fres]
        groups=[*groups,*group_forces,*group_speeds,self.ax_pressure_angle,self.ax_sliding_velocity]
        #groups.extend(*self.group_forces,*self.group_speeds,*self.ax_pressure_angle,*self.ax_sliding_velocity)
        for i,item in enumerate(groups):
            item.set_visible(show)
        self.showLines(False)


    def showLines(self,show:bool=False):
        group=self.get_lines()
        group.extend((self.ax_pressure_angle,self.ax_sliding_velocity))
        for i,item in enumerate(group):
            item.set_visible(show)

    def showSpeeds(self,show:bool=False):
        group_speeds=[ self.vt1_arrow,self.vt2_arrow,self.vg_arrow]
        group_speeds.append(self.ax_sliding_velocity)
        for i,item in enumerate(group_speeds):
            item.set_visible(show)
        for v in (self.vt1_arrow,self.vt2_arrow):
            v.set_visible(False)
    def showForces(self,show:bool=False):
        group_forces=[self.Fn_vec,self.Fr_vec,self.arrow_Fres]
        group_forces.append(self.ax_sliding_velocity)
        for i,item in enumerate(group_forces):
            item.set_visible(show)
    
    def showSlider(self,show:bool=False):
        self.vis_slider=show
        slider=self.slider
        if slider:
            if show:
                slider.ax.set_position([0.2, 0.05, 0.3, 0.02])
            else:
                #slider.valtext.set_visible(False)
                slider.ax.set_position([0.2, 1.1, 0.3, 0.02])
            #slider.ax.remove()  # Remove the slider axis from the figure
            
        else:
            return
    #endregion

    #region Animation stuff
    def animate(self,val,framenumber=None):
        start = time.perf_counter()
        self.slider.set_val(val)
        #self.update_plot(val)
        end = time.perf_counter()
        print(f"Angle {val:.2f}: Execution time: {end - start:.7f} seconds")
        # if framenumber:
        #     self.title.set_text("frame{self.framenumber}")
        # print(f"animate {val:.1f}°")

    def zoomIn(self,val:float=1):
        ax=self.ax
        center=np.array((0,self.gear.r2))
        height=self.gear.r1*2
        width=height*19.2/10.8
        xmin,xmax,ymin,ymax=self.limits0
        xmin_zoom,xmax_zoom=center[0]-height/2,center[0]-height/2+width
        ymin_zoom,ymax_zoom=center[1]-height/2,center[1]+height/2
        xmin_new=val*xmin_zoom+xmin*(1-val)
        ymin_new=val*ymin_zoom+ymin*(1-val)
        xmax_new=val*xmax_zoom+xmax*(1-val)
        ymax_new=val*ymax_zoom+ymax*(1-val)
        ax.set_xlim(xmin_new,xmax_new)
        ax.set_ylim(ymin_new,ymax_new)

    def Animation1(self,save:bool=False,**kwargs):  
        '''Animation1 full meshing'''
        print("Animation1 started")
        self.showSlider()
        frames_per_anim=kwargs.get("frames_per_anim",self.frames_per_anim)
        fps=kwargs.get("fps",self.fps)
        self.title.set_text("Eccentric Gear Mesh")
        self.showAll(False)
        angles1=np.linspace(self.slider.valmin+self.slider.valstep,self.slider.valmin*-1,num=frames_per_anim,endpoint=True)
        interval=1000/fps
        #print(f"interval:{interval} frames={frames_per_anim} fps={fps} duration={frames_per_anim/(60*fps)}min\t len={len(angles1)}")
        anim=FuncAnimation(self.fig,self.update_plot,interval=interval,frames=angles1,repeat=False)
        if save:
            writervideo = mpl.animation.FFMpegWriter(fps=fps)
            anim.save('animation1.mp4', writer=writervideo)
            print("animation1.mp4 saved")
    
    def Animation11(self,save:bool=False,**kwargs):  
        '''Animation11 full meshing'''
        frames_per_anim=kwargs.get("frames_per_anim",self.frames_per_anim)
        fps=kwargs.get("fps",self.fps)
        self.title.set_text("Eccentric Gear Mesh")
        self.showAll(False)
        angles1=np.linspace(self.slider.valmin+self.slider.valstep,self.slider.valmin*-1,num=frames_per_anim,endpoint=True)
        # Set up the writer
        writer = mpl.animation.FFMpegWriter(fps=fps)
        print("animation11")
        # Create and save the animation
        with writer.saving(self.fig, "animation11.mp4", dpi=100):
            for i in angles1:  # 100 frames
                self.slider.set_val(i)
                writer.grab_frame()
        print("animation11.mp4 saved")

    def Animation2(self,save:bool=False,**kwargs):  
        '''Animation2 lines and pressure angle Plot'''
        print("Animation2 started")
        frames_per_anim=kwargs.get("frames_per_anim",self.frames_per_anim)
        fps=kwargs.get("fps",self.fps)
        self.showAll(False)
        self.showLines(True)
        self.showForces(False)
        self.showSpeeds(False)
        self.ax_sliding_velocity.set_visible(False)
        angles=np.linspace(self.slider.valmin+self.slider.valstep*0,self.slider.valmax,num=frames_per_anim)
        anim=FuncAnimation(self.fig,self.update_plot,interval=1000/fps,frames=angles,repeat=False)
        if save:
            writervideo = mpl.animation.FFMpegWriter(fps=fps)
            anim.save('animation2.mp4', writer=writervideo)
            print("animation2.mp4 saved")
    def Animation3(self,save:bool=False,**kwargs):  
        '''Animation3 sliding velocity $v_g$'''
        print("Animation3 started")
        self.showSlider()
        frames_per_anim=kwargs.get("frames_per_anim",self.frames_per_anim)
        fps=kwargs.get("fps",self.fps)
        angles=np.linspace(self.slider.valmin+self.slider.valstep*0,self.slider.valmax,num=frames_per_anim)
        self.title.set_text("sliding velocity $v_g$")
        self.showLines(False)
        self.showForces(False)
        self.showSpeeds(True)
        self.zoomIn(1)
        anim=FuncAnimation(self.fig,self.update_plot,interval=1000/fps,frames=angles,repeat=False)
        if save:
            writervideo = mpl.animation.FFMpegWriter(fps=fps)
            anim.save('animation3.mp4', writer=writervideo)
            print("animation3.mp4 saved")

    def Animation4(self,save:bool=False,**kwargs):  
        '''Animation4 forces'''
        print("Animation4 started")
        self.showSlider()
        frames_per_anim=kwargs.get("frames_per_anim",self.frames_per_anim)
        fps=kwargs.get("fps",self.fps)
        angles=np.linspace(self.slider.valmin+self.slider.valstep*0,self.slider.valmax,num=frames_per_anim)
        self.title.set_text("Forces")
        
        self.showLines(False)
        self.showSpeeds(False)
        self.showForces(True)
        self.ax_pressure_angle.set_visible(True)
        self.ax_sliding_velocity.set_visible(False)
        self.zoomIn(1)
        anim=FuncAnimation(self.fig,self.update_plot,interval=1000/fps,frames=angles,repeat=False)
        if save:
            writervideo = mpl.animation.FFMpegWriter(fps=fps)
            anim.save('animation4.mp4', writer=writervideo)
            print("animation4.mp4 saved")
    #endregion

if __name__ == "__main__":
    # gear_pair = InvoluteGear(
    #     z1=32,
    #     z2=30,
    #     a=62,
    #     alpha=20*pi/180,
    #     hap1=1,hap2=1,
    #     hfp1=1.25,hfp2=1.25,
    #     x1=0,x2=0
    # )
    gear_pair = InvoluteGear(
        z1=16,
        z2=-32,
        m=1,
        alpha=25*pi/180,
        hap1=1,hap2=1,
        hfp1=1.25,hfp2=1.25,
        x1=0,x2=0
    )
    ###epsilon_alpha=2
    # gear_pair = InvoluteGear(
    #     z1=23,
    #     z2=40,
    #     m=1,
    #     alpha=20*pi/180,
    #     hap1=1.25,hap2=1.25,
    #     hfp1=1.4,hfp2=1.4,
    #     x1=0,x2=0
    # )
    gear_pair.print_Parameter()

    gear_meshing_plot = GearMeshingPlot(gear_pair)
    #gear_meshing_plot.gear.rF1=16.1
    gear_meshing_plot.setup_plot()
