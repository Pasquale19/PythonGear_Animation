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
import time
#Styles
mpl.rcParams['lines.linewidth'] = 2 #documentation https://matplotlib.org/stable/api/matplotlib_configuration_api.html
plt.rcParams['axes.spines.top']=False
plt.rcParams['axes.spines.right']=False

# Define your class here
class GearMeshingPlot:
    def __init__(self, gear:EccentricCycloidGear):
        self.gear = gear
        self.fig = None
        self.ax = None
        self.ax_pressure_angle = None
        self.ax_sliding_velocity = None
        self.slider = None
        self.vis_slider=True
        self.angle_text = None
        self.zeta0=gear.zetaA


        # Plot elements to be updated
        self.polyPatch_Gear1 = None
        self.polyPatch_Gear2 = None
        self.Fn_vec = None
        self.Fr_vec = None
        self.arrow_F2=None
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
        self.group_forces=[ self.Fn_vec,self.Fr_vec,self.arrow_F2]
        self.group_speeds=[ self.vt1_arrow,self.vt2_arrow,self.vg_arrow]

        #style elements
        self.limits0=None
        self.arrowprops=dict(head_length=1.5, head_width=1.2, lw=1,length_includes_head=True)
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

        height_main_plot=gear.a+gear.da1/2+5
        width=height_main_plot*19.2/10.8
        width_main_plot=gear.da2+2

        self.ax.set_xlim(-width_main_plot/2, -width_main_plot/2+width)
        lb_y=-5
        self.ax.set_ylim(lb_y,lb_y+height_main_plot)
        margin=0.1
        width_subPlot=(1-width_main_plot/width-margin)
        height_subPlot=1/3
        xleft=width_main_plot/width+margin
        self.ax_pressure_angle=self.ax.inset_axes([xleft,1-(margin+height_subPlot),width_subPlot,height_subPlot])
        self.ax_sliding_velocity=self.ax.inset_axes([xleft,margin,width_subPlot,height_subPlot])
        # self.ax = self.fig.add_subplot(gs[:, 0])
        # self.ax_pressure_angle = self.fig.add_subplot(gs[0, 1])
        # self.ax_sliding_velocity = self.fig.add_subplot(gs[1, 1])


        zeta0=self.zeta0
        # Initial plot configuration
        self.plotPoints(ax=self.ax,gear=self.gear)
        self.setup_slider()
        self.plot_initial_gear_profiles()
        self.setup_pressure_angle_plot(ax=self.ax_pressure_angle,gear=self.gear,zeta=zeta0)
        self.setup_sliding_velocity_plot(ax=self.ax_sliding_velocity,gear=self.gear,zeta=zeta0)
        
        self.init_angle(zeta0)
        #self.plotExtras()
        # Connect event handlers
        self.fig.canvas.mpl_connect('key_press_event', self.key_press)

        # Style the plot
        
        ax=self.ax
        ax.set_position([0.05, 0.05, 0.9, 0.9])
        ax.set_xlim(-width_main_plot/2, -width_main_plot/2+width)
        lb_y=-5
        ax.set_ylim(lb_y,lb_y+height_main_plot)
        self.limits0=(*ax.get_xlim(),*ax.get_ylim())
        self.ax.axis('off')
        ax.set_aspect('equal', 'box')
        self.title=ax.text((width_main_plot/2)/width, 1, 'Eccentric Cycloid Gear',transform=ax.transAxes,horizontalalignment='center',
            verticalalignment='top',fontsize=16,
            bbox=dict(facecolor='white', edgecolor='white', alpha=1))
        parameter_list:str=f"$z_1$={gear.z1}" "\t\t" f"$z_2$={gear.z2}" "\n" f"$r_1$={gear.r1:.2f}" "\t" f"$e$={gear.e:.2f}" "\n" r"s$\lambda$" f"={gear.lambda_}""\t\t" r"$\beta$=0°" "\n" r"$\varphi_{j1}$" f"={gear.phi_j1*180/pi:.1f}°" "\t" "$c*$" f"={gear.c_star:.1f}"
        self.parameter_list=ax.text(0.02, 0.05, parameter_list,transform=ax.transAxes,horizontalalignment='left',
            verticalalignment='bottom',fontsize=10, animated=False,
            bbox=dict(facecolor='yellow', edgecolor='red', alpha=0.8))
        ax.legend(loc='upper left',bbox_to_anchor=(0.0,1),framealpha=0, prop={'size': 10})
        hideOnClick(self.fig,self.ax)
        #self.createAnimation(frames_per_anim=10,fps=1)
        
        
        plt.show()
        #plt.close()
    
    def plotExtras(self):
        angle=pi/self.gear.z1+self.gear.zetaA*self.gear.i
        F10=np.array((0,self.gear.a))+np.array((sin(angle),-cos(angle)))*self.gear.qF1
        self._F10=F10
        print(f"F10={F10}")
        self.circle_F1=Circle(F10,radius=self.gear.rF1,linestyle=":",edgecolor="blue",facecolor="none",linewidth=2,alpha=1,label="F1")
        self.ax.add_patch(self.circle_F1)

        K0=self.gear.p_pOA(self.gear.zetaA)
        self._K0=K0
        self.marker_K,=self.ax.plot([K0[0]],[K0[1]],c="yellow",label=f"head",marker='o')
        self.slider.on_changed(self.updateExtras)


    def updateExtras(self,zeta:float):
        """zeta in rad"""
        zeta=zeta*pi/180
        from Utilities.rotatePoint import rotate_point
        rotation_angle2=zeta-self.gear.zetaA
        kappa=rotation_angle2*self.gear.i
        #kappa=zeta*self.gear.i
        #print(f"updateExtras {zeta*180/pi:.2f}°")
        rotation = Affine2D().rotate_around(0, self.gear.a, kappa)
        rotation2 = Affine2D().rotate_around(0, 0, -rotation_angle2)

        F1s=rotate_point(self._F10,kappa,center=(0,self.gear.a),use_radians=True)
        Ks=rotate_point(self._K0,-rotation_angle2,center=(0,0),use_radians=True)
        distance=-F1s+Ks
        dist=np.hypot(distance[0],distance[1])
        self.circle_F1.center=F1s[0],F1s[1]
        self.marker_K.set_data([Ks[0]],[Ks[1]])
        # self.circle_F1.set_transform(rotation + self.ax.transData)
        # self.marker_K.set_transform(rotation2 + self.ax.transData)
        info:str=f"K with distance {dist:.2f} vs rF1={self.gear.rF1:.2f}"
        print(info)
        self.marker_K.set_label(info)
        self.title.set_text(info)

    def plotPoints(self,ax,gear,fontsize:int=12,color='black'):
        A=gear.A
        if A:
            Ax,Ay=A
            ax.annotate('A',(Ax+0.3,Ay+1),color=color, fontsize=fontsize)
            line_A=ax.plot([Ax],[Ay],c=color,label=f"A  start of meshing",marker='o')
        ax.plot([0],[gear.r2],c=color,marker='o',label="C pitch point")
        ax.text(0-0.01,gear.r2-2,'C',c=color, fontsize=fontsize,ha='center',va='top')
        Ex,Ey=gear.E
        ax.plot([Ex],[Ey],c=color,label=f"E end of meshing",marker='o')
        ax.text(Ex-1,Ey-1.6,'E',c=color, fontsize=fontsize,ha='right',va='top')

        pointLabels="$O_1$","$O_2$"
        ax.scatter(x=0,y=gear.a,c='b')
        ax.text(0,gear.a+1,pointLabels[0],c='b', fontsize=fontsize,clip_on = True)
        ax.scatter(x=0,y=0,c='g')
        ax.text(0,0+1,pointLabels[1],c='g', fontsize=fontsize,clip_on = True)

    def plot_initial_gear_profiles(self):
        # Plot initial points and gear profiles
        self.plot_reference_circles()
        self.plot_gear_polygons()
        self.plot_contact_path()
        self.initialize_vectors()
        self.add_tangent_and_normal_lines()

        self.angle_text = self.ax.text(0.05, 0.95, '', transform=self.ax.transAxes, verticalalignment='top')
        self.angle_text.set_visible(False)

    def plot_reference_circles(self):
        # Plot reference circles
        r2_circle = Circle((0, 0), self.gear.r2, linestyle='-.', color='g', fill=False)
        r1_circle = Circle((0, self.gear.a), self.gear.r1, linestyle='-.', color='b', fill=False)
        rInter_circle = Circle((0, self.gear.a), self.gear.rInter, linestyle='-.', color='m', fill=False,label="$r_{inter}$"f"={self.gear.rInter:.2f}")
        da1_circle = Circle((0, self.gear.a), self.gear.da1 / 2, linestyle=':', color='b', fill=False)
        da2_circle = Circle((0, 0), self.gear.ra2, linestyle=':', color='g', fill=False)
        for circle in [r2_circle, r1_circle, da1_circle, da2_circle,rInter_circle]:
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
        cycloidal_gear2 = self.gear.ShapelyCycloidal2(ax=self.ax)
        #cycloidal_gear2 = self.gear.ShapelyCycloidal3(ax=self.ax)
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
        self.PC,=self.ax.plot(*p_poa,'yo')

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
        tangent=ax.arrow(x=P1_t[0],y=P1_t[1],dx=dt[0],dy=dt[1],**props,clip_on=True)
        normal=ax.arrow(x=P1_n[0],y=P1_n[1],dx=dn[0],dy=dn[1],**props,clip_on=True)
        text_normal=ax.text(*P1_n,"flank normal",va="bottom",ha="left",fontsize=8,clip_on=True)
        text_tangent=ax.text(*P1_t,"flank tangent",va="top",ha="left",fontsize=8,rotation=0,clip_on=True)
        dct=np.array((cos(0),sin(0)))*length
        P1_ct=np.array((0,gear.r2))-2*dct/3
        self.common_tangent=ax.arrow(x=P1_ct[0],y=P1_ct[1],dx=dct[0],dy=dct[1],**props,clip_on=True)
        self.text_common_tangent=ax.text(*P1_ct,"common tangent",va="bottom",ha="left",fontsize=8,clip_on=True)
        return tangent,normal,text_tangent,text_normal

    def initialize_vectors(self):
        zeta = self.gear.zetaA
        self.Fn_vec, self.Fr_vec = self.init_force_vector(zeta)
        self.vt1_arrow, self.vt2_arrow, self.vg_arrow = self.init_speed_vector(zeta)

    def init_angle(self,zeta:float):
        zeta0=zeta
        ax=self.ax
        xt,yt=(0,self.gear.r2)
        radius=self.angle_radius
        self.text_alpha=ax.text(xt+radius*0.9,yt-1,r"$\alpha$",va="top",ha="right",fontsize=12,color="black")
        alpha=90+self.gear.xi(zeta0)*180/pi
        print(f"alpha0={alpha:.1f}")
        #self.arc_alpha=Arc((xt,yt),width=width*2,height=width*2,theta1=-alpha,theta2=0,color="black",label=rf"$\alpha$={alpha:.1f}")
        xy=calculate_arcPolygon(center=(xt,yt),radius=radius,start_angle=0,end_angle=-alpha)
        #self.arc_alpha,=ax.plot(x,y,color="black")
        self.arc_alpha=Polygon(xy, closed=True,fill="grey",facecolor="grey",edgecolor="black",alpha=0.3)
        ax.add_patch(self.arc_alpha)

        #init zeta


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
        print("Fr not correct at init")
        Fr_arrow=ax.arrow(xc-dr[0],yc-dr[1],*dr,label="$F_r$",color='r',**default_arrow_props)
        dF2=-dn-dr
        self.arrow_F2=ax.arrow(xc-dF2[0],yc-dF2[1],*dF2,label="$F_{2}$",color='magenta',**default_arrow_props)
        return Fn_arrow,Fr_arrow

    def init_speed_vector(self,zeta,**kwargs):
        #allgemeines verzahnungsgesetz: Kraft muss immer durch den Wälzpunkt verlaufeb
        default_arrow_props = self.arrowprops
        gear=self.gear
        ax=self.ax
        w_length=kwargs.pop('w_length',1)
        vt1,vt2,vg=gear.v_gear(zeta,w_length)
        xc,yc=gear.p_pOA(zeta)
        # vt1_arrow=ax.arrow(xc-vt1[0],yc-vt1[1],*vt1,label="$v_{t1}$",color='gold',**default_arrow_props)
        vt1_arrow=ax.arrow(xc-vt1[0],yc-vt1[1],*vt1,color='gold',**default_arrow_props)
        vt2_arrow=ax.arrow(xc-vt2[0],yc-vt2[1],*vt2,color='salmon',**default_arrow_props)
        vt1_arrow.set_visible(False)
        vt2_arrow.set_visible(False)
        vg_arrow=ax.arrow(xc-vg[0],yc-vg[1],*vg,label="$v_g$",color='goldenrod',**default_arrow_props)
        return vt1_arrow,vt2_arrow,vg_arrow

    def add_tangent_and_normal_lines(self):
        zeta = self.zeta0
        self.tangent, self.normal, self.text_tangent, self.text_normal = self.add_lines(zeta)

    def setup_slider(self):
        zeta = self.zeta0*180/pi
        ax_slider = plt.axes([0.2, 0.05, 0.3, 0.02]) # [left, bottom, width, height]
        self.slider = Slider(ax_slider, 'Rotation', self.gear.zetaA*180/pi, self.gear.zetaE*180/pi, valinit=zeta,valstep=0.5)
        self.slider.on_changed(self.update_plot)



    #region update functions
    def update_plot(self, zeta:float):
        artists = []
        zeta=zeta/180*pi
        kappa = zeta * self.gear.i
        alpha_deg=90+self.gear.xi(zeta)*180/pi
        alpha_rad=alpha_deg*pi/180
        rotation = Affine2D().rotate_around(0, self.gear.a, kappa)
        self.polyPatch_Gear1.set_transform(rotation + self.ax.transData)
        artists.append(self.polyPatch_Gear1)
        rotation2 = Affine2D().rotate_around(0, 0, -zeta)
        self.polyPatch_Gear2.set_transform(rotation2 + self.ax.transData)
        artists.append(self.polyPatch_Gear2)

        # Update contact point
        p_pOA = self.gear.p_pOA(zeta)
        self.PC.set_data([p_pOA[0]], [p_pOA[1]])
        artists.append(self.PC)
        #self.PC.set_label(f"C = {zeta * 180 / np.pi:.1f}°")
        #self.angle_text.set_text(f"ζ={zeta * 180 / np.pi:.1f}°")
        artists.append(self.angle_text)
        # Update lines and vectors
        artists.extend(self.update_lines(zeta,alpha=alpha_rad,p_pOA=p_pOA))
        artists.extend(self.update_vectors(zeta,alpha=alpha_rad,p_pOA=p_pOA))
        artists.extend(self.update_pressureAnglePlot(zeta,alpha=alpha_deg))
        artists.extend(self.update_SlidingVelocityPlot(zeta))
        artist_alpha,=self.updateAlpha(alpha=float(alpha_deg))
        #self.updateExtras(zeta)
        artists.append(artist_alpha)
        #self.fig.savefig(f"export/frame{self.framenumber:06d}.png", format='png', dpi=120)
        self.framenumber+=1
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

    def update_vectors(self,zeta:float,alpha:float,p_pOA,**kwargs):
        gear=self.gear
        vt1_arrow,vt2_arrow,vg_arrow=self.vt1_arrow,self.vt2_arrow,self.vg_arrow
        if zeta==0:
            zeta=0.001
        w_length=kwargs.pop('w_length',1)
        vt1,vt2,vg=gear.v_gear(zeta,w_length)
        xc,yc=p_pOA

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

        if Fn_arrow:
            dn=np.array((cos(-alpha),sin(-alpha)))*F1_length
            Fn_arrow.set_data(x=xc-dn[0],y=yc-dn[1],dx=dn[0],dy=dn[1])
        if Fr_arrow:
            unit_vector = -vg / np.linalg.norm(vg)
            dr=unit_vector*abs(F1_length)*mu
            Fr_arrow.set_data(x=xc-dr[0],y=yc-dr[1],dx=dr[0],dy=dr[1])
        d2=-dn-dr
        self.arrow_F2.set_data(x=xc-d2[0],y=yc-d2[1],dx=d2[0],dy=d2[1])
        artists=[self.arrow_F2,self.vg_arrow,self.vt1_arrow,self.vt2_arrow,self.Fn_vec,self.Fr_vec]
        return artists

    def updateAlpha(self,alpha:float):
        angle_alpha=self.arc_alpha
        radius=self.angle_radius
        xy=calculate_arcPolygon(center=(0,self.gear.r2),radius=radius,start_angle=0,end_angle=-alpha)
        angle_alpha.set_xy(xy)
        return angle_alpha,

    def update_lines(self,zeta: float,alpha:float,p_pOA,**kwargs):
        length=kwargs.get("length",90)
        C=p_pOA
        s_alpha,c_alpha=sin(alpha),cos(alpha)
        dt=np.array((-s_alpha,-c_alpha))*length*-1
        dn=np.array((-c_alpha,s_alpha))*length
        P1_t=C-dt/2
        P1_n=C-dn/2
        
        self.tangent.set_data(x=P1_t[0],y=P1_t[1],dx=dt[0],dy=dt[1])
        self.normal.set_data(x=P1_n[0],y=P1_n[1],dx=dn[0],dy=dn[1])
        self.text_normal.set_position(P1_n)
        self.text_tangent.set_position(P1_t)
        artists=[self.tangent,self.normal,self.text_normal,self.text_tangent]
        return artists
    #endregion
   
    def key_press(self, event,**kwargs):
        if event==None:
            key=key=kwargs.get("key",'2')
        else: key=event.key
        
        if key == '1':
            vis_v = not self.ax_sliding_velocity.get_visible()
            self.showSpeeds(vis_v)
            if vis_v:
                self.title.set_text("sliding velocity $v_g$")   
        if key == '0':
            F_vec=self.Fn_vec,self.Fr_vec
            vis_v = not F_vec[0].get_visible()
            self.showForces(vis_v)
            if vis_v:
                self.title.set_text("Forces")            
        if key == '2':
            lines=self.tangent,self.normal,self.text_tangent,self.text_normal,self.common_tangent,self.text_common_tangent,self.text_alpha,self.arc_alpha
            vis_v = not lines[0].get_visible()
            self.showLines(vis_v)
        if key=="4":
            self.showAll(False)
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
    def setup_pressure_angle_plot(self,ax,gear,zeta:float=0):
        # gear.A()
        # gear.E()
        zetaA=gear.zetaA
        zetaE=gear.zetaE
        angles=np.linspace(zetaA,zetaE,num=400)
        
        #xis=np.array([(gear.xi(a))*180/pi for a in angles])
        alphas=np.array([90+(gear.xi(a))*180/pi for a in angles])
        
        # x_coords,y_coords=gear.calculate_path_of_contact(start=zetaA,end=zetaE,num_points=len(angles))
        # lengths = np.zeros(len(angles))
        # for i in range(1, len(angles)):
        #     lengths[i] = lengths[i-1] + np.sqrt((x_coords[i] - x_coords[i-1])**2 + (y_coords[i] - y_coords[i-1])**2)
        ax.plot(angles*180/pi,alphas,label="α_t",color="black", animated=False)
        ax.set_xlabel("Rotation ζ in °")
        ax.set_ylabel(r"Pressure angle $\alpha$ in °")
        self.text_pressureAngle=ax.text(0.5,0.9,r"pressure angle $\alpha$",
            horizontalalignment='center', 
            verticalalignment='top', 
            transform=ax.transAxes)
        ax.grid(axis='y')
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        #alpha_zeta=gear.alpha_t(zeta)*180/pi

        
        self.marker_pressureAngle,=ax.plot(angles[0]*180/pi,alphas[0],marker='o',color="black",linestyle="",linewidth=4)

        # y_max=max(alphas)+5
        # xmin,xmax=min(angles*180/pi),max(angles*180/pi)+3
        # ax.set_ylim(0,y_max)
        # ax.set_xlim(xmin,xmax)
        # arrowprops=dict(head_length=1.5, head_width=1.2, linewidth=1,length_includes_head=True, clip_on=False,color="black")
        # ax.arrow(x=xmin,y=0,dx=xmax-xmin,dy=0,**arrowprops)
        # ax.arrow(x=xmin,y=0,dx=0,dy=y_max,**arrowprops)
        
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
        group_forces=[self.Fn_vec,self.Fr_vec,self.arrow_F2]
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
        group_forces=[self.Fn_vec,self.Fr_vec,self.arrow_F2]
        group_forces.append(self.ax_sliding_velocity)
        for i,item in enumerate(group_forces):
            item.set_visible(show)
    
    def showSlider(self,show:bool=False):
        self.vis_slider=show
        if slider:=self.slider:
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

    def createAnimation(self,frames_per_anim:int=100,fps:int=10):
        #self.key_press(None,key="1")
        #self.slider.set_visible(False)
        writervideo = mpl.animation.FFMpegWriter(fps=fps) 
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
    gear_pair = EccentricCycloidGear(
        z1=3,
        z2=6,
        a=100,
        rA_star=1.0,
        st_star=1.0,
        phi_j1=0.1*pi/180,
        phi_As=49*pi/180,
        phi_Ae=170*pi/180,
        lambda_=0.97,
        c_star=0.125
    )
    # gear_pair = EccentricCycloidGear(
    #     z1=3,
    #     z2=6,
    #     a=100,
    #     rA_star=1.0,
    #     st_star=1.0,
    #     phi_j1=1*pi/180,
    #     phi_As=60*pi/180,
    #     phi_Ae=170*pi/180,
    #     lambda_=0.999
    # )
    gear_meshing_plot = GearMeshingPlot(gear_pair)
    #gear_meshing_plot.gear.rF1=16.1
    gear_meshing_plot.setup_plot()
