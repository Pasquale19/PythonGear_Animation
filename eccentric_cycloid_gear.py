from dataclasses import dataclass, field
from math import pi, cos, sin, acos, atan2, sqrt,atan,asin
from typing import Tuple, Dict
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from shapely.geometry import Point, Polygon
from shapely import affinity as affinity
from intersection import find_intersection
from Utilities.Converter import convert_to_rack
#from deprecated import deprecated

@dataclass
class EccentricCycloidGear:
    """
    Represents an eccentric cycloidal gear pair.

    This class encapsulates the geometry and calculations for an eccentric cycloidal gear pair,
    including methods to generate profiles and calculate the path of contact.

    Attributes:
        z1 (int): Number of teeth of arc gear
        z2 (int): Number of teeth of cycloid gear
        a (float): Center distance
        kwargs (Dict): Additional parameters for gear configuration
    """

    z1: int
    z2: int
    a: float
    kwargs: Dict = field(default_factory=dict)

    def __init__(self, z1: int, z2: int, a: float, **kwargs):
        """
        Initializes the gear pair with given parameters and calculates derived values.

        Args:
            z1 (int): Number of teeth of arc gear
            z2 (int): Number of teeth of cycloid gear
            a (float): Center distance
            **kwargs: Additional parameters for gear configuration
        """
        self.z1 = z1
        self.z2 = z2
        self.a = a

        self.m = self.a / (self.z1 / 2 + self.z2 / 2)
        self.d1 = self.m * self.z1
        self.d2 = self.m * self.z2
        self.r1 = self.d1 / 2
        self.r2 = self.d2 / 2
        
        self.i = self.z2 / self.z1
        self.rw1 = self.r1
        self.rw2 = self.r2

        self.rA_star = kwargs.get('rA_star', 1.0)
        self.c_star = kwargs.get('c_star', 0.1)
        self.st_star = kwargs.get('st_star', 1.0)
        self.phi_j1 = kwargs.get('phi_j1', 0.0)
        self.phi_As = kwargs.get('phi_As', 60*pi/180)
        self.phi_Ae = kwargs.get('phi_Ae', 170*pi/180)
        self.lambda_ = kwargs.get('lambda_', 0.97)
        self.delta = kwargs.get('delta', 0.1) #equidistant distance
        self.e = self.r1*self.lambda_
        self.calculate_arc_gear_geometry()
        self.calculate_cycloid_gear_geometry()
        self.calculate_A()
        self.calculate_E()
        phi_As=self.phi_As
        A_start=np.array((sin(0), -cos(0))) * self.e + np.array((-sin(phi_As), cos(phi_As))) * self.rA
        self.rInter=np.hypot(A_start[0],A_start[1])
        
    
    def calculate_A(self,**kwargs):
        """
        calculate start point of mesh and set self.A and self.zetaA
        """
        steps=kwargs.get('steps',250)
        def Kreis_ra2(t):
            return -sin(t)*self.ra2,cos(t)*self.ra2
        t_range_circle=np.array((pi/(self.z2*4),3/4*pi/self.z2))
        t_range_contact=np.array((-0,-2*pi/self.z2))
        
        A=find_intersection(Kreis_ra2,self.p_pOA,t_range=t_range_circle,t_range2=t_range_contact,precision=1,steps=steps,return_t=True)
        if A:
            self.A=A[0]
            self.zetaA=A[2]
            print(f"zeta A={self.zetaA*180/pi:.1f}")
            print(f"Intersection A found at approximately: {A}")
        else:
            print(f"No intersection for A found in the range({t_range_contact*180/pi})")

    def calculate_E(self,**kwargs):
        """
        end point of mesh
        """
        steps=kwargs.get('steps',250)
        def Kreis_ra1(t):
            return sin(t)*self.ra1,-cos(t)*self.ra1+self.a
        
        t_range_circle=np.array((0,1*pi/self.z2))
        t_range_contact=np.array((-0,-pi/(self.z2*2)))
        
        E=find_intersection(Kreis_ra1,self.p_pOA,t_range=t_range_circle,t_range2=t_range_contact,precision=0.5,steps=steps,return_t=True)

        if E:
            self.E=E[0]
            self.zetaE=E[2]
            print(f"Intersection E found at approximately: {E}")
        else:
            print(f"No intersection for E found in the range({t_range_contact*180/pi})")
        return E[0]
    

    
    def calculate_arc_gear_geometry(self):
        """
        Calculates the geometry parameters for the arc gear.
        """
        e=self.e
        self.c=self.c_star*self.m
        rA=self.rA_star * e*sqrt(2-2*cos(pi/(2*self.z1)))
        self.rA = rA
        self.phi_rA = 2 * acos((2 * self.e**2 - self.rA**2) / (2 * self.e**2))-self.st_star*pi/(self.z1*1) #richtig
        self.phi_s1 = self.phi_rA  - self.phi_j1
        self.phi_s1 = self.phi_rA  + self.phi_j1
        self.phi_tooth=pi/(self.z1)*self.st_star #angle of one tooth
        
        #phne Backlash phi_j1 funktioniert nicht
        q_F1 = e * sin(self.phi_As) / sin(pi -(pi/self.z1) -self.phi_As-self.phi_rA/2)
        self.rF1 = sqrt(q_F1**2 + self.e**2 - 2 * q_F1 * self.e * cos(pi/self.z1+self.phi_rA/2)) - self.rA
        #mit Backlash phi_j1
        q_F1 = e * sin(self.phi_As) / sin(pi -(pi/self.z1) -self.phi_As-self.phi_s1/2)
        self.rF1 = sqrt(q_F1**2 + self.e**2 - 2 * q_F1 * self.e * cos(pi/self.z1+self.phi_s1/2)) - self.rA


        self.qF1=q_F1
        self.df1 = 2 * (q_F1 - self.rF1)
        self.rf1 = self.df1 / 2
        self.da1 = 2 * (self.e - self.rA * cos(self.phi_Ae))
        self.ra1 = self.da1 / 2
        
        self.r_start=sqrt(e**2+rA**2-2*e*rA*cos(self.phi_As))

    def calculate_cycloid_gear_geometry(self):
        """
        Calculates the geometry parameters for the cycloid gear.
        """
        #self.df2 = 2 * (self.a - (self.rA+self.e)-self.c)
        self.df2 = 2 * (self.a - self.ra1-self.c)
        self.rf2 = self.df2 / 2
        self.da2 = 2 * (self.a - self.rf1-self.c)
        self.ra2 = self.da2 / 2

    def calculateRequired_rF1(self,**kwargs):
        num=kwargs.get("num",1000)
        angle=pi/self.z1+self.zetaA*self.i
        F10=np.array((0,self.a))+np.array((sin(angle),-cos(angle)))*self.qF1
        K0=self.p_pOA(self.zetaA)
        angles=np.linspace(self.zetaA,self.zetaE,num=num,endpoint=False)
        angles2=-self.zetaA+angles
        angles1=angles2*self.i
        from Utilities.rotatePoint import rotate_point
        F1s=rotate_point(F10,angles1,center=(0,self.a),use_radians=True)
        Ks=rotate_point(K0,-angles2,center=(0,0),use_radians=True)
        ax=kwargs.get("ax",False)
        distance=-F1s+Ks
        d=np.hypot(distance[:,0],distance[:,1])
        if ax:
            ax.plot(F1s[:,0],F1s[:,1],label="F1s",marker='o')
            for i, p in enumerate(F1s):
                c_F1=patches.Circle(p,radius=self.rF1,fill=False,color=(1/num*i,1/num*i,1/num*i),linestyle=":",alpha=1)
                ax.add_patch(c_F1)
            for i,di in enumerate(d):
                print(f"distance is {di:.2f} vs rF1 {self.rF1:.2f}")
            ax.plot(Ks[:,0],Ks[:,1],label="Ks",marker='o')
            #ax.plot(d[:,0],d[:,1],label="d")
        rF1_required=max(d)
        print(f"current rF1={self.rF1:.2f} vs needed {rF1_required:.2f} ")


    def get_arc_profile(self,rotation_angle=0, num_points: int = 100,full_gear:bool=True,**kwargs) -> Tuple[Tuple[float, ...], Tuple[float, ...]]:
        """
        Generates the profile of the arc gear with mathematical formulas

        Args:
            num_points (int): Number of points to generate for the profile.

        Returns:
            Tuple[Tuple[float, ...], Tuple[float, ...]]: x and y coordinates of the arc profile.
        """
        x_coords=[]
        y_coords=[]
        
        
        phi_As=self.phi_As
        phi_Ae=self.phi_Ae
        angles_rightFlank=np.linspace(phi_As,phi_Ae,num=int(num_points/4))
        phi_s1=kwargs.get("phi_s1",self.phi_s1)
        center_x=kwargs.get("center_x",0)
        center_y=kwargs.get("center_y",self.a)
        #right_flank
        for a in angles_rightFlank:
            theta=rotation_angle-phi_s1/2 #winkel des rihtungsvektors
            a2=-a+theta
            x,y=np.array((sin(theta),-cos(theta)))*self.e+np.array((0,self.a))+np.array((-sin(a2),cos(a2)))*self.rA
            x_coords.append(x)
            y_coords.append(y)
                
        x_top,y_top=np.array((sin(rotation_angle),-cos(rotation_angle)))*(self.ra1)+np.array((0,self.a))
        x_coords.append(x_top)
        y_coords.append(y_top)
        #left_flank
        for a in angles_rightFlank[::-1]: #reverse order
            theta=rotation_angle+phi_s1/2 #winkel des rihtungsvektors
            a2=a+theta
            x,y=np.array((sin(theta),-cos(theta)))*self.e+np.array((0,self.a))+np.array((-sin(a2),cos(a2)))*self.rA
            x_coords.append(x)
            y_coords.append(y)
        
        #left foot
        phi_Fs=pi-pi/self.z1-phi_As
        angles_foot=np.linspace(phi_Fs,-phi_Fs,num=int(num_points/2),endpoint=False)
        for a in angles_foot:
            theta=rotation_angle-pi/self.z1 #winkel des rihtungsvektors
            a2=-a+theta
            x,y=np.array((sin(theta),-cos(theta)))*self.qF1+np.array((0,self.a))+np.array((-sin(a2),cos(a2)))*self.rF1
            x_coords.append(x)
            y_coords.append(y)
        #return x_coords, y_coords
        xy=np.array([np.array(x_coords),np.array(y_coords)])
        translated_points=-np.array([[0],[self.a]])+xy
        xy_gear=np.empty((0, 2))
        empty_list=[]
        for i in range(0,self.z1,1):
            angle_rad=2*pi/self.z1*i
            rotation_matrix = np.array([[np.cos(angle_rad), -np.sin(angle_rad)],
                                [np.sin(angle_rad),  np.cos(angle_rad)]])
            rotated_points = translated_points.T @ rotation_matrix
            empty_list.append(rotated_points)
            #xy_gear=np.append(xy_gear,rotated_points,axis=None)
        xy_gear=np.vstack(empty_list)
        if full_gear:
            return xy_gear[:,0]+center_x,xy_gear[:,1]+center_y
        else:
            return x_coords+center_x, y_coords+center_y
        
    
        #left_foot



    def get_arc_profile2(self, rotation_angle: float = 0, num_points: int = 100, full_gear: bool = True,**kwargs) -> Tuple[np.ndarray, np.ndarray]:
        '''faster method to calculate the profile of the first gear'''
        phi_As, phi_Ae = self.phi_As, self.phi_Ae
        angles_right_flank = np.linspace(phi_As, phi_Ae, num=max(1, num_points // 4))
        phi_s1=kwargs.get("phi_s1",self.phi_s1)
        center_x=kwargs.get("center_x",0)
        center_y=kwargs.get("center_y",self.a)
        #print(f"{phi_s1=}")
        def calculate_coords(angle: float, theta: float) -> Tuple[float, float]:
            a2 = -angle + theta
            return np.array((sin(theta), -cos(theta))) * self.e + np.array((0, self.a)) + np.array((-sin(a2), cos(a2))) * self.rA

        
        # Right flank
        coords = [calculate_coords(a, rotation_angle - phi_s1 / 2) for a in angles_right_flank]

        # Top point
        coords.append(np.array((sin(rotation_angle), -cos(rotation_angle))) * self.ra1 + np.array((0, self.a)))

        # Left flank
        coords.extend(calculate_coords(-a, rotation_angle + phi_s1 / 2) for a in reversed(angles_right_flank))

        # Left foot
        phi_Fs = pi - pi / self.z1 - phi_As
        num=max(1, num_points // 2)
        step=(2*phi_Fs)/(num)
        angles_foot = np.linspace(phi_Fs-step, -phi_Fs, num=num, endpoint=False)
        theta = rotation_angle - pi / self.z1
        coords.extend(np.array((np.sin(theta), -np.cos(theta))) * self.qF1 + np.array((0, self.a)) + 
                    np.array((-np.sin(-a + theta), np.cos(-a + theta))) * self.rF1 for a in angles_foot)

        xy = np.array(coords).T - np.array([[0], [self.a]])
       
        if not full_gear:
            return xy[0]+center_x , xy[1]+ center_y

        rotation_matrices = [
            np.array([[np.cos(angle), -np.sin(angle)],
                    [np.sin(angle), np.cos(angle)]])
            for angle in np.linspace(0, 2 * np.pi, self.z1, endpoint=False)
        ]

        xy_gear = np.vstack([xy.T @ matrix for matrix in rotation_matrices])
        return xy_gear[:, 0]+center_x, xy_gear[:, 1] + center_y


    def get_arc_profile3(self, rotation_angle: float = 0, num_points: int = 100, full_gear: bool = True,**kwargs) -> Tuple[np.ndarray, np.ndarray]:
        '''faster method to calculate the profile of the first gear'''
        phi_As, phi_Ae = self.phi_As, self.phi_Ae
        angles_right_flank = np.linspace(phi_As, phi_Ae, num=max(1, num_points // 4))
        phi_s1=kwargs.get("phi_s1",self.phi_s1)
        center_x=kwargs.get("center_x",0)
        center_y=kwargs.get("center_y",self.a)
        #print(f"{phi_s1=}")
        def calculate_coords(angle: float, theta: float) -> Tuple[float, float]:
            a2 = -angle + theta
            return np.array((sin(theta), -cos(theta))) * self.e + np.array((0, self.a)) + np.array((-sin(a2), cos(a2))) * self.rA
        phi_Fs = pi - pi / self.z1 - phi_As
        num=max(1, num_points // 2)
        angles_footl = np.linspace(0, -phi_Fs, num=num, endpoint=False)
        theta2 = rotation_angle + pi / self.z1
        coords2=[np.array((np.sin(theta2), -np.cos(theta2))) * self.qF1 + np.array((0, self.a)) + 
                    np.array((-np.sin(-a + theta2), np.cos(-a + theta2))) * self.rF1 for a in angles_footl]

        coords=coords2
        
        # Right flank
        coords.extend(calculate_coords(a, rotation_angle - phi_s1 / 2) for a in angles_right_flank)

        # Top point
        coords.append(np.array((sin(rotation_angle), -cos(rotation_angle))) * self.ra1 + np.array((0, self.a)))

        # Left flank
        coords.extend(calculate_coords(-a, rotation_angle + phi_s1 / 2) for a in reversed(angles_right_flank))

        # Left foot
        
        
        step=(2*phi_Fs)/(num)
        angles_foot = np.linspace(phi_Fs-step, -0, num=num, endpoint=False)
        theta = rotation_angle - pi / self.z1
        coords.extend(np.array((np.sin(theta), -np.cos(theta))) * self.qF1 + np.array((0, self.a)) + 
                    np.array((-np.sin(-a + theta), np.cos(-a + theta))) * self.rF1 for a in angles_foot)
        
        

        xy = np.array(coords).T - np.array([[0], [self.a]])
       
        if not full_gear:
            return xy[0]+center_x , xy[1]+ center_y

        rotation_matrices = [
            np.array([[np.cos(angle), -np.sin(angle)],
                    [np.sin(angle), np.cos(angle)]])
            for angle in np.linspace(0, 2 * np.pi, self.z1, endpoint=False)
        ]

        xy_gear = np.vstack([xy.T @ matrix for matrix in rotation_matrices])
        return xy_gear[:, 0]+center_x, xy_gear[:, 1] + center_y

    def get_rackProfile(self,num=1000):
        
        x,y=self.get_arc_profile2(rotation_angle=pi,num_points=num,full_gear=True)
        xn,yn=convert_to_rack(x,y,self.r1,center_x=0,center_y=self.a)
        return xn,yn
    
    def calculateFlankParameter(self,num_points:int=1000,epsilon=0.0001):
        x,y=self.get_rackProfile(num=num_points)
        arr=np.zeros(shape=(3,num_points-1))
        p0=np.array((x[0],y[0]))
        for i in range(1,num_points,1):
            p1=np.array((x[i],y[i]))
            d=p0-p1
            d=d/np.linalg.norm(d)
            alpha=np.arctan2(d[1],d[0])
            arr[1,i-1]=alpha
            n=(0-p1[1])/cos(alpha)
            if round(abs(alpha)-pi/2,4)==0:
                print(f"{p1=}")
            if p1[1]**2<epsilon:
                xr=p1[0]
                print(f"smaller epsilon {p1[1]=}")
            else:
                xr=p1[0]+n*-sin(alpha)
            phi=xr/self.r1
            arr[0,i-1]=phi
            #arr[0,i-1]=xr
            arr[2,i-1]=n
            print(f"phi={phi*180/pi:.1f}°\t xr={xr:.1f}\t alpha={alpha*180/pi:.5f} \t n={n:.1f}")
            p0=p1
        return arr
    
    def plotFlankParameter(self):
        from Utilities.Converter import calculateFlankParameter
        xr,yr=self.get_rackProfile(num=num_points)
        plt.rcParams['lines.markersize'] = 1
        num_points=100
        arr=self.calculateFlankParameter(num_points=num_points)
        
        
        xg,yg=self.get_arc_profile2(rotation_angle=pi,num_points=num_points,full_gear=False)
        xr,yr=convert_to_rack(xg,yg,self.r1,center_x=0,center_y=self.a)
        alpha=arr[1,:]
        n=arr[2,:]
        x=arr[0,:]
        fig,(ax1,ax2)=plt.subplots(2)
        
        ax1.set_aspect("equal")
        c=patches.Circle(xy=(0,0),radius=self.r1,label="r1",fill=False,linestyle=":")
        ax1.add_patch(c)
        ax1.hlines(0,xmin=-30,xmax=30)
        ax1.grid()
        scatter = ax1.scatter(xr, yr, c=range(len(xr)), cmap='viridis',label="rack")
        plt.colorbar(scatter,label='index')
        #ax1.plot(xr,yr,label=r"$rack$")
        ax1.scatter(xg, yg-self.a, c=range(len(xg)), cmap='viridis',label="gear")
        #ax1.plot(xg,yg-self.a,label=r"$gear$")
        ax1.legend()
        ax2.scatter(x, n, c=range(len(n)), cmap='viridis',label="n")
        #ax2.plot(x,n,label="n")
        ax2.scatter(x, alpha*180/pi, c=range(len(x)), cmap='viridis',label=r"$\alpha$")
        #ax2.plot(x,alpha*180/pi,label=)
        ax2.legend()
        
        plt.show()

            

    def get_cycloid_profile(self, num_points: int = 100) -> Tuple[Tuple[float, ...], Tuple[float, ...]]:
        """
        Generates the profile of the cycloid gear stars with -1.5*pi/z2

        Args:
            num_points (int): Number of points to generate for the profile.

        Returns:
            Tuple[Tuple[float, ...], Tuple[float, ...]]: x and y coordinates of the cycloid profile.
        """
        theta_range = np.linspace(0, 2*pi, num_points)
        x_coords, y_coords = self.get_cycloid_profile_points(theta_range)
    
        return x_coords, y_coords
        x_coords = []
        y_coords = []
        for theta in theta_range:
            x, y = self.get_cycloid_profile_point(theta)
            x_coords.append(x)
            y_coords.append(y)
        return tuple(x_coords), tuple(y_coords)
    
    def get_cycloid_profile_point(self, zeta: float) -> Tuple[float, float]:
        """
        p_pC Calculates a single point on the cycloid profile p_pC

        Args:
            theta (float): Angle parameter for the cycloid curve.

        Returns:
            Tuple[float, float]: x and y coordinates of the point on the cycloid profile.
        """
      
        a=self.a
        e=self.e
        i=self.i
        kappa=i*zeta
        xi=self.xi(zeta)
        xy=a*np.array((-sin(zeta),cos(zeta)))-e*np.array((-sin(zeta+kappa),cos(zeta+kappa)))-self.rA*np.array((-sin(-xi+zeta),cos(-xi+zeta)))
        return xy

    def xi(self,zeta:float=0):
        """
        contact angle ξ
        """
        lam=self.lambda_
        i=self.i
        return np.arctan(lam*np.sin(i*zeta)/(1-lam*np.cos(i*zeta)))
    
    def alpha_t(self,zeta:float=0):
        """
        pressure angle alpha
        """
        return pi/2-self.xi(zeta)

    def list_parameters(self) -> str:
        """
        Returns a string containing all parameters of the gear pair, each on a new line.

        Returns:
            str: A string with all parameters and their values.
        """
        params = [
            f"d1 = {self.d1:.1f}",
            f"r1 = {self.r1:.1f}",
            f"e = {self.e:.1f}",
            f"qF1 = {self.qF1:.1f}",
            f"rw1 = {self.rw1:.1f}",
            f"rw2 = {self.rw2:.1f}",
            f"rA_star = {self.rA_star:.1f}",
            f"st_star = {self.st_star:.1f}",
            f"phi_j1 = {self.phi_j1*180/pi:.1f}",
            f"phi_As = {self.phi_As*180/pi:.1f}",
            f"phi_Ae = {self.phi_Ae*180/pi:.1f}",
            f"lambda_ = {self.lambda_:.2f}",
            f"phi_tooth = {self.phi_tooth*180/pi:.1f}",
            f"rA = {self.rA:.1f}",
            f"phi_rA = {self.phi_rA*180/pi:.1f}",
            f"phi_s1 = {self.phi_s1*180/pi:.1f}",
            f"rF1 = {self.rF1:.1f}",
            f"df1 = {self.df1:.1f}",
            f"rf1 = {self.rf1:.1f}",
            f"da1 = {self.da1:.1f}",
            f"ra1 = {self.ra1:.1f}",
            f"df2 = {self.df2:.1f}"
        ]
        return "\n".join(params)    
    
    def p_pOA(self,zeta:float):
        '''contact point'''
        a=self.a
        xi=self.xi(zeta)
        e=self.e
        x=a*0-e*-sin(self.i*zeta)-self.rA*-sin(-xi)
        y=a*1-e*cos(self.i*zeta)-self.rA*cos(-xi)
        return np.array((x,y))
    
    def p_OA(self,zeta:float):
        """Punktvektor A"""
        p_OA=self.a*np.array((-sin(zeta),cos(zeta)))-self.e*np.array(((-sin(zeta+self.i*zeta),cos(zeta+self.i*zeta))))
        return p_OA
        
    def calculate_path_of_contact(self, num_points: int = 1000,**kwargs) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculates the path of contact between the arc and cycloid gears using Shapely.

        Args:
            num_points (int): Number of points to check for contact.

        Returns:
            Tuple[np.ndarray, np.ndarray]: x and y coordinates of the path of contact.
        """
        default_start=self.zetaA if self.zetaA else -pi/(self.z1*2)
        default_end=self.zetaE if self.zetaE else pi/(self.z1*2)
        start=kwargs.get('start',default_start)
        end=kwargs.get('end',default_end)
        angles=np.linspace(start,end,num=num_points,endpoint=True)
        e=self.e
        a=self.a
        xp=[]
        yp=[]
        
        for zeta in angles:
            xi=self.xi(zeta)
            x=a*0-e*-sin(self.i*zeta)-self.rA*-sin(-xi)
            y=a*1-e*cos(self.i*zeta)-self.rA*cos(-xi)
            xp.append(x)
            yp.append(y)
        return xp,yp
    
    def calculate_path_of_contact2(self, num_points: int = 1000, **kwargs) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculates the path of contact between the arc and cycloid gears using Shapely.

        Args:
            num_points (int): Number of points to check for contact.

        Returns:
            Tuple[np.ndarray, np.ndarray]: x and y coordinates of the path of contact.
        """
        default_start = self.zetaA if self.zetaA is not None else -np.pi / (self.z1 * 2)
        default_end = self.zetaE if self.zetaE is not None else np.pi / (self.z1 * 2)
        start = kwargs.get('start', default_start)
        end = kwargs.get('end', default_end)
        
        angles = np.linspace(start, end, num=num_points, endpoint=True)
        
        # Vectorize calculations
        xi_values = self.xi(angles)
        sin_i_zeta = np.sin(self.i * angles)
        cos_i_zeta = np.cos(self.i * angles)
        sin_xi = np.sin(-xi_values)
        cos_xi = np.cos(-xi_values)
        
        # Calculate x and y coordinates in a vectorized manner
        x = -self.e * sin_i_zeta + self.rA * sin_xi
        y = self.a - self.e * cos_i_zeta - self.rA * cos_xi
        
        return x, y

    def v_gear(self,zeta:float,w1:float=1):
        p_poa=self.p_pOA(zeta)
        a=self.a
        v1=w1*np.array((-p_poa[1]+a,p_poa[0]))
        v2=w1/self.i*np.array((p_poa[1],-p_poa[0]))
        alpha_t=self.alpha_t(zeta)
        vg=np.dot(v1-v2,np.array((sin(alpha_t),cos(alpha_t))))
        vt1=np.dot(v1,np.array((-sin(alpha_t),cos(alpha_t))))
        vt2=np.dot(v2,np.array((-sin(alpha_t),cos(alpha_t))))
        vt1_vec=vt1*np.array((-sin(alpha_t),cos(alpha_t)))
        vt2_vec=vt2*np.array((-sin(alpha_t),cos(alpha_t)))
        vg=(vt1-vt2)*np.array((-sin(alpha_t),cos(alpha_t)))
       
        return vt1_vec,vt2_vec,-vg

    def Kg(self,zeta:float):
        p_poa=self.p_pOA(zeta)
        a=self.a
        w1=1
        v1=w1*np.array((-p_poa[1]+a,p_poa[0]))
        v2=w1/self.i*np.array((p_poa[1],-p_poa[0]))
        alpha_t=self.alpha_t(zeta)
        vg=np.dot(v1-v2,np.array((sin(alpha_t),cos(alpha_t))))
        vt1=np.dot(v1,np.array((-sin(alpha_t),cos(alpha_t))))
        vt2=np.dot(v2,np.array((-sin(alpha_t),cos(alpha_t))))
        vt1_vec=vt1*np.array((-sin(alpha_t),cos(alpha_t)))
        vt2_vec=vt2*np.array((-sin(alpha_t),cos(alpha_t)))
        vg=(vt1-vt2)*np.array((-sin(alpha_t),cos(alpha_t)))
        Kg=-(vt1-vt2)/(self.r1*w1)
        return Kg


    def ShapelyGear(self,**kwargs):
        """
        Calculates the polygon of a gear with arc-shaped teeth.

        Parameters:
        r (float): Pitch radius of the gear
        z (int): Number of teethqF1
        center_x (float): X-coordinate of the gear center (default is 0)
        center_y (float): Y-coordinate of the gear center (default is 0)
        sector (tuple): Sector of the gear to be considered (default is (pi/2, 2*pi+pi/2))
        Rr (float): Radius of the tooth arc (default is r*sin(pi/(z*2)))

        Returns:
        shapely.geometry.Polygon: Polygon representing the gear
        """       
        r=z=kwargs.get('r',self.r1)
        z=kwargs.get('z',self.z1)
        #m=2*r/z
        da=self.da1
        center_x=kwargs.get('center_x',0)
        center_y=kwargs.get('center_y',0)
        initial_angle=kwargs.get('initial_angle',pi/2)
        rA= kwargs.get('rA',0)
        if rA <=0:
            rA=r*sin(pi/(z*2))
        rA=self.rA
        angles = np.linspace(initial_angle%(2*pi), 2*np.pi+initial_angle%(2*pi), z, endpoint=False)
        
        e=self.e
        xp = [e * np.cos(angle) + center_x for angle in angles]
        yp = [e * np.sin(angle) + center_y for angle in angles]

        #gear=Point(0,0).buffer(r_circle)
        gear=Point(center_x,center_y).buffer(e*0.9)
        print("Startkreis zu groß")
        angles_foot=[a+pi/z for a in angles]
        rF=self.rF1
        qF1=self.qF1
        
        
        #erstellen der Zahnköpfe
        for a in angles:
            a1=a-self.phi_s1
            circle_a=Point(cos(a1)*e+center_x,sin(a1)*e+center_y).buffer(rA)
            gear=gear.union(circle_a)
            a2=a+self.phi_s1
            circle_a=Point(cos(a2)*e+center_x,sin(a2)*e+center_y).buffer(rA)
            gear=gear.union(circle_a)
        
        #erstellen der Zahnfüße
        for a in angles_foot:
            circle_f=Point(cos(a)*qF1+center_x,sin(a)*qF1+center_y).buffer(rF)
            gear=gear.difference(circle_f)

        
        
        # for i, (x, y) in enumerate(list(zip(xp,yp))):
        #     circle_a = Point(x, y).buffer(rA)#Kopf
        #     gear=gear.union(circle_a)
        r_f=self.rf1
        #gear=gear.union(Point(center_x,center_y).buffer(r_f))
        
            
        # Cut out tips
        circle=Point(center_x,center_y).buffer(da/2)
        gear=gear.intersection(circle)
        return gear
    
    def ShapelyCycloidal(self,**kwargs):
        """
        Calculates the polygon of a gear with arc-shaped teeth.

        Parameters:
        r (float): Pitch radius of the gear
        z (int): Number of teethqF1
        center_x (float): X-coordinate of the gear center (default is 0)
        center_y (float): Y-coordinate of the gear center (default is 0)
        sector (tuple): Sector of the gear to be considered (default is (pi/2, 2*pi+pi/2))
        Rr (float): Radius of the tooth arc (default is r*sin(pi/(z*2)))

        Returns:
        shapely.geometry.Polygon: Polygon representing the gear
        """
        ra2=self.ra2
        a=self.a
        base=Point(0,0).buffer(ra2)
        
        x1,y1=self.get_arc_profile2(0)
        cutter=Polygon(list(zip(x1,y1)))
        cutter=affinity.translate(cutter,xoff=0,yoff=-a)
        num_points=z=kwargs.get('num_points',360)
        angles=np.linspace(0,2*pi,num=num_points,endpoint=False)        
        for theta in angles:
            rotation_angle=theta+theta*self.i
            cutter2=affinity.rotate(cutter,angle=rotation_angle,origin=(0,0),use_radians=True)
            cutter2=affinity.translate(cutter2,xoff=-sin(theta)*a,yoff=cos(theta)*a)
            base=base.difference(cutter2)
        return base
    
    def ShapelyCycloidal2(self,**kwargs):
        """
        Calculates the polygon of a gear with arc-shaped teeth.

        Parameters:
        r (float): Pitch radius of the gear
        z (int): Number of teethqF1
        center_x (float): X-coordinate of the gear center (default is 0)
        center_y (float): Y-coordinate of the gear center (default is 0)
        sector (tuple): Sector of the gear to be considered (default is (pi/2, 2*pi+pi/2))
        Rr (float): Radius of the tooth arc (default is r*sin(pi/(z*2)))
        num (int): number of points to calculate default=6000
        Returns:
        shapely.geometry.Polygon: Polygon representing the gear
        """
        num=kwargs.get("num",6000)
        theta_range = np.linspace(-1.3*pi/self.z2, 1.3*pi/self.z2, num=num//self.z2,endpoint=False)
        x_coords = [0]
        y_coords = [self.a]
        for theta in theta_range:
            x, y = self.get_cycloid_profile_point(theta)
            x_coords.append(x)
            y_coords.append(y)

        base=Point(0,0).buffer(self.ra2)
        cutter=Polygon(list(zip(x_coords,y_coords)))
        
        #return cutter
        angles=np.linspace(0,2*pi,num=int(self.z2),endpoint=False)        
        for theta in angles:
            cutter2=affinity.rotate(cutter,angle=theta,origin=(0,0),use_radians=True)
            base=base.difference(cutter2)
        return base

    def ShapelyCycloidal3(self,**kwargs):
        """
        Calculates the polygon of a gear with arc-shaped teeth.

        Parameters:
        r (float): Pitch radius of the gear
        z (int): Number of teethqF1
        center_x (float): X-coordinate of the gear center (default is 0)
        center_y (float): Y-coordinate of the gear center (default is 0)
        sector (tuple): Sector of the gear to be considered (default is (pi/2, 2*pi+pi/2))
        Rr (float): Radius of the tooth arc (default is r*sin(pi/(z*2)))
        num (int): number of points to calculate default=6000
        Returns:
        shapely.geometry.Polygon: Polygon representing the gear
        """
        num=kwargs.get("num",6000)
        ax=kwargs.get("ax",False)
        theta_range = np.linspace(-1.3*pi/self.z2, 1.3*pi/self.z2, num=num//self.z2,endpoint=False)
        x_coords = [0]
        y_coords = [self.a]
        for theta in theta_range:
            x, y = self.get_cycloid_profile_point(theta)
            x_coords.append(x)
            y_coords.append(y)

        base=Point(0,0).buffer(self.ra2)
        cutter=Polygon(list(zip(x_coords,y_coords)))
        
        #return cutter
        angles=np.linspace(0,2*pi,num=int(self.z2),endpoint=False)        
        for theta in angles:
            cutter2=affinity.rotate(cutter,angle=theta,origin=(0,0),use_radians=True)
            #base=base.difference(cutter2)
        
        arc_x, arc_y = self.get_arc_profile2(rotation_angle=0, num_points=1000, full_gear=True,phi_s1=0,center_x=0,center_y=0)
        cutter=Polygon(list(zip(arc_x, arc_y )))
        angles=np.linspace(0,2*pi,num=num,endpoint=False)
        a=self.a 
        for theta in angles:
            rotation=theta*(1+self.i)
            cutter2=affinity.rotate(cutter,angle=rotation,origin=(0,0),use_radians=True)
            x,y=np.array((-sin(theta),cos(theta)))*a
            cutter2=affinity.translate(cutter2,xoff=x,yoff=y)
            base=base.difference(cutter2)
        if ax:
            ax.plot(arc_x,arc_y,label="cuter")
            x,y=cutter2.exterior.xy
            ax.plot(x,y,label="cuter trans")
        return base

    def ShapelyCycloidal(self,**kwargs):
        """
        Calculates the polygon of a gear with arc-shaped teeth.

        Parameters:
        r (float): Pitch radius of the gear
        z (int): Number of teethqF1
        center_x (float): X-coordinate of the gear center (default is 0)
        center_y (float): Y-coordinate of the gear center (default is 0)
        sector (tuple): Sector of the gear to be considered (default is (pi/2, 2*pi+pi/2))
        Rr (float): Radius of the tooth arc (default is r*sin(pi/(z*2)))

        Returns:
        shapely.geometry.Polygon: Polygon representing the gear
        """
        ra2=self.ra2
        a=self.a
        base=Point(0,0).buffer(ra2)
        
        x1,y1=self.get_arc_profile2(0)
        cutter=Polygon(list(zip(x1,y1)))
        cutter=affinity.translate(cutter,xoff=0,yoff=-a)
        num_points=z=kwargs.get('num_points',360)
        angles=np.linspace(0,2*pi,num=num_points,endpoint=False)        
        for theta in angles:
            rotation_angle=theta+theta*self.i
            cutter2=affinity.rotate(cutter,angle=rotation_angle,origin=(0,0),use_radians=True)
            cutter2=affinity.translate(cutter2,xoff=-sin(theta)*a,yoff=cos(theta)*a)
            base=base.difference(cutter2)
        return base
    
if __name__ == "__main__":
    # Create an instance of EccentricCycloidGear
    from HideOnClick2 import hideOnClick
    gear_pair = EccentricCycloidGear(
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
    
    
    fig, (ax, ax1) = plt.subplots(1, 2)
    ax.set_aspect("equal")
    gear_pair.rF1=16.4
    gear_pair.calculateRequired_rF1(ax=ax,num=10)
    plt.show()
    # vs=gear_pair.v_gear(-5*pi/180)
    # for v in vs:
    #     print(f"{v}")
    # #gear_pair.A()
    # #gear_pair.E()
    # zetaA=gear_pair.zetaA
    # zetaE=gear_pair.zetaE

    # angles=np.linspace(-pi/(gear_pair.z1*2)*0,pi/(gear_pair.z1*2),num=1000)
    # angles=np.linspace(-zetaA,-zetaE,num=1000)
    
    # xis=np.array([(gear_pair.xi(a))*180/pi for a in angles])
    # alphas=np.array([(pi/2-gear_pair.xi(a))*180/pi for a in angles])
    
    # x_coords,y_coords=gear_pair.calculate_path_of_contact(start=zetaA,end=zetaE,num_points=len(angles))
    # lengths = np.zeros(len(angles))
    # for i in range(1, len(angles)):
    #     lengths[i] = lengths[i-1] + np.sqrt((x_coords[i] - x_coords[i-1])**2 + (y_coords[i] - y_coords[i-1])**2)
    # # Create subplots
    # 
    # ax.plot(angles*180/pi,xis,c='r',label="ξ")
    # ax.plot(angles*180/pi,alphas,label="α_t")
    # ax.set_xlabel("angle ξ [°]")
    # plt.title("Transverse pressure angle[°]")
    # ax.legend()
    # ax.grid()
    # ax1.plot(lengths,alphas,label="α_t")
    # ax1.set_xlabel("Path of contact [mm]")
    # ax1.axhline(y=48, color='black', linestyle=':')
    # ax1.axhline(y=27, color='black', linestyle=':')
    # ax1.grid()
    # plt.show()
