from dataclasses import dataclass, field
from math import pi, cos, sin, acos, atan2, sqrt,atan,asin
from typing import Tuple, Dict
import numpy as np
from matplotlib import pyplot as plt
from shapely.geometry import Point, Polygon
from shapely import affinity as affinity
from intersection import find_intersection
#import matplotlib.axes as ax
#from deprecated import deprecated
from matplotlib.patches import Circle, FancyArrow,Arc,Polygon
from eccentric_cycloid_gear import EccentricCycloidGear


class EccentricCycloidGear_Test(EccentricCycloidGear):
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

    def __init__(self, z1: int, z2: int, a: float, **kwargs):
        """
        Initializes the gear pair with given parameters and calculates derived values.

        Args:
            z1 (int): Number of teeth of arc gear
            z2 (int): Number of teeth of cycloid gear
            a (float): Center distance
            **kwargs: Additional parameters for gear configuration
        """
        super().__init__(z1, z2, a, **kwargs)
        # self.z1 = z1
        # self.z2 = z2
        # self.a = a

        # self.m = self.a / (self.z1 / 2 + self.z2 / 2)
        # self.d1 = self.m * self.z1
        # self.d2 = self.m * self.z2
        # self.r1 = self.d1 / 2
        # self.r2 = self.d2 / 2
        
        # self.i = self.z2 / self.z1
        # self.rw1 = self.r1
        # self.rw2 = self.r2

        # self.rA_star = kwargs.get('rA_star', 1.0)
        # self.c_star = kwargs.get('c_star', 0.1)
        # self.st_star = kwargs.get('st_star', 1.0)
        # self.phi_j1 = kwargs.get('phi_j1', 0.0)
        # self.phi_As = kwargs.get('phi_As', 60*pi/180)
        # self.phi_Ae = kwargs.get('phi_Ae', 170*pi/180)
        # self.lambda_ = kwargs.get('lambda_', 0.97)
        # self.delta = kwargs.get('delta', 0.1) #equidistant distance
        # self.e = self.r1*self.lambda_
        # self.calculate_arc_gear_geometry()
        # self.calculate_cycloid_gear_geometry()
        # self.calculate_A()
        # self.calculate_E()
    
    
    
    def ShapelyGear(self,ax,**kwargs):
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
    

        #gear=Point(0,0).buffer(r_circle)
        gear=Point(center_x,center_y).buffer(e*0.9)
        print("Startkreis zu groß")
        angles_foot=[a+pi/z for a in angles]
        rF=self.rF1
        qF1=self.qF1
        
        ce=Circle((0,0),radius=e, linestyle='-.', color='g', fill=False,label=f"e{e:.2f}")
        ax.add_patch(ce)
        c2=Circle((0,0),radius=qF1, linestyle='-.', color='g', fill=False,label=f"qF1{qF1:.2f}")
        ax.add_patch(c2)
        #erstellen der Zahnköpfe
        for i,a in enumerate(angles):
            a1=a-self.phi_s1
            A1=cos(a1)*e+center_x,sin(a1)*e+center_y
            circle_a=Point(A1[0],A1[1]).buffer(rA)
            ax.scatter(A1[0],A1[1],color="black",label=f"A1{i}")
            gear=gear.union(circle_a)
            a2=a+self.phi_s1
            A2=cos(a2)*e+center_x,sin(a2)*e+center_y
            circle_a=Point(A2[0],A2[1]).buffer(rA)
            gear=gear.union(circle_a)
            ax.scatter(A1[0],A1[1],color="black",label=f"A2{i}")
            c1=Circle(A1,radius=rA, linestyle='-.', color='b', fill=False,label=f"c1{i}")
            ax.add_patch(c1)
            c2=Circle(A2,radius=rA, linestyle='-.', color='g', fill=False,label=f"c2{i}")
            ax.add_patch(c2)
        
        #erstellen der Zahnfüße
        for i,a in enumerate(angles_foot):
            F1=cos(a)*qF1+center_x,sin(a)*qF1+center_y
            circle_f=Point(F1[0],F1[1]).buffer(rF)
            gear=gear.difference(circle_f)
            c2=Circle(F1,radius=rF, linestyle='-.', color='g', fill=False,label=f"cf{i}")
            ax.add_patch(c2)

        
        cr=Circle((0,0),radius=self.rInter, linestyle=':', color='m', fill=False,label=f"rInter{self.rInter:.1f}",linewidth=2)
        ax.add_patch(cr)
        # for i, (x, y) in enumerate(list(zip(xp,yp))):
        #     circle_a = Point(x, y).buffer(rA)#Kopf
        #     gear=gear.union(circle_a)
        r_f=self.rf1
        #gear=gear.union(Point(center_x,center_y).buffer(r_f))
        
            
        # Cut out tips
        circle=Point(center_x,center_y).buffer(da/2)
        gear=gear.intersection(circle)
        return gear
    


    
if __name__ == "__main__":
    # Create an instance of EccentricCycloidGear
    gear_pair = EccentricCycloidGear_Test(
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

    #gear_pair.A()
    #gear_pair.E()
    zetaA=gear_pair.zetaA
    zetaE=gear_pair.zetaE

    angles=np.linspace(-pi/(gear_pair.z1*2)*0,pi/(gear_pair.z1*2),num=1000)
    angles=np.linspace(-zetaA,-zetaE,num=1000)
    
    xis=np.array([(gear_pair.xi(a))*180/pi for a in angles])
    alphas=np.array([(pi/2-gear_pair.xi(a))*180/pi for a in angles])
    
    x_coords,y_coords=gear_pair.calculate_path_of_contact(start=zetaA,end=zetaE,num_points=len(angles))
    lengths = np.zeros(len(angles))
    for i in range(1, len(angles)):
        lengths[i] = lengths[i-1] + np.sqrt((x_coords[i] - x_coords[i-1])**2 + (y_coords[i] - y_coords[i-1])**2)
    # Create subplots
    from HideOnClick2 import hideOnClick
    fig, ax = plt.subplots()
    g=gear_pair.ShapelyGear(ax)
    x,y=g.exterior.xy
    ax.plot(x,y,c='r',label="gear")
    plt.title("Transverse pressure angle[°]")
    ax.legend()
    hideOnClick(fig,ax)
    ax.grid()
    plt.show()
