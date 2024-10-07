import matplotlib.pyplot as plt
import numpy as np
from math import pi,cos,sin,tan,atan2,atan,sqrt,acos
from shapely import Polygon
from matplotlib.patches import Circle

def Gear2D(z:int,m, alpha=20*pi/180,x=0,**kwargs):
    '''berechnet ein Evolventenverzahntes Zahnrad
    rotation_angle=0 -> zahnluecke bei 0'''
    delta=tan(alpha)-alpha
    hap=kwargs.get("hap",1)
    hfp=kwargs.get("hfp",1.25)
    inner=kwargs.get("inner",False)
    anz=kwargs.get("num",50)
    rotation_angle=kwargs.get("rotation_angle",0*pi/180)
    if z<0:
        inner=True
        z=abs(z)
    
    M=kwargs.get("M",np.array((0,0)))
    ha=m*(hap+x)
    hf=(hfp-x)*m
    d=m*z
    da=m*z+2*ha
    df=d-2*hf
    delta=tan(alpha)-alpha
    
    jt=kwargs.get("j") if kwargs.get("j") is not None else 0
    phi0=pi/z/2-delta+jt/d
    if inner:
        df=d-2*ha
        da=d+2*hf
        phi0=pi/z/2-delta-jt/d
        #print(f"inner={df=} {da=}")

    d_b=cos(alpha)*m*z

    xp=[]
    yp=[]
    
    
    

    #print(f'phi0={phi0*180/pi:.1f}° and delta={delta*180/pi:.1f}')
    def inv(t,rb=1,t0=0):
        return (rb*(t*sin(t+t0)+cos(t+t0)),rb*(sin(t+t0)-t*cos(t+t0)))

    xp=[]
    yp=[]
    xp.append(1*df/2)
    yp.append(0)
    xp.append(cos(phi0)*df/2)
    yp.append(sin(abs(phi0))*df/2)
    alpha_a=acos(d_b/da)

    diameter=np.linspace(d_b,da,num=int(anz/2))
    for d in diameter:
        alpha_i=acos(d_b/d )


    for i in range(anz):
        phi=1/anz*i
        xi,yi=inv(phi,d_b/2,phi0)
        ri=sqrt(xi**2+yi**2)
        if ri<=da/2 and ri>=df/2:
            xp.append(xi)
            yp.append(yi)
        # else:
        #     print(f' db={d_b:.1f}  df={df:.1f} and da={da:.1f} and di={ri*2:.1f}')


    import copy
    xp2=copy.deepcopy(xp)
    yp2=copy.deepcopy(yp)
    xp2.reverse()
    yp2.reverse()
    yp2 = [ x * -1 for x in yp2 ]
    xp2.extend(xp)
    yp2.extend(yp)
    xp2.append(cos(pi/(z))*da/2)
    yp2.append(sin(pi/(z))*da/2)

    debug=kwargs.get("debug",True)
    if debug:
        print(f"da={da:.1f}")
        print(f"db={d_b:.1f}")
        print(f"df={df:.1f}")

    x_gear=[]
    y_gear=[]
    for i in range(abs(z)):
        angle=2*pi/z*i

        for idx,x_i in enumerate(xp2, 0):
            y_i=yp2[idx]
            x_gear.append(cos(angle)*x_i-sin(angle)*y_i)
            y_gear.append(sin(angle)*x_i+cos(angle)*y_i)
    x_gear.append(x_gear[0])
    y_gear.append(y_gear[0])
    gear_arr=np.column_stack((x_gear,y_gear))
    gear=gear_arr@np.array([[cos(rotation_angle),-sin(rotation_angle)],[sin(rotation_angle),cos(rotation_angle)]]).T
    gear= np.column_stack((gear[:,0]+M[0],gear[:,1]+M[1]))
    return gear[:,0],gear[:,1]
    return gear


def involute_coordinates(t, r, a):
    x = r * (np.cos(t + a) + t * np.sin(t + a))
    y = r * (np.sin(t + a) - t * np.cos(t + a))
    return x, y


def InvoluteGear2D(z: int, m, alpha=20 * pi / 180, x=0, **kwargs):
    """
    Calculates the 2D coordinates of an involute gear profile.
    
    Parameters:
    z (int): Number of teeth (negative values for internal gears).
    m (float): Module (size of the gear teeth).
    alpha (float): Pressure angle in radians. Default is 20°.
    x (float): Profile shift coefficient. Default is 0.
    kwargs (optional):
        hap (float): Addendum coefficient (height of the teeth above the pitch circle). Default is 1.
        hfp (float): Dedendum coefficient (height of the teeth below the pitch circle). Default is 1.25.
        inner (bool): Indicates if the gear is internal. Default is False.
        num (int): Number of points used to define the involute curve. Default is 50.
        rotation_angle (float): Angle to rotate the gear in radians. Default is 0.
        j (float): Shift angle coefficient for internal gears. Default is 0.
        M (np.array): Translation vector for the gear center. Default is (0, 0).
        debug (bool): Enables debugging outputs. Default is False.
    
    Returns:
    tuple: (x, y) coordinates of the gear profile after rotation and translation.
    """
    hap = kwargs.get("hap", 1)
    hfp = kwargs.get("hfp", 1.25)
    inner = kwargs.get("inner", False)
    anz = kwargs.get("num", 50)
    rotation_angle = kwargs.get("rotation_angle", 0)
    
    if z < 0:
        inner = True
        z = abs(z)

    M = kwargs.get("M", np.array((0, 0)))
    ha = m * (hap + x)
    hf = (hfp - x) * m
    d = m * z
    da = d + 2 * ha
    df = d - 2 * hf
    delta = tan(alpha) - alpha
    jt = kwargs.get("j", 0)

    phi0 = pi / z / 2 - delta + jt / d
    if inner:
        df = d - 2 * ha
        da = d + 2 * hf
        phi0 = pi / z / 2 - delta - jt / d

    db = cos(alpha) * m * z  # Base circle diameter

    # Initialize gear coordinate lists
    x_gear, y_gear, xp, yp = [], [], [], []

    # Helper function for involute calculation
    def inv(a):
        return tan(a) - a

    print(f"db={db:.1f} \t df={df:.1f} \t da={da:.1f} ")

    # Calculate the involute angle for a given diameter
    def phi_y(dy):
        a_y = np.arccos(db / dy)
        print(f"ay={a_y*180/pi:.1f}")
        sy = dy * ((pi + 4 * x * tan(alpha)) / (2 *abs(z)) + inv(alpha) - inv(a_y))
        print(f"sy={sy:.1f}\t ry={dy/2:.2f} x={x}")
        return sy / dy

    # Generate involute points
    diameter = np.linspace(max(db, df), da, num=int(anz / 2), endpoint=True)
    phis = [phi_y(d) for d in diameter]

    # First point on the dedendum circle
    xp.append(cos(phis[0]) * df / 2)
    yp.append(sin(phis[0]) * df / 2)
    
    for i, phi in enumerate(phis):
        #print(f"phi_{i}={phi:.2f}")
        xi, yi = np.array((cos(phi), sin(phi))) * diameter[i] / 2
        xp.append(xi)
        yp.append(yi)

    # Mirror the involute curve for the other side of the tooth
    xp2, yp2 = xp[::-1], [-y for y in yp[::-1]]

    # Construct the full tooth profile
    x_gear.extend(xp)
    y_gear.extend(yp)
    x_gear.append(cos(pi / z * 0) * da / 2) #append head
    y_gear.append(sin(pi / z * 0) * da / 2)
    x_gear.extend(xp2)
    y_gear.extend(yp2)
    import copy
    xp2=copy.deepcopy(x_gear)
    yp2=copy.deepcopy(y_gear)
    xp2.append(cos(-pi/(z))*df/2)
    yp2.append(sin(-pi/(z))*df/2)

    # Handle multiple teeth by rotating the first tooth profile
    for i in range(1, abs(z)*1):
        angle = 2 * pi / z * i
        for x_i, y_i in zip(xp2, yp2):
            x_gear.append(cos(angle) * x_i + sin(angle) * y_i)
            y_gear.append(-sin(angle) * x_i + cos(angle) * y_i)

    # Close the profile by adding the first point at the end
    x_gear.append(x_gear[0])
    y_gear.append(y_gear[0])

    # Convert to numpy arrays and apply rotation and translation
    gear_arr = np.column_stack((x_gear, y_gear))
    rotation_matrix = np.array([[cos(rotation_angle), -sin(rotation_angle)], 
                                [sin(rotation_angle), cos(rotation_angle)]])
    gear_rotated = gear_arr @ rotation_matrix.T
    gear_translated = gear_rotated + M

    print(f"rotation_agle={rotation_angle*180/pi:.1f}°")
    # Return the final coordinates
    return gear_translated[:, 0], gear_translated[:, 1]




def InvoluteGear2D2(z: int, m, alpha=20 * pi / 180, x=0, **kwargs):
    """
    Calculates the 2D coordinates of an involute gear profile.
    
    Parameters:
    z (int): Number of teeth (negative values for internal gears).
    m (float): Module (size of the gear teeth).
    alpha (float): Pressure angle in radians. Default is 20°.
    x (float): Profile shift coefficient. Default is 0.
    kwargs (optional):
        hap (float): Addendum coefficient (height of the teeth above the pitch circle). Default is 1.
        hfp (float): Dedendum coefficient (height of the teeth below the pitch circle). Default is 1.25.
        inner (bool): Indicates if the gear is internal. Default is False.
        num (int): Number of points used to define the involute curve. Default is 50.
        rotation_angle (float): Angle to rotate the gear in radians. Default is 0.
        j (float): Shift angle coefficient for internal gears. Default is 0.
        M (np.array): Translation vector for the gear center. Default is (0, 0).
        debug (bool): Enables debugging outputs. Default is False.
    
    Returns:
    tuple: (x, y) coordinates of the gear profile after rotation and translation.
    """
    hap = kwargs.get("hap", 1)
    hfp = kwargs.get("hfp", 1.25)
    inner = kwargs.get("inner", False)
    anz = kwargs.get("num", 50)
    rotation_angle = kwargs.get("rotation_angle", 0)
    print(f"rotation_angle={rotation_angle*180/pi:.1f}°")
    s=np.sign(z)
    if z < 0:
        inner = True
        z = abs(z)

    M = kwargs.get("M", np.array((0, 0)))
    ha = m * (hap + x)
    hf = (hfp - x) * m
    d = m * z
    da = d + 2 * ha
    df = d - 2 * hf
    delta = tan(alpha) - alpha
    jt = kwargs.get("j", 0)

    phi0 = pi / z / 2 - delta + jt / d
    if inner:
        da=m*(z+2*x+2*hfp)
        df=m*(z+2*x-2*hap)
    else:
        da=m*(z+2*x+2*hap)
        df=m*(z+2*x-2*hfp)
    print(f"da={da:.2f} df={df:.2f}")
    db = cos(alpha) * m * z  # Base circle diameter

    # Initialize gear coordinate lists
    x_gear, y_gear, xp, yp = [], [], [], []

    # Helper function for involute calculation
    def inv(a):
        return tan(a) - a

    print(f"rb={db/2:.1f} \t rf={df/2:.2f} \t ra={da/2:.1f} ra2={m*(z/2+x-hfp):.2f}")

    # Calculate the involute angle for a given diameter
    s0=m*(pi/2+2*x*tan(alpha)*s)
    print(f"s0={s0}")
    def phi_y(dy):
        a_y = np.arccos(db / dy)
        #print(f"ay={a_y*180/pi:.1f}")
        sy = dy * (s0/d + inv(alpha) - inv(a_y))
        #print(f"sy={sy:.1f}\t ry={dy/2:.2f} x={x}")
        return sy / dy-jt*s

    # Generate involute points
    if inner:
        diameter = np.linspace(min(db, df), da, num=int(anz / 2), endpoint=True)
    else:
        diameter = np.linspace(max(db, df), da, num=int(anz / 2), endpoint=True)
    diameter = np.linspace(max(db, df), da, num=int(anz / 2), endpoint=True)
    phis = [phi_y(d) for d in diameter]

    s=np.sign(z)
    # First point on the dedendum circle
    print(f"df={df:.1f}")
    xp.append(cos(phis[0]) * df/2 )
    yp.append(sin(phis[0]) * df/2 )
    
    for i, phi in enumerate(phis):
        #print(f"phi_{i}={phi*180/pi:.2f} \t r={diameter[i]/2:.2f}")
        xi, yi = np.array((cos(phi), sin(phi))) * diameter[i] / 2
        xp.append(xi)
        yp.append(yi)


    # Mirror the involute curve for the other side of the tooth
    xp2, yp2 = xp[::-1], [-y for y in yp[::-1]]

    ax=kwargs.get("ax",False)

    if ax:
        ax.plot(xp,yp,label="involute1",color="red")
        ax.plot(xp2,yp2,label="involute2",color="magenta")

    # Construct the full tooth profile
    x_gear.extend(xp)
    y_gear.extend(yp)
    x_gear.append(cos(pi / z * 0) * da / 2) #append head
    y_gear.append(sin(pi / z * 0) * da / 2)
    x_gear.extend(xp2)
    y_gear.extend(yp2)
    import copy
    xp2=copy.deepcopy(x_gear)
    yp2=copy.deepcopy(y_gear)
    xp2.append(cos(-pi/(z))*df/2)
    yp2.append(sin(-pi/(z))*df/2)

    if ax:
        ax.plot(x_gear,y_gear,label="xy gear",color="magenta",marker='o',markersize=2,linestyle="",alpha=0.1)

    # Handle multiple teeth by rotating the first tooth profile
    for i in range(1, abs(z)*1):
        angle = 2 * pi / z * i
        for x_i, y_i in zip(xp2, yp2):
            x_gear.append(cos(angle) * x_i + sin(angle) * y_i)
            y_gear.append(-sin(angle) * x_i + cos(angle) * y_i)

    if ax:
        ax.plot(x_gear,y_gear,label="xy gear",color="blue",marker='o',markersize=2,linestyle="",lw=0.1,alpha=0.1)


    # Close the profile by adding the first point at the end
    x_gear.append(x_gear[0])
    y_gear.append(y_gear[0])

    #return np.array(x_gear)+M[0],np.array(y_gear+M[1])
    # Convert to numpy arrays and apply rotation and translation
    gear_arr = np.column_stack((x_gear, y_gear))
    rotation_matrix = np.array([[cos(rotation_angle), -sin(rotation_angle)], 
                                [sin(rotation_angle), cos(rotation_angle)]])
    #gear_rotated = gear_arr @ rotation_matrix.T
    gear_rotated = gear_arr
    gear_translated = gear_rotated + M

    print(f"rotation_agle={rotation_angle*180/pi:.1f}°")
    # Return the final coordinates
    return gear_translated[:, 0], gear_translated[:, 1]


def InvoluteGear2D3(z: int, m, alpha=20 * pi / 180, x=0, **kwargs):
    """
    Calculates the 2D coordinates of an involute gear profile, works for positive and negativ z
    
    Parameters:
    z (int): Number of teeth (negative values for internal gears).
    m (float): Module (size of the gear teeth).
    alpha (float): Pressure angle in radians. Default is 20°.
    x (float): Profile shift coefficient. Default is 0.
    kwargs (optional):
        hap (float): Addendum coefficient (height of the teeth above the pitch circle). Default is 1.
        hfp (float): Dedendum coefficient (height of the teeth below the pitch circle). Default is 1.25.
        inner (bool): Indicates if the gear is internal. Default is False.
        num (int): Number of points used to define the involute curve. Default is 50.
        rotation_angle (float): Angle to rotate the gear in radians. Default is 0.
        j (float): Shift angle coefficient for internal gears. Default is 0.
        M (np.array): Translation vector for the gear center. Default is (0, 0).
        debug (bool): Enables debugging outputs. Default is False.
    
    Returns:
    tuple: (x, y) coordinates of the gear profile after rotation and translation.
    """
    hap = kwargs.get("hap", 1)
    hfp = kwargs.get("hfp", 1.25)
    anz = kwargs.get("num", 50)
    rotation_angle = kwargs.get("rotation_angle", 0)
    print(f"rotation_angle={rotation_angle*180/pi:.1f}°")
    s=np.sign(z)


    M = kwargs.get("M", np.array((0, 0)))
    ra=m*(z/2+x+hap)
    rf=m*(z/2+x-hfp)
    d = m * z
    da = ra*2
    df = rf*2
    jt = kwargs.get("j", 0)
    db = cos(alpha) * m * z  # Base circle diameter

    # Initialize gear coordinate lists
    x_gear, y_gear, xp, yp = [], [], [], []

    # Helper function for involute calculation
    def inv(a):
        return tan(a) - a

    #print(f"rb={db/2:.1f} \t rf={df/2:.2f} \t ra={da/2:.1f} ra2={m*(z/2+x-hfp):.2f}")

    # Calculate the involute angle for a given diameter
    s0=m*(pi/2+2*x*tan(alpha))
    def phi_y(dy):
        a_y = np.arccos(db / dy)     
        sy = dy * (s0/d + inv(alpha) - inv(a_y))
        return sy / dy-jt*s

    # Generate involute points
    if z<0:
        diameter = np.linspace(min(df,db), da, num=int(anz / 2), endpoint=True)
    else:
        diameter = np.linspace(max(df,db), da, num=int(anz / 2), endpoint=True)
    phis = [phi_y(d) for d in diameter]
    
    # First point on the dedendum circle
    xp.append(cos(phis[0]) * df/2 )
    yp.append(sin(phis[0]) * df/2 )
    
    for i, phi in enumerate(phis):
        xi, yi = np.array((cos(phi), sin(phi))) * diameter[i] / 2
        xp.append(xi)
        yp.append(yi)


    # Mirror the involute curve for the other side of the tooth
    xp2, yp2 = xp[::-1], [-y for y in yp[::-1]]

    ax=kwargs.get("ax",False)

    if ax:
        ax.plot(xp,np.array(yp)+M[1],label="involute1",color="red")
        ax.plot(xp2,np.array(yp2),label="involute2",color="magenta")
        ax.plot(cos(pi / z * 0) * da / 2,sin(pi / z * 0) * da / 2,label="Kopfpunkt",marker='o',markersize=3,color="red")
        ax.plot(cos(pi / z ) * df / 2,sin(pi / z) * df / 2,label="mittlerer Fußpunkt",marker='o',markersize=3,color="red")
        FP1=np.array([cos(phis[0]) * df/2,sin(phis[0]) * df/2])
        FP2=np.array([cos(-phis[0]) * min(abs(df),abs(da))/2,sin(-phis[0]) * min(abs(df),abs(da))/2])
        KP1=np.array([cos(phis[-1]) * da/2,sin(phis[-1]) * da/2])
        ax.plot(*KP1 ,label=f"KP1= {KP1}",marker='o',markersize=6,color="green")
        ax.annotate('$K_1$',KP1,color="green", fontsize=8)
        #FP1=np.array([cos(phis[0]) * min(abs(df),abs(da))/2,sin(phis[0]) * min(abs(df),abs(da))/2])
        #FP2=np.array([cos(-phis[0]) * min(abs(df),abs(da))/2,sin(-phis[0]) * min(abs(df),abs(da))/2])
        ax.plot(FP1[0],FP1[1],label=f"Fußpunkt1= {FP1}",marker='o',markersize=6,color="green")
        #ax.plot(cos(-phis[0]) * df/2,sin(-phis[0]) * df/2 ,label=f"Fußpunkt2= {FP2}",marker='o',markersize=6,color="green")
        ax.annotate('F1',FP1,color="green", fontsize=8)
        #ax.annotate('F2',FP2,color="green", fontsize=8)
        circle_f=Circle(M,radius=abs(rf),label=f"$r_f$={rf:.2f}",linestyle='-', color='b', fill=False,alpha=0.5)
        circle_a=Circle(M,radius=abs(ra),label=f"$r_a$={ra:.2f}",linestyle='-', color='b', fill=False,alpha=0.5)
        circle_b=Circle(M,radius=abs(db/2),label=f"$r_b$={db/2:.2f}",linestyle=':', color='b', fill=False,alpha=0.5)
        circle_r=Circle(M,radius=abs(d/2),label=f"$r$={d/2:.2f}",linestyle='-.', color='b', fill=False,alpha=0.5)
        for c in [circle_a,circle_b,circle_f,circle_r]:
            ax.add_patch(c)

    # Construct the full tooth profile
    x_gear.extend(xp)
    y_gear.extend(yp)
    x_gear.append(cos(pi / z * 0) * da / 2) #append head
    y_gear.append(sin(pi / z * 0) * da / 2)
    x_gear.extend(xp2)
    y_gear.extend(yp2)
    import copy
    xp2=copy.deepcopy(x_gear)
    yp2=copy.deepcopy(y_gear)
    xp2.append(cos(-pi/(z))*df/2)
    yp2.append(sin(-pi/(z))*df/2)

    if ax:
        ax.plot(x_gear,y_gear,label="xy gear",color="magenta",marker='o',markersize=2,linestyle="",alpha=0.1)

    # Handle multiple teeth by rotating the first tooth profile
    for i in range(1, abs(z)*1):
        angle = 2 * pi / z * i
        for x_i, y_i in zip(xp2, yp2):
            x_gear.append(cos(angle) * x_i + sin(angle) * y_i)
            y_gear.append(-sin(angle) * x_i + cos(angle) * y_i)

    if ax:
        ax.plot(x_gear,y_gear,label="xy gear",color="blue",marker='o',markersize=2,linestyle="",lw=0.1,alpha=0.1)


    # Close the profile by adding the first point at the end
    x_gear.append(x_gear[0])
    y_gear.append(y_gear[0])

    # Convert to numpy arrays and apply rotation and translation
    gear_arr = np.column_stack((x_gear, y_gear))
    rotation_matrix = np.array([[cos(rotation_angle), -sin(rotation_angle)], 
                                [sin(rotation_angle), cos(rotation_angle)]])
    gear_rotated = gear_arr @ rotation_matrix.T
    #gear_rotated = gear_arr
    gear_translated = gear_rotated + M

    return gear_translated[:, 0], gear_translated[:, 1]



def InvoluteGear_poly(z:int,m, alpha=20*pi/180,x=0,beta=0*pi/180,M=np.array((0,0)),inner:bool=False,**kwargs)->Polygon:
    xg,yg=Gear2D(z,m,alpha,x,beta,M,inner,kwargs)
    coords=coords = list(zip(xg, yg))
    return Polygon(coords)




if __name__ =="__main__":
    hap=1
    hfp=1.25
    m=1
    z=23
    da=m*(z+hap*2)
    df=m*(z-hfp*2)
    alpha=20*pi/180
    db=z*m*cos(alpha)
    M=np.array((23,12))
    x,y=Gear2D(z,m,M=M,hap=hap,hfp=hfp,num=23,rotation_angle=0*pi/180)
    fig,ax=plt.subplots()
    ax.plot(x,y,marker=".", markersize=1,linestyle="")
    x2,y2=InvoluteGear2D(z,m,M=M,hap=hap,hfp=hfp,num=23,rotation_angle=0*pi/180)
    ax.plot(x2,y2,marker=".", markersize=1,linestyle="-",c="magenta")

    z3=-31
    x3,y3=InvoluteGear2D(z3,m,M=M,hap=hap,hfp=hfp,num=23,rotation_angle=1*pi/z3)
    ax.plot(x3,y3,marker=".", markersize=1,linestyle="-",c="orange",label=f"z={z3}")

    ax.hlines(M[1],M[0]-da/2,M[0]+da/2)
    ax.vlines(M[0],M[1]-da/2,M[1]+da/2)
    ax.set_aspect("equal")
    #ax.grid(True)

    ra_circle = Circle(M, da/2,linestyle='-.', color='g', fill=False)
    rf_circle = Circle(M, df/2, linestyle='-.', color='g', fill=False)
    rb_circle = Circle(M, db/2, linestyle='-.', color='g', fill=False)
    for circle in [rf_circle, ra_circle,rb_circle]:
        ax.add_patch(circle)
    ax.legend()
    plt.show()