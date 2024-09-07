import numpy as np
from math import pi,sin,cos,tan,acos,asin,atan


def calculate_polar_coordinates(x,y,center_x=0,center_y=0):
    x0,y0=x-center_x,y-center_y
    phi=np.arctan2(y0,x0)-pi/2
    r=np.hypot(x0,y0)
    return phi,r

def convert_to_rack(x,y,radius,center_x=0,center_y=0):
    phi,r=calculate_polar_coordinates(x,y,center_x=center_x,center_y=center_y)
    xr=np.sin(phi)*radius
    yr=r-radius
    return xr,yr

import numpy as np
import matplotlib.pyplot as plt

def plot_polar_coordinates(r, theta, title="Polar Coordinate Plot"):
    """
    Plot polar coordinates using Matplotlib.
    
    Parameters:
    r (array-like): Radial coordinates
    theta (array-like): Angular coordinates in radians
    title (str): Title of the plot (default: "Polar Coordinate Plot")
    """
    # Create a new figure
    fig = plt.figure(figsize=(10, 10))
    
    # Add a polar subplot
    ax = fig.add_subplot(111, projection='polar')
    
    # Plot the polar coordinates
    ax.plot(theta, r)
    
    # Set the title
    ax.set_title(title, fontsize=16)
    
    # Set the angular labels to degrees
    ax.set_thetagrids(np.degrees(ax.get_xticks()))
    
    # Add gridlines
    ax.grid(True)
    
    # Show the plot
    plt.show()


def calculateFlankParameter(x,y,radius=33,epsilon=0.0001):
        num_points=len(x)
        arr=np.zeros(shape=(3,num_points-1))
        for i in range(1,num_points-1,1):
            p0=np.array((x[i-1],y[i-1]))
            p_i=np.array((x[i],y[i]))
            p2=np.array((x[i+1],y[i+1]))
            p_i1=(p_i+p0)/2
            p_i2=(p_i+p2)/2
            d=(-p_i1+p_i2)
            d=d/np.linalg.norm(d)
            alpha=np.arctan2(d[1],d[0])
            #alpha=np.arctan2(d[0],d[1])
            arr[1,i-1]=alpha
            n=(p_i[1])/cos(alpha) #solte stimmen
            if round(abs(alpha)-pi/2,4)==0:
                print(f"{p_i=}")
            
            if p_i[1]**2<epsilon:
                xr=p_i[0]
                print(f"smaller epsilon {p_i[1]=}")
            else:
                xr=p_i[0]+n*sin(alpha)
            phi=xr/radius
            arr[0,i-1]=phi
            arr[0,i-1]=xr
            arr[2,i-1]=n
            print(f"phi={phi*180/pi:.1f}°\t xr={xr:.1f}\t alpha={alpha*180/pi:.1f} \t n={n:.1f} \t d={d} \t p_i1={p_i1} \t p_i2={p_i2}")
            
        return arr

if __name__ =="__main__":
    import matplotlib.pyplot  as plt
    from GearRack import GearRack2D
    theta=np.linspace(pi/4,3*pi/4,num=10,endpoint=False)

    rack_x,rack_y=GearRack2D(modul=1)
    x_rack=np.linspace(rack_x[0],rack_x[-1],40,endpoint=False)
    y_rack=np.interp(x_rack,xp=rack_x,fp=rack_y,right=rack_y[-1],left=rack_y[0])
    info=f"ymax={max(y_rack):.2f}""\t"f"ymin={min(y_rack):.2f}"
    arr=calculateFlankParameter(x_rack,y_rack,radius=10)
    x1=arr[0,:]
    alpha=arr[1,:]
    n=arr[2,:]
    info2=f"nmax={max(n):.2f}" "\t" f"nmin={min(n):.2f}"
    fig,(ax,ax2)=plt.subplots(2)
    ax.set_aspect("equal")
    ax.set_title(info+ "\n"+ info2)
    ax.plot(rack_x,rack_y)
    ax.plot(x_rack,y_rack,marker='o',markersize=2,linestyle="",label="interp rack",color="black")
    x,y=np.cos(theta)*100,np.sin(theta)*100
    #xr,yr=convert_to_rack(x,y,100)
    ax2.plot(x1,alpha*180/pi,label=r'$\alpha$',marker="o",markersize="3")
    ax2.set_ylabel("$alpha_P$ in °")
    ax2.set_ylabel(r"$\alpha_P$ in ° ($\alpha$=" f"{alpha[0]*180/pi}")
    ax2.plot(x1,n,label="n",color="red",marker="o",markersize="3")
    # ax.plot(xr,yr,marker='o')
    # 
    # ax.plot(x,y,marker='o')
    ax.legend()
    ax2.legend()
    import sys
    import os
    up1 = os.path.abspath('..')
    up1 = os.path.abspath('')
    sys.path.insert(0, up1)
    # for p in sys.path:
    #     print(p)
    
    from HideOnClick2 import hideOnClick
    hideOnClick(fig,ax)
    
    plt.show()