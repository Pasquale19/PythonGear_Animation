import numpy as np
from math import pi,cos,sin,tan,tau
from typing import Tuple, Dict


def GearRack2D(modul: float=1, pressure_angle_alpha: float=20*pi/180, profile_shift_x: float=0, 
                                    center_x: float=0, center_y: float=0, number_of_teeth: int=1) -> tuple[np.ndarray, np.ndarray]:
    pitch = modul * np.pi
    addendum = (1 + profile_shift_x) * modul
    dedendum = (1.25 - profile_shift_x) * modul
    pressure_angle_rad = pressure_angle_alpha
    
    x_coords = [center_x-pi*modul/2-(number_of_teeth-1)*modul*pi/2]
    y_coords = [center_y - dedendum]

    for i in range(number_of_teeth):
        start_x = center_x + i * pitch -pitch/4
        
        # Left side of tooth
        x_coords.extend([start_x- dedendum * np.tan(pressure_angle_rad), start_x + addendum * np.tan(pressure_angle_rad)])
        y_coords.extend([center_y - dedendum, center_y + addendum])
        
        # Top of tooth
        x_coords.append(start_x + pitch/2 - addendum * np.tan(pressure_angle_rad))
        y_coords.append(center_y + addendum)
        
        # Right side of tooth
        x_coords.append(start_x + pitch/2 + dedendum * np.tan(pressure_angle_rad)+(number_of_teeth-1)*modul*pi/2)
        y_coords.append(center_y - dedendum)
# Closing the polygon
    x_coords.append(center_x+modul*pi/2)
    y_coords.append(y_coords[-1])

    return np.array(x_coords), np.array(y_coords)

def KreisbogenRack(ra: float = 0, num_points: int = 100, **kwargs) -> Tuple[np.ndarray, np.ndarray]:
    '''method to calculate the profile of an gear rack with circular teeth'''
    F1=np.array((-ra*2,0))
    F2=np.array((ra*2,0))
    center_x=kwargs.get("center_x",0)
    center_y=kwargs.get("center_y",0)
    def calculate_coords(angles: float,P=np.array((0,0))) -> Tuple[float, float]:
        
        return P+ [np.array((np.cos(a ), np.sin(a ))) * ra for a in angles]
    coords=[]
    num=max(1, num_points // 2)
    angles_footl = np.linspace(-pi/2, 0, num=num//4, endpoint=False)
    coords.extend(calculate_coords(angles_footl,F1))
    angles_head = np.linspace(pi,0, num=num//2, endpoint=False)
    coords.extend(calculate_coords(angles_head,np.array((0,0))))
    angles_footr= np.linspace(pi, 3/4*tau, num=num//4, endpoint=False)
    coords.extend(calculate_coords(angles_footr,F2))
    coords = np.array(coords) + np.array((center_x, center_y))
    return coords[:, 0], coords[:, 1]


def plot_kreisbogenrack(ra_values, num_points=100):
    """
    Function to plot the KreisbogenRack profile for different `ra` values.
    
    Parameters:
    - ra_values (list): List of `ra` values to plot.
    - num_points (int): Number of points to calculate for each profile.
    """
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    
    for ra in ra_values:
        x, y = KreisbogenRack(ra=ra, num_points=num_points)
        plt.plot(x, y, label=f'ra = {ra}')
    
    plt.title("KreisbogenRack Profile")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend()
    plt.grid(True)
    plt.axis('equal')
    plt.show()

if __name__ == "__main__":
    ra_values = [0.5, 1.0, 1.5]
    plot_kreisbogenrack(ra_values, num_points=100)
