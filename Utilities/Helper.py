import numpy as np

def calculate_arc(center, radius, start_angle, end_angle, num_points=100):
    '''calculates the coordinates of an arc, the angles must be defined in deg'''
    theta = np.linspace(np.deg2rad(start_angle), np.deg2rad(end_angle), num_points)
    x = center[0] + radius * np.cos(theta)
    y = center[1] + radius * np.sin(theta)
    return x, y


def calculate_arcPolygon(center, radius, start_angle, end_angle, num_points=100,use_radians=False):
    '''Calculates the coordinates of an arc, the angles must be defined in degrees'''
    # Convert angles from degrees to radians
    if use_radians:
        theta = np.linspace(start_angle, end_angle, num_points)
    else:
        theta = np.linspace(np.deg2rad(start_angle), np.deg2rad(end_angle), num_points)
    
    # Calculate the arc points directly
    xy = np.array(center) + np.array([np.cos(theta), np.sin(theta)]).T * radius
    
    # Append the center point as the last point
    xy = np.vstack([xy, center])
    
    return xy