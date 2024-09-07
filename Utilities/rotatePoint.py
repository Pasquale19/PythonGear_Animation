import numpy as np
from typing import Tuple

def rotate_point(points, angle, center,use_radians=True):
    """
    Rotate a point around a given center by a specified angle.

    Parameters:
    point (np.ndarray): A 2D array representing the coordinates of the point (shape: (2,)).
    angle (float): The rotation angle in degrees.
    center (np.ndarray): A 2D array representing the center of rotation (shape: (2,)).

    Returns:
    np.ndarray: The new coordinates of the rotated point.
    """
    # Convert angle from degrees to radians
    if use_radians:
        angle_rad=angle
    else:
        angle_rad = np.radians(angle)

    # Create the rotation matrix
    rotation_matrix = np.array([[np.cos(angle_rad), -np.sin(angle_rad)],
                                 [np.sin(angle_rad),  np.cos(angle_rad)]])

    # Translate point to origin (subtract center)
    if isinstance(center,Tuple):
        center=np.array(center)

    print(f"points={points} \n {points.shape}")
    # Ensure points is a 2D array
    points = np.atleast_2d(points)
    print(f"points2D={points} \n {points.shape}")

    translated_points = points - center

    # Rotate the points
    rotated_points = translated_points @ rotation_matrix.T  # Use transpose

    # Translate back to the original center
    new_points = rotated_points + center

     # Check if input is a single point or multiple points
    single_point =  new_points.shape[0] == 1 
    if single_point:
        return new_points[0]
    else:
        return new_points.reshape(-1, 2)  # Ensure output shape is (N, 2)
    

if __name__ =="__main__":
    # Example usage with a single point
    single_point = np.array([2, 3])
    angle = 90  # degrees
    center = np.array([1, 1])

    rotated_single_point = rotate_point(single_point, angle, center,use_radians=False)
    print("Rotated Single Point:", rotated_single_point)
    print("Shape of rotated single point:", rotated_single_point.shape)

    # Example usage with multiple points
    multiple_points = np.array([[2, 3], [4, 5], [1, 0]])
    rotated_multiple_points = rotate_point(multiple_points, angle, center,use_radians=False)
    print("Rotated Multiple Points:\n", rotated_multiple_points)
    print("Shape of rotated multiple points:", rotated_multiple_points.shape)