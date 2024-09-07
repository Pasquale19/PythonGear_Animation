import numpy as np
import matplotlib.pyplot as plt

# Constants
num_teeth = 3
base_radius = 10.0
addendum = 2.0
dedendum = 2.5
pressure_angle_base = 20  # degrees
angle_increment = 2 * np.pi / num_teeth

# Function to calculate the pressure angle as a function of the rotation angle
def pressure_angle(theta):
    # Example: linear increase in pressure angle with respect to the rotation angle
    return pressure_angle_base + 5 * np.sin(3 * theta)  # degrees

# Function to calculate the tooth profile
def calculate_tooth_profile():
    profile_points = []
    
    for i in range(num_teeth):
        theta_start = i * angle_increment
        theta_end = theta_start + angle_increment

        # Arc for the dedendum (gear foot)
        theta = np.linspace(theta_start, theta_start + angle_increment / 2, 100)
        x_dedendum = (base_radius - dedendum) * np.cos(theta)
        y_dedendum = (base_radius - dedendum) * np.sin(theta)
        
        # Arc for the addendum (tooth head)
        theta = np.linspace(theta_start + angle_increment / 2, theta_end, 100)
        x_addendum = (base_radius + addendum) * np.cos(theta)
        y_addendum = (base_radius + addendum) * np.sin(theta)
        
        profile_points.append((x_dedendum, y_dedendum))
        profile_points.append((x_addendum, y_addendum))
    
    return profile_points

# Generate data
theta_values = np.linspace(0, 2 * np.pi, 1000)
pressure_angles = pressure_angle(theta_values)
gear_profile = calculate_tooth_profile()

# Plotting
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Plot the pressure angle
axs[0].plot(theta_values, pressure_angles, label="Pressure Angle")
axs[0].set_title("Pressure Angle Plot")
axs[0].set_xlabel("Rotation Angle (radians)")
axs[0].set_ylabel("Pressure Angle (degrees)")
axs[0].grid(True)
axs[0].legend()

# Plot the gear profile
for x, y in gear_profile:
    axs[1].plot(x, y, 'b')

# Add base circle for reference
circle = plt.Circle((0, 0), base_radius, color='r', fill=False, linestyle='--')
axs[1].add_artist(circle)

axs[1].set_title("Gear Tooth Profile")
axs[1].set_aspect('equal')
axs[1].grid(True)

plt.show()
