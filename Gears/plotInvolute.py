import numpy as np
import matplotlib.pyplot as plt
from math import pi

# Define the function
def inv(t, rb=1, t0=0):
    x = rb * (t * np.sin(t + t0) + np.cos(t + t0))
    y = rb * (np.sin(t + t0) - t * np.cos(t + t0))
    return x, y

# Set the range for t
fig,ax=plt.subplots()
t = np.linspace(0, np.pi, 1000)
rb=-1
# Get the x and y values from the function
x, y = inv(t,rb=rb)

# Plot the function
ax.plot(x, y, label="inv(t)")

x2, y2 = inv(t,rb=-rb,t0=pi)
ax.plot(x2,y2,label="involut pi")
# Plot the circle with radius rb
circle = plt.Circle((0, 0), radius=abs(rb), color='r', fill=False, linestyle='--', label='Circle with radius rb')

# Add the circle to the plot
plt.gca().add_patch(circle)

# Set equal aspect ratio for the plot to ensure the circle is not distorted
plt.gca().set_aspect('equal', adjustable='box')

# Add labels and title
plt.xlabel('x')
plt.ylabel('y')
plt.title('Plot of inv(t) and a Circle with radius rb')
plt.legend()

# Show the plot
plt.grid(True)
plt.show()
