from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Parameters
ringColor       = (1,0,0)   # Color of the Silicon
backgroundColor = (1,1,1)   # Color of not Silicon

numberRings = 3
ringRadius  = 5             # Ring radius [µm]
waveguideWidth = 0.5        # Width of the waveguide [µm]
couplingProfile = [0.1, 0.01, 0.001]

couplingApoPlot = False

# Methods
def drawRing(center=(0,0), width=0.5, radius=5):

    incircle = plt.Circle((center[0], center[1]), radius=radius-width/2, facecolor=backgroundColor)
    outcircle = plt.Circle((center[0], center[1]), radius=radius+width/2, facecolor=ringColor)
    ax.add_patch(outcircle)
    ax.add_patch(incircle)


ax=plt.gca()

# Rings
import numpy as np
for ring in np.linspace(0,numberRings-1,numberRings):
    drawRing(center=(0, ring*( ringRadius*2 + waveguideWidth)))

# Bus
rect = patches.Rectangle((-(ringRadius+waveguideWidth),-(ringRadius+2*waveguideWidth)), 2*(ringRadius+waveguideWidth), waveguideWidth, facecolor=ringColor)
ax.add_patch(rect)

rect = patches.Rectangle((-(ringRadius+waveguideWidth), (numberRings-1)*(ringRadius*2 + waveguideWidth)), 2*(ringRadius+waveguideWidth), waveguideWidth, facecolor=ringColor)
ax.add_patch(rect)

plt.axis('scaled')
plt.show()

# Plot coupling apodizing profile
# Make this plot share the y axis of the geometry plot
if couplingApoPlot:
    plt.figure()
    plt.plot([1,2,3], couplingProfile)
    plt.xlabel('Ring number [-]')
    plt.ylabel('Coupling coefficient [-]')
    plt.show()