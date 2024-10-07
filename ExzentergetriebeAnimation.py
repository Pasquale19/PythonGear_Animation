import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Circle, FancyArrow,Arc,Polygon
from matplotlib.transforms import Affine2D
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider
import numpy as np
from Gears.BaseGear2 import InvoluteGear
from math import pi,sin,cos,asin,acos,atan,tan,atan2,hypot,tau
from HideOnClick2 import hideOnClick
from Utilities.Helper import calculate_arcPolygon
from matplotlib.animation import FuncAnimation

from GearAnimation import GearMeshingPlot


if __name__ == "__main__":
    gear_pair = InvoluteGear(
        z1=50,
        z2=-51,
        m=2.15,
        alpha=35*pi/180,
        hap1=0.4,hap2=0.4,
        hfp1=0.5,hfp2=0.5,
        x1=0.2,x2=-0.45
    )

    gear_pair.print_Parameter()

    gear_meshing_plot = GearMeshingPlot(gear_pair)
    gear_meshing_plot.setup_plot()