from FlagComplex import FlagComplex
from Examples import construct_from_ratio
from ProjGeometryUtility import connecting_line, transform_four_points
import numpy as np
import drawSvg as draw
from DrawingUtility import drawcirc
from FlagTesselator import FlagTesselator

fcomplex = FlagComplex()

other_points = [np.array([1, 1, 1]), np.array([1, -1, 1]), np.array([-1, -1, 1]), np.array([-1, 1, 1])]
p0 = np.array([0, 1, 1])
q0 = other_points[0]

p1 = np.array([1, 0, 1])
q1 = other_points[1]

p2 = np.array([0, -1, 1])
q2 = other_points[2]

p3 = np.array([-1, 0, 1])
q3 = other_points[3]



fcomplex.add_flag(p0,q0)
fcomplex.add_flag(p1,q1)
fcomplex.add_flag(p2,q2)
fcomplex.add_flag(p3,q3)

fcomplex.set_projection_plane(np.array([0,0,1]))

import copy
fcomplex1 = copy.deepcopy(fcomplex)
ftesselator1 = FlagTesselator(fcomplex1, steps = 4)

import copy

d = fcomplex.initialize_canvas()

fcomplex.create_triangulation()
fcomplex.draw_complex()

ftesselator = FlagTesselator(fcomplex, steps = 5)
d1 = copy.deepcopy(d)
ftesselator.print_image(old = d1, with_tiles = False)
fcomplex.visualize(d1,with_helper_lines=True, with_inner_triangles = False, with_label= True)
#d1.savePng("ef_4.png")
d1