from flagcomplex import FlagComplex
from .DrawingUtility import drawl, draw_flag, drawtri, draw_polygon, drawcirc
from flagcomplex.ProjGeometryUtility import line_intersection, connecting_line
from flagcomplex.EuklGeometryUtility import rotate_vectors
import numpy as np


class FlagComplexSvg(FlagComplex):

    def __init__(self):
        super().__init__()

        self.image_bottom_line = None
        self.image_top_line = None
        self.image_left_line = None
        self.image_right_line = None

    def initialize_canvas(self, width=800, height = 450, origin=(-400, -225)):
        """
        Sets up a drawSvg object that can be used as a canvas for drawing.

        :param width:
        :param height:
        :param origin:
        :return:
        """
        rotation_matrix = rotate_vectors(np.array([0, 0, 1]), self.projection_plane)

        self.image_bottom_line = connecting_line(np.array([origin[0]/100, origin[1]/100, 1]), np.array([(origin[0] + width)/100, origin[1]/100, 1]))
        self.image_top_line = connecting_line(np.array([origin[0]/100, (origin[1] + height)/100, 1]), np.array([(origin[0] + width)/100, (origin[1] + height)/100, 1]))
        self.image_left_line = connecting_line(np.array([origin[0]/100, (origin[1])/100, 1]), np.array([(origin[0])/100, (origin[1] + height)/100, 1]))
        self.image_right_line = connecting_line(np.array([(origin[0] + width)/100, (origin[1])/100, 1]), np.array([(origin[0] + width)/100, (origin[1] + height)/100, 1]))

        self.image_bottom_line = np.matmul(rotation_matrix, self.image_bottom_line)
        self.image_top_line = np.matmul(rotation_matrix, self.image_top_line)
        self.image_left_line = np.matmul(rotation_matrix, self.image_left_line)
        self.image_rigth_line = np.matmul(rotation_matrix, self.image_right_line)

        import drawSvg as draw
        d = draw.Drawing(width, height, origin=origin)
        # White background for the pictures
        drawcirc(d, [0, 0], 8)
        # Render size for the pictures (resolution in terms of image width)
        d.setRenderSize(w=2560)

        return d

    def visualize_polygon(self, d, polygon, col):
        """
        Projects the points of the polygon to two dimensions and prints them to the draw object.

        :param d: a DrawSvg object
        :param polygon: a list of points (3dim numpy arrays)
        :param col: a color string
        :return:
        """
        drawpoly = [self.get_two_dimensional_point(p) for p in polygon]
        draw_polygon(d, [drawpoly], col=[col], wop=True)

    def visualize_connecting_line(self, d, point1, point2, col):
        """
        Draws the line between two points. It is possible that one of the two points is at infinity, but the line is
        still visible. The function takes care of this case by calculating a third point on the line (intersection
        points with the image frame.

        :param d:
        :param point1: a 3dim numpy array
        :param point2: a 3dim numpy array
        :param col: a color string
        :return:
        """
        p1_infinite = False
        p2_infinite = False
        drawp1 = None
        drawp2 = None
        try:
            drawp1 = self.get_two_dimensional_point(point1)
        except AssertionError:
            p1_infinite = True
        try:
            drawp2 = self.get_two_dimensional_point(point2)
        except AssertionError:
            p2_infinite = True
        if not p1_infinite:
            if p2_infinite:
                try:
                    drawp1 =self.get_two_dimensional_point(line_intersection(self.image_bottom_line, connecting_line(point1, point2)))
                    drawp2 = self.get_two_dimensional_point(line_intersection(self.image_top_line, connecting_line(point1, point2)))
                except AssertionError:
                    drawp1 = self.get_two_dimensional_point(
                        line_intersection(self.image_left_line, connecting_line(point1, point2)))
                    drawp2 = self.get_two_dimensional_point(
                        line_intersection(self.image_right_line, connecting_line(point1, point2)))
                drawl(d, drawp2, drawp1, col=col)
            else:
                drawl(d, drawp1, drawp2, col=col)
        else:
            if p2_infinite:
                print("Warning! The line between " + str(point1) + " and " + str(point2) + " is a line at infinity.")
            else:
                try:
                    drawp1 =self.get_two_dimensional_point(line_intersection(self.image_bottom_line, connecting_line(point1, point2)))
                    drawp2 = self.get_two_dimensional_point(line_intersection(self.image_top_line, connecting_line(point1, point2)))
                except AssertionError:
                    drawp1 = self.get_two_dimensional_point(
                        line_intersection(self.image_left_line, connecting_line(point1, point2)))
                    drawp2 = self.get_two_dimensional_point(
                        line_intersection(self.image_right_line, connecting_line(point1, point2)))
                drawl(d, drawp2, drawp1, col=col)

    def visualize(self, d, col=dict(), with_inner_triangles = False, with_middle_triangles = True,
                  with_helper_lines=False, with_ellipse = False, with_label=False, ellipse_flags = None,
                  with_flags = True):
        """
        Visualizes the flag complex. Possible keys for defining colors are currently:
        inner, middle, ellipse, helper, flags

        :param d: a drawSvg.drawing object for visualizing points, lines etc.
        :param col: a dictionary of strings defining colors of different objects.
        :param with_inner_triangles: visualizes the inner triangles u_i (compare [WZ17] for notation)
        :param with_middle_triangles: visualizes the middle triangle p_i
        :param with_helper_lines: a boolean value, visualize the helper lines defining the inner triangle
        :param with_label: a boolean value, set True for labeling the points self.ps[i] with numbers
        :return:
        """
        self.d = d
        self.get_projected_ps()
        self.get_projected_qs()

        if "middle" not in col:
            col["middle"] = "blue"
        if "inner" not in col:
            col["inner"] = "red"
        if "ellipse" not in col:
            col["ellipse"] = "green"
        if "helper" not in col:
            col["helper"] = "grey"
        if "flags" not in col:
            col["flags"] = "black"

        if with_ellipse:
            flag1 = self.get_flag(ellipse_flags[0])
            flag2 = self.get_flag(ellipse_flags[1])
            point = self.ps[ellipse_flags[2]]
            ellipse = self.get_conic_through_flags_and_point(flag1, flag2, point, resolution=32)
            self.visualize_polygon(d, ellipse, col=col["ellipse"])

        # Draw the helper lines which define the inner triangle
        if with_helper_lines:
            for triangle in self.triangles:
                tps = self.get_middle_triangle(triangle)
                tqs = self.get_outer_triangle(triangle)
                for i in range(3):
                    self.visualize_connecting_line(d, tps[i], tqs[i], col=col["helper"])

        # Draw the inner triangles (u_1, u_2, u_3)
        if with_inner_triangles:
            for triangle in self.triangles:
                us = self.get_inner_triangle(triangle)
                drawus = [self.get_two_dimensional_point(us[i]) for i in range(3)]
                drawtri(d, drawus, col=col["inner"])

        # Draw the points p and the flags (i.e. the lines that the points are on).
        if with_flags:
            for i in range(self.n):
                if with_label:
                    label = "p"+str(i)
                else:
                    label = ""
                draw_flag(d, self.drawps[i], self.drawqs[i], col=col["flags"], label=label)
                #drawpt(d, self.drawqs[i], col=col)

        # Draw the triangle connecting the points p
        if with_middle_triangles:
            for triangle in self.triangles:
                points = self.get_middle_triangle(triangle)
                points = [self.get_two_dimensional_point(points[i]) for i in range(3)]
                drawtri(d, points, col=col["middle"], pts=True)

        return d