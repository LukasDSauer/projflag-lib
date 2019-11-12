from flagcomplex import FlagTesselator
from .DrawingUtility import draw_path


class FlagTesselatorSvg(FlagTesselator):
    def __init__(self, flagcomplex, steps=3):
        super().__init__(self, flagcomplex, steps=3)

    def print_image(self, color="#006666", with_hull=True, with_initial=True, with_tiles=False, without_points=True,
                    fill_convex=False, old=None,
                    data=None, image_offset=0):
        """
        This function performs the tesselation and draws it. It can also draw on existing data.

        :param color: the color used for drawing the convex set
        :param with_hull: a boolean value. Do you want to draw the complex hull?
        :param with_initial: a boolean value. Do you want to highlight the initial polygon?
        :param without_points: a boolean value. Do you want to highlight the points of the polygons?
        :param old: an existing drawSvg object that you want to draw on
        :param data: existing data for the initial polygon, for the convex hull and for the tiles
        :param image_offset:
        :return:
        """
        if data is not None:
            ini, hull, tiles = data[0], data[1], data[2]
        else:
            ini, hull, tiles = self.generate_tesselation()
        if old is not None:
            d = old
        else:
            d = self.generate_canvas(hull, image_offset)
        if with_tiles:
            for word in self.words:
                draw_path(d, tiles[word], col=color, wop=without_points, fill=0.1, width=0.2)
        if with_initial:
            draw_path(d, ini, col=color, wop=without_points, fill=1, width=0)
        if with_hull:
            if fill_convex:
                fill = 0.2
            else:
                fill = 0
            draw_path(d, hull, col=color, wop=without_points, fill=fill, width=0.4)
        # # d.setPixelScale(2)  # Set number of pixels pert[0][0] geometry unit

        return d