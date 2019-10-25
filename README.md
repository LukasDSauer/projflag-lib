# FlagComplex

This is a project to visualize tuples of flags in the real projective plane. A flag is simply a point lying on a projective line.

## Project structure
* `FlagComplex` contains all the structure necessary for saving a tuple of flags, using "three-dimensional" projective coordinates. It also contains transformations (eruption flow, bulge flow and shear flow) and visualization functions for projecting the points to two-dimensional coordinates.

* `FlagTesselator`: A positive tuple of flags defines a polygon. By iteratively reflecting these polygon, a convex set can be generated. This class generates such a convex set.

* `tutorials/` contains Jupyter notebooks with commented examples of the code functionality. The `temp/` subfolder is necessary for generating image files on the way.

* `DrawingUtility` contains helper functions for drawing SVG pictures in Jupyter.

* `Model3DMaker` contains helper functions to generate 3D models.

* `EuklGeometryUtility` and `ProjGeometryUtility` contain helper functions for calculations.

* `Examples` contains some helper functions to generate exemplary flag tuples: random tuples and tuples with a certain triple ratio.
