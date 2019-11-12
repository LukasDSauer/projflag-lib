# flagcomplex

This is a project to perform geometric transformations tuples of flags in the real projective plane. A flag is simply a point lying on a projective line.

## Project structure
* `flagcomplex` contains the package's source files. 
  * `FlagComplex` contains the FlagComplex class. It implements all the structure necessary for saving a tuple of flags, using "three-dimensional" projective coordinates. It also contains transformations (eruption flow, bulge flow and shear flow) and visualization functions for projecting the points to two-dimensional coordinates.
  * `FlagTesselator`: A positive tuple of flags defines a polygon. By iteratively reflecting these polygon, a convex set can be generated. This class generates such a convex set.
  * `EuklGeometryUtility` and `ProjGeometryUtility` contain helper functions for calculations.
  * `AutoComplex` contains some helper functions to generate exemplary flag tuples: random tuples and tuples with a certain triple ratio.
