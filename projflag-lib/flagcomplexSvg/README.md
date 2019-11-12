# flagcomplexSvg

This is a project to visualize tuples of flags in the real projective plane. A flag is simply a point lying on a projective line.
It uses `drawSvg` for the purpose of visualization and the `flagcomplex` package for the geometric calculations.
It can nicely be used together with Jupyter notebooks.

## Project structure

* `flagcomplex` contains the package's source files. 
  * `FlagComplexSvg` contains the FlagComplexSvg class. It is a child of the `flagcomplex.FlagComplex` class and implements
   functions necessary for visualizing SVG images of the tuples of flags.
  * `DrawingUtility` contains helper functions for drawing SVG pictures with drawSvg.
  * `Model3DMaker` contains helper functions to generate 3D models.
