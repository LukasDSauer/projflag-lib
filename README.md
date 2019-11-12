# projflag-lib

This library is a collection of several packages for transforming and visualizing tuples of flags in the real
projective plane. For now, each package has to be installed separately. In the future, it
would be nice to install all of them in a bundle.

## Library structure
* `flagcomplex/` contains the flagcomplex package. It is a package for representing the projective coordinates of a flag.
Furthermore it implements several transformations, such as the eruption flow, the bulge flow and the shear flow.
* `flagcomplexSvg/` is a package intended for visualizing tuples of flags implemented in the flagcomplex package.
* `examples/` contains Jupyter notebooks with commented examples of the code functionality. The `temp/` subfolder is necessary for generating image files on the way.
