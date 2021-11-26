# projflag-lib

This library is a collection of several packages for transforming and visualizing tuples of flags in the real
projective plane. For now, each package has to be installed separately. In the future, it
would be nice to install all of them in a bundle.

This library was developed by Lukas D. Sauer at the Heidelberg Institute of Theoretical Studies in 2019. A special thanks
goes to Sven Grützmacher, whose Master's thesis code contributed to this library.

## Setup

If you already have Python installed on your machine, the flagcomplex and the flagcomplexSvg package can be installed with pip.

1. Open a command prompt of your choice. (Optionally [activate a virtual environment](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment).)
2. Type `pip install -e git+https://github.com/LukasDSauer/projflag-lib.git#egg=flagcomplex&subdirectory=projflag-lib/flagcomplex` for installing the current version.

## Library structure
* `projflag-lib/flagcomplex/` contains the flagcomplex package. It is a package for representing the projective coordinates of a flag.
Furthermore it implements several transformations, such as the eruption flow, the bulge flow and the shear flow.
* `projflag-lib/flagcomplexSvg/` is a package intended for visualizing tuples of flags implemented in the flagcomplex package.
* `examples/` contains Jupyter notebooks with commented examples of the code functionality. The `temp/` subfolder is necessary for generating image files on the way

## Versioning of library releases
This library consists of two packages, flagcomplex and flagcomplexSvg. Both of the packages have independent version numbering in the format `1.1.1`.
The library releases are tagged `v1.2.3` or `v1.2.3a`, where ...
* ... `v` means "version",
* ... `1.2.3` is the version number of the flagcomplex package, and
* ... `a` is a lower case letter that is appended in case the flagcomplexSvg package is updated without the flagcomplex package being updated.

An exemplary versioning history could be `v1.0.0` (flagcomplex 1.0.0, flagcomplexSvg 1.0.0), `v1.0.1` (flagcomplex 1.0.1, flagcomplexSvg 1.0.0), `v1.0.1a` (flagcomplex 1.0.1, flagcomplexSvg 1.0.1) and finally `v2.0.0` (flagcomplex 2.0.0, flagcomplexSvg 1.0.1).
