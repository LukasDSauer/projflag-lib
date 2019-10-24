# In a programming environment, you may not need the DrawingUtilities and the tutorials. This bash
# script generates a raw version of the project without the unnecessary files.

mkdir flagcomplex-raw/

cp FlagComplex.py flagcomplex-raw/
cp EuklGeometryUtility.py flagcomplex-raw/
cp ProjGeometryUtility.py flagcomplex-raw/
cp LinAlgebraUtility.py flagcomplex-raw/

printf "Successfully generated a raw flagcomplex project."