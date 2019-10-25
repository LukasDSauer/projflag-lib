# In a programming environment, you may not needus the tutorials. This bash
# script generates a raw version of the project without the unnecessary files.

mkdir flagcomplex/

cp FlagComplex.py flagcomplex/
cp EuklGeometryUtility.py flagcomplex/
cp ProjGeometryUtility.py flagcomplex/
cp DrawingUtility.py flagcomplex/

printf "Successfully generated a raw flagcomplex project."