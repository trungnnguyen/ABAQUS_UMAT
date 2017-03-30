set -x

 pgf90 -Mpreprocess  -Minfo=ccff -DNOArrayVersion  -DNOScalarVersion -DArray_Scalar_Version -tp=p7 -r8 -ta=host umat_PGI.f Test_Umat.for