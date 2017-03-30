set -x

 pgf90 -Mpreprocess -fast  -DNOArrayVersion  -DNOScalarVersion -DArray_Scalar_Version  -Minfo=accel -r8 -tp=p7 -ta=host  umat_PGI.f Test_Umat.for