set -x

 pgf90 -Mpreprocess -DGPU -r8 -fast -DNOArrayVersion  -DNOScalarVersion -DArray_Scalar_Version  -Minfo=accel -ta=tesla:pinned,fastmath -tp=p7 umat_PGI.f Test_Umat.for 