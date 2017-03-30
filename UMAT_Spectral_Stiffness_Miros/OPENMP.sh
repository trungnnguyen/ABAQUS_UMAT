set -x

 pgf90 -Mpreprocess -DOMP -DNOArrayVersion  -DNOScalarVersion -DArray_Scalar_Version -fast -Minfo=accel,mp -r8 -tp=p7 -mp -ta=host umat_PGI.f Test_Umat.for