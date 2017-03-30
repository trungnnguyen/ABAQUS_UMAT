set -x

<<<<<<< .merge_file_LMh3dV
 pgf90 -Mpreprocess -DGPU -r8 -fast -DNOArrayVersion  -DNOScalarVersion -DArray_Scalar_Version -Minfo=accel -ta=multicore -tp=p7 umat_PGI.f Test_Umat.for 
=======
 pgf90 -Mpreprocess -DGPU -fast -r8 -Minfo=accel -ta=multicore -tp=p7 umat_PGI.f Test_Umat.for
>>>>>>> .merge_file_9hUJ9T
