set -x

 pgf90 -Mpreprocess -DGPU -r8 -fast -Minfo=accel -ta=tesla -tp=p7  umat.f Test_Umat.for