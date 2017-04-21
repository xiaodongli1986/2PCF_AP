
gfortran -o LSS_test LSS_main_LSS_ximu.f90 LSS_main_program.f90

gfortran -o LSS_MCMC LSS_main_LSS_ximu.f90 LSS_main_ChisqLike.f90 $delm
