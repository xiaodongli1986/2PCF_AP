#gfortran -o LSS_test LSS_main_LSS_ximu.f90 LSS_main_program.f90

#gfortran -o LSS_MCMC LSS_main_LSS_ximu.f90 LSS_main_ChisqLike.f90 $delm


#lib = lib/libde.a
EXE1 = AP_test
EXE2 = AP_MCMC
EXE3 = AP_grid
EXE4 = AP_MCMC_wbinned
EXE5 = AP_Rhct
EXE6 = AP_grid_wcdm
EXE7 = AP_coupled_de

OBJS = AP_tools.o

F90C = gfortran  
F90FLAGS = 

default: $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7)

AP_test.o : $(OBJS) 
AP_MCMC.o : $(OBJS) 
AP_grid.o : $(OBJS) 
AP_MCMC_wbinned.o : $(OBJS) 
AP_Rhct.o : $(OBJS) 
AP_coupled_de.o : $(OBJS) 
AP_grid_wcdm.o : $(OBJS)
#AP_test: $(OBJS) AP_test.o
#AP_MCMC: $(OBJS) AP_MCMC.o


$(EXE1): $(OBJS) AP_test.o
	$(F90C)  -o $(EXE1) $(OBJS) AP_test.o 
$(EXE2): $(OBJS) AP_MCMC.o
	$(F90C)  -o $(EXE2) $(OBJS) AP_MCMC.o -I/home/xiaodongli/software/delib/mods -lde
$(EXE3): $(OBJS) AP_grid.o
	$(F90C)  -o $(EXE3) $(OBJS) AP_grid.o -I/home/xiaodongli/software/delib/mods -lde
$(EXE4): $(OBJS) AP_MCMC_wbinned.o
	$(F90C)  -o $(EXE4) $(OBJS) AP_MCMC_wbinned.o -I/home/xiaodongli/software/delib/mods -lde
$(EXE5): $(OBJS) AP_Rhct.o
	$(F90C)  -o $(EXE5) $(OBJS) AP_Rhct.o -I/home/xiaodongli/software/delib/mods -lde
$(EXE6): $(OBJS) AP_grid_wcdm.o
	$(F90C)  -o $(EXE6) $(OBJS) AP_grid_wcdm.o -I/home/xiaodongli/software/delib/mods -lde
$(EXE7): $(OBJS) AP_coupled_de.o
	$(F90C)  -o $(EXE7) $(OBJS) AP_coupled_de.o -I/home/xiaodongli/software/delib/mods -lde



%.o: %.f90
	$(F90C) -c $*.f90 -lde -I/home/xiaodongli/software/delib/mods 
clean :
	rm -rf *.o *.mod *.a *.so  $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7)
	
