#Macro definition
FC = ifort
FFLAGS =  -O3 -qopenmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
LDFLAGS =  $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
           -Wl,--start-group \
           $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
           $(MKLROOT)/lib/intel64/libmkl_core.a \
           $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
           -Wl,--end-group -lpthread -lm
OBJDPL = nn.o dip.o maind.o orient.o
OBJSOC = nn.o soc.o mains.o orient.o
#end of Macro definition
testd.x: $(OBJDPL)
	$(FC) $(FFLAGS) $(OBJDPL) -o testd.x $(LDFLAGS)

tests.x: $(OBJSOC)
	$(FC) $(FFLAGS) $(OBJSOC) -o tests.x $(LDFLAGS)

clean:
	rm -f *.o *.mod a.out *.x *.exe

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
