#Makefile

FF     = gfortran
FFLAGs = -O3

all: repmoms.o 
	$(FF) $(FFLAGS) repmoms.o -o repmoms $(shell lhapdf-config --ldflags)

%.o: %.f
	            $(FF) $(FFLAGS) -c $<

clean: 
	rm -rf *.o repmoms *~
