FC=mpif90
###FC=ifort
###FCFLAGS= -openmp -shared-intel -mcmodel=medium -fast
###FCFLAGS= -openmp -mcmodel=large -i-dynamic -ipo -parallel
###FCFLAGS=-fast
FCFLAGS=-O3 
OBJS=input.o global.o main.o init.o output.o creator.o cont.o smear.o setrho.o setv.o smooth.o fields.o move.o curnt.o dran2.o

.SUFFIXES: .f90 .o

em2: $(OBJS)
	$(FC) -o $@ $(FCFLAGS) $(OBJS)

.f90.o:
	$(FC) -o $@ $(FCFLAGS) -c $<

input.o   : input.f90
global.o  : global.f90
main.o    : main.f90
init.o    : init.f90
output.o  : output.f90
creator.o : creator.f90
cont.o    : cont.f90
smear.o   : smear.f90
#waveint.o : waveint.f90
setrho.o  : setrho.f90
setv.o    : setv.f90
smooth.o  : smooth.f90
fields.o  : fields.f90
move.o    : move.f90
curnt.o   : curnt.f90
dran2.o   : dran2.f90
#fft.o     : fft.f90

.PHONY    : clean
clean     :
	rm -f em2 *.mod $(OBJS)