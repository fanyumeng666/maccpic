# ccp1d -- 1D electrostatic PIC/MCC (Fortran 90 + MPI)
#   make            optimised executable  ./ccp1d      (objects in obj/)
#   make debug      checked executable    ./ccp1d_dbg  (objects in obj_dbg/)
#   make check      build and run the unit tests in tests/ with the debug flags
#   make clean
FC      = mpif90
OPT     = -O3 -march=native -fimplicit-none -ffree-line-length-none
DBG     = -g -O0 -ffree-line-length-none -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -Wextra \
          -finit-real=snan -finit-integer=-99999999 -fimplicit-none -Wno-unused-dummy-argument -Wno-integer-division
MPIRUN  = mpirun --allow-run-as-root --oversubscribe

SRC  = input.f90 global.f90 random.f90 poisson.f90 deposit.f90 push.f90 version.f90 diag.f90 init.f90 main.f90
OBJ  = $(patsubst %.f90,obj/%.o,$(SRC))
OBJD = $(patsubst %.f90,obj_dbg/%.o,$(SRC))
# objects the tests link against (everything except main)
LIBD = $(filter-out obj_dbg/main.o,$(OBJD))

TESTS = tests/poisson/test_poisson tests/deposit/test_deposit tests/random/test_random tests/langmuir/langmuir
RNGTESTS = tests/random/test_tails tests/random/test_seeds tests/random/test_corr

.PHONY: all debug check check-rng clean version

all: ccp1d
debug: ccp1d_dbg

ccp1d: $(OBJ)
	$(FC) $(OPT) -o $@ $^
ccp1d_dbg: $(OBJD)
	$(FC) $(DBG) -o $@ $^

obj/%.o: %.f90 | obj
	$(FC) $(OPT) -J obj -c $< -o $@
obj_dbg/%.o: %.f90 | obj_dbg
	$(FC) $(DBG) -J obj_dbg -c $< -o $@
obj obj_dbg:
	mkdir -p $@

# version.f90 carries the git commit into info.txt; regenerated only when it changes
version:
	@c=$$(git rev-parse --short HEAD 2>/dev/null || echo unknown); \
	 git diff --quiet 2>/dev/null || c="$$c-dirty"; \
	 printf 'module version\n  implicit none\n  character(len=*), parameter :: git_commit = "%s"\nend module version\n' "$$c" > version.f90.tmp; \
	 if cmp -s version.f90.tmp version.f90; then rm -f version.f90.tmp; else mv version.f90.tmp version.f90; fi
version.f90: version

# module dependencies (both object directories)
define DEPS
$(1)/global.o:  $(1)/input.o
$(1)/deposit.o: $(1)/global.o
$(1)/push.o:    $(1)/global.o
$(1)/diag.o:    $(1)/global.o $(1)/random.o $(1)/poisson.o $(1)/deposit.o $(1)/push.o $(1)/version.o
$(1)/init.o:    $(1)/global.o $(1)/random.o $(1)/poisson.o $(1)/deposit.o $(1)/diag.o
$(1)/main.o:    $(1)/global.o $(1)/random.o $(1)/poisson.o $(1)/deposit.o $(1)/push.o $(1)/diag.o $(1)/init.o
endef
$(eval $(call DEPS,obj))
$(eval $(call DEPS,obj_dbg))

# ---- tests ----
tests/%: tests/%.f90 $(LIBD)
	$(FC) $(DBG) -I obj_dbg -o $@ $< $(LIBD)

check: ccp1d_dbg $(TESTS)
	@fail=0; \
	for t in tests/poisson/test_poisson tests/deposit/test_deposit tests/random/test_random; do \
	  echo "==== $$t"; $(MPIRUN) -np 1 ./$$t || fail=1; \
	done; \
	echo "==== tests/deposit/test_deposit (2 ranks)"; $(MPIRUN) -np 2 ./tests/deposit/test_deposit || fail=1; \
	echo "==== tests/langmuir/langmuir"; $(MPIRUN) -np 1 ./tests/langmuir/langmuir || fail=1; \
	echo "==== tests/restart/run.sh"; sh tests/restart/run.sh || fail=1; \
	if [ $$fail -eq 0 ]; then echo "ALL TESTS PASSED"; else echo "SOME TESTS FAILED"; exit 1; fi

# slower statistical tests of the generator in the loading pattern (~1 min, optimised build)
tests/random/test_tails tests/random/test_seeds tests/random/test_corr: tests/random/%: tests/random/%.f90 obj/random.o
	$(FC) $(OPT) -I obj -o $@ $< obj/random.o
check-rng: $(RNGTESTS)
	@for t in $(RNGTESTS); do echo "==== $$t"; ./$$t; done

clean:
	rm -rf obj obj_dbg ccp1d ccp1d_dbg version.f90 $(TESTS) $(RNGTESTS)
