FC = mpif90
FFLAGS=

SRCS	= main.f90
OBJS    = main.o
EXEC    = exe
MODS	= amsta01probleme.mod amsta01sparse.mod amsta01maillage.mod amsta01solveur.mod
MODOBJS = amsta01probleme.o amsta01sparse.o amsta01maillage.o amsta01solveur.o

CURRENTMOD = amsta01maillage

all: $(EXEC)

amsta01probleme.o: amsta01maillage.mod amsta01sparse.mod
amsta01solveur.o: amsta01probleme.o amsta01maillage.mod amsta01sparse.mod

main.o: $(MODS)

%.o %.mod: %.f90
	$(FC) $(FFLAGS) -c $<
	touch $*.o $*.mod

main.o: main.f90
	$(FC) $(FFLAGS) -c $<

exe: $(MODOBJS) main.o
	$(FC) $(LDFLAGS) -o $@ $^

run: exe
	./$<

clean :
	rm -f $(OBJS) $(MODOBJS) $(EXEC) $(MODS)
