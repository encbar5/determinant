IDIR=src
ODIR=obj
BIN=det.out

# mpifort -std=f95 gemv.f -lscalapack -lblas
# module load OpenMPI
# module load ScaLAPACK
# module load PETSc

CC=mpifort
CFLAGS=-I$(IDIR)
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
LIBS=-lscalapack -lblas

SRC = $(wildcard $(IDIR)/*.f)
OBJ = $(ODIR)/gemv.o #$(patsubst $(IDIR)/%.f,$(ODIR)/%.o,$(SRC))


$(ODIR)/%.o: $(IDIR)/%.f
	$(CC) -o $@ $< $(CFLAGS)

$(BIN): $(OBJ)
	$(CC) -o $@ $^ $(LFLAGS) $(LIBS)

.PHONY: clean
.PHONY: test
test:
	echo $(DEPS); echo $(OBJ); echo $(SRC)

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ $(BIN) *.mod
