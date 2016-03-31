NVCC      = nvcc
CFLAGS    = -O3
LDFLAGS   =

OBJ     = time_it.o
HDRS    = time_it.h

main: Assignment3.o $(OBJ)
	$(NVCC) $(CFLAGS) Assignment3.o $(OBJ) $(LDFLAGS) -o Assignment3

print: Assgn3.o $(OBJ)
	$(NVCC) $(CFLAGS) -arch sm_20 Assgn3.o $(OBJ) $(LDFLAGS) -o Assgn3
a: Assgn3.o $(OBJ)
	$(NVCC) $(CFLAGS) Assgn3.o $(OBJ) $(LDFLAGS) -o Assgn3
clean:
	rm Assignment3 Assignment3.o Assgn3 Assgn3.o


.SUFFIXES: .c .cu .o
.c.o: $(HDRS)
	$(NVCC) -c $(CFLAGS) $<
.cu.o: $(HDRS)
	$(NVCC) -c $(CFLAGS) $<

