NVCC      = nvcc
CFLAGS    = -O3
LDFLAGS   =

OBJ     = time_it.o
HDRS    = time_it.h

q1: Assgn3_q1.o $(OBJ)
	$(NVCC) $(CFLAGS) Assgn3_q1.o $(OBJ) $(LDFLAGS) -o Assgn3_q1
q3: Assgn3_q3.o $(OBJ)
	$(NVCC) $(CFLAGS) Assgn3_q3.o $(OBJ) $(LDFLAGS) -o Assgn3_q3

clean:
	rm Assignment3 Assignment3.o Assgn3 Assgn3.o Assgn3_q1 Assgn3_q1.o


.SUFFIXES: .c .cu .o
.c.o: $(HDRS)
	$(NVCC) -c $(CFLAGS) $<
.cu.o: $(HDRS)
	$(NVCC) -c $(CFLAGS) $<

