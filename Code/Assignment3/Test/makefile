NVCC      = nvcc
CFLAGS    = -O3 
LDFLAGS   =

OBJ     = time_it.o
HDRS    = time_it.h

perc: perc.o $(OBJ)
	$(NVCC) $(CFLAGS) perc.o $(OBJ) $(LDFLAGS) -o perc

examples: examples.o
	$(NVCC) $(CFLAGS) examples.o $(LDFLAGS) -o examples

timing: timing.o $(OBJ)
	$(NVCC) $(CFLAGS) timing.o $(OBJ) $(LDFLAGS) -o timing


.SUFFIXES: .c .cu .o
.c.o: $(HDRS)
	$(NVCC) -c $(CFLAGS) $<
.cu.o: $(HDRS)
	$(NVCC) -c $(CFLAGS) $<
