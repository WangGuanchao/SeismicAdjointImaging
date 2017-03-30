#-----------------------------------------------------------------------------
# Makefile
#-----------------------------------------------------------------------------

include ../make.inc

#-----------------------------------------------------------------------------
# Programs
#-----------------------------------------------------------------------------

PROGS = $(wildcard *.c)
OBJES = $(patsubst %.c, %.o, $(PROGS))
FUNCS = $(patsubst %.c, %  , $(PROGS))

all : $(OBJES) $(FUNCS)

$(FUNCS) : $(patsubst %, %.o, $@) $(include)
	$(MPICC) $(CFLAGS) $(patsubst %, %.o, $@) $(libary) -o ../bin/$@ -I$(include)

%.o : %.c
	$(MPICC) -c $(CFLAGS) $< -o $@ -I$(include)

clean:
	rm *.a -f
	rm *.o -f
