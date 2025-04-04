AR = ar r
CC = gfortran
RM = rm -f
CFLAGS = -fPIC -D CERNLIB_DOUBLE

ARCHIVE = ../libwclbes.a
OBJS    = wclbes.o clogam64.o r8dp.o r7dp.o r6dp.o r5dp.o r4dp.o r3dp.o r2dp.o r1dp.o cdigam64.o mtlset.o mtlprt.o lenocc.o abend.o

.PHONY: default all clean

%.o: %.F
	$(CC) -c -o $@ $^ $(CFLAGS)

all: $(ARCHIVE)

$(ARCHIVE): $(OBJS)
	$(RM) $@
	$(AR) $@ $(OBJS)
clean: 
	$(RM) *.o core *~

very-clean:
	$(RM) *.o core *~ $(ARCHIVE)
