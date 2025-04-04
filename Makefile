#IDIR =../include
TARGET = ckin
CC = gcc
CFLAGS = -Wall -O2 -DCOMPILE_DIR=\"$(PWD)\" -DHAVE_WCLBES
LIBS = -lm ./libwclbes.a -lgfortran

.PHONY: default all clean

default: $(TARGET)
all: default

ckin:
	$(CC) *.c -o $@ $(CFLAGS) $(LIBS)

clean:
	\rm -f $(TARGET)
