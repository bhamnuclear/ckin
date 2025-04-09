CC = gcc
WCLBES = NO
CERNLIB = NO
CFLAGS = -Wall -lm -O2 -pedantic -DCOMPILE_DIR=\"$(PWD)\"


all: check

clean:
	rm -f ckin check.c check.out
check:
	@echo "int main() { return 0; }" > check.c
	@$(CC) check.c -lmathlib -lkernlib -o check.out > /dev/null 2>&1; \
	if [ $$? -eq 0 ]; then \
		echo "Compiling with cernlib"; \
		$(MAKE) ckin_cernlib; \
	else \
		$(CC) check.c ./libwclbes.a -o check.out > /dev/null 2>&1; \
		if [ $$? -eq 0 ]; then \
			echo "Compiling with wclbes"; \
			$(MAKE) ckin_wclbes; \
		else \
			echo "Compiling"; \
			$(MAKE) ckin; \
		fi \
	fi
	@rm -f check.c check.out



ckin :
	$(CC) -o ckin ckin.c $(CFLAGS)

ckin_wclbes :
	$(CC) -o ckin ckin.c $(CFLAGS) -DHAVE_WCLBES libwclbes.a -lgfortran

ckin_cernlib :
	$(CC) -o ckin ckin.c $(CFLAGS) -DHAVE_WCLBES -lmathlib -lkernlib -lgfortran
