CC = gcc
WCLBES = NO
CERNLIB = NO
CFLAGS = -Wall -lm -O2 -pedantic -DCOMPILE_DIR=\"$(PWD)\"


all: check

clean:
	\rm -f ckin
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



ckin: ckin.c
	$(CC) $< -o $@ $(CFLAGS)

ckin_wclbes: ckin.c
	$(CC) $< -o $@ $(CFLAGS) -DHAVE_WCLBES libwclbes.a -lgfortran; \
        \rm -f ckin; \
        ln -s ckin_wclbes ckin

ckin_cernlib: ckin.c
	$(CC) $< -o $@ $(CFLAGS) -DHAVE_WCLBES -lmathlib -lkernlib -lgfortran; \
        \rm -f ckin; \
        ln -s ckin_cernlib ckin
