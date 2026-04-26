CC = gcc
CFLAGS = -O3 -march=native -Wall
LDFLAGS = -lgmp -lm

all: qs_scale

qs_scale: qs_scale.c
	$(CC) $(CFLAGS) qs_scale.c -o qs_scale $(LDFLAGS)

.PHONY: clean run-tests

clean:
	rm -f qs_scale

run-tests: qs_scale
	@echo "=== QS-EDEN scale tests ==="
	./qs_scale 50 800 30000
	./qs_scale 60 2000 100000
	./qs_scale 70 5000 300000
	./qs_scale 80 12000 800000

# install dependency on Debian/Ubuntu:
#   sudo apt install libgmp-dev
