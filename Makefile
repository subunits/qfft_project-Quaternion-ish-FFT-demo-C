# Makefile for qfft project
CC = gcc
CFLAGS = -O2 -std=c99 -Wall
LDFLAGS =
SRC = src/qfft.c
OUT = qfft

# To build with FFTW, run: make WITH_FFTW=1
ifdef WITH_FFTW
    CFLAGS += -DWITH_FFTW
    LDFLAGS += -lfftw3 -lm
endif

all: $(OUT)

$(OUT): $(SRC)
	$(CC) $(CFLAGS) -o $(OUT) $(SRC) $(LDFLAGS)

clean:
	rm -f $(OUT) *.o
