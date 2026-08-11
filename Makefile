# Makefile for qfft project (all files in root)
CC = gcc
CFLAGS = -O2 -std=c99 -D_GNU_SOURCE -Wall
LDFLAGS = -lm
SRC = qfft.c
OUT = qfft

# Target bundle directory structure
APP_NAME = QFFT.app
APP_BIN_DIR = $(APP_NAME)/Contents/MacOS

# To build with FFTW, run: make WITH_FFTW=1
ifdef WITH_FFTW
    CFLAGS += -DWITH_FFTW
    LDFLAGS += -lfftw3
endif

all: $(OUT) bundle

$(OUT): $(SRC)
	$(CC) $(CFLAGS) -o $(OUT) $(SRC) $(LDFLAGS)

bundle: $(OUT)
	@mkdir -p $(APP_BIN_DIR)
	cp $(OUT) $(APP_BIN_DIR)/
	@echo "Packaged $(OUT) into $(APP_BIN_DIR)/"

clean:
	rm -f $(OUT) *.o
	rm -f $(APP_BIN_DIR)/$(OUT)
