# QFFT Project — Quaternion-ish Fourier Transform in C

This project implements a **per-channel 2D Fourier transform** for RGB images, packaged both as a command-line tool and a Mac OS drag-and-drop `.app` for Apple Silicon. While the project is inspired by quaternion Fourier transforms, it currently **does not implement full quaternion arithmetic**; each color channel is transformed independently.

---

## Project Contents

- `src/qfft.c` — C implementation:
  - Naive 2D DFT per channel.
  - Optional FFTW-based 2D FFT (`make WITH_FFTW=1`).
- `Makefile` — build instructions.
- `test/generate_test_ppm.py` — generates 128×128 `.ppm` test images.
- `QFFT.app/` — drag-and-drop app bundle for Apple Silicon(zip file).
- `README.md` — this documentation.
- `LICENSE` — MIT license.

---

## Build & Run (Command-Line)

### Prerequisites
- C compiler (`gcc` or `clang`).
- Optional: FFTW3 development libraries.

### Build naive DFT
```bash
make
