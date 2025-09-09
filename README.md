# qfft_project — Quaternion-ish FFT demo (C)

This repository contains a small C demonstration that packs RGB into quaternion-like components (scalar part = 0) and computes **per-channel 2D transforms**. It's intentionally simple so you can extend it into a proper Quaternion Fourier Transform later.

## Project contents
- `src/qfft.c` — C implementation:
  - Naive 2D DFT implementation (always available).
  - Optional FFTW3-enabled 2D FFT (compile with `WITH_FFTW=1` in the Makefile if FFTW is installed).
- `Makefile` — builds the project. Use `make` or `make WITH_FFTW=1`.
- `test/generate_test_ppm.py` — Python script to produce a sample `test.ppm` (100×100 gradient).
- `LICENSE` — MIT license.

## Roadmap (complete)
This roadmap shows practical next steps from prototype → production.

**Phase 1 — Prototype (done here)**
1. Pack R,G,B into imaginary quaternion components and compute per-channel DFT (naive).  
2. Export magnitude images per channel for inspection.

**Phase 2 — Performance & correctness**
1. Add FFT-based implementation (FFTW) to accelerate transforms (implemented conditionally in this repo).  
2. Add separable real-to-complex optimizations and zero-centering (shift DC to center of image).  
3. Add inverse transform and round-trip validation tools.

**Phase 3 — Quaternion Fourier Transform (QFT)**
1. Choose a QFT definition (several exist; decide axis/order and handle noncommutativity). Common choices:
   - Side-by-side complexification: treat quaternion as pair of complex numbers and apply 2D complex FFTs.
   - Full quaternion exponential-based QFT (requires careful ordering of exponentials).
2. Implement quaternion arithmetic (or include a small quaternion helper library).
3. Implement forward and inverse QFT, test with synthetic vector fields and polarization examples.

**Phase 4 — Applications & UX**
1. Polarization-aware Fourier optics demos (store Jones vectors or Stokes parameters).
2. Quaternion holography prototypes (encode multi-channel holograms).
3. Provide Python bindings / visualization scripts (matplotlib) for easy inspection.

**Phase 5 — Packaging & Distribution**
1. Add unit tests & CI (GitHub Actions).
2. Create Dockerfile for reproducible builds (include FFTW).
3. Publish examples and documentation.

## Build & Run
### Prerequisites
- A C compiler (gcc/clang).
- For FFTW support: install FFTW3 development libraries (platform-specific).

### Build (naive DFT)
```
make
```

### Build (FFTW, if available)
```
make WITH_FFTW=1
```

### Run
Generate a test image:
```
python3 test/generate_test_ppm.py test/test.ppm
```
Run the program:
```
./qfft test/test.ppm out
```
Outputs:
- `out_R.pgm`, `out_G.pgm`, `out_B.pgm` — frequency magnitude images (PGM).

## Notes & Tips
- The naive DFT is O(N^4) and will be slow for images larger than ~128×128. Use FFTW for larger images.
- The output magnitude is log-scaled for visualization. To view, convert to PNG with ImageMagick:
  ```
  convert out_R.pgm out_R.png
  ```
- This project is a starting point. If you'd like, I can:
  - Implement a chosen QFT definition in C.
  - Add tests, inverse transform, and visualization scripts.
  - Provide a prebuilt binary (Linux x86_64) — but tell me your target platform.

