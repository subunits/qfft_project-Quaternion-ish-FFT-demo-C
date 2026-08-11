# QFFT — Quaternion Fourier Transform in C

A pure, zero-dependency C implementation of the 2D Quaternion Discrete Fourier Transform (QDFT) and its exact Inverse (IQDFT), treating RGB pixel data as pure quaternions in 4D hypercomplex space.

---

## Overview

Unlike standard Fourier transforms that process color channels independently, QFFT treats RGB pixel data as **pure quaternions in 4D space (ℍ)**:

```
q(m,n) = 0 + r·i + g·j + b·k
```

Using a configurable unit pure quaternion axis **μ** (defaulting to `(i + j + k) / √3`, where μ² = −1), the forward and inverse transforms apply a **two-sided exponential kernel**:

```
F(u,v) = Σ_m Σ_n  e^{μ·α} · f(m,n) · e^{μ·β}
```

This two-sided structure is essential: because quaternion multiplication is non-commutative, applying kernels on both sides couples the RGB channels in a way a scalar DFT — or a single-sided kernel — cannot. The axis μ is a runtime parameter, allowing different transform behaviours without recompilation.

---

## Features

- **True Two-Sided QFFT** — Left and right exponential kernels are applied separately, preserving the non-commutative channel coupling that defines a correct quaternion transform.
- **Configurable μ axis** — `PureUnitQuat` is passed into both `forward_qfft` and `inverse_qfft`; swap in any unit pure quaternion at the call site.
- **Exact Inverse Reconstruction** — Round-trip RGB error is at floating-point machine epsilon (~10⁻¹³); scalar `w` drift is verified to remain near zero.
- **Robust Memory Handling** — All `malloc` calls are checked; failures print a diagnostic and exit cleanly.
- **Output Clamping** — Reconstructed pixel channels are clamped to `[0, 255]` before any downstream use.
- **Zero External Dependencies** — Built entirely with standard C libraries (`math.h`, `stdlib.h`, `stdio.h`).
- **Cross-Platform Build** — Compiles with `gcc` or `clang` on Linux, macOS (including Apple Silicon via native toolchain), and other POSIX systems. Optional FFTW linkage available for future FFT acceleration.
- **App Bundle Packaging** — Automated macOS `.app` directory structure via `make bundle`.

---

## Build & Run

### Prerequisites

- C compiler (`gcc` or `clang`)
- `make`

### Compilation

```bash
make
```

To build with FFTW support (for future acceleration):

```bash
make WITH_FFTW=1
```

To clean build artifacts:

```bash
make clean
```

### Execution

```bash
./qfft
```

Or from the app bundle:

```bash
./QFFT.app/Contents/MacOS/qfft
```

### Generate a test image

```bash
mkdir -p test
python3 generate_test_ppm.py   # writes test/test.ppm (128×128 RGB gradient)
```

---

## Repository Structure

| Path | Description |
|------|-------------|
| `qfft.c` | Core implementation: quaternion arithmetic, two-sided `forward_qfft`, `inverse_qfft`, pixel clamping, reconstruction validation, and test runner. |
| `Makefile` | Build rules and macOS `.app` bundle assembler. Optional `WITH_FFTW=1` flag. |
| `generate_test_ppm.py` | Python script that generates a 128×128 P6 PPM test image with RGB gradients and a radial blue pattern. |
| `QFFT.app/` | Standalone macOS application bundle directory. |
| `README.md` | Project documentation. |
| `LICENSE` | MIT License. |

---

## Implementation Notes

This is an **O(M²·N²) reference implementation** — correct and exact, but not intended for large images. For production use, the inner loops should be replaced with a row-column decomposition and butterfly FFT, reducing complexity to O(MN log MN).

The `w` (scalar) component of every input pixel is set to zero, consistent with the pure-quaternion image model. After an inverse transform, `w` should return to near-zero; the validation routine checks this explicitly.

---

## License

[MIT License](LICENSE)
