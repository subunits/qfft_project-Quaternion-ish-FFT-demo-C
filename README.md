# QFFT — Quaternion Fourier Transform in C

A pure, zero-dependency C implementation of the 2D Quaternion Discrete Fourier Transform (QDFT) and Inverse Transform (IQDFT), tailored for full 4D spatial color processing and packaged for Apple Silicon macOS systems.

---

## Overview

Unlike standard Fourier transforms that process color channels independently, QFFT treats RGB pixel data as **pure quaternions in 4D space (ℍ)**:

```
q(m,n) = 0 + r·i + g·j + b·k
```

Using a unit pure quaternion axis **μ = (i + j + k) / √3** (where μ² = −1), the forward and inverse transforms apply full 4D hypercomplex arithmetic to preserve inter-channel phase relationships across spatial domain transformations.

---

## Features

- **Full 4D Quaternion Engine** — Implements hypercomplex quaternion multiplication, forward transform, and exact inverse reconstruction.
- **Zero External Dependencies** — Built entirely with standard C libraries (`math.h`, `stdlib.h`, `stdio.h`).
- **Apple Silicon Optimized** — Native ARM compilation flag configuration with low-overhead memory allocation routines.
- **App Bundle Packaging** — Automated target directory structure setup for native macOS `.app` distribution.

---

## Build & Run

### Prerequisites

- C Compiler (`gcc` or `clang`)
- `make` utility

### Compilation

Build the executable and package the macOS app bundle:

```bash
make
```

To clean build artifacts:

```bash
make clean
```

### Execution

Run the binary directly from the terminal:

```bash
./qfft
```

Or execute from the packaged bundle:

```bash
./QFFT.app/Contents/MacOS/qfft
```

---

## Repository Structure

| Path | Description |
|------|-------------|
| `qfft.c` | Core C source code containing quaternion matrix routines, `forward_qfft`, `inverse_qfft`, and test runner. |
| `Makefile` | Build automation rules and `.app` bundle structure assembler. |
| `QFFT.app/` | Standalone macOS application bundle directory. |
| `README.md` | Project documentation. |
| `LICENSE` | MIT License. |

---

## License

[MIT License](LICENSE)
