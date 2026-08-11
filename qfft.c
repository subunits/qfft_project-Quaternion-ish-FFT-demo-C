#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------------------------------------------------------------------
 * Quaternion structure
 *   w  – scalar (real) part
 *   x  – i component  (Red channel)
 *   y  – j component  (Green channel)
 *   z  – k component  (Blue channel)
 * -------------------------------------------------------------------------*/
typedef struct {
    double w;
    double x;
    double y;
    double z;
} Quaternion;

/* Pure unit quaternion used as the transform axis μ.
 * The caller is responsible for ensuring |μ| = 1.
 * Both the left and right exponential kernels share the same μ,
 * which is the standard single-axis 2-D QFFT convention. */
typedef struct {
    double x;
    double y;
    double z;
} PureUnitQuat;

/* Default axis: μ = (i + j + k) / √3 */
static PureUnitQuat default_mu(void) {
    double inv = 1.0 / sqrt(3.0);
    PureUnitQuat mu = { inv, inv, inv };
    return mu;
}

/* ---------------------------------------------------------------------------
 * Memory helpers
 * -------------------------------------------------------------------------*/
Quaternion **allocate_qmatrix(int width, int height) {
    Quaternion **mat = (Quaternion **)malloc(height * sizeof(Quaternion *));
    if (!mat) {
        fprintf(stderr, "allocate_qmatrix: outer malloc failed\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < height; i++) {
        mat[i] = (Quaternion *)malloc(width * sizeof(Quaternion));
        if (!mat[i]) {
            fprintf(stderr, "allocate_qmatrix: inner malloc failed at row %d\n", i);
            exit(EXIT_FAILURE);
        }
    }
    return mat;
}

void free_qmatrix(Quaternion **mat, int height) {
    for (int i = 0; i < height; i++) free(mat[i]);
    free(mat);
}

/* ---------------------------------------------------------------------------
 * Quaternion arithmetic
 * -------------------------------------------------------------------------*/

/* Hamilton product: q1 * q2 */
static Quaternion q_mul(Quaternion q1, Quaternion q2) {
    Quaternion r;
    r.w = q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z;
    r.x = q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y;
    r.y = q1.w*q2.y - q1.x*q2.z + q1.y*q2.w + q1.z*q2.x;
    r.z = q1.w*q2.z + q1.x*q2.y - q1.y*q2.x + q1.z*q2.w;
    return r;
}

/* Quaternion addition */
static Quaternion q_add(Quaternion a, Quaternion b) {
    Quaternion r = { a.w+b.w, a.x+b.x, a.y+b.y, a.z+b.z };
    return r;
}

/* Scale all components */
static Quaternion q_scale(Quaternion q, double s) {
    Quaternion r = { q.w*s, q.x*s, q.y*s, q.z*s };
    return r;
}

/* Build e^{μ·θ} = cos(θ) + μ·sin(θ) */
static Quaternion q_exp_mu(PureUnitQuat mu, double theta) {
    double c = cos(theta);
    double s = sin(theta);
    Quaternion r = { c, mu.x*s, mu.y*s, mu.z*s };
    return r;
}

/* Clamp a double to [lo, hi] */
static double clamp(double v, double lo, double hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

/* ---------------------------------------------------------------------------
 * 2-D Quaternion Fourier Transform (two-sided, single-axis μ)
 *
 *   F(u,v) = Σ_m Σ_n  e^{μ·α_{um}} · f(m,n) · e^{μ·β_{vn}}
 *
 *   where  α_{um} = -2π·u·m / M
 *          β_{vn} = -2π·v·n / N
 *
 * Splitting the exponent across rows (left) and columns (right) is the
 * standard approach for a separable two-sided QFFT.  It preserves the
 * non-commutative structure that couples the RGB channels in a way a
 * scalar DFT cannot.
 *
 * NOTE: O(M²·N²) — intentionally a reference/verification implementation.
 *       Do not use on large images without replacing the inner loops with
 *       a row-column decomposition + butterfly FFT.
 * -------------------------------------------------------------------------*/
void forward_qfft(Quaternion **in, Quaternion **out,
                  int width, int height, PureUnitQuat mu) {
    for (int u = 0; u < height; u++) {
        for (int v = 0; v < width; v++) {
            Quaternion sum = {0.0, 0.0, 0.0, 0.0};

            for (int m = 0; m < height; m++) {
                /* Left kernel: e^{μ · (-2π·u·m/M)} */
                double left_angle = -2.0 * M_PI * (double)u * m / height;
                Quaternion left_k = q_exp_mu(mu, left_angle);

                for (int n = 0; n < width; n++) {
                    /* Right kernel: e^{μ · (-2π·v·n/N)} */
                    double right_angle = -2.0 * M_PI * (double)v * n / width;
                    Quaternion right_k = q_exp_mu(mu, right_angle);

                    /* Two-sided product: left · f(m,n) · right */
                    Quaternion prod = q_mul(q_mul(left_k, in[m][n]), right_k);
                    sum = q_add(sum, prod);
                }
            }
            out[u][v] = sum;
        }
    }
}

/* ---------------------------------------------------------------------------
 * Inverse 2-D QFFT
 *
 *   f(m,n) = (1/MN) · Σ_u Σ_v  e^{μ·α'_{um}} · F(u,v) · e^{μ·β'_{vn}}
 *
 *   where  α'_{um} = +2π·u·m / M  (positive angle)
 *          β'_{vn} = +2π·v·n / N
 * -------------------------------------------------------------------------*/
void inverse_qfft(Quaternion **freq, Quaternion **spatial,
                  int width, int height, PureUnitQuat mu) {
    double scale = 1.0 / ((double)width * height);

    for (int m = 0; m < height; m++) {
        for (int n = 0; n < width; n++) {
            Quaternion sum = {0.0, 0.0, 0.0, 0.0};

            for (int u = 0; u < height; u++) {
                double left_angle = 2.0 * M_PI * (double)u * m / height;
                Quaternion left_k = q_exp_mu(mu, left_angle);

                for (int v = 0; v < width; v++) {
                    double right_angle = 2.0 * M_PI * (double)v * n / width;
                    Quaternion right_k = q_exp_mu(mu, right_angle);

                    Quaternion prod = q_mul(q_mul(left_k, freq[u][v]), right_k);
                    sum = q_add(sum, prod);
                }
            }
            spatial[m][n] = q_scale(sum, scale);
        }
    }
}

/* ---------------------------------------------------------------------------
 * Clamp all pixel channels of a reconstructed matrix to [0, 255].
 * Floating-point accumulation can push values slightly out of range.
 * -------------------------------------------------------------------------*/
void clamp_pixels(Quaternion **mat, int width, int height) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            /* w should remain near 0; clamp it too for safety */
            mat[i][j].w = clamp(mat[i][j].w, 0.0, 255.0);
            mat[i][j].x = clamp(mat[i][j].x, 0.0, 255.0);
            mat[i][j].y = clamp(mat[i][j].y, 0.0, 255.0);
            mat[i][j].z = clamp(mat[i][j].z, 0.0, 255.0);
        }
    }
}

/* ---------------------------------------------------------------------------
 * Validation helpers
 * -------------------------------------------------------------------------*/

/* Compute max absolute error across all pixels and all RGB channels.
 * Also checks that the scalar (w) channel returned to near-zero. */
static void print_reconstruction_error(Quaternion **original,
                                       Quaternion **reconstructed,
                                       int width, int height) {
    double max_rgb_err = 0.0;
    double max_w_err   = 0.0;

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            max_rgb_err = fmax(max_rgb_err,
                fabs(original[i][j].x - reconstructed[i][j].x));
            max_rgb_err = fmax(max_rgb_err,
                fabs(original[i][j].y - reconstructed[i][j].y));
            max_rgb_err = fmax(max_rgb_err,
                fabs(original[i][j].z - reconstructed[i][j].z));
            /* w was 0 in all input pixels; verify it stayed near 0 */
            max_w_err = fmax(max_w_err, fabs(reconstructed[i][j].w));
        }
    }
    printf("Max RGB reconstruction error : %e\n", max_rgb_err);
    printf("Max scalar (w) drift         : %e  (expect ~0)\n", max_w_err);
}

/* ---------------------------------------------------------------------------
 * main
 * -------------------------------------------------------------------------*/
int main(void) {
    const int WIDTH  = 8;
    const int HEIGHT = 8;

    /* Pixel value bounds used during initialisation.
     * With (i+j)*R_SCALE: max = (7+7)*15 = 210  ✓ within [0,255]
     * With  i*C_SCALE   : max =  7*20   = 140  ✓
     * With  j*C_SCALE   : max =  7*20   = 140  ✓  */
    const double R_SCALE = 15.0;
    const double C_SCALE = 20.0;

    PureUnitQuat mu = default_mu();   /* μ = (i+j+k)/√3 — change here if needed */

    Quaternion **input        = allocate_qmatrix(WIDTH, HEIGHT);
    Quaternion **freq         = allocate_qmatrix(WIDTH, HEIGHT);
    Quaternion **reconstructed = allocate_qmatrix(WIDTH, HEIGHT);

    /* Initialise spatial colour data; scalar part w = 0 throughout */
    for (int i = 0; i < HEIGHT; i++) {
        for (int j = 0; j < WIDTH; j++) {
            input[i][j].w = 0.0;
            input[i][j].x = (i + j) * R_SCALE;  /* Red   */
            input[i][j].y =  i      * C_SCALE;  /* Green */
            input[i][j].z =  j      * C_SCALE;  /* Blue  */
        }
    }

    printf("Executing Forward QFFT (two-sided, μ = (i+j+k)/√3)...\n");
    forward_qfft(input, freq, WIDTH, HEIGHT, mu);

    printf("Executing Inverse QFFT...\n");
    inverse_qfft(freq, reconstructed, WIDTH, HEIGHT, mu);

    /* Clamp pixel values to [0,255] before any downstream use */
    clamp_pixels(reconstructed, WIDTH, HEIGHT);

    /* --- Spot-check a single pixel --- */
    printf("\nSample pixel [2][2]:\n");
    printf("  Original      RGB : (%.2f, %.2f, %.2f)\n",
           input[2][2].x, input[2][2].y, input[2][2].z);
    printf("  Reconstructed RGB : (%.2f, %.2f, %.2f)\n",
           reconstructed[2][2].x, reconstructed[2][2].y, reconstructed[2][2].z);

    /* --- Full-matrix error report --- */
    printf("\nFull-matrix validation:\n");
    print_reconstruction_error(input, reconstructed, WIDTH, HEIGHT);

    free_qmatrix(input,         HEIGHT);
    free_qmatrix(freq,          HEIGHT);
    free_qmatrix(reconstructed, HEIGHT);

    return EXIT_SUCCESS;
}
