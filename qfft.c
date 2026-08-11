#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Quaternion structure holding full 4D frequency/spatial domain state
typedef struct {
    double w; // Real / scalar
    double x; // Pure i (Red channel)
    double y; // Pure j (Green channel)
    double z; // Pure k (Blue channel)
} Quaternion;

// Allocate 2D matrix of Quaternions
Quaternion** allocate_qmatrix(int width, int height) {
    Quaternion** mat = (Quaternion**)malloc(height * sizeof(Quaternion*));
    for (int i = 0; i < height; i++) {
        mat[i] = (Quaternion*)malloc(width * sizeof(Quaternion));
    }
    return mat;
}

// Free 2D matrix of Quaternions
void free_qmatrix(Quaternion** mat, int height) {
    for (int i = 0; i < height; i++) {
        free(mat[i]);
    }
    free(mat);
}

// Quaternion multiplication: q1 * q2
Quaternion q_mul(Quaternion q1, Quaternion q2) {
    Quaternion res;
    res.w = q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z;
    res.x = q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y;
    res.y = q1.w*q2.y - q1.x*q2.z + q1.y*q2.w + q1.z*q2.x;
    res.z = q1.w*q2.z + q1.x*q2.y - q1.y*q2.x + q1.z*q2.w;
    return res;
}

// Forward 2D QFFT (Spatial Domain -> Quaternion Frequency Domain)
void forward_qfft(Quaternion** in, Quaternion** out, int width, int height) {
    double mu_x = 1.0 / sqrt(3.0);
    double mu_y = 1.0 / sqrt(3.0);
    double mu_z = 1.0 / sqrt(3.0);

    for (int u = 0; u < height; u++) {
        for (int v = 0; v < width; v++) {
            Quaternion sum = {0.0, 0.0, 0.0, 0.0};

            for (int m = 0; m < height; m++) {
                for (int n = 0; n < width; n++) {
                    double angle = -2.0 * M_PI * (((double)u * m / height) + ((double)v * n / width));
                    double cos_a = cos(angle);
                    double sin_a = sin(angle);

                    // Exponential kernel: e^(mu * angle)
                    Quaternion kernel = {
                        cos_a,
                        mu_x * sin_a,
                        mu_y * sin_a,
                        mu_z * sin_a
                    };

                    Quaternion prod = q_mul(in[m][n], kernel);
                    sum.w += prod.w;
                    sum.x += prod.x;
                    sum.y += prod.y;
                    sum.z += prod.z;
                }
            }
            out[u][v] = sum;
        }
    }
}

// Inverse 2D QFFT (Quaternion Frequency Domain -> Spatial Image)
void inverse_qfft(Quaternion** freq, Quaternion** spatial, int width, int height) {
    double mu_x = 1.0 / sqrt(3.0);
    double mu_y = 1.0 / sqrt(3.0);
    double mu_z = 1.0 / sqrt(3.0);
    double scale = 1.0 / (width * height);

    for (int m = 0; m < height; m++) {
        for (int n = 0; n < width; n++) {
            Quaternion sum = {0.0, 0.0, 0.0, 0.0};

            for (int u = 0; u < height; u++) {
                for (int v = 0; v < width; v++) {
                    // Positive angle for inverse transform
                    double angle = 2.0 * M_PI * (((double)u * m / height) + ((double)v * n / width));
                    double cos_a = cos(angle);
                    double sin_a = sin(angle);

                    Quaternion inv_kernel = {
                        cos_a,
                        mu_x * sin_a,
                        mu_y * sin_a,
                        mu_z * sin_a
                    };

                    Quaternion prod = q_mul(freq[u][v], inv_kernel);
                    sum.w += prod.w;
                    sum.x += prod.x;
                    sum.y += prod.y;
                    sum.z += prod.z;
                }
            }

            // Apply 1/(M*N) scaling normalization
            spatial[m][n].w = sum.w * scale;
            spatial[m][n].x = sum.x * scale;
            spatial[m][n].y = sum.y * scale;
            spatial[m][n].z = sum.z * scale;
        }
    }
}

int main(int argc, char** argv) {
    int width = 8;
    int height = 8;

    Quaternion** input = allocate_qmatrix(width, height);
    Quaternion** freq = allocate_qmatrix(width, height);
    Quaternion** reconstructed = allocate_qmatrix(width, height);

    // Initialize spatial color data (R=x, G=y, B=z)
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            input[i][j].w = 0.0;
            input[i][j].x = (i + j) * 15.0; // Red
            input[i][j].y = (i * 20.0);      // Green
            input[i][j].z = (j * 20.0);      // Blue
        }
    }

    printf("Executing Forward QFFT...\n");
    forward_qfft(input, freq, width, height);

    printf("Executing Inverse QFFT...\n");
    inverse_qfft(freq, reconstructed, width, height);

    printf("\nSample Pixel Check (Original vs Reconstructed):\n");
    printf("Original RGB:      (%.1f, %.1f, %.1f)\n", input[2][2].x, input[2][2].y, input[2][2].z);
    printf("Reconstructed RGB: (%.1f, %.1f, %.1f)\n", reconstructed[2][2].x, reconstructed[2][2].y, reconstructed[2][2].z);

    free_qmatrix(input, height);
    free_qmatrix(freq, height);
    free_qmatrix(reconstructed, height);

    return 0;
}
