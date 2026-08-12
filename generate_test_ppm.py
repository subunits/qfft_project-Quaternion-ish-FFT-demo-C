#!/usr/bin/env python3
# generate_test_ppm.py
# Creates test/test.ppm (P6) with simple RGB gradients and a radial blue pattern.
import os
import math

w, h = 128, 128
out_dir  = "test"
out_path = os.path.join(out_dir, "test.ppm")

# Create output directory if it doesn't exist
os.makedirs(out_dir, exist_ok=True)

with open(out_path, "wb") as f:
    f.write(b"P6\n%d %d\n255\n" % (w, h))
    cx, cy   = w // 2, h // 2
    max_dist = math.hypot(cx, cy)
    for y in range(h):
        for x in range(w):
            r = int(255 * x / (w - 1))
            g = int(255 * y / (h - 1))
            d = math.hypot(x - cx, y - cy)
            b = int(255 * max(0.0, 1.0 - d / max_dist))
            f.write(bytes((r, g, b)))

print("Wrote", out_path)
