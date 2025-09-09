#!/usr/bin/env python3
# generate_test_ppm.py
# Creates test/test.ppm (P6) with simple RGB gradients and some circles
import sys, math
w,h = 128,128
out_path = "test/test.ppm"
with open(out_path, "wb") as f:
    f.write(b"P6\n%d %d\n255\n" % (w,h))
    for y in range(h):
        for x in range(w):
            r = int(255 * x / (w-1))
            g = int(255 * y / (h-1))
            # add a radial pattern for blue
            cx,cy = w//2, h//2
            d = math.hypot(x-cx, y-cy)
            b = int(255 * max(0, 1 - d/(math.hypot(cx,cy))))
            f.write(bytes((r,g,b)))
print("Wrote", out_path)
