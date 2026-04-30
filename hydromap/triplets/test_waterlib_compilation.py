#!/usr/bin/env python3

import numpy as np
import sys

print("\nTesting waterlib compilation...")

try:
    import waterlib as wl
except ImportError:
    sys.exit(
        "\nERROR: Could not import waterlib.\n"
        "Please compile it first:\n\n"
        "    python setup.py build_ext --inplace\n\n"
        "Make sure waterlib.c exists and is compiled into a .so file.\n"
    )

# make tiny fake data (one oxygen, a few waters nearby)
subPos = np.array([[0.0, 0.0, 0.0]])           # one central atom
Pos    = np.array([
    [0.5, 0.0, 0.0],
    [0.0, 0.5, 0.0],
    [0.0, 0.0, 0.5],
    [1.5, 0.0, 0.0],                           # outside cutoff
])
BoxDims = np.array([10.0, 10.0, 10.0])

lowCut  = 0.0
highCut = 1.0

try:
    angles, counts = wl.triplet_angles(subPos, Pos, BoxDims, lowCut, highCut)
except Exception as e:
    sys.exit(
        f"\nERROR: waterlib compiled, but triplet_angles call failed.\n"
        f"Reason: {e}\n"
    )

# test expected behavior
if len(angles) == 3:
    print("✅ waterlib imported and triplet_angles worked correctly!")
    print(f"Computed angles (degrees): {angles}")
else:
    sys.exit(
        "\nERROR: waterlib triplet_angles returned unexpected result.\n"
        "Please recompile waterlib or check your build settings.\n"
    )
