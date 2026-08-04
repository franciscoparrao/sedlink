#!/usr/bin/env python3
"""Independent reference implementation of the Index of Connectivity.

Implements IC = log10(D_up / D_dn) following Borselli et al. (2008) and
Cavalli et al. (2013) with NumPy, matching sedlink's documented numeric
conventions:

- priority-flood-epsilon pit filling (Barnes 2014) in float32, raising
  to ``nextafter(spill)`` with (elevation, index) min-heap ordering;
- D8 steepest-descent flow directions by drop/distance, cardinals
  checked first (SurtGIS convention 1=E..8=SE, CCW from East);
- slope from central differences on the original DEM, stored as float32;
- slope gradient tan(theta) clamped to [min_slope, 1.0];
- D_up = Wbar * Sbar * sqrt(A) with upslope means over the contributing
  area (cell count including self);
- D_dn = sum d_i / (W_i * S_i) along the D8 path to the first stream
  cell; stream cells get IC = +clamp; unreachable cells are NaN.

This is a cross-implementation check written from the papers' equations,
NOT SedInConnect itself; true SedInConnect parity remains a separate task.

Usage: python3 tools/reference_ic.py
Regenerates the fixture files under crates/core/tests/data/.
"""

import heapq
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "crates", "core", "tests", "data")

ROWS, COLS = 40, 40
CELLSIZE = 10.0
THRESHOLD = 40.0
MIN_SLOPE = 0.005
MAX_SLOPE = 1.0
CLAMP = 10.0
MIN_WEIGHT = 0.001

# SurtGIS D8 convention: 1=E, 2=NE, 3=N, 4=NW, 5=W, 6=SW, 7=S, 8=SE.
D8_OFFSETS = [(0, 1), (-1, 1), (-1, 0), (-1, -1),
              (0, -1), (1, -1), (1, 0), (1, 1)]
D8_DISTANCE = [1.0, np.sqrt(2.0), 1.0, np.sqrt(2.0),
               1.0, np.sqrt(2.0), 1.0, np.sqrt(2.0)]
# Cardinals first (E, N, W, S), then diagonals (NE, NW, SW, SE).
CHECK_ORDER = [0, 2, 4, 6, 1, 3, 5, 7]


def make_dem():
    """Deterministic synthetic DEM: south slope + valley + fractional
    noise + an interior depression + a NoData patch."""
    z = np.zeros((ROWS, COLS), dtype=np.float32)
    for r in range(ROWS):
        for c in range(COLS):
            base = 100.0 - r * 0.9 - abs(c - 20) * 0.35
            noise = ((r * 7 + c * 13) % 5) * 0.07
            z[r, c] = np.float32(base + noise)
    # Interior depression (fractional depth).
    z[14:18, 18:22] -= np.float32(2.3)
    # NoData patch.
    z[5:9, 5:9] = np.nan
    return z


def make_weight():
    """Deterministic weighting factor in [0.2, 0.8]."""
    w = np.zeros((ROWS, COLS), dtype=np.float64)
    for r in range(ROWS):
        for c in range(COLS):
            w[r, c] = 0.2 + 0.6 * (((r * 3 + c * 5) % 7) / 6.0)
    return w


def priority_fill(z, solid):
    filled = z.copy()  # float32
    visited = np.zeros((ROWS, COLS), dtype=bool)
    heap = []
    for c in range(COLS):
        for r in (0, ROWS - 1):
            if not solid[r, c] and not visited[r, c]:
                visited[r, c] = True
                heapq.heappush(heap, (float(filled[r, c]), r * COLS + c))
    for r in range(1, ROWS - 1):
        for c in (0, COLS - 1):
            if not solid[r, c] and not visited[r, c]:
                visited[r, c] = True
                heapq.heappush(heap, (float(filled[r, c]), r * COLS + c))

    nbrs = [(-1, 0), (1, 0), (0, -1), (0, 1),
            (-1, -1), (-1, 1), (1, -1), (1, 1)]
    while heap:
        spill, idx = heapq.heappop(heap)
        spill32 = np.float32(spill)
        r, c = divmod(idx, COLS)
        for dr, dc in nbrs:
            nr, nc = r + dr, c + dc
            if not (0 <= nr < ROWS and 0 <= nc < COLS):
                continue
            if visited[nr, nc] or solid[nr, nc]:
                continue
            visited[nr, nc] = True
            if filled[nr, nc] <= spill32:
                filled[nr, nc] = np.nextafter(
                    spill32, np.float32(np.inf), dtype=np.float32)
            heapq.heappush(heap, (float(filled[nr, nc]), nr * COLS + nc))
    return filled


def flow_dir(filled, solid):
    fdir = np.zeros((ROWS, COLS), dtype=np.uint8)
    for r in range(ROWS):
        for c in range(COLS):
            if solid[r, c]:
                continue
            here = float(filled[r, c])
            best_dir, best_grad = 0, 0.0
            for k in CHECK_ORDER:
                dr, dc = D8_OFFSETS[k]
                nr, nc = r + dr, c + dc
                if not (0 <= nr < ROWS and 0 <= nc < COLS):
                    continue
                if solid[nr, nc]:
                    continue
                grad = (here - float(filled[nr, nc])) / D8_DISTANCE[k]
                if grad > best_grad:
                    best_grad, best_dir = grad, k + 1
            fdir[r, c] = best_dir
    return fdir


def downstream_of(fdir):
    ds = np.full(ROWS * COLS, -1, dtype=np.int64)
    for r in range(ROWS):
        for c in range(COLS):
            d = fdir[r, c]
            if d == 0:
                continue
            dr, dc = D8_OFFSETS[d - 1]
            ds[r * COLS + c] = (r + dr) * COLS + (c + dc)
    return ds


def accumulate(ds, solid_flat, local):
    """Topological upslope accumulation (including self)."""
    n = ROWS * COLS
    total = local.astype(np.float64).copy()
    pending = np.zeros(n, dtype=np.int64)
    for i in range(n):
        if solid_flat[i] or ds[i] < 0:
            continue
        pending[ds[i]] += 1
    stack = [i for i in range(n) if not solid_flat[i] and pending[i] == 0]
    while stack:
        i = stack.pop()
        j = ds[i]
        if j >= 0:
            total[j] += total[i]
            pending[j] -= 1
            if pending[j] == 0:
                stack.append(j)
    return total


def slope_radians(z, solid):
    """Central differences on the original DEM, one-sided at edges,
    stored as float32 (matching sedlink)."""
    slope = np.zeros((ROWS, COLS), dtype=np.float32)

    def sample(r, c):
        if not (0 <= r < ROWS and 0 <= c < COLS) or solid[r, c]:
            return None
        return float(z[r, c])

    for r in range(ROWS):
        for c in range(COLS):
            if solid[r, c]:
                continue
            here = sample(r, c)
            pairs = []
            for (a, b) in (((r, c - 1), (r, c + 1)), ((r - 1, c), (r + 1, c))):
                va, vb = sample(*a), sample(*b)
                if va is not None and vb is not None:
                    g = (vb - va) / (2.0 * CELLSIZE)
                elif va is not None:
                    g = (here - va) / CELLSIZE
                elif vb is not None:
                    g = (vb - here) / CELLSIZE
                else:
                    g = 0.0
                pairs.append(g)
            mag = np.sqrt(pairs[0] ** 2 + pairs[1] ** 2)
            slope[r, c] = np.float32(min(max(np.arctan(mag), 0.0), np.pi / 2))
    return slope


def compute_ic(z, weight):
    solid = ~np.isfinite(z)
    solid_flat = solid.ravel()
    filled = priority_fill(z, solid)
    fdir = flow_dir(filled, solid)
    ds = downstream_of(fdir)
    acc = accumulate(ds, solid_flat, np.ones(ROWS * COLS))
    slope = slope_radians(z, solid)

    s_local = np.clip(np.tan(slope.astype(np.float64)),
                      MIN_SLOPE, MAX_SLOPE).ravel()
    w_local = np.maximum(weight, MIN_WEIGHT).ravel()

    is_stream = acc >= THRESHOLD

    # D_dn: memoized downstream trace to the first stream cell.
    n = ROWS * COLS
    d_dn = np.full(n, np.nan)
    d_dn[is_stream] = 0.0
    fdir_flat = fdir.ravel()

    def resolve(i):
        chain = []
        cur = i
        while np.isnan(d_dn[cur]):
            d = fdir_flat[cur]
            if d == 0 or ds[cur] < 0:
                return  # unreachable; leave NaN
            chain.append(cur)
            cur = ds[cur]
        base = d_dn[cur]
        for j in reversed(chain):
            step = D8_DISTANCE[fdir_flat[j] - 1] * CELLSIZE
            base = base + step / (w_local[j] * s_local[j])
            d_dn[j] = base

    for i in range(n):
        if not solid_flat[i] and np.isnan(d_dn[i]):
            resolve(i)

    # Guard against unreachable chains being re-walked endlessly: resolve
    # leaves NaN, and re-walking is harmless for fixture generation.

    # D_up with upslope means.
    sum_w = accumulate(ds, solid_flat, w_local)
    sum_s = accumulate(ds, solid_flat, s_local)
    w_mean = sum_w / acc
    s_mean = sum_s / acc
    d_up = w_mean * s_mean * np.sqrt(acc * CELLSIZE * CELLSIZE)

    ic = np.full(n, np.nan)
    for i in range(n):
        if solid_flat[i]:
            continue
        if np.isnan(d_dn[i]):
            continue
        if d_dn[i] < 1e-10:
            ic[i] = CLAMP
        else:
            ic[i] = min(max(np.log10(d_up[i] / d_dn[i]), -CLAMP), CLAMP)

    d_up[solid_flat] = np.nan
    d_dn[solid_flat] = np.nan
    acc = acc.copy()

    return (acc.reshape(ROWS, COLS), d_up.reshape(ROWS, COLS),
            d_dn.reshape(ROWS, COLS), ic.reshape(ROWS, COLS))


def write_grid(path, arr, fmt):
    with open(path, "w", encoding="ascii") as f:
        f.write(f"{ROWS} {COLS} {CELLSIZE}\n")
        for r in range(arr.shape[0]):
            f.write(" ".join(fmt % v for v in arr[r]) + "\n")


def main():
    os.makedirs(DATA, exist_ok=True)
    z = make_dem()
    w = make_weight()
    acc, d_up, d_dn, ic = compute_ic(z, w)

    write_grid(os.path.join(DATA, "parity_dem.txt"), z, "%.9g")
    write_grid(os.path.join(DATA, "parity_weight.txt"), w, "%.17g")
    write_grid(os.path.join(DATA, "parity_acc.txt"), acc, "%.17g")
    write_grid(os.path.join(DATA, "parity_d_up.txt"), d_up, "%.17g")
    write_grid(os.path.join(DATA, "parity_d_dn.txt"), d_dn, "%.17g")
    write_grid(os.path.join(DATA, "parity_ic.txt"), ic, "%.17g")

    valid = ic[np.isfinite(ic)]
    print(f"fixture written to {os.path.normpath(DATA)}")
    print(f"IC: valid={valid.size} min={valid.min():.4f} "
          f"max={valid.max():.4f} nan={np.isnan(ic).sum()}")
    print(f"streams: {(acc >= THRESHOLD).sum()} cells")


if __name__ == "__main__":
    main()
