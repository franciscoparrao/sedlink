//! Flat-surface drainage resolution (Garbrecht & Martz 1997; Barnes,
//! Lehman & Mulla 2014).
//!
//! Priority-flood filling levels depressions to their exact spill
//! elevation, and natural terrain contains plains and terraces: in both
//! cases cells have no strictly lower neighbour and steepest descent
//! stalls. This module assigns drainage across those flats with the
//! two-gradient method:
//!
//! 1. **Towards lower terrain**: a multi-source BFS from the flat's
//!    outlets (cells with a strictly lower neighbour, or on the grid
//!    boundary) labels each flat cell with its hop distance to an
//!    outlet. Descending this field always terminates at an outlet, so
//!    it is the primary key.
//! 2. **Away from higher terrain**: a BFS from cells adjacent to higher
//!    terrain breaks ties so parallel artificial paths hug the flat
//!    centre instead of the high edge (Garbrecht & Martz's second
//!    gradient).
//!
//! Both fields are BFS distances from seed sets, hence independent of
//! traversal order — the resolution is deterministic and reproducible
//! across implementations.
//!
//! Flats with no outlet at all (fully enclosed by `NoData` or higher
//! terrain) keep their cells as pits.

// Flat membership is *exact* elevation equality by design: cells raised
// by priority-flood share the spill value bit-for-bit, and natural flats
// are equal by definition. A tolerance would merge distinct terraces.
#![allow(clippy::float_cmp)]

use crate::network::D8_OFFSETS;

/// Sentinel: cell is not part of a resolvable flat (or unreached).
const UNSET: u32 = u32::MAX;

/// BFS offsets for the cells of flat regions.
pub(crate) struct FlatOffsets {
    /// Hop distance to the nearest flat outlet (UNSET off-flat or when
    /// the flat has no outlet).
    to_lower: Vec<u32>,
    /// Hop distance from the nearest higher-terrain edge (UNSET when
    /// the flat has no high edge or off-flat).
    from_higher: Vec<u32>,
    rows: usize,
    cols: usize,
}

impl FlatOffsets {
    /// Compute flat offsets for a (filled) DEM. Returns `None` when the
    /// DEM has no flat cells needing resolution.
    pub(crate) fn compute(z: &[f32], solid: &[bool], rows: usize, cols: usize) -> Option<Self> {
        let n = rows * cols;

        let neighbours = |idx: usize| {
            let r = (idx / cols) as isize;
            let c = (idx % cols) as isize;
            D8_OFFSETS.iter().filter_map(move |&(dr, dc)| {
                let (nr, nc) = (r + dr, c + dc);
                if nr >= 0 && nc >= 0 && (nr as usize) < rows && (nc as usize) < cols {
                    Some(nr as usize * cols + nc as usize)
                } else {
                    None
                }
            })
        };

        // Cells needing resolution: valid, no strictly lower valid
        // neighbour. (True pits after filling are flats or genuinely
        // closed cells.)
        let mut needs = vec![false; n];
        let mut any = false;
        for idx in 0..n {
            if solid[idx] {
                continue;
            }
            let lower = neighbours(idx).any(|nb| !solid[nb] && z[nb] < z[idx]);
            if !lower {
                needs[idx] = true;
                any = true;
            }
        }
        if !any {
            return None;
        }

        // Flat membership: expand from `needs` cells across equal
        // elevation, so outlet cells of the same flat (which do have a
        // lower neighbour) join the region.
        let mut in_flat = vec![false; n];
        let mut queue: Vec<usize> = (0..n).filter(|&i| needs[i]).collect();
        for &i in &queue {
            in_flat[i] = true;
        }
        while let Some(idx) = queue.pop() {
            for nb in neighbours(idx) {
                if !solid[nb] && !in_flat[nb] && z[nb] == z[idx] {
                    in_flat[nb] = true;
                    queue.push(nb);
                }
            }
        }

        // Seeds.
        let mut to_lower = vec![UNSET; n];
        let mut from_higher = vec![UNSET; n];
        let mut low_q = std::collections::VecDeque::new();
        let mut high_q = std::collections::VecDeque::new();
        for idx in 0..n {
            if !in_flat[idx] {
                continue;
            }
            let on_boundary = {
                let (r, c) = (idx / cols, idx % cols);
                r == 0 || c == 0 || r == rows - 1 || c == cols - 1
            };
            let has_lower = neighbours(idx).any(|nb| !solid[nb] && z[nb] < z[idx]);
            let has_higher = neighbours(idx).any(|nb| !solid[nb] && z[nb] > z[idx]);
            if has_lower || on_boundary {
                to_lower[idx] = 0;
                low_q.push_back(idx);
            }
            if has_higher {
                from_higher[idx] = 0;
                high_q.push_back(idx);
            }
        }

        // Multi-source BFS (order-independent distances).
        let bfs = |mut q: std::collections::VecDeque<usize>, field: &mut Vec<u32>| {
            while let Some(idx) = q.pop_front() {
                let d = field[idx];
                for nb in neighbours(idx) {
                    if in_flat[nb] && z[nb] == z[idx] && field[nb] == UNSET {
                        field[nb] = d + 1;
                        q.push_back(nb);
                    }
                }
            }
        };
        bfs(low_q, &mut to_lower);
        bfs(high_q, &mut from_higher);

        Some(Self {
            to_lower,
            from_higher,
            rows,
            cols,
        })
    }

    /// D8 direction (1–8, `SurtGIS` convention) draining a flat cell, or
    /// `None` when the cell is not on a resolvable flat.
    ///
    /// Chooses, among equal-elevation flat neighbours strictly closer to
    /// an outlet, the one farthest from higher terrain; remaining ties
    /// break cardinal-first.
    pub(crate) fn direction(&self, idx: usize, z: &[f32], solid: &[bool]) -> Option<u8> {
        // Cardinals first (indices into D8_OFFSETS), then diagonals.
        const ORDER: [usize; 8] = [0, 2, 4, 6, 1, 3, 5, 7];

        let own = self.to_lower[idx];
        if own == UNSET || own == 0 {
            // Not resolvable, or an outlet cell (drains normally or off-grid).
            return None;
        }

        let r = (idx / self.cols) as isize;
        let c = (idx % self.cols) as isize;
        let mut best: Option<(u32, i64, u8)> = None; // (to_lower, -from_higher, dir)

        for &k in &ORDER {
            let (dr, dc) = D8_OFFSETS[k];
            let (nr, nc) = (r + dr, c + dc);
            if nr < 0 || nc < 0 || nr as usize >= self.rows || nc as usize >= self.cols {
                continue;
            }
            let nb = nr as usize * self.cols + nc as usize;
            if solid[nb] || z[nb] != z[idx] || self.to_lower[nb] == UNSET {
                continue;
            }
            if self.to_lower[nb] >= own {
                continue;
            }
            let fh = self.from_higher[nb];
            // Prefer smaller to_lower; tie → larger from_higher (farther
            // from high terrain; UNSET counts as farthest); tie → first
            // in cardinal-first order (strict < keeps the earlier hit).
            let key = (
                self.to_lower[nb],
                -i64::from(if fh == UNSET { u32::MAX - 1 } else { fh }),
            );
            if best.is_none_or(|(bl, bh, _)| (key.0, key.1) < (bl, bh)) {
                best = Some((key.0, key.1, (k + 1) as u8));
            }
        }

        best.map(|(_, _, d)| d)
    }
}
