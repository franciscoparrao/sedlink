//! Row-parallel execution helper.
//!
//! Per-cell kernels (flow directions, slope, D∞ facets) are
//! embarrassingly parallel by row. With the `parallel` feature (default)
//! rows are processed on the rayon thread pool; without it the same
//! closure runs sequentially, so results are identical either way.

/// Apply `f(row_index, row_slice)` to each row of a row-major grid.
pub(crate) fn for_each_row<T: Send, F>(data: &mut [T], cols: usize, f: F)
where
    F: Fn(usize, &mut [T]) + Sync + Send,
{
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        data.par_chunks_mut(cols)
            .enumerate()
            .for_each(|(r, row)| f(r, row));
    }
    #[cfg(not(feature = "parallel"))]
    {
        data.chunks_mut(cols)
            .enumerate()
            .for_each(|(r, row)| f(r, row));
    }
}
