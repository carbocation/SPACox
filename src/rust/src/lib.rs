use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;
use std::thread;

const SUM_BLOCK_SIZE: usize = 256;

fn calculate_point(
    residuals: &[f64],
    zero_count: usize,
    total_count: usize,
    residual_min: f64,
    residual_max: f64,
    t: f64,
) -> (f64, f64, f64) {
    let mut shift = if t >= 0.0 {
        t * residual_max
    } else {
        t * residual_min
    };
    if zero_count != 0 {
        shift = shift.max(0.0);
    }

    let zero_weight = (zero_count as f64) * (-shift).exp();
    let mut weight_sum = zero_weight;
    let mut weighted_sum = 0.0;
    let mut weighted_square_sum = 0.0;

    for block in residuals.chunks(SUM_BLOCK_SIZE) {
        let mut block_weight = 0.0;
        let mut block_weighted = 0.0;
        let mut block_squared = 0.0;
        for &residual in block {
            let weight = (t * residual - shift).exp();
            block_weight += weight;
            block_weighted += weight * residual;
            block_squared += weight * residual * residual;
        }
        weight_sum += block_weight;
        weighted_sum += block_weighted;
        weighted_square_sum += block_squared;
    }

    let weighted_mean = weighted_sum / weight_sum;
    let weighted_second = weighted_square_sum / weight_sum;
    let mut weighted_variance = weighted_second - weighted_mean * weighted_mean;
    let variance_scale = weighted_second
        .max(weighted_mean * weighted_mean)
        .max(f64::MIN_POSITIVE);

    if !weighted_variance.is_finite() || weighted_variance < f64::EPSILON.sqrt() * variance_scale {
        let mut variance_numerator = zero_weight * weighted_mean * weighted_mean;
        for block in residuals.chunks(SUM_BLOCK_SIZE) {
            let mut block_variance = 0.0;
            for &residual in block {
                let centered = residual - weighted_mean;
                let weight = (t * residual - shift).exp();
                block_variance += weight * centered * centered;
            }
            variance_numerator += block_variance;
        }
        weighted_variance = variance_numerator / weight_sum;
    }

    (
        shift + weight_sum.ln() - (total_count as f64).ln(),
        weighted_mean,
        weighted_variance,
    )
}

fn calculate_chunk(
    residuals: &[f64],
    zero_count: usize,
    total_count: usize,
    residual_min: f64,
    residual_max: f64,
    t_values: &[f64],
    k0: &mut [f64],
    k1: &mut [f64],
    k2: &mut [f64],
) {
    for index in 0..t_values.len() {
        let (point_k0, point_k1, point_k2) = calculate_point(
            residuals,
            zero_count,
            total_count,
            residual_min,
            residual_max,
            t_values[index],
        );
        k0[index] = point_k0;
        k1[index] = point_k1;
        k2[index] = point_k2;
    }
}

fn calculate_all(
    residuals: &[f64],
    zero_count: usize,
    t_values: &[f64],
    requested_threads: usize,
    k0: &mut [f64],
    k1: &mut [f64],
    k2: &mut [f64],
) -> Result<(), ()> {
    if residuals.is_empty() || t_values.is_empty() {
        return Err(());
    }

    if residuals.iter().any(|value| !value.is_finite())
        || t_values.iter().any(|value| !value.is_finite())
    {
        return Err(());
    }

    let residual_min = residuals.iter().copied().fold(f64::INFINITY, f64::min);
    let residual_max = residuals.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let total_count = residuals.len().checked_add(zero_count).ok_or(())?;

    let available_threads = thread::available_parallelism()
        .map(|count| count.get())
        .unwrap_or(1);
    let worker_count = if requested_threads == 0 {
        available_threads
    } else {
        requested_threads.min(available_threads)
    }
    .min(t_values.len())
    .max(1);

    if worker_count == 1 {
        calculate_chunk(
            residuals,
            zero_count,
            total_count,
            residual_min,
            residual_max,
            t_values,
            k0,
            k1,
            k2,
        );
    } else {
        let chunk_size = (t_values.len() + worker_count - 1) / worker_count;
        thread::scope(|scope| {
            for (((t_chunk, k0_chunk), k1_chunk), k2_chunk) in t_values
                .chunks(chunk_size)
                .zip(k0.chunks_mut(chunk_size))
                .zip(k1.chunks_mut(chunk_size))
                .zip(k2.chunks_mut(chunk_size))
            {
                scope.spawn(move || {
                    calculate_chunk(
                        residuals,
                        zero_count,
                        total_count,
                        residual_min,
                        residual_max,
                        t_chunk,
                        k0_chunk,
                        k1_chunk,
                        k2_chunk,
                    );
                });
            }
        });
    }

    if k0.iter().any(|value| !value.is_finite())
        || k1.iter().any(|value| !value.is_finite())
        || k2.iter().any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(());
    }

    Ok(())
}

#[no_mangle]
pub unsafe extern "C" fn spacox_cgf_compute(
    residuals: *const f64,
    residual_count: usize,
    zero_count: usize,
    t_values: *const f64,
    t_count: usize,
    thread_count: usize,
    k0: *mut f64,
    k1: *mut f64,
    k2: *mut f64,
) -> i32 {
    if residuals.is_null()
        || t_values.is_null()
        || k0.is_null()
        || k1.is_null()
        || k2.is_null()
        || residual_count == 0
        || t_count == 0
    {
        return 1;
    }

    let result = catch_unwind(AssertUnwindSafe(|| {
        let residuals = slice::from_raw_parts(residuals, residual_count);
        let t_values = slice::from_raw_parts(t_values, t_count);
        let k0 = slice::from_raw_parts_mut(k0, t_count);
        let k1 = slice::from_raw_parts_mut(k1, t_count);
        let k2 = slice::from_raw_parts_mut(k2, t_count);

        calculate_all(residuals, zero_count, t_values, thread_count, k0, k1, k2)
    }));

    match result {
        Ok(Ok(())) => 0,
        Ok(Err(())) => 2,
        Err(_) => 3,
    }
}

#[cfg(test)]
mod tests {
    use super::calculate_all;

    #[test]
    fn overflow_tail_is_finite() {
        let residuals = [-8.0, 0.5, 1.0];
        let t_values = [-100.0, 0.0, 100.0];
        let mut k0 = [0.0; 3];
        let mut k1 = [0.0; 3];
        let mut k2 = [0.0; 3];

        calculate_all(&residuals, 1, &t_values, 1, &mut k0, &mut k1, &mut k2).unwrap();

        assert!(k0.iter().all(|value| value.is_finite()));
        assert!(k1.iter().all(|value| value.is_finite()));
        assert!(k2.iter().all(|value| value.is_finite() && *value >= 0.0));
        assert!((k0[0] - (800.0 - 4.0_f64.ln())).abs() < 1e-12);
        assert!((k1[0] + 8.0).abs() < 1e-12);
    }

    #[test]
    fn thread_count_does_not_change_results() {
        let residuals: Vec<f64> = (1..=1000)
            .map(|index| ((index as f64) / 37.0).sin())
            .collect();
        let t_values: Vec<f64> = (-50..=50).map(|value| value as f64 / 5.0).collect();
        let mut single_k0 = vec![0.0; t_values.len()];
        let mut single_k1 = vec![0.0; t_values.len()];
        let mut single_k2 = vec![0.0; t_values.len()];
        let mut parallel_k0 = vec![0.0; t_values.len()];
        let mut parallel_k1 = vec![0.0; t_values.len()];
        let mut parallel_k2 = vec![0.0; t_values.len()];

        calculate_all(
            &residuals,
            4,
            &t_values,
            1,
            &mut single_k0,
            &mut single_k1,
            &mut single_k2,
        )
        .unwrap();
        calculate_all(
            &residuals,
            4,
            &t_values,
            4,
            &mut parallel_k0,
            &mut parallel_k1,
            &mut parallel_k2,
        )
        .unwrap();

        assert_eq!(single_k0, parallel_k0);
        assert_eq!(single_k1, parallel_k1);
        assert_eq!(single_k2, parallel_k2);
    }
}
