#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#include <limits.h>
#include <stdint.h>
#include <stddef.h>

extern int spacox_cgf_compute(const double *residuals,
                              size_t residual_count,
                              size_t zero_count,
                              const double *t_values,
                              size_t t_count,
                              size_t thread_count,
                              size_t use_simd,
                              double *k0,
                              double *k1,
                              double *k2);
extern int spacox_cgf_simd_available(void);

/* Returning to the C layer between batches lets R process user interrupts.
 * This mirrors the R backend's work budget: at most 256 points and roughly
 * five million residual/point operations between interrupt checks. */
#define CGF_MAX_BATCH_SIZE 256
#define CGF_BATCH_WORK_BUDGET 5000000

SEXP C_spacox_cgf_simd_available_call(void)
{
    return Rf_ScalarLogical(spacox_cgf_simd_available() != 0);
}

SEXP C_spacox_empirical_cgf_call(SEXP residuals,
                                  SEXP t_values,
                                  SEXP zero_count,
                                  SEXP thread_count,
                                  SEXP use_simd)
{
    if (TYPEOF(residuals) != REALSXP)
        Rf_error("residuals must be a double vector");
    if (TYPEOF(t_values) != REALSXP)
        Rf_error("t_values must be a double vector");

    R_xlen_t residual_count_x = XLENGTH(residuals);
    R_xlen_t t_count_x = XLENGTH(t_values);
    if (residual_count_x == 0)
        Rf_error("residuals must not be empty");
    if (t_count_x == 0)
        Rf_error("t_values must not be empty");
    if (t_count_x > INT_MAX)
        Rf_error("too many empirical CGF grid points");

    double zero_count_value = Rf_asReal(zero_count);
    if (!R_FINITE(zero_count_value) ||
        zero_count_value < 0 ||
        zero_count_value > (double) SIZE_MAX ||
        zero_count_value != (double) ((size_t) zero_count_value))
        Rf_error("zero_count must be a non-negative whole number");

    int thread_count_value = Rf_asInteger(thread_count);
    if (thread_count_value == NA_INTEGER || thread_count_value < 0)
        Rf_error("thread_count must be a non-negative integer");

    int use_simd_value = Rf_asLogical(use_simd);
    if (use_simd_value == NA_LOGICAL)
        Rf_error("use_simd must be TRUE or FALSE");

    SEXP output = PROTECT(Rf_allocMatrix(REALSXP, (int) t_count_x, 3));
    double *output_data = REAL(output);
    size_t residual_count_size = (size_t) residual_count_x;
    size_t t_count_size = (size_t) t_count_x;
    size_t batch_capacity = CGF_BATCH_WORK_BUDGET / residual_count_size;
    if (batch_capacity < 1)
        batch_capacity = 1;
    if (batch_capacity > CGF_MAX_BATCH_SIZE)
        batch_capacity = CGF_MAX_BATCH_SIZE;

    R_CheckUserInterrupt();
    for (size_t offset = 0; offset < t_count_size; offset += batch_capacity) {
        size_t remaining = t_count_size - offset;
        size_t batch_size = remaining < batch_capacity
                                ? remaining
                                : batch_capacity;

        int status = spacox_cgf_compute(REAL(residuals),
                                        residual_count_size,
                                        (size_t) zero_count_value,
                                        REAL(t_values) + offset,
                                        batch_size,
                                        (size_t) thread_count_value,
                                        (size_t) use_simd_value,
                                        output_data + offset,
                                        output_data + t_count_size + offset,
                                        output_data + 2 * t_count_size + offset);

        if (status != 0) {
            UNPROTECT(1);
            if (status == 2)
                Rf_error("Rust CGF backend received a non-finite input");
            if (status == 3)
                Rf_error("Rust CGF backend panicked while calculating the CGF");
            Rf_error("Rust CGF backend received an invalid input");
        }

        /* This is safe because no Rust frames or worker threads are active
         * after spacox_cgf_compute returns. */
        R_CheckUserInterrupt();
    }

    UNPROTECT(1);
    return output;
}
