#ifndef ALGORITHM_RUNNER_H
#define ALGORITHM_RUNNER_H

#include <stddef.h>

#include "accel_sample.h"

/*
 * Runs all currently implemented step-detection algorithms
 * against the supplied recording.
 *
 * Returns:
 *   0  on success
 *  -1  on failure
 */
int run_all_algorithms(
    const accel_sample_t *samples,
    size_t sample_count,
    double sample_rate_hz,
    const char *debug_output_path);


#endif /* ALGORITHM_RUNNER_H */