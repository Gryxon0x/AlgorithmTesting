#ifndef ALGORITHM_RUNNER_H
#define ALGORITHM_RUNNER_H

#include <stddef.h>

#include "accel_sample.h"

typedef struct
{
    double sample_rate_hz;

    size_t hop_size;

    /*
     * Frequency band in which we search
     * for walking-related periodicity.
     */
    double search_low_hz;
    double search_high_hz;

} frequency_detector_config_t;


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