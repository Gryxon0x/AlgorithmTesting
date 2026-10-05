#ifndef PEAK_DETECTOR_H
#define PEAK_DETECTOR_H

#include <stdbool.h>
#include <stddef.h>

#include "accel_sample.h"

typedef struct
{
    double low_hz;
    double high_hz;
    double min_event_interval_s;
    double min_prominence_g;
    bool use_min_height;
    double min_height_g;
    unsigned filter_order;
} peak_detector_config_t;

typedef struct
{
    size_t sample_index;
    double filtered_peak_g;
    double prominence_g;
} peak_event_t;

typedef struct
{
    peak_event_t *events;
    size_t event_capacity;
    size_t event_count;

    /* Caller-owned array with one element per input sample. */
    double *filtered_signal_g;
    size_t filtered_capacity;
} peak_detection_result_t;

void peak_detector_default_config(peak_detector_config_t *config);

/*
 * Offline reference implementation of the Python/SciPy detector.
 *
 * This intentionally uses the full recording because scipy.signal.sosfiltfilt
 * is a zero-phase forward/backward filter and peak prominence is also naturally
 * an offline operation. The function is plain C and independent of CSV/OS APIs,
 * but it is not the final streaming MCU implementation.
 *
 * Returns 0 on success, non-zero on invalid input/allocation failure.
 */
int peak_detector_detect(
    const accel_sample_t *samples,
    size_t sample_count,
    double sample_rate_hz,
    const peak_detector_config_t *config,
    peak_detection_result_t *result);

#endif /* PEAK_DETECTOR_H */
