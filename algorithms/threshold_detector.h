#ifndef THRESHOLD_DETECTOR_H
#define THRESHOLD_DETECTOR_H

#include <stdbool.h>
#include <stdint.h>

#include "accel_sample.h"

typedef struct
{
    double threshold_g;
    uint32_t min_event_interval_ms;
} threshold_detector_config_t;

typedef struct
{
    threshold_detector_config_t config;

    bool has_previous_sample;
    double previous_magnitude_g;

    bool has_last_event;
    uint32_t last_event_timestamp_ms;

    uint32_t event_count;
} threshold_detector_t;


void threshold_detector_default_config(
    threshold_detector_config_t *config);

void threshold_detector_init(
    threshold_detector_t *detector,
    const threshold_detector_config_t *config);

/*
 * Process one accelerometer sample.
 *
 * Returns true when a new threshold-crossing event is detected.
 *
 * If magnitude_g is non-NULL, the calculated acceleration magnitude
 * is written there for debugging/logging.
 */
bool threshold_detector_process(
    threshold_detector_t *detector,
    const accel_sample_t *sample,
    double *magnitude_g);

uint32_t threshold_detector_get_count(
    const threshold_detector_t *detector);

#endif /* THRESHOLD_DETECTOR_H */