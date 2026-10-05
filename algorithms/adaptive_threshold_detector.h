#ifndef ADAPTIVE_THRESHOLD_DETECTOR_H
#define ADAPTIVE_THRESHOLD_DETECTOR_H

#include <stdbool.h>
#include <stdint.h>

#include "accel_sample.h"

typedef struct
{
    double baseline_alpha;
    double activity_alpha;

    double sensitivity;

    double min_threshold_g;
    double max_threshold_g;

    uint32_t min_event_interval_ms;

} adaptive_threshold_config_t;


typedef struct
{
    double magnitude_g;
    double baseline_g;
    double activity_g;
    double threshold_g;

    bool event_detected;

} adaptive_threshold_output_t;


typedef struct
{
    adaptive_threshold_config_t config;

    bool initialized;

    double baseline_g;
    double activity_g;
    double threshold_g;

    double previous_magnitude_g;
    double previous_threshold_g;

    bool has_last_event;
    uint32_t last_event_timestamp_ms;

    uint32_t event_count;

} adaptive_threshold_detector_t;


void adaptive_threshold_detector_default_config(
    adaptive_threshold_config_t *config);


void adaptive_threshold_detector_init(
    adaptive_threshold_detector_t *detector,
    const adaptive_threshold_config_t *config);


bool adaptive_threshold_detector_process(
    adaptive_threshold_detector_t *detector,
    const accel_sample_t *sample,
    adaptive_threshold_output_t *output);


uint32_t adaptive_threshold_detector_get_count(
    const adaptive_threshold_detector_t *detector);


#endif /* ADAPTIVE_THRESHOLD_DETECTOR_H */