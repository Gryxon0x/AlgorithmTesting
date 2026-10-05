#include "threshold_detector.h"

#include <math.h>
#include <stddef.h>


void threshold_detector_default_config(
    threshold_detector_config_t *config)
{
    if (config == NULL)
    {
        return;
    }

    config->threshold_g = 1.10;
    config->min_event_interval_ms = 280U;
}


void threshold_detector_init(
    threshold_detector_t *detector,
    const threshold_detector_config_t *config)
{
    if ((detector == NULL) || (config == NULL))
    {
        return;
    }

    detector->config = *config;

    detector->has_previous_sample = false;
    detector->previous_magnitude_g = 0.0;

    detector->has_last_event = false;
    detector->last_event_timestamp_ms = 0U;

    detector->event_count = 0U;
}


bool threshold_detector_process(
    threshold_detector_t *detector,
    const accel_sample_t *sample,
    double *magnitude_g)
{
    double magnitude;
    bool threshold_crossed = false;
    bool interval_ok = true;
    bool event_detected = false;

    if ((detector == NULL) || (sample == NULL))
    {
        return false;
    }

    magnitude = sqrt(
        sample->x_g * sample->x_g +
        sample->y_g * sample->y_g +
        sample->z_g * sample->z_g);

    if (magnitude_g != NULL)
    {
        *magnitude_g = magnitude;
    }

    if (detector->has_previous_sample)
    {
        threshold_crossed =
            (detector->previous_magnitude_g < detector->config.threshold_g) &&
            (magnitude >= detector->config.threshold_g);

        if (detector->has_last_event)
        {
            const uint32_t elapsed_ms =
                sample->timestamp_ms -
                detector->last_event_timestamp_ms;

            interval_ok =
                elapsed_ms >= detector->config.min_event_interval_ms;
        }

        if (threshold_crossed && interval_ok)
        {
            detector->event_count++;
            detector->last_event_timestamp_ms =
                sample->timestamp_ms;
            detector->has_last_event = true;

            event_detected = true;
        }
    }

    detector->previous_magnitude_g = magnitude;
    detector->has_previous_sample = true;

    return event_detected;
}


uint32_t threshold_detector_get_count(
    const threshold_detector_t *detector)
{
    if (detector == NULL)
    {
        return 0U;
    }

    return detector->event_count;
}