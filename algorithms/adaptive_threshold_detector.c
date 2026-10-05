#include "adaptive_threshold_detector.h"

#include <math.h>
#include <stddef.h>


static double clamp_double(
    double value,
    double minimum,
    double maximum)
{
    if (value < minimum)
    {
        return minimum;
    }

    if (value > maximum)
    {
        return maximum;
    }

    return value;
}


static double calculate_magnitude_g(
    const accel_sample_t *sample)
{
    return sqrt(
        sample->x_g * sample->x_g +
        sample->y_g * sample->y_g +
        sample->z_g * sample->z_g);
}


void adaptive_threshold_detector_default_config(
    adaptive_threshold_config_t *config)
{
    if (config == NULL)
    {
        return;
    }

    config->baseline_alpha = 0.004;
    config->activity_alpha = 0.010;

    config->sensitivity = 1.0;

    config->min_threshold_g = 1.05;
    config->max_threshold_g = 1.35;

    config->min_event_interval_ms = 280U;
}


void adaptive_threshold_detector_init(
    adaptive_threshold_detector_t *detector,
    const adaptive_threshold_config_t *config)
{
    if ((detector == NULL) || (config == NULL))
    {
        return;
    }

    detector->config = *config;

    detector->initialized = false;

    detector->baseline_g = 0.0;
    detector->activity_g = 0.0;
    detector->threshold_g = config->min_threshold_g;

    detector->previous_magnitude_g = 0.0;
    detector->previous_threshold_g =
        config->min_threshold_g;

    detector->has_last_event = false;
    detector->last_event_timestamp_ms = 0U;

    detector->event_count = 0U;
}


bool adaptive_threshold_detector_process(
    adaptive_threshold_detector_t *detector,
    const accel_sample_t *sample,
    adaptive_threshold_output_t *output)
{
    double magnitude_g;
    double threshold_for_this_sample;
    double baseline_for_this_sample;
    double activity_for_this_sample;
    double deviation_g;

    bool crossed_threshold;
    bool interval_ok = true;
    bool event_detected = false;

    if ((detector == NULL) || (sample == NULL))
    {
        return false;
    }

    magnitude_g = calculate_magnitude_g(sample);

    if (!detector->initialized)
    {
        detector->baseline_g = magnitude_g;
        detector->activity_g = 0.0;

        detector->threshold_g =
            clamp_double(
                detector->baseline_g,
                detector->config.min_threshold_g,
                detector->config.max_threshold_g);

        detector->previous_magnitude_g =
            magnitude_g;

        detector->previous_threshold_g =
            detector->threshold_g;

        detector->initialized = true;

        if (output != NULL)
        {
            output->magnitude_g = magnitude_g;
            output->baseline_g =
                detector->baseline_g;
            output->activity_g =
                detector->activity_g;
            output->threshold_g =
                detector->threshold_g;
            output->event_detected = false;
        }

        return false;
    }

    threshold_for_this_sample =
        detector->threshold_g;

    baseline_for_this_sample =
        detector->baseline_g;

    activity_for_this_sample =
        detector->activity_g;

    crossed_threshold =
        (detector->previous_magnitude_g <
         detector->previous_threshold_g) &&
        (magnitude_g >= threshold_for_this_sample);

    if (detector->has_last_event)
    {
        uint32_t elapsed_ms;

        elapsed_ms =
            sample->timestamp_ms -
            detector->last_event_timestamp_ms;

        interval_ok =
            elapsed_ms >=
            detector->config.min_event_interval_ms;
    }

    if (crossed_threshold && interval_ok)
    {
        detector->event_count++;

        detector->last_event_timestamp_ms =
            sample->timestamp_ms;

        detector->has_last_event = true;

        event_detected = true;
    }

    deviation_g =
        fabs(magnitude_g - detector->baseline_g);

    detector->baseline_g +=
        detector->config.baseline_alpha *
        (magnitude_g - detector->baseline_g);

    detector->activity_g +=
        detector->config.activity_alpha *
        (deviation_g - detector->activity_g);

    detector->threshold_g =
        clamp_double(
            detector->baseline_g +
            detector->config.sensitivity *
            detector->activity_g,
            detector->config.min_threshold_g,
            detector->config.max_threshold_g);

    detector->previous_magnitude_g =
        magnitude_g;

    detector->previous_threshold_g =
        threshold_for_this_sample;

    if (output != NULL)
    {
        output->magnitude_g =
            magnitude_g;

        output->baseline_g =
            baseline_for_this_sample;

        output->activity_g =
            activity_for_this_sample;

        output->threshold_g =
            threshold_for_this_sample;

        output->event_detected =
            event_detected;
    }

    return event_detected;
}


uint32_t adaptive_threshold_detector_get_count(
    const adaptive_threshold_detector_t *detector)
{
    if (detector == NULL)
    {
        return 0U;
    }

    return detector->event_count;
}