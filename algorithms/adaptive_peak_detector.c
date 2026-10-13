#include "adaptive_peak_detector.h"

#include <math.h>
#include <stddef.h>


#define TWO_PI 6.28318530717958647692


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


static double filter_sample(
    adaptive_peak_detector_t *detector,
    double input_g)
{
    double hp_output;
    double lp_output;

    /*
     * Initialize the filter from the first input sample.
     */
    if (!detector->filter_initialized)
    {
        detector->hp_previous_input_g = input_g;
        detector->hp_previous_output_g = 0.0;
        detector->lp_previous_output_g = 0.0;

        detector->filter_initialized = true;

        return 0.0;
    }

    /*
     * First-order high-pass:
     *
     * y[n] = alpha *
     *        (y[n-1] + x[n] - x[n-1])
     */
    hp_output =
        detector->hp_alpha *
        (
            detector->hp_previous_output_g +
            input_g -
            detector->hp_previous_input_g
        );

    detector->hp_previous_input_g = input_g;
    detector->hp_previous_output_g = hp_output;

    /*
     * First-order low-pass:
     *
     * y[n] =
     * y[n-1] + alpha * (x[n] - y[n-1])
     */
    lp_output =
        detector->lp_previous_output_g +
        detector->lp_alpha *
        (
            hp_output -
            detector->lp_previous_output_g
        );

    detector->lp_previous_output_g = lp_output;

    return lp_output;
}


void adaptive_peak_detector_default_config(
    adaptive_peak_config_t *config)
{
    if (config == NULL)
    {
        return;
    }

    config->sample_rate_hz = 50.0;

    config->low_cutoff_hz = 0.9;
    config->high_cutoff_hz = 8.0;

    /*
     * Starting adaptation rates.
     *
     * These will later be tuned using the longer
     * chronological recordings.
     */
    config->baseline_alpha = 0.010;
    config->activity_alpha = 0.020;

    config->height_multiplier = 1.5;
    config->prominence_multiplier = 1.0;

    config->min_peak_height_g = 0.030;
    config->max_peak_height_g = 0.500;

    config->min_prominence_g = 0.020;
    config->max_prominence_g = 0.500;

    /*
     * At 50 Hz, 120 ms is approximately six samples.
     */
    config->confirmation_window_ms = 120U;

    config->min_event_interval_ms = 280U;
}


int adaptive_peak_detector_init(
    adaptive_peak_detector_t *detector,
    const adaptive_peak_config_t *config)
{
    double dt_s;
    double hp_rc_s;
    double lp_rc_s;

    if ((detector == NULL) || (config == NULL))
    {
        return -1;
    }

    if (!(config->sample_rate_hz > 0.0))
    {
        return -1;
    }

    if (!(config->low_cutoff_hz > 0.0))
    {
        return -1;
    }

    if (!(config->high_cutoff_hz >
          config->low_cutoff_hz))
    {
        return -1;
    }

    if (!(config->high_cutoff_hz <
          (config->sample_rate_hz / 2.0)))
    {
        return -1;
    }

    if ((config->baseline_alpha <= 0.0) ||
        (config->baseline_alpha > 1.0))
    {
        return -1;
    }

    if ((config->activity_alpha <= 0.0) ||
        (config->activity_alpha > 1.0))
    {
        return -1;
    }

    if (config->min_peak_height_g >
        config->max_peak_height_g)
    {
        return -1;
    }

    if (config->min_prominence_g >
        config->max_prominence_g)
    {
        return -1;
    }

    detector->config = *config;

    /*
     * Calculate first-order RC filter coefficients.
     */
    dt_s = 1.0 / config->sample_rate_hz;

    hp_rc_s =
        1.0 /
        (TWO_PI * config->low_cutoff_hz);

    lp_rc_s =
        1.0 /
        (TWO_PI * config->high_cutoff_hz);

    detector->hp_alpha =
        hp_rc_s /
        (hp_rc_s + dt_s);

    detector->lp_alpha =
        dt_s /
        (lp_rc_s + dt_s);

    detector->filter_initialized = false;

    detector->hp_previous_input_g = 0.0;
    detector->hp_previous_output_g = 0.0;
    detector->lp_previous_output_g = 0.0;

    detector->adaptive_initialized = false;

    detector->baseline_g = 0.0;
    detector->activity_g = 0.0;

    detector->required_height_g =
        config->min_peak_height_g;

    detector->required_prominence_g =
        config->min_prominence_g;

    detector->history_count = 0U;

    detector->prev2_filtered_g = 0.0;
    detector->prev1_filtered_g = 0.0;

    detector->prev1_timestamp_ms = 0U;

    detector->prev1_required_height_g =
        config->min_peak_height_g;

    detector->prev1_required_prominence_g =
        config->min_prominence_g;

    detector->prev1_baseline_g = 0.0;
    detector->prev1_activity_g = 0.0;

    detector->left_valley_initialized = false;
    detector->left_valley_g = 0.0;

    detector->candidate_active = false;

    detector->candidate_peak_g = 0.0;
    detector->candidate_timestamp_ms = 0U;

    detector->candidate_left_valley_g = 0.0;
    detector->candidate_right_valley_g = 0.0;

    detector->candidate_required_height_g =
        config->min_peak_height_g;

    detector->candidate_required_prominence_g =
        config->min_prominence_g;

    detector->candidate_baseline_g = 0.0;
    detector->candidate_activity_g = 0.0;

    detector->has_last_event = false;
    detector->last_event_timestamp_ms = 0U;

    detector->event_count = 0U;

    return 0;
}


bool adaptive_peak_detector_process(
    adaptive_peak_detector_t *detector,
    const accel_sample_t *sample,
    adaptive_peak_output_t *output)
{
    double magnitude_g;
    double filtered_g;

    double baseline_for_this_sample;
    double activity_for_this_sample;

    double required_height_for_this_sample;
    double required_prominence_for_this_sample;

    double deviation_g;

    bool event_detected = false;
    uint32_t event_timestamp_ms = 0U;

    double confirmed_prominence_g = 0.0;

    if ((detector == NULL) || (sample == NULL))
    {
        return false;
    }

    magnitude_g = calculate_magnitude_g(sample);

    filtered_g =
        filter_sample(detector, magnitude_g);

    /*
     * Initialize adaptive statistics.
     */
    if (!detector->adaptive_initialized)
    {
        detector->baseline_g = filtered_g;
        detector->activity_g = 0.0;

        detector->required_height_g =
            detector->config.min_peak_height_g;

        detector->required_prominence_g =
            detector->config.min_prominence_g;

        detector->adaptive_initialized = true;
    }

    /*
     * Use adaptive state from previous samples
     * for the current classification.
     */
    baseline_for_this_sample =
        detector->baseline_g;

    activity_for_this_sample =
        detector->activity_g;

    required_height_for_this_sample =
        detector->required_height_g;

    required_prominence_for_this_sample =
        detector->required_prominence_g;

    /*
     * If a candidate peak is active, continue
     * looking for its right-side valley.
     */
    if (detector->candidate_active)
    {
        uint32_t elapsed_ms;

        if (filtered_g <
            detector->candidate_right_valley_g)
        {
            detector->candidate_right_valley_g =
                filtered_g;
        }

        elapsed_ms =
            sample->timestamp_ms -
            detector->candidate_timestamp_ms;

        if (elapsed_ms >=
            detector->config.confirmation_window_ms)
        {
            double reference_valley_g;
            bool height_ok;
            bool prominence_ok;
            bool interval_ok = true;

            /*
             * Local prominence approximation:
             *
             * peak - max(left valley, right valley)
             */
            reference_valley_g =
                detector->candidate_left_valley_g;

            if (detector->candidate_right_valley_g >
                reference_valley_g)
            {
                reference_valley_g =
                    detector->candidate_right_valley_g;
            }

            confirmed_prominence_g =
                detector->candidate_peak_g -
                reference_valley_g;

            height_ok =
                detector->candidate_peak_g >=
                detector->candidate_required_height_g;

            prominence_ok =
                confirmed_prominence_g >=
                detector->candidate_required_prominence_g;

            if (detector->has_last_event)
            {
                const uint32_t event_interval_ms =
                    detector->candidate_timestamp_ms -
                    detector->last_event_timestamp_ms;

                interval_ok =
                    event_interval_ms >=
                    detector->config.min_event_interval_ms;
            }

            if (height_ok &&
                prominence_ok &&
                interval_ok)
            {
                detector->event_count++;

                detector->last_event_timestamp_ms =
                    detector->candidate_timestamp_ms;

                detector->has_last_event = true;

                event_detected = true;

                event_timestamp_ms =
                    detector->candidate_timestamp_ms;
            }

            /*
             * The right valley becomes the starting
             * left valley for the next candidate.
             */
            detector->left_valley_g =
                detector->candidate_right_valley_g;

            detector->left_valley_initialized = true;

            detector->candidate_active = false;
        }
    }
    else
    {
        /*
         * Detect a local maximum with one-sample delay:
         *
         * prev1 > prev2
         * prev1 >= current
         */
        if ((detector->history_count >= 2U) &&
            (detector->prev1_filtered_g >
             detector->prev2_filtered_g) &&
            (detector->prev1_filtered_g >=
             filtered_g))
        {
            detector->candidate_active = true;

            detector->candidate_peak_g =
                detector->prev1_filtered_g;

            detector->candidate_timestamp_ms =
                detector->prev1_timestamp_ms;

            detector->candidate_required_height_g =
                detector->prev1_required_height_g;

            detector->candidate_required_prominence_g =
                detector->prev1_required_prominence_g;

            detector->candidate_baseline_g =
                detector->prev1_baseline_g;

            detector->candidate_activity_g =
                detector->prev1_activity_g;

            if (detector->left_valley_initialized)
            {
                detector->candidate_left_valley_g =
                    detector->left_valley_g;
            }
            else
            {
                detector->candidate_left_valley_g =
                    detector->prev2_filtered_g;
            }

            /*
             * Current sample is already on the
             * right side of the detected peak.
             */
            detector->candidate_right_valley_g =
                filtered_g;
        }
        else
        {
            /*
             * Continue tracking the lowest point
             * before the next candidate peak.
             */
            if (!detector->left_valley_initialized)
            {
                detector->left_valley_g =
                    filtered_g;

                detector->left_valley_initialized =
                    true;
            }
            else if (filtered_g <
                     detector->left_valley_g)
            {
                detector->left_valley_g =
                    filtered_g;
            }
        }
    }

    /*
     * Update adaptive statistics after making the
     * classification decision for the current sample.
     */
    deviation_g =
        fabs(
            filtered_g -
            detector->baseline_g);

    detector->baseline_g +=
        detector->config.baseline_alpha *
        (
            filtered_g -
            detector->baseline_g
        );

    detector->activity_g +=
        detector->config.activity_alpha *
        (
            deviation_g -
            detector->activity_g
        );

    detector->required_height_g =
        clamp_double(
            detector->baseline_g +
            detector->config.height_multiplier *
            detector->activity_g,

            detector->config.min_peak_height_g,
            detector->config.max_peak_height_g);

    detector->required_prominence_g =
        clamp_double(
            detector->config.prominence_multiplier *
            detector->activity_g,

            detector->config.min_prominence_g,
            detector->config.max_prominence_g);

    /*
     * Store current sample so that on the next call
     * it can become the center of a local maximum.
     */
    if (detector->history_count == 0U)
    {
        detector->prev1_filtered_g =
            filtered_g;

        detector->prev1_timestamp_ms =
            sample->timestamp_ms;

        detector->prev1_required_height_g =
            required_height_for_this_sample;

        detector->prev1_required_prominence_g =
            required_prominence_for_this_sample;

        detector->prev1_baseline_g =
            baseline_for_this_sample;

        detector->prev1_activity_g =
            activity_for_this_sample;

        detector->history_count = 1U;
    }
    else
    {
        detector->prev2_filtered_g =
            detector->prev1_filtered_g;

        detector->prev1_filtered_g =
            filtered_g;

        detector->prev1_timestamp_ms =
            sample->timestamp_ms;

        detector->prev1_required_height_g =
            required_height_for_this_sample;

        detector->prev1_required_prominence_g =
            required_prominence_for_this_sample;

        detector->prev1_baseline_g =
            baseline_for_this_sample;

        detector->prev1_activity_g =
            activity_for_this_sample;

        if (detector->history_count < 2U)
        {
            detector->history_count++;
        }
    }

    if (output != NULL)
    {
        output->magnitude_g =
            magnitude_g;

        output->filtered_g =
            filtered_g;

        output->baseline_g =
            baseline_for_this_sample;

        output->activity_g =
            activity_for_this_sample;

        output->required_height_g =
            required_height_for_this_sample;

        output->required_prominence_g =
            required_prominence_for_this_sample;

        output->candidate_active =
            detector->candidate_active;

        output->candidate_peak_g =
            detector->candidate_peak_g;

        output->candidate_prominence_g =
            confirmed_prominence_g;

        output->event_detected =
            event_detected;

        output->event_timestamp_ms =
            event_timestamp_ms;
    }

    return event_detected;
}


uint32_t adaptive_peak_detector_get_count(
    const adaptive_peak_detector_t *detector)
{
    if (detector == NULL)
    {
        return 0U;
    }

    return detector->event_count;
}