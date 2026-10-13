#include "algorithm_runner.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "peak_detector.h"
#include "threshold_detector.h"
#include "adaptive_threshold_detector.h"
#include "adaptive_peak_detector.h"
#include "frequency_detector.h"


static int write_peak_debug_csv(
    const char *path,
    const accel_sample_t *samples,
    size_t count,
    const peak_detection_result_t *result)
{
    FILE *file;
    size_t i;
    size_t event_cursor = 0U;

    if (path == NULL)
    {
        return 0;
    }

    file = fopen(path, "w");

    if (file == NULL)
    {
        return -1;
    }

    fprintf(
        file,
        "sample_index,"
        "t_ms,"
        "phone_time_est_ms,"
        "magnitude_g,"
        "filtered_signal_g,"
        "is_peak,"
        "prominence_g\n");

    for (i = 0U; i < count; ++i)
    {
        const double magnitude =
            sqrt(
                samples[i].x_g * samples[i].x_g +
                samples[i].y_g * samples[i].y_g +
                samples[i].z_g * samples[i].z_g);

        int is_peak = 0;
        double prominence = 0.0;

        if ((event_cursor < result->event_count) &&
            (result->events[event_cursor].sample_index == i))
        {
            is_peak = 1;

            prominence =
                result->events[event_cursor].prominence_g;

            ++event_cursor;
        }

        fprintf(
            file,
            "%zu,%u,%.3f,%.9f,%.9f,%d,%.9f\n",
            i,
            samples[i].timestamp_ms,
            samples[i].phone_time_est_ms,
            magnitude,
            result->filtered_signal_g[i],
            is_peak,
            prominence);
    }

    fclose(file);

    return 0;
}


static int run_reference_peak_detector(
    const accel_sample_t *samples,
    size_t sample_count,
    double sample_rate_hz,
    const char *debug_output_path)
{
    peak_detector_config_t config;
    peak_detection_result_t result;
    size_t i;

    result.events =
        (peak_event_t *)calloc(
            sample_count,
            sizeof(*result.events));

    result.filtered_signal_g =
        (double *)calloc(
            sample_count,
            sizeof(*result.filtered_signal_g));

    result.event_capacity = sample_count;
    result.filtered_capacity = sample_count;
    result.event_count = 0U;

    if ((result.events == NULL) ||
        (result.filtered_signal_g == NULL))
    {
        fprintf(
            stderr,
            "Could not allocate peak-detector output.\n");

        free(result.events);
        free(result.filtered_signal_g);

        return -1;
    }

    peak_detector_default_config(&config);

    if (peak_detector_detect(
            samples,
            sample_count,
            sample_rate_hz,
            &config,
            &result) != 0)
    {
        fprintf(
            stderr,
            "Reference peak detector failed.\n");

        free(result.events);
        free(result.filtered_signal_g);

        return -1;
    }

    printf("\nReference peak detector\n");
    printf("-----------------------\n");

    printf(
        "Band-pass: %.3f - %.3f Hz\n",
        config.low_hz,
        config.high_hz);

    printf(
        "Minimum interval: %.3f s\n",
        config.min_event_interval_s);

    printf(
        "Minimum prominence: %.3f g\n",
        config.min_prominence_g);

    printf(
        "Detected peaks: %zu\n\n",
        result.event_count);

    printf(
        "#\tindex\tt_ms\tfiltered_g\tprominence_g\n");

    for (i = 0U; i < result.event_count; ++i)
    {
        const peak_event_t *event =
            &result.events[i];

        printf(
            "%zu\t%zu\t%u\t%.9f\t%.9f\n",
            i + 1U,
            event->sample_index,
            samples[event->sample_index].timestamp_ms,
            event->filtered_peak_g,
            event->prominence_g);
    }

    if (debug_output_path != NULL)
    {
        if (write_peak_debug_csv(
                debug_output_path,
                samples,
                sample_count,
                &result) != 0)
        {
            fprintf(
                stderr,
                "Could not write debug CSV '%s'.\n",
                debug_output_path);

            free(result.events);
            free(result.filtered_signal_g);

            return -1;
        }

        printf(
            "\nPeak debug CSV: %s\n",
            debug_output_path);
    }

    free(result.events);
    free(result.filtered_signal_g);

    return 0;
}


static void run_fixed_threshold_detector(
    const accel_sample_t *samples,
    size_t sample_count)
{
    threshold_detector_config_t config;
    threshold_detector_t detector;
    size_t i;

    threshold_detector_default_config(&config);

    threshold_detector_init(
        &detector,
        &config);

    printf("\nFixed threshold detector\n");
    printf("------------------------\n");

    printf(
        "Threshold: %.3f g\n",
        config.threshold_g);

    printf(
        "Minimum interval: %u ms\n\n",
        config.min_event_interval_ms);

    printf(
        "#\tindex\tt_ms\tmagnitude_g\n");

    for (i = 0U; i < sample_count; ++i)
    {
        double magnitude_g;

        if (threshold_detector_process(
                &detector,
                &samples[i],
                &magnitude_g))
        {
            printf(
                "%u\t%zu\t%u\t%.6f\n",
                threshold_detector_get_count(&detector),
                i,
                samples[i].timestamp_ms,
                magnitude_g);
        }
    }

    printf(
        "\nDetected threshold events: %u\n",
        threshold_detector_get_count(&detector));
}


static void run_adaptive_threshold_detector(
    const accel_sample_t *samples,
    size_t sample_count)
{
    adaptive_threshold_config_t config;
    adaptive_threshold_detector_t detector;
    size_t i;

    adaptive_threshold_detector_default_config(
        &config);

    adaptive_threshold_detector_init(
        &detector,
        &config);

    printf("\nAdaptive threshold detector\n");
    printf("---------------------------\n");

    printf(
        "Baseline alpha: %.4f\n",
        config.baseline_alpha);

    printf(
        "Activity alpha: %.4f\n",
        config.activity_alpha);

    printf(
        "Sensitivity: %.3f\n",
        config.sensitivity);

    printf(
        "Threshold bounds: %.3f - %.3f g\n",
        config.min_threshold_g,
        config.max_threshold_g);

    printf(
        "Minimum interval: %u ms\n\n",
        config.min_event_interval_ms);

    printf(
        "#\tindex\tt_ms\t"
        "magnitude_g\tthreshold_g\t"
        "baseline_g\tactivity_g\n");

    for (i = 0U; i < sample_count; ++i)
    {
        adaptive_threshold_output_t output;

        if (adaptive_threshold_detector_process(
                &detector,
                &samples[i],
                &output))
        {
            printf(
                "%u\t%zu\t%u\t"
                "%.6f\t%.6f\t%.6f\t%.6f\n",

                adaptive_threshold_detector_get_count(
                    &detector),

                i,
                samples[i].timestamp_ms,

                output.magnitude_g,
                output.threshold_g,
                output.baseline_g,
                output.activity_g);
        }
    }

    printf(
        "\nDetected adaptive-threshold events: %u\n",
        adaptive_threshold_detector_get_count(
            &detector));
}


static int run_adaptive_peak_detector(
    const accel_sample_t *samples,
    size_t sample_count,
    double sample_rate_hz)
{
    adaptive_peak_config_t config;
    adaptive_peak_detector_t detector;
    size_t i;

    adaptive_peak_detector_default_config(
        &config);

    /*
     * Use the measured sample rate rather than
     * assuming exactly 50 Hz.
     */
    config.sample_rate_hz =
        sample_rate_hz;

    if (adaptive_peak_detector_init(
            &detector,
            &config) != 0)
    {
        fprintf(
            stderr,
            "Adaptive peak detector initialization failed.\n");

        return -1;
    }

    printf("\nAdaptive peak detector\n");
    printf("----------------------\n");

    printf(
        "Sample rate: %.3f Hz\n",
        config.sample_rate_hz);

    printf(
        "Band-pass: %.3f - %.3f Hz\n",
        config.low_cutoff_hz,
        config.high_cutoff_hz);

    printf(
        "Baseline alpha: %.4f\n",
        config.baseline_alpha);

    printf(
        "Activity alpha: %.4f\n",
        config.activity_alpha);

    printf(
        "Height multiplier: %.3f\n",
        config.height_multiplier);

    printf(
        "Prominence multiplier: %.3f\n",
        config.prominence_multiplier);

    printf(
        "Confirmation window: %u ms\n",
        config.confirmation_window_ms);

    printf(
        "Minimum interval: %u ms\n\n",
        config.min_event_interval_ms);

    printf(
        "#\tconfirm_index\tpeak_t_ms\t"
        "peak_g\tprominence_g\n");

    for (i = 0U; i < sample_count; ++i)
    {
        adaptive_peak_output_t output;

        if (adaptive_peak_detector_process(
                &detector,
                &samples[i],
                &output))
        {
            printf(
                "%u\t%zu\t%u\t%.6f\t%.6f\n",

                adaptive_peak_detector_get_count(
                    &detector),

                i,

                output.event_timestamp_ms,
                output.candidate_peak_g,
                output.candidate_prominence_g);
        }
    }

    printf(
        "\nDetected adaptive-peak events: %u\n",
        adaptive_peak_detector_get_count(
            &detector));

    return 0;
}

static int run_frequency_analysis(
    const accel_sample_t *samples,
    size_t sample_count,
    double sample_rate_hz)
{
    frequency_detector_config_t config;
    frequency_detector_t detector;

    size_t i;

    frequency_detector_default_config(
        &config);

    /*
     * Use the actual sample rate measured
     * from the recording.
     */
    config.sample_rate_hz =
        sample_rate_hz;

    if (frequency_detector_init(
            &detector,
            &config) != 0)
    {
        fprintf(
            stderr,
            "Frequency detector initialization failed.\n");

        return -1;
    }

    printf("\nFrequency analysis\n");
    printf("------------------\n");

    printf(
        "FFT size: %u samples\n",
        FREQUENCY_FFT_SIZE);

    printf(
        "Hop size: %zu samples\n",
        config.hop_size);

    printf(
        "Sample rate: %.3f Hz\n",
        config.sample_rate_hz);

    printf(
        "Frequency resolution: %.3f Hz/bin\n",
        config.sample_rate_hz /
        (double)FREQUENCY_FFT_SIZE);

    printf(
        "Search range: %.3f - %.3f Hz\n\n",
        config.search_low_hz,
        config.search_high_hz);

    printf(
        "#\tstart_ms\tend_ms\t"
        "mean_g\t"
        "dominant_bin\t"
        "dominant_hz\t"
        "dominant_power\n");

    for (i = 0U;
         i < sample_count;
         ++i)
    {
        frequency_spectrum_output_t spectrum_output;

        /*
         * Stage 1:
         * Magnitude -> mean removal -> Hann -> FFT.
         */
        if (frequency_detector_process(
                &detector,
                &samples[i],
                &spectrum_output))
        {
            frequency_dominant_result_t dominant;

            /*
             * Stage 2:
             * Find strongest component in the
             * configured walking-frequency band.
             */
            if (frequency_detector_find_dominant(
                    &detector,
                    &dominant) != 0)
            {
                fprintf(
                    stderr,
                    "Dominant-frequency analysis failed "
                    "for window %u.\n",
                    spectrum_output.window_number);

                return -1;
            }

            printf(
                "%u\t%u\t%u\t"
                "%.6f\t"
                "%zu\t"
                "%.6f\t"
                "%.9f\n",

                spectrum_output.window_number,

                spectrum_output.window_start_ms,
                spectrum_output.window_end_ms,

                spectrum_output.mean_magnitude_g,

                dominant.dominant_bin,
                dominant.dominant_frequency_hz,
                dominant.dominant_power);
        }
    }

    printf(
        "\nGenerated frequency windows: %u\n",
        detector.window_count);

    return 0;
}


int run_all_algorithms(
    const accel_sample_t *samples,
    size_t sample_count,
    double sample_rate_hz,
    const char *debug_output_path)
{
    if ((samples == NULL) ||
        (sample_count == 0U) ||
        !(sample_rate_hz > 0.0))
    {
        return -1;
    }

    printf("\n========================================\n");
    printf("STEP DETECTION TEST\n");
    printf("========================================\n");

    printf(
        "Samples: %zu\n",
        sample_count);

    printf(
        "Sample rate: %.3f Hz\n",
        sample_rate_hz);

    /*
     * Heuristic algorithms.
     */
    run_fixed_threshold_detector(
        samples,
        sample_count);

    if (run_reference_peak_detector(
            samples,
            sample_count,
            sample_rate_hz,
            debug_output_path) != 0)
    {
        return -1;
    }

    /*
     * Adaptive algorithms.
     */
    run_adaptive_threshold_detector(
        samples,
        sample_count);

    if (run_adaptive_peak_detector(
            samples,
            sample_count,
            sample_rate_hz) != 0)
    {
        return -1;
    }

    /*
        Signal Based Algorithms
    */

    if (run_frequency_detector_analysis(
        samples,
        sample_count,
        sample_rate_hz) != 0)
    {
        return -1;
    }

    printf("\n========================================\n");

    return 0;
}