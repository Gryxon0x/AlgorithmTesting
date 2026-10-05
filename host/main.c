#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "csv_reader.h"
#include "peak_detector.h"
#include "threshold_detector.h"
#include "adaptive_threshold_detector.h"

static int compare_double(const void *left, const void *right)
{
    const double a = *(const double *)left;
    const double b = *(const double *)right;
    return (a > b) - (a < b);
}

static double estimate_sample_rate_hz(const accel_sample_t *samples, size_t count)
{
    double *dt;
    double median_dt_ms;
    size_t i;

    if (count < 2U)
    {
        return 0.0;
    }

    dt = (double *)malloc((count - 1U) * sizeof(*dt));
    if (dt == NULL)
    {
        return 0.0;
    }

    for (i = 0U; i + 1U < count; ++i)
    {
        dt[i] = samples[i + 1U].phone_time_est_ms - samples[i].phone_time_est_ms;
        if (!(dt[i] > 0.0))
        {
            free(dt);
            return 0.0;
        }
    }

    qsort(dt, count - 1U, sizeof(*dt), compare_double);

    if (((count - 1U) % 2U) != 0U)
    {
        median_dt_ms = dt[(count - 1U) / 2U];
    }
    else
    {
        const size_t upper = (count - 1U) / 2U;
        median_dt_ms = 0.5 * (dt[upper - 1U] + dt[upper]);
    }

    free(dt);
    return 1000.0 / median_dt_ms;
}

static int write_debug_csv(
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

    fprintf(file,
            "sample_index,t_ms,phone_time_est_ms,magnitude_g,filtered_signal_g,is_peak,prominence_g\n");

    for (i = 0U; i < count; ++i)
    {
        const double magnitude = sqrt(samples[i].x_g * samples[i].x_g +
                                      samples[i].y_g * samples[i].y_g +
                                      samples[i].z_g * samples[i].z_g);
        int is_peak = 0;
        double prominence = 0.0;

        if ((event_cursor < result->event_count) &&
            (result->events[event_cursor].sample_index == i))
        {
            is_peak = 1;
            prominence = result->events[event_cursor].prominence_g;
            ++event_cursor;
        }

        fprintf(file, "%zu,%u,%.3f,%.9f,%.9f,%d,%.9f\n",
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

int main(int argc, char *argv[])
{
    const char *csv_path;
    const char *debug_output_path = NULL;
    csv_reader_t reader;
    csv_read_status_t status;
    accel_sample_t *samples = NULL;
    size_t sample_count = 0U;
    size_t sample_capacity = 0U;
    double sample_rate_hz;
    peak_detector_config_t config;
    peak_detection_result_t result;
    size_t i;

    if ((argc < 2) || (argc > 3))
    {
        fprintf(stderr, "Usage: %s <csv-file> [debug-output.csv]\n", argv[0]);
        fprintf(stderr, "Example: %s data/example_session.csv output/peak_debug.csv\n", argv[0]);
        return EXIT_FAILURE;
    }

    csv_path = argv[1];
    if (argc == 3)
    {
        debug_output_path = argv[2];
    }

    if (csv_reader_open(&reader, csv_path) != 0)
    {
        fprintf(stderr, "CSV error: %s\n", csv_reader_error(&reader));
        return EXIT_FAILURE;
    }

    while (1)
    {
        accel_sample_t sample;
        status = csv_reader_next(&reader, &sample);

        if (status == CSV_READ_EOF)
        {
            break;
        }
        if (status == CSV_READ_ERROR)
        {
            fprintf(stderr, "CSV error: %s\n", csv_reader_error(&reader));
            csv_reader_close(&reader);
            free(samples);
            return EXIT_FAILURE;
        }

        if (sample_count == sample_capacity)
        {
            const size_t new_capacity = (sample_capacity == 0U) ? 512U : sample_capacity * 2U;
            accel_sample_t *grown = (accel_sample_t *)realloc(
                samples, new_capacity * sizeof(*samples));
            if (grown == NULL)
            {
                fprintf(stderr, "Out of memory while loading samples.\n");
                csv_reader_close(&reader);
                free(samples);
                return EXIT_FAILURE;
            }
            samples = grown;
            sample_capacity = new_capacity;
        }

        samples[sample_count++] = sample;
    }

    csv_reader_close(&reader);

    if (sample_count < 16U)
    {
        fprintf(stderr, "Recording is too short for the zero-phase filter.\n");
        free(samples);
        return EXIT_FAILURE;
    }

    sample_rate_hz = estimate_sample_rate_hz(samples, sample_count);
    if (!(sample_rate_hz > 0.0))
    {
        fprintf(stderr, "Could not determine a valid sample rate.\n");
        free(samples);
        return EXIT_FAILURE;
    }

    result.events = (peak_event_t *)calloc(sample_count, sizeof(*result.events));
    result.filtered_signal_g = (double *)calloc(sample_count, sizeof(*result.filtered_signal_g));
    result.event_capacity = sample_count;
    result.filtered_capacity = sample_count;
    result.event_count = 0U;

    if ((result.events == NULL) || (result.filtered_signal_g == NULL))
    {
        fprintf(stderr, "Out of memory while allocating detector output.\n");
        free(result.events);
        free(result.filtered_signal_g);
        free(samples);
        return EXIT_FAILURE;
    }

    peak_detector_default_config(&config);

    if (peak_detector_detect(samples, sample_count, sample_rate_hz, &config, &result) != 0)
    {
        fprintf(stderr, "Peak detector failed.\n");
        free(result.events);
        free(result.filtered_signal_g);
        free(samples);
        return EXIT_FAILURE;
    }

    threshold_detector_config_t threshold_config;
    threshold_detector_t threshold_detector;

    threshold_detector_default_config(&threshold_config);

    threshold_detector_init(
        &threshold_detector,
        &threshold_config);

    printf("\nThreshold detector\n");
    printf("------------------\n");
    printf(
        "Threshold: %.3f g\n",
        threshold_config.threshold_g);

    printf(
        "Minimum interval: %u ms\n\n",
        threshold_config.min_event_interval_ms);

    printf("#\tindex\tt_ms\tmagnitude_g\n");

    for (i = 0U; i < sample_count; ++i)
    {
        double magnitude_g;

        if (threshold_detector_process(
                &threshold_detector,
                &samples[i],
                &magnitude_g))
        {
            printf(
                "%u\t%zu\t%u\t%.6f\n",
                threshold_detector_get_count(
                    &threshold_detector),
                i,
                samples[i].timestamp_ms,
                magnitude_g);
        }
    }

    printf(
        "\nDetected threshold events: %u\n",
        threshold_detector_get_count(
            &threshold_detector));

    adaptive_threshold_config_t adaptive_config;
    adaptive_threshold_detector_t adaptive_detector;

    adaptive_threshold_detector_default_config(
        &adaptive_config);

    adaptive_threshold_detector_init(
        &adaptive_detector,
        &adaptive_config);

    printf("\nAdaptive threshold detector\n");
    printf("---------------------------\n");

    printf(
        "Baseline alpha: %.4f\n",
        adaptive_config.baseline_alpha);

    printf(
        "Activity alpha: %.4f\n",
        adaptive_config.activity_alpha);

    printf(
        "Sensitivity: %.3f\n",
        adaptive_config.sensitivity);

    printf(
        "Threshold bounds: %.3f - %.3f g\n",
        adaptive_config.min_threshold_g,
        adaptive_config.max_threshold_g);

    printf(
        "Minimum interval: %u ms\n\n",
        adaptive_config.min_event_interval_ms);

    printf(
        "#\tindex\tt_ms\tmagnitude_g\t"
        "threshold_g\tbaseline_g\tactivity_g\n");

    for (i = 0U; i < sample_count; ++i)
    {
        adaptive_threshold_output_t output;

        if (adaptive_threshold_detector_process(
                &adaptive_detector,
                &samples[i],
                &output))
        {
            printf(
                "%u\t%zu\t%u\t"
                "%.6f\t%.6f\t%.6f\t%.6f\n",

                adaptive_threshold_detector_get_count(
                    &adaptive_detector),

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
            &adaptive_detector));

    printf("File: %s\n", csv_path);
    printf("Samples: %zu\n", sample_count);
    printf("Sample rate: %.3f Hz\n", sample_rate_hz);
    printf("Config: %.3f-%.3f Hz, min interval %.3f s, min prominence %.3f g\n",
           config.low_hz, config.high_hz,
           config.min_event_interval_s, config.min_prominence_g);
    printf("Detected peaks: %zu\n\n", result.event_count);
    printf("#\tindex\tt_ms\tfiltered_g\tprominence_g\n");

    for (i = 0U; i < result.event_count; ++i)
    {
        const peak_event_t *event = &result.events[i];
        printf("%zu\t%zu\t%u\t%.9f\t%.9f\n",
               i + 1U,
               event->sample_index,
               samples[event->sample_index].timestamp_ms,
               event->filtered_peak_g,
               event->prominence_g);
    }

    if (debug_output_path != NULL)
    {
        if (write_debug_csv(debug_output_path, samples, sample_count, &result) != 0)
        {
            fprintf(stderr, "Could not write debug CSV '%s'.\n", debug_output_path);
            free(result.events);
            free(result.filtered_signal_g);
            free(samples);
            return EXIT_FAILURE;
        }
        printf("\nDebug CSV: %s\n", debug_output_path);
    }

    free(result.events);
    free(result.filtered_signal_g);
    free(samples);
    return EXIT_SUCCESS;
}
