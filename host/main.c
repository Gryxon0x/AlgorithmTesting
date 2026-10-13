#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "csv_reader.h"
#include "algorithm_runner.h"

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

    printf("File: %s\n", csv_path);

    if (run_all_algorithms(
        samples,
        sample_count,
        sample_rate_hz,
        debug_output_path) != 0)
    {
        fprintf(
            stderr,
            "One or more algorithms failed.\n");

        free(samples);

        return EXIT_FAILURE;
    }

    free(samples);

    return EXIT_SUCCESS;
}
