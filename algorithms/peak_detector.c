#include "peak_detector.h"

#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SOS_SECTIONS 2U
#define SOS_COEFFS 6U

typedef struct
{
    size_t index;
    double height;
    int accepted;
    int visited;
} candidate_peak_t;

static double magnitude_g(const accel_sample_t *sample)
{
    return sqrt(sample->x_g * sample->x_g +
                sample->y_g * sample->y_g +
                sample->z_g * sample->z_g);
}

static int design_butterworth_bandpass_order2(
    double sample_rate_hz,
    double low_hz,
    double high_hz,
    double sos[SOS_SECTIONS][SOS_COEFFS])
{
    const double inv_sqrt2 = 0.70710678118654752440;
    const double warped_low = 4.0 * tan(M_PI * low_hz / sample_rate_hz);
    const double warped_high = 4.0 * tan(M_PI * high_hz / sample_rate_hz);
    const double bandwidth = warped_high - warped_low;
    const double center = sqrt(warped_low * warped_high);
    const double complex prototype[2] = {
        -inv_sqrt2 + inv_sqrt2 * I,
        -inv_sqrt2 - inv_sqrt2 * I
    };
    double complex analog_poles[4];
    double complex digital_poles[4];
    double complex pole_product = 1.0 + 0.0 * I;
    double complex positive_imag[2];
    size_t positive_count = 0U;
    size_t i;
    double gain;
    double complex low_angle_pole;
    double complex high_angle_pole;

    if (!(sample_rate_hz > 0.0) || !(low_hz > 0.0) ||
        !(low_hz < high_hz) || !(high_hz < sample_rate_hz / 2.0))
    {
        return -1;
    }

    for (i = 0U; i < 2U; ++i)
    {
        const double complex scaled = prototype[i] * bandwidth / 2.0;
        const double complex root = csqrt(scaled * scaled - center * center);
        analog_poles[2U * i] = scaled + root;
        analog_poles[2U * i + 1U] = scaled - root;
    }

    for (i = 0U; i < 4U; ++i)
    {
        pole_product *= (4.0 - analog_poles[i]);
        digital_poles[i] = (4.0 + analog_poles[i]) / (4.0 - analog_poles[i]);

        if ((cimag(digital_poles[i]) > 0.0) && (positive_count < 2U))
        {
            positive_imag[positive_count++] = digital_poles[i];
        }
    }

    if (positive_count != 2U)
    {
        return -1;
    }

    gain = bandwidth * bandwidth * creal(16.0 / pole_product);

    if (fabs(carg(positive_imag[0])) < fabs(carg(positive_imag[1])))
    {
        low_angle_pole = positive_imag[0];
        high_angle_pole = positive_imag[1];
    }
    else
    {
        low_angle_pole = positive_imag[1];
        high_angle_pole = positive_imag[0];
    }

    /* Matches SciPy zpk2sos pairing for this 2nd-order band-pass design. */
    sos[0][0] = gain;
    sos[0][1] = 2.0 * gain;
    sos[0][2] = gain;
    sos[0][3] = 1.0;
    sos[0][4] = -2.0 * creal(high_angle_pole);
    sos[0][5] = creal(high_angle_pole * conj(high_angle_pole));

    sos[1][0] = 1.0;
    sos[1][1] = -2.0;
    sos[1][2] = 1.0;
    sos[1][3] = 1.0;
    sos[1][4] = -2.0 * creal(low_angle_pole);
    sos[1][5] = creal(low_angle_pole * conj(low_angle_pole));

    return 0;
}

static void biquad_steady_state_zi(const double row[SOS_COEFFS], double zi[2])
{
    const double b0 = row[0];
    const double b1 = row[1];
    const double b2 = row[2];
    const double a1 = row[4];
    const double a2 = row[5];
    const double B1 = b1 - a1 * b0;
    const double B2 = b2 - a2 * b0;

    zi[0] = (B1 + B2) / (1.0 + a1 + a2);
    zi[1] = (1.0 + a1) * zi[0] - B1;
}

static void sos_steady_state_zi(
    const double sos[SOS_SECTIONS][SOS_COEFFS],
    double zi[SOS_SECTIONS][2])
{
    double scale = 1.0;
    size_t section;

    for (section = 0U; section < SOS_SECTIONS; ++section)
    {
        double section_zi[2];
        const double numerator_sum = sos[section][0] + sos[section][1] + sos[section][2];
        const double denominator_sum = sos[section][3] + sos[section][4] + sos[section][5];

        biquad_steady_state_zi(sos[section], section_zi);
        zi[section][0] = scale * section_zi[0];
        zi[section][1] = scale * section_zi[1];
        scale *= numerator_sum / denominator_sum;
    }
}

static void sos_filter_in_place(
    const double sos[SOS_SECTIONS][SOS_COEFFS],
    const double zi[SOS_SECTIONS][2],
    double initial_scale,
    double *values,
    size_t count)
{
    double state[SOS_SECTIONS][2];
    size_t section;
    size_t i;

    for (section = 0U; section < SOS_SECTIONS; ++section)
    {
        state[section][0] = zi[section][0] * initial_scale;
        state[section][1] = zi[section][1] * initial_scale;
    }

    for (i = 0U; i < count; ++i)
    {
        double x = values[i];

        for (section = 0U; section < SOS_SECTIONS; ++section)
        {
            const double b0 = sos[section][0];
            const double b1 = sos[section][1];
            const double b2 = sos[section][2];
            const double a1 = sos[section][4];
            const double a2 = sos[section][5];
            const double y = b0 * x + state[section][0];

            state[section][0] = b1 * x - a1 * y + state[section][1];
            state[section][1] = b2 * x - a2 * y;
            x = y;
        }

        values[i] = x;
    }
}

static void reverse_array(double *values, size_t count)
{
    size_t i;
    for (i = 0U; i < count / 2U; ++i)
    {
        const double tmp = values[i];
        values[i] = values[count - 1U - i];
        values[count - 1U - i] = tmp;
    }
}

static int sosfiltfilt_like_scipy(
    const double *signal,
    size_t count,
    const double sos[SOS_SECTIONS][SOS_COEFFS],
    double *filtered)
{
    /* SciPy default: 3 * (2 * n_sections + 1) = 15 for these two SOS. */
    const size_t edge = 15U;
    const size_t extended_count = count + 2U * edge;
    double zi[SOS_SECTIONS][2];
    double *extended;
    size_t i;

    if ((signal == NULL) || (filtered == NULL) || (count <= edge))
    {
        return -1;
    }

    extended = (double *)malloc(extended_count * sizeof(*extended));
    if (extended == NULL)
    {
        return -1;
    }

    /* Odd extension, equivalent to scipy.signal.sosfiltfilt(..., padtype="odd"). */
    for (i = 0U; i < edge; ++i)
    {
        extended[i] = 2.0 * signal[0] - signal[edge - i];
    }
    memcpy(extended + edge, signal, count * sizeof(*signal));
    for (i = 0U; i < edge; ++i)
    {
        extended[edge + count + i] = 2.0 * signal[count - 1U] - signal[count - 2U - i];
    }

    sos_steady_state_zi(sos, zi);

    sos_filter_in_place(sos, zi, extended[0], extended, extended_count);
    reverse_array(extended, extended_count);
    sos_filter_in_place(sos, zi, extended[0], extended, extended_count);
    reverse_array(extended, extended_count);

    memcpy(filtered, extended + edge, count * sizeof(*filtered));
    free(extended);
    return 0;
}

static double peak_prominence(const double *signal, size_t count, size_t peak)
{
    const double peak_height = signal[peak];
    double left_min = peak_height;
    double right_min = peak_height;
    size_t i;

    i = peak;
    while (i > 0U)
    {
        --i;
        if (signal[i] > peak_height)
        {
            break;
        }
        if (signal[i] < left_min)
        {
            left_min = signal[i];
        }
    }

    i = peak + 1U;
    while (i < count)
    {
        if (signal[i] > peak_height)
        {
            break;
        }
        if (signal[i] < right_min)
        {
            right_min = signal[i];
        }
        ++i;
    }

    return peak_height - ((left_min > right_min) ? left_min : right_min);
}

static void sort_events_by_index(peak_event_t *events, size_t count)
{
    size_t i;
    for (i = 1U; i < count; ++i)
    {
        const peak_event_t value = events[i];
        size_t j = i;
        while ((j > 0U) && (events[j - 1U].sample_index > value.sample_index))
        {
            events[j] = events[j - 1U];
            --j;
        }
        events[j] = value;
    }
}

void peak_detector_default_config(peak_detector_config_t *config)
{
    if (config == NULL)
    {
        return;
    }

    config->low_hz = 0.9;
    config->high_hz = 8.0;
    config->min_event_interval_s = 0.28;
    config->min_prominence_g = 0.02;
    config->use_min_height = false;
    config->min_height_g = 0.0;
    config->filter_order = 2U;
}

int peak_detector_detect(
    const accel_sample_t *samples,
    size_t sample_count,
    double sample_rate_hz,
    const peak_detector_config_t *config,
    peak_detection_result_t *result)
{
    double sos[SOS_SECTIONS][SOS_COEFFS];
    double *magnitude = NULL;
    candidate_peak_t *candidates = NULL;
    size_t candidate_count = 0U;
    size_t accepted_count = 0U;
    size_t min_distance_samples;
    size_t i;

    if ((samples == NULL) || (config == NULL) || (result == NULL) ||
        (result->events == NULL) || (result->filtered_signal_g == NULL) ||
        (sample_count < 3U) || (result->event_capacity < sample_count) ||
        (result->filtered_capacity < sample_count) ||
        !(sample_rate_hz > 0.0) ||
        !(config->low_hz > 0.0) || !(config->low_hz < config->high_hz) ||
        !(config->high_hz < sample_rate_hz / 2.0) ||
        !(config->min_event_interval_s > 0.0) ||
        !(config->min_prominence_g >= 0.0) ||
        (config->filter_order != 2U))
    {
        return -1;
    }

    result->event_count = 0U;

    magnitude = (double *)malloc(sample_count * sizeof(*magnitude));
    candidates = (candidate_peak_t *)calloc(sample_count, sizeof(*candidates));
    if ((magnitude == NULL) || (candidates == NULL))
    {
        free(magnitude);
        free(candidates);
        return -1;
    }

    for (i = 0U; i < sample_count; ++i)
    {
        magnitude[i] = magnitude_g(&samples[i]);
    }

    if (design_butterworth_bandpass_order2(
            sample_rate_hz, config->low_hz, config->high_hz, sos) != 0 ||
        sosfiltfilt_like_scipy(
            magnitude, sample_count, sos, result->filtered_signal_g) != 0)
    {
        free(magnitude);
        free(candidates);
        return -1;
    }

    /* Local maxima, including flat plateaus using the midpoint like SciPy. */
    i = 1U;
    while (i + 1U < sample_count)
    {
        if (result->filtered_signal_g[i - 1U] < result->filtered_signal_g[i])
        {
            size_t right = i;
            while ((right + 1U < sample_count) &&
                   (result->filtered_signal_g[right + 1U] == result->filtered_signal_g[i]))
            {
                ++right;
            }

            if ((right + 1U < sample_count) &&
                (result->filtered_signal_g[right] > result->filtered_signal_g[right + 1U]))
            {
                const size_t peak = (i + right) / 2U;
                const double height = result->filtered_signal_g[peak];

                if (!config->use_min_height || (height >= config->min_height_g))
                {
                    candidates[candidate_count].index = peak;
                    candidates[candidate_count].height = height;
                    ++candidate_count;
                }
            }
            i = right + 1U;
        }
        else
        {
            ++i;
        }
    }

    /* Python round() uses round-to-nearest-even; nearbyint matches that default. */
    min_distance_samples = (size_t)nearbyint(config->min_event_interval_s * sample_rate_hz);
    if (min_distance_samples < 1U)
    {
        min_distance_samples = 1U;
    }

    /* SciPy evaluates distance before prominence and gives higher peaks priority. */
    for (i = 0U; i < candidate_count; ++i)
    {
        size_t j;
        size_t best = (size_t)-1;
        double best_height = -INFINITY;
        int keep = 1;

        for (j = 0U; j < candidate_count; ++j)
        {
            if (!candidates[j].visited && candidates[j].height > best_height)
            {
                best = j;
                best_height = candidates[j].height;
            }
        }

        if (best == (size_t)-1)
        {
            break;
        }

        candidates[best].visited = 1;

        for (j = 0U; j < candidate_count; ++j)
        {
            if (candidates[j].accepted)
            {
                const size_t a = candidates[best].index;
                const size_t b = candidates[j].index;
                const size_t distance = (a > b) ? (a - b) : (b - a);
                if (distance < min_distance_samples)
                {
                    keep = 0;
                    break;
                }
            }
        }

        if (keep)
        {
            candidates[best].accepted = 1;
        }
    }

    for (i = 0U; i < candidate_count; ++i)
    {
        if (candidates[i].accepted)
        {
            const size_t peak = candidates[i].index;
            const double prominence = peak_prominence(
                result->filtered_signal_g, sample_count, peak);

            if (prominence >= config->min_prominence_g)
            {
                result->events[accepted_count].sample_index = peak;
                result->events[accepted_count].filtered_peak_g = result->filtered_signal_g[peak];
                result->events[accepted_count].prominence_g = prominence;
                ++accepted_count;
            }
        }
    }

    sort_events_by_index(result->events, accepted_count);
    result->event_count = accepted_count;

    free(magnitude);
    free(candidates);
    return 0;
}
