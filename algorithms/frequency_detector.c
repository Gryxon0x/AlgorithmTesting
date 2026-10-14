#include "frequency_detector.h"

#include <math.h>
#include <stddef.h>
#include <string.h>


#define TWO_PI 6.28318530717958647692


typedef struct
{
    double real;
    double imag;

} complex_sample_t;


static double calculate_magnitude_g(
    const accel_sample_t *sample)
{
    return sqrt(
        sample->x_g * sample->x_g +
        sample->y_g * sample->y_g +
        sample->z_g * sample->z_g);
}


static void fft_128(
    complex_sample_t *data)
{
    size_t i;
    size_t j = 0U;
    size_t length;

    /*
     * Bit-reversal permutation.
     */
    for (i = 1U; i < FREQUENCY_FFT_SIZE; ++i)
    {
        size_t bit =
            FREQUENCY_FFT_SIZE >> 1U;

        while ((j & bit) != 0U)
        {
            j ^= bit;
            bit >>= 1U;
        }

        j ^= bit;

        if (i < j)
        {
            complex_sample_t temporary =
                data[i];

            data[i] = data[j];
            data[j] = temporary;
        }
    }

    /*
     * Iterative radix-2 Cooley-Tukey FFT.
     */
    for (length = 2U;
         length <= FREQUENCY_FFT_SIZE;
         length <<= 1U)
    {
        const double angle =
            -TWO_PI / (double)length;

        const double step_real =
            cos(angle);

        const double step_imag =
            sin(angle);

        size_t block_start;

        for (block_start = 0U;
             block_start < FREQUENCY_FFT_SIZE;
             block_start += length)
        {
            double w_real = 1.0;
            double w_imag = 0.0;

            const size_t half_length =
                length / 2U;

            size_t k;

            for (k = 0U;
                 k < half_length;
                 ++k)
            {
                const size_t even_index =
                    block_start + k;

                const size_t odd_index =
                    even_index + half_length;

                const complex_sample_t even =
                    data[even_index];

                const complex_sample_t odd =
                    data[odd_index];

                const double odd_real =
                    odd.real * w_real -
                    odd.imag * w_imag;

                const double odd_imag =
                    odd.real * w_imag +
                    odd.imag * w_real;

                data[even_index].real =
                    even.real + odd_real;

                data[even_index].imag =
                    even.imag + odd_imag;

                data[odd_index].real =
                    even.real - odd_real;

                data[odd_index].imag =
                    even.imag - odd_imag;

                {
                    const double next_w_real =
                        w_real * step_real -
                        w_imag * step_imag;

                    const double next_w_imag =
                        w_real * step_imag +
                        w_imag * step_real;

                    w_real = next_w_real;
                    w_imag = next_w_imag;
                }
            }
        }
    }
}


static void calculate_spectrum(
    frequency_detector_t *detector,
    frequency_spectrum_output_t *output)
{
    complex_sample_t fft_data[FREQUENCY_FFT_SIZE];

    double mean_magnitude = 0.0;

    size_t n;
    size_t bin;

    /*
     * Mean does not depend on circular-buffer order.
     */
    for (n = 0U;
         n < FREQUENCY_FFT_SIZE;
         ++n)
    {
        mean_magnitude +=
            detector->magnitude_buffer[n];
    }

    mean_magnitude /=
        (double)FREQUENCY_FFT_SIZE;

    /*
     * detector->write_index points to the oldest
     * sample once the circular buffer is full.
     *
     * Copy samples into chronological order,
     * remove DC, then apply the Hann window.
     */
    for (n = 0U;
         n < FREQUENCY_FFT_SIZE;
         ++n)
    {
        const size_t buffer_index =
            (detector->write_index + n) %
            FREQUENCY_FFT_SIZE;

        const double centered =
            detector->magnitude_buffer[buffer_index] -
            mean_magnitude;

        const double hann =
            0.5 *
            (
                1.0 -
                cos(
                    TWO_PI *
                    (double)n /
                    (double)(FREQUENCY_FFT_SIZE - 1U)
                )
            );

        fft_data[n].real =
            centered * hann;

        fft_data[n].imag = 0.0;
    }

    fft_128(fft_data);

    /*
     * Convert complex FFT values into a one-sided
     * power spectrum.
     */
    for (bin = 0U;
         bin < FREQUENCY_SPECTRUM_BINS;
         ++bin)
    {
        double power =
            fft_data[bin].real *
            fft_data[bin].real +
            fft_data[bin].imag *
            fft_data[bin].imag;

        power /=
            (double)(
                FREQUENCY_FFT_SIZE *
                FREQUENCY_FFT_SIZE);

        /*
         * Double non-DC/non-Nyquist bins because
         * we're keeping only the positive half of
         * the spectrum.
         */
        if ((bin != 0U) &&
            (bin != FREQUENCY_FFT_SIZE / 2U))
        {
            power *= 2.0;
        }

        detector->power_spectrum[bin] =
            power;
    }

    detector->window_count++;

    if (output != NULL)
    {
        const size_t oldest_index =
            detector->write_index;

        const size_t newest_index =
            (
                detector->write_index +
                FREQUENCY_FFT_SIZE -
                1U
            ) %
            FREQUENCY_FFT_SIZE;

        output->window_number =
            detector->window_count;

        output->window_start_ms =
            detector->timestamp_buffer[
                oldest_index];

        output->window_end_ms =
            detector->timestamp_buffer[
                newest_index];

        output->mean_magnitude_g =
            mean_magnitude;
    }
}


void frequency_detector_default_config(
    frequency_detector_config_t *config)
{
    if (config == NULL)
    {
        return;
    }

    config->sample_rate_hz = 50.0;

    /*
     * 75% overlap:
     *
     * FFT size = 128
     * hop      = 32
     */
    config->hop_size = 32U;

    /*
     * Initial walking-related search range.
     */
    config->search_low_hz = 0.5;
    config->search_high_hz = 4.0;

    config->quality_neighbor_bins = 1U;

    config->min_spectral_quality = 0.55;
    config->min_band_power = 0.000001;

    config->step_frequency_multiplier = 1.0;
}


int frequency_detector_init(
    frequency_detector_t *detector,
    const frequency_detector_config_t *config)
{
    if ((detector == NULL) ||
        (config == NULL))
    {
        return -1;
    }

    if (!(config->sample_rate_hz > 0.0))
    {
        return -1;
    }

    if ((config->hop_size == 0U) ||
        (config->hop_size >
         FREQUENCY_FFT_SIZE))
    {
        return -1;
    }

    if (!(config->search_low_hz >= 0.0))
    {
        return -1;
    }

    if (!(config->search_high_hz >
        config->search_low_hz))
    {
        return -1;
    }

    if (!(config->search_high_hz <=
        config->sample_rate_hz / 2.0))
    {
        return -1;
    }

    if ((config->min_spectral_quality < 0.0) ||
    (config->min_spectral_quality > 1.0))
    {
        return -1;
    }

    if (config->min_band_power < 0.0)
    {
        return -1;
    }

    if (!(config->step_frequency_multiplier > 0.0))
    {
        return -1;
    }

    memset(
        detector,
        0,
        sizeof(*detector));

    detector->config = *config;

    return 0;
}


bool frequency_detector_process(
    frequency_detector_t *detector,
    const accel_sample_t *sample,
    frequency_spectrum_output_t *output)
{
    double magnitude_g;

    if ((detector == NULL) ||
        (sample == NULL))
    {
        return false;
    }

    magnitude_g =
        calculate_magnitude_g(sample);

    detector->magnitude_buffer[
        detector->write_index] =
        magnitude_g;

    detector->timestamp_buffer[
        detector->write_index] =
        sample->timestamp_ms;

    detector->write_index =
        (
            detector->write_index + 1U
        ) %
        FREQUENCY_FFT_SIZE;

    /*
     * Initial window fill.
     */
    if (detector->sample_count <
        FREQUENCY_FFT_SIZE)
    {
        detector->sample_count++;

        /*
         * First spectrum as soon as all 128
         * samples have arrived.
         */
        if (detector->sample_count ==
            FREQUENCY_FFT_SIZE)
        {
            detector->samples_since_spectrum = 0U;

            calculate_spectrum(
                detector,
                output);

            return true;
        }

        return false;
    }

    /*
     * Buffer already full.
     *
     * Calculate another FFT once hop_size new
     * samples have arrived.
     */
    detector->samples_since_spectrum++;

    if (detector->samples_since_spectrum >=
        detector->config.hop_size)
    {
        detector->samples_since_spectrum = 0U;

        calculate_spectrum(
            detector,
            output);

        return true;
    }

    return false;
}


const double *frequency_detector_get_power_spectrum(
    const frequency_detector_t *detector)
{
    if (detector == NULL)
    {
        return NULL;
    }

    return detector->power_spectrum;
}


double frequency_detector_get_bin_frequency_hz(
    const frequency_detector_t *detector,
    size_t bin)
{
    if ((detector == NULL) ||
        (bin >= FREQUENCY_SPECTRUM_BINS))
    {
        return -1.0;
    }

    return
        (double)bin *
        detector->config.sample_rate_hz /
        (double)FREQUENCY_FFT_SIZE;
}

int frequency_detector_find_dominant(
    const frequency_detector_t *detector,
    frequency_dominant_result_t *result)
{
    size_t bin;

    bool found = false;

    size_t best_bin = 0U;
    double best_power = 0.0;

    if ((detector == NULL) ||
        (result == NULL))
    {
        return -1;
    }

    /*
     * Search only the configured frequency band.
     */
    for (bin = 0U;
         bin < FREQUENCY_SPECTRUM_BINS;
         ++bin)
    {
        const double frequency_hz =
            frequency_detector_get_bin_frequency_hz(
                detector,
                bin);

        const double power =
            detector->power_spectrum[bin];

        if (frequency_hz <
            detector->config.search_low_hz)
        {
            continue;
        }

        if (frequency_hz >
            detector->config.search_high_hz)
        {
            break;
        }

        if ((!found) ||
            (power > best_power))
        {
            best_bin = bin;
            best_power = power;

            found = true;
        }
    }

    if (!found)
    {
        return -1;
    }

    result->dominant_bin =
        best_bin;

    result->dominant_frequency_hz =
        frequency_detector_get_bin_frequency_hz(
            detector,
            best_bin);

    result->dominant_power =
        best_power;

    return 0;
}

int frequency_detector_evaluate_quality(
    const frequency_detector_t *detector,
    const frequency_dominant_result_t *dominant,
    frequency_quality_result_t *result)
{
    size_t bin;

    size_t region_start;
    size_t region_end;

    double total_band_power = 0.0;
    double dominant_region_power = 0.0;

    if ((detector == NULL) ||
        (dominant == NULL) ||
        (result == NULL))
    {
        return -1;
    }

    if (dominant->dominant_bin >=
        FREQUENCY_SPECTRUM_BINS)
    {
        return -1;
    }

    /*
     * Determine dominant-frequency neighborhood.
     */
    if (dominant->dominant_bin >
        detector->config.quality_neighbor_bins)
    {
        region_start =
            dominant->dominant_bin -
            detector->config.quality_neighbor_bins;
    }
    else
    {
        region_start = 0U;
    }

    region_end =
        dominant->dominant_bin +
        detector->config.quality_neighbor_bins;

    if (region_end >=
        FREQUENCY_SPECTRUM_BINS)
    {
        region_end =
            FREQUENCY_SPECTRUM_BINS - 1U;
    }

    /*
     * Sum all power inside the configured
     * walking-frequency search band.
     */
    for (bin = 0U;
         bin < FREQUENCY_SPECTRUM_BINS;
         ++bin)
    {
        const double frequency_hz =
            frequency_detector_get_bin_frequency_hz(
                detector,
                bin);

        const double power =
            detector->power_spectrum[bin];

        /*
         * Ignore frequencies below the
         * analysis band.
         */
        if (frequency_hz <
            detector->config.search_low_hz)
        {
            continue;
        }

        /*
         * Frequencies increase with bin number,
         * so we can stop once we pass the band.
         */
        if (frequency_hz >
            detector->config.search_high_hz)
        {
            break;
        }

        total_band_power += power;

        /*
         * Only include dominant-region power
         * if the bin is also inside the configured
         * search band.
         */
        if ((bin >= region_start) &&
            (bin <= region_end))
        {
            dominant_region_power +=
                power;
        }
    }

    result->total_band_power =
        total_band_power;

    result->dominant_region_power =
        dominant_region_power;

    /*
     * Avoid division by zero.
     */
    if (total_band_power > 0.0)
    {
        result->spectral_quality =
            dominant_region_power /
            total_band_power;
    }
    else
    {
        result->spectral_quality = 0.0;
    }

    /*
     * Require BOTH:
     *
     * 1. enough overall signal energy
     * 2. enough spectral concentration
     */
    result->periodic_activity_detected =
        (
            total_band_power >=
            detector->config.min_band_power
        ) &&
        (
            result->spectral_quality >=
            detector->config.min_spectral_quality
        );

    return 0;
}

int frequency_detector_update_cadence(
    frequency_detector_t *detector,
    const frequency_spectrum_output_t *spectrum,
    const frequency_dominant_result_t *dominant,
    const frequency_quality_result_t *quality,
    frequency_cadence_result_t *result)
{
    double integration_duration_s;
    double step_frequency_hz;
    double cadence_spm;

    double steps_added = 0.0;

    if ((detector == NULL) ||
        (spectrum == NULL) ||
        (dominant == NULL) ||
        (quality == NULL) ||
        (result == NULL))
    {
        return -1;
    }

    /*
     * Determine how much NEW time this FFT
     * result represents.
     *
     * The first FFT covers the complete initial
     * 128-sample window.
     *
     * Later FFTs overlap previous windows, so
     * only the time since the previous FFT result
     * may be integrated.
     */
    if (!detector->integration_initialized)
    {
        integration_duration_s =
            (double)FREQUENCY_FFT_SIZE /
            detector->config.sample_rate_hz;

        detector->integration_initialized = true;
    }
    else
    {
        const uint32_t elapsed_ms =
            spectrum->window_end_ms -
            detector->previous_window_end_ms;

        integration_duration_s =
            (double)elapsed_ms /
            1000.0;
    }

    detector->previous_window_end_ms =
        spectrum->window_end_ms;

    /*
     * Convert dominant spectral frequency into
     * estimated step frequency.
     */
    step_frequency_hz =
        dominant->dominant_frequency_hz *
        detector->config.step_frequency_multiplier;

    cadence_spm =
        step_frequency_hz * 60.0;

    /*
     * Only periodic windows contribute steps.
     */
    if (quality->periodic_activity_detected)
    {
        steps_added =
            step_frequency_hz *
            integration_duration_s;

        detector->fractional_steps +=
            steps_added;

        /*
         * Convert accumulated fractional steps
         * into whole events while preserving the
         * remainder.
         */
        while (detector->fractional_steps >= 1.0)
        {
            detector->estimated_step_count++;

            detector->fractional_steps -= 1.0;
        }
    }

    result->cadence_valid =
        quality->periodic_activity_detected;

    result->estimated_step_frequency_hz =
        step_frequency_hz;

    result->estimated_cadence_spm =
        cadence_spm;

    result->integration_duration_s =
        integration_duration_s;

    result->estimated_steps_added =
        steps_added;

    result->fractional_steps =
        detector->fractional_steps;

    result->estimated_step_count =
        detector->estimated_step_count;

    return 0;
}