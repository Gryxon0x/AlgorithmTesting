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