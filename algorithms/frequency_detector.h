#ifndef FREQUENCY_DETECTOR_H
#define FREQUENCY_DETECTOR_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "accel_sample.h"


#define FREQUENCY_FFT_SIZE       128U
#define FREQUENCY_SPECTRUM_BINS  (FREQUENCY_FFT_SIZE / 2U + 1U)


typedef struct
{
    double sample_rate_hz;

    /*
     * Number of new samples between consecutive FFTs.
     *
     * 32 samples at 50 Hz = 0.64 seconds.
     */
    size_t hop_size;

} frequency_detector_config_t;


typedef struct
{
    uint32_t window_number;

    uint32_t window_start_ms;
    uint32_t window_end_ms;

    double mean_magnitude_g;

} frequency_spectrum_output_t;


typedef struct
{
    frequency_detector_config_t config;

    /*
     * Circular buffer containing the latest 128
     * acceleration-magnitude samples.
     */
    double magnitude_buffer[FREQUENCY_FFT_SIZE];
    uint32_t timestamp_buffer[FREQUENCY_FFT_SIZE];

    size_t write_index;
    size_t sample_count;
    size_t samples_since_spectrum;

    /*
     * One-sided spectrum:
     *
     * bin 0 ... bin 64
     */
    double power_spectrum[FREQUENCY_SPECTRUM_BINS];

    uint32_t window_count;

} frequency_detector_t;


void frequency_detector_default_config(
    frequency_detector_config_t *config);


int frequency_detector_init(
    frequency_detector_t *detector,
    const frequency_detector_config_t *config);


/*
 * Process one accelerometer sample.
 *
 * Returns true whenever a new FFT spectrum has
 * been calculated.
 */
bool frequency_detector_process(
    frequency_detector_t *detector,
    const accel_sample_t *sample,
    frequency_spectrum_output_t *output);


/*
 * Pointer to the most recently calculated
 * one-sided power spectrum.
 */
const double *frequency_detector_get_power_spectrum(
    const frequency_detector_t *detector);


/*
 * Convert an FFT-bin index to Hz.
 */
double frequency_detector_get_bin_frequency_hz(
    const frequency_detector_t *detector,
    size_t bin);


#endif /* FREQUENCY_DETECTOR_H */