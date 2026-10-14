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

    /*
     * Frequency band in which we search
     * for walking-related periodicity.
     */
    double search_low_hz;
    double search_high_hz;

    /*
     * Number of FFT bins on either side of the
     * dominant bin included in its local power.
     *
     * 1 means:
     * dominant_bin - 1
     * dominant_bin
     * dominant_bin + 1
     */
    size_t quality_neighbor_bins;

    /*
     * Minimum fraction of walking-band power
     * that must be concentrated around the
     * dominant frequency.
     */
    double min_spectral_quality;

    /*
     * Prevent very-low-energy/noisy windows from
     * being classified as periodic merely because
     * one bin happens to dominate.
     */
    double min_band_power;

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

typedef struct
{
    size_t dominant_bin;

    double dominant_frequency_hz;
    double dominant_power;

} frequency_dominant_result_t;

typedef struct 
{
    double total_band_power;
    double dominant_region_power;
    double spectral_quality;
    bool periodic_activity_detected;
} frequency_quality_result_t;


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

/*
 * Stage 2:
 * Find the strongest FFT component inside
 * the configured analysis-frequency range.
 *
 * Returns:
 *   0  success
 *  -1  error
 */
int frequency_detector_find_dominant(
    const frequency_detector_t *detector,
    frequency_dominant_result_t *result);

/*
 * Stage 3:
 *
 * Evaluate how strongly the spectrum is
 * concentrated around the dominant frequency.
 *
 * Returns:
 *   0  success
 *  -1  error
 */
int frequency_detector_evaluate_quality(
    const frequency_detector_t *detector,
    const frequency_dominant_result_t *dominant,
    frequency_quality_result_t *result);


#endif /* FREQUENCY_DETECTOR_H */