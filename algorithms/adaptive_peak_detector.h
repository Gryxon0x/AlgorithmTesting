#ifndef ADAPTIVE_PEAK_DETECTOR_H
#define ADAPTIVE_PEAK_DETECTOR_H

#include <stdbool.h>
#include <stdint.h>

#include "accel_sample.h"


typedef struct
{
    /* Expected sample rate. */
    double sample_rate_hz;

    /*
     * Causal band-pass approximation:
     * first-order high-pass followed by first-order low-pass.
     */
    double low_cutoff_hz;
    double high_cutoff_hz;

    /*
     * Adaptive statistics:
     *
     * baseline = EMA(filtered signal)
     * activity = EMA(abs(filtered - baseline))
     */
    double baseline_alpha;
    double activity_alpha;

    /*
     * Adaptive peak requirements:
     *
     * required height =
     *     baseline + height_multiplier * activity
     *
     * required prominence =
     *     prominence_multiplier * activity
     */
    double height_multiplier;
    double prominence_multiplier;

    double min_peak_height_g;
    double max_peak_height_g;

    double min_prominence_g;
    double max_prominence_g;

    /*
     * Time allowed after a local maximum to find
     * a right-side valley for prominence estimation.
     */
    uint32_t confirmation_window_ms;

    /*
     * Minimum interval between accepted steps.
     */
    uint32_t min_event_interval_ms;

} adaptive_peak_config_t;


typedef struct
{
    double magnitude_g;
    double filtered_g;

    double baseline_g;
    double activity_g;

    double required_height_g;
    double required_prominence_g;

    bool candidate_active;

    double candidate_peak_g;
    double candidate_prominence_g;

    bool event_detected;
    uint32_t event_timestamp_ms;

} adaptive_peak_output_t;


typedef struct
{
    adaptive_peak_config_t config;

    /*
     * Filter coefficients/state.
     */
    double hp_alpha;
    double lp_alpha;

    bool filter_initialized;

    double hp_previous_input_g;
    double hp_previous_output_g;
    double lp_previous_output_g;

    /*
     * Adaptive statistics.
     */
    bool adaptive_initialized;

    double baseline_g;
    double activity_g;

    double required_height_g;
    double required_prominence_g;

    /*
     * Previous filtered samples used for
     * local maximum detection.
     */
    uint32_t history_count;

    double prev2_filtered_g;
    double prev1_filtered_g;

    uint32_t prev1_timestamp_ms;

    double prev1_required_height_g;
    double prev1_required_prominence_g;

    double prev1_baseline_g;
    double prev1_activity_g;

    /*
     * Valley tracking.
     */
    bool left_valley_initialized;
    double left_valley_g;

    /*
     * Peak awaiting prominence confirmation.
     */
    bool candidate_active;

    double candidate_peak_g;
    uint32_t candidate_timestamp_ms;

    double candidate_left_valley_g;
    double candidate_right_valley_g;

    double candidate_required_height_g;
    double candidate_required_prominence_g;

    double candidate_baseline_g;
    double candidate_activity_g;

    /*
     * Accepted events.
     */
    bool has_last_event;
    uint32_t last_event_timestamp_ms;

    uint32_t event_count;

} adaptive_peak_detector_t;


void adaptive_peak_detector_default_config(
    adaptive_peak_config_t *config);


int adaptive_peak_detector_init(
    adaptive_peak_detector_t *detector,
    const adaptive_peak_config_t *config);


bool adaptive_peak_detector_process(
    adaptive_peak_detector_t *detector,
    const accel_sample_t *sample,
    adaptive_peak_output_t *output);


uint32_t adaptive_peak_detector_get_count(
    const adaptive_peak_detector_t *detector);


#endif /* ADAPTIVE_PEAK_DETECTOR_H */