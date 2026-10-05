#include "dummy_detector.h"

void dummy_detector_init(dummy_detector_t *detector)
{
    detector->samples_processed = 0U;
    detector->step_count = 0U;
}

bool dummy_detector_process(dummy_detector_t *detector, const accel_sample_t *sample)
{
    (void)sample;
    detector->samples_processed++;

    /*
     * Placeholder detector: no step-detection logic yet.
     * Returning true will later mean "a step was detected for this sample".
     */
    return false;
}

uint32_t dummy_detector_get_step_count(const dummy_detector_t *detector)
{
    return detector->step_count;
}
