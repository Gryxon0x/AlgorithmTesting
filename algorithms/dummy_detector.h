#ifndef DUMMY_DETECTOR_H
#define DUMMY_DETECTOR_H

#include <stdbool.h>
#include <stdint.h>

#include "accel_sample.h"

typedef struct
{
    uint32_t samples_processed;
    uint32_t step_count;
} dummy_detector_t;

void dummy_detector_init(dummy_detector_t *detector);
bool dummy_detector_process(dummy_detector_t *detector, const accel_sample_t *sample);
uint32_t dummy_detector_get_step_count(const dummy_detector_t *detector);

#endif /* DUMMY_DETECTOR_H */
