#ifndef ACCEL_SAMPLE_H
#define ACCEL_SAMPLE_H

#include <stdint.h>

typedef struct
{
    uint32_t timestamp_ms;
    double phone_time_est_ms;

    int16_t x_raw;
    int16_t y_raw;
    int16_t z_raw;

    /* Calibrated acceleration from ax_mg/ay_mg/az_mg, converted to g. */
    double x_g;
    double y_g;
    double z_g;
} accel_sample_t;

#endif /* ACCEL_SAMPLE_H */
