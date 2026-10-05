#ifndef CSV_READER_H
#define CSV_READER_H

#include <stddef.h>
#include <stdio.h>

#include "accel_sample.h"

#define CSV_READER_LINE_BUFFER_SIZE 2048U
#define CSV_READER_MAX_COLUMNS 64U

typedef enum
{
    CSV_READ_OK = 0,
    CSV_READ_EOF,
    CSV_READ_ERROR
} csv_read_status_t;

typedef struct
{
    FILE *file;
    char line[CSV_READER_LINE_BUFFER_SIZE];
    size_t line_number;

    int timestamp_column;
    int phone_time_column;
    int x_raw_column;
    int y_raw_column;
    int z_raw_column;
    int x_mg_column;
    int y_mg_column;
    int z_mg_column;

    char error_message[256];
} csv_reader_t;

int csv_reader_open(csv_reader_t *reader, const char *path);
csv_read_status_t csv_reader_next(csv_reader_t *reader, accel_sample_t *sample);
void csv_reader_close(csv_reader_t *reader);
const char *csv_reader_error(const csv_reader_t *reader);

#endif /* CSV_READER_H */
