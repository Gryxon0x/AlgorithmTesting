#include "csv_reader.h"

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static size_t split_csv_line(char *line, char *fields[], size_t max_fields)
{
    size_t count = 0U;
    char *cursor = line;

    if ((line == NULL) || (fields == NULL) || (max_fields == 0U))
    {
        return 0U;
    }

    while ((count < max_fields) && (*cursor != '\0'))
    {
        char *start = cursor;
        int in_quotes = 0;

        if (*cursor == '"')
        {
            in_quotes = 1;
            start = ++cursor;
        }

        fields[count++] = start;

        while (*cursor != '\0')
        {
            if (in_quotes)
            {
                if ((*cursor == '"') && (cursor[1] == '"'))
                {
                    cursor += 2;
                    continue;
                }
                if (*cursor == '"')
                {
                    in_quotes = 0;
                    *cursor = '\0';
                    cursor++;
                    if (*cursor == ',')
                    {
                        *cursor = '\0';
                        cursor++;
                    }
                    break;
                }
            }
            else if (*cursor == ',')
            {
                *cursor = '\0';
                cursor++;
                break;
            }
            else if ((*cursor == '\r') || (*cursor == '\n'))
            {
                *cursor = '\0';
                break;
            }

            cursor++;
        }
    }

    return count;
}

static int find_column(char *fields[], size_t field_count, const char *name)
{
    size_t i;
    for (i = 0U; i < field_count; ++i)
    {
        if (strcmp(fields[i], name) == 0)
        {
            return (int)i;
        }
    }
    return -1;
}

static int parse_u32(const char *text, uint32_t *value)
{
    char *end = NULL;
    unsigned long parsed;

    errno = 0;
    parsed = strtoul(text, &end, 10);
    if ((errno != 0) || (end == text) || (*end != '\0') || (parsed > UINT32_MAX))
    {
        return -1;
    }

    *value = (uint32_t)parsed;
    return 0;
}

static int parse_i16(const char *text, int16_t *value)
{
    char *end = NULL;
    long parsed;

    errno = 0;
    parsed = strtol(text, &end, 10);
    if ((errno != 0) || (end == text) || (*end != '\0') ||
        (parsed < INT16_MIN) || (parsed > INT16_MAX))
    {
        return -1;
    }

    *value = (int16_t)parsed;
    return 0;
}

static int parse_double(const char *text, double *value)
{
    char *end = NULL;
    double parsed;

    errno = 0;
    parsed = strtod(text, &end);
    if ((errno != 0) || (end == text) || (*end != '\0'))
    {
        return -1;
    }

    *value = parsed;
    return 0;
}

int csv_reader_open(csv_reader_t *reader, const char *path)
{
    char *fields[CSV_READER_MAX_COLUMNS];
    size_t field_count;

    if ((reader == NULL) || (path == NULL))
    {
        return -1;
    }

    memset(reader, 0, sizeof(*reader));
    reader->timestamp_column = -1;
    reader->phone_time_column = -1;
    reader->x_raw_column = -1;
    reader->y_raw_column = -1;
    reader->z_raw_column = -1;
    reader->x_mg_column = -1;
    reader->y_mg_column = -1;
    reader->z_mg_column = -1;

    reader->file = fopen(path, "r");
    if (reader->file == NULL)
    {
        (void)snprintf(reader->error_message, sizeof(reader->error_message),
                       "Could not open '%s'.", path);
        return -1;
    }

    if (fgets(reader->line, (int)sizeof(reader->line), reader->file) == NULL)
    {
        (void)snprintf(reader->error_message, sizeof(reader->error_message),
                       "CSV file is empty or unreadable.");
        csv_reader_close(reader);
        return -1;
    }

    reader->line_number = 1U;
    field_count = split_csv_line(reader->line, fields, CSV_READER_MAX_COLUMNS);

    reader->timestamp_column = find_column(fields, field_count, "t_ms");
    reader->phone_time_column = find_column(fields, field_count, "phone_time_est_ms");
    reader->x_raw_column = find_column(fields, field_count, "ax_raw");
    reader->y_raw_column = find_column(fields, field_count, "ay_raw");
    reader->z_raw_column = find_column(fields, field_count, "az_raw");
    reader->x_mg_column = find_column(fields, field_count, "ax_mg");
    reader->y_mg_column = find_column(fields, field_count, "ay_mg");
    reader->z_mg_column = find_column(fields, field_count, "az_mg");

    if ((reader->timestamp_column < 0) || (reader->phone_time_column < 0) ||
        (reader->x_raw_column < 0) || (reader->y_raw_column < 0) ||
        (reader->z_raw_column < 0) || (reader->x_mg_column < 0) ||
        (reader->y_mg_column < 0) || (reader->z_mg_column < 0))
    {
        (void)snprintf(reader->error_message, sizeof(reader->error_message),
                       "CSV must contain t_ms, phone_time_est_ms, ax_raw, ay_raw, "
                       "az_raw, ax_mg, ay_mg and az_mg columns.");
        csv_reader_close(reader);
        return -1;
    }

    return 0;
}

csv_read_status_t csv_reader_next(csv_reader_t *reader, accel_sample_t *sample)
{
    char *fields[CSV_READER_MAX_COLUMNS];
    size_t field_count;
    int max_required_column;
    double x_mg;
    double y_mg;
    double z_mg;

    if ((reader == NULL) || (reader->file == NULL) || (sample == NULL))
    {
        return CSV_READ_ERROR;
    }

    if (fgets(reader->line, (int)sizeof(reader->line), reader->file) == NULL)
    {
        if (feof(reader->file) != 0)
        {
            return CSV_READ_EOF;
        }

        (void)snprintf(reader->error_message, sizeof(reader->error_message),
                       "Read error near line %zu.", reader->line_number + 1U);
        return CSV_READ_ERROR;
    }

    reader->line_number++;
    field_count = split_csv_line(reader->line, fields, CSV_READER_MAX_COLUMNS);

    max_required_column = reader->timestamp_column;
#define MAX_COL(col) do { if ((col) > max_required_column) max_required_column = (col); } while (0)
    MAX_COL(reader->phone_time_column);
    MAX_COL(reader->x_raw_column);
    MAX_COL(reader->y_raw_column);
    MAX_COL(reader->z_raw_column);
    MAX_COL(reader->x_mg_column);
    MAX_COL(reader->y_mg_column);
    MAX_COL(reader->z_mg_column);
#undef MAX_COL

    if (field_count <= (size_t)max_required_column)
    {
        (void)snprintf(reader->error_message, sizeof(reader->error_message),
                       "Line %zu has too few columns.", reader->line_number);
        return CSV_READ_ERROR;
    }

    if ((parse_u32(fields[reader->timestamp_column], &sample->timestamp_ms) != 0) ||
        (parse_double(fields[reader->phone_time_column], &sample->phone_time_est_ms) != 0) ||
        (parse_i16(fields[reader->x_raw_column], &sample->x_raw) != 0) ||
        (parse_i16(fields[reader->y_raw_column], &sample->y_raw) != 0) ||
        (parse_i16(fields[reader->z_raw_column], &sample->z_raw) != 0) ||
        (parse_double(fields[reader->x_mg_column], &x_mg) != 0) ||
        (parse_double(fields[reader->y_mg_column], &y_mg) != 0) ||
        (parse_double(fields[reader->z_mg_column], &z_mg) != 0))
    {
        (void)snprintf(reader->error_message, sizeof(reader->error_message),
                       "Invalid numeric value on line %zu.", reader->line_number);
        return CSV_READ_ERROR;
    }

    sample->x_g = x_mg / 1000.0;
    sample->y_g = y_mg / 1000.0;
    sample->z_g = z_mg / 1000.0;

    return CSV_READ_OK;
}

void csv_reader_close(csv_reader_t *reader)
{
    if ((reader != NULL) && (reader->file != NULL))
    {
        fclose(reader->file);
        reader->file = NULL;
    }
}

const char *csv_reader_error(const csv_reader_t *reader)
{
    return (reader != NULL) ? reader->error_message : "Unknown CSV reader error.";
}
