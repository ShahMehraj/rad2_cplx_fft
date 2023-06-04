/**********************************************************************************
 * @file wav_reader.h
 * @brief WAV file reader utility functions
 *
 * This header file provides utility functions for reading audio data from WAV files.
 * It includes a structure definition for the WAV file header and a function to read
 * the audio samples from a WAV file into a float array.
 *
 * The `read_wav` function reads the audio samples from a given WAV file and returns
 * the number of samples and a pointer to a dynamically allocated float array that
 * holds the real part of the samples. The imaginary part of the audio samples is
 * assumed to be 0, as the WAV file format primarily stores real audio data.
 *
 *
 * Usage example:
 * `````````````````````````````````````````````````````````````````````````````````
 * uint32_t num_samples;
 * float* ptr_real;
 * read_wav("audio.wav", &num_samples, &ptr_real);
 * // Process the audio data...
 * free(ptr_real);  // Don't forget to free the allocated memory

 * `````````````````````````````````````````````````````````````````````````````````
 * @date 2023-06-04
***********************************************************************************/


#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char riff[4];
    uint32_t overall_size;
    char wave[4];
    char fmt_chunk_marker[4];
    uint32_t length_of_fmt;
    uint16_t format_type;
    uint16_t channels;
    uint32_t sample_rate;
    uint32_t byterate;
    uint16_t block_align;
    uint16_t bits_per_sample;
    char data_chunk_header[4];
    uint32_t data_size;
} Header;

void read_wav(const char* file_name, uint32_t* num_samples, float** ptr_real) {
    FILE* file = fopen(file_name, "rb");
    if (file == NULL) {
        perror("Failed to open file");
        *num_samples = 0;
        *ptr_real = NULL;
        return;
    }

    Header header;
    fread(&header, sizeof(Header), 1, file);

    *num_samples = header.data_size / (header.channels * (header.bits_per_sample / 8));
    *ptr_real = (float*)malloc(sizeof(float) * (*num_samples));

    for (uint32_t i = 0; i < *num_samples; i++) {
        for (uint16_t j = 0; j < header.channels; j++) {
            int16_t data_in_channel = 0;
            fread(&data_in_channel, sizeof(int16_t), 1, file);
            (*ptr_real)[(i * header.channels + j)] = (float)data_in_channel;
        }
    }

    fclose(file);
}
