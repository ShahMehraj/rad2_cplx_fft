/**
 * @file main.c
 * @brief WAV file processing using FFT
 *
 * This program reads audio samples from a WAV file, performs Fast Fourier Transform (FFT) on the samples,
 * and prints the real and imaginary parts of the FFT result.
 *
 * The program processes the audio samples in blocks of size `BLOCK_SIZE`. Processing in larger block sizes
 * may lead to unresponsiveness, so it is recommended to keep `BLOCK_SIZE` as 1024 or lower.
 *
 * For the remaining samples that are not a multiple of `BLOCK_SIZE`, the program finds the largest power of 4
 * that is less than or equal to the remaining samples and performs FFT on that size.
 *
 *
 * Note: This implementation assumes the availability of the `read_wav` function from the "wave.h" library,
 * as well as the `impeghd_rad2_cplx_fft` function from the "fft.h" library. Make sure to include the correct
 * header files and provide the necessary library files when compiling the program.
 *
 * @date 2023-06-04
 **/

#include <stdio.h>
#include <time.h>
#include "wave.h"
#include "impeghd_type_def.h"
#include "fft.h"

#define BLOCK_SIZE 1024

int main() {
    const char* file_name = "1khz_Sine_44_1khz.wav";
    UWORD32 num_samples;
    FLOAT32* ptr_real;
    read_wav(file_name, &num_samples, &ptr_real);
    //num_samples = (num_samples / BLOCK_SIZE) * BLOCK_SIZE;
    FLOAT32* ptr_imag = calloc(num_samples, sizeof(FLOAT32));
    FLOAT32* ptr_scratch = malloc(sizeof(FLOAT32) * (4 * 1024));

    UWORD32 num_of_iterations = num_samples / BLOCK_SIZE;

    for (UWORD32 i = 0; i < num_of_iterations; i++) {
        
        UWORD32 start_index = i * BLOCK_SIZE;
        UWORD32 end_index = (i + 1) * BLOCK_SIZE;
        if(end_index > num_samples) break;

        impeghd_rad2_cplx_fft(
          &ptr_real[start_index], 
          &ptr_imag[start_index], 
          BLOCK_SIZE, 
          ptr_scratch);
    }
    /**
     The following block finds the largest
     power of 4 less than or equal to the
     remaining_samples and computes their fft
    **/
    UWORD32 remaining_samples = num_samples % BLOCK_SIZE;
    if (remaining_samples > 0) {
        UWORD32 start_index = num_samples - remaining_samples;
        UWORD32 end_index = num_samples;

        UWORD32 power_of_4 = 1;
        while (power_of_4 * 4 <= remaining_samples) {
            power_of_4 *= 4;
        }

        impeghd_rad2_cplx_fft(
            &ptr_real[start_index],
            &ptr_imag[start_index],
            power_of_4,
            ptr_scratch
        );
    }
    for(UWORD32 i = 0; i < num_samples; i++){
      printf("%f %c %f i\n", ptr_real[i], (ptr_imag[i] < 0) ? '-' : '+', fabs(ptr_imag[i]));
    }

    free(ptr_real);
    free(ptr_imag);
    free(ptr_scratch); 

    return 0;
}
