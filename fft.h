#define _USE_MATH_DEFINES
#include <math.h>
#include "impeghd_type_def.h"
#include "ia_core_coder_constants.h"
#include "impeghd_fft_ifft_rom.h"
#include "impeghd_fft_ifft.h"
#include "impeghd_intrinsics_flt.h"
#include "ia_core_coder_basic_ops32.h"
VOID impeghd_rad2_cplx_fft(FLOAT32 *ptr_real, FLOAT32 *ptr_imag, WORD32 n_points,
                           FLOAT32 *ptr_scratch)
{
  WORD32 i, j, k, n_stages, h2;
  FLOAT32 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
  WORD32 del, nodespacing, in_loop_cnt;
  WORD32 not_power_4;
  WORD32 dig_rev_shift;
  WORD32 m_points = n_points;
  FLOAT32 *ptr_x = ptr_scratch;
  FLOAT32 *y = ptr_scratch + 2048;
  FLOAT32 *ptr_y = y;
  const FLOAT32 *ptr_w;

  dig_rev_shift = norm32(m_points) + 1 - 16;
  n_stages = 30 - norm32(m_points);
  not_power_4 = n_stages & 1;

  n_stages = n_stages >> 1;

  ptr_w = ia_fft_twiddle_table_float;
  for (i = 0; i < n_points; i++)
  {
    ptr_x[2 * i] = ptr_real[i];
    ptr_x[2 * i + 1] = ptr_imag[i];
  }

  for (i = 0; i < n_points; i += 4)
  {
    FLOAT32 *inp = ptr_x;
    FLOAT32 tmk;

    DIG_REV(i, dig_rev_shift, h2);
    if (not_power_4)
    {
      h2 += 1;
      h2 &= ~1;
    }
    inp += (h2);

    x0r = *inp;
    x0i = *(inp + 1);
    inp += (n_points >> 1);

    x1r = *inp;
    x1i = *(inp + 1);
    inp += (n_points >> 1);

    x2r = *inp;
    x2i = *(inp + 1);
    inp += (n_points >> 1);

    x3r = *inp;
    x3i = *(inp + 1);

    x0r = ia_add_flt(x0r, x2r);
    x0i = ia_add_flt(x0i, x2i);

    tmk = ia_sub_flt(x0r, x2r);
    x2r = ia_sub_flt(tmk, x2r);
    tmk = ia_sub_flt(x0i, x2i);
    x2i = ia_sub_flt(tmk, x2i);

    x1r = ia_add_flt(x1r, x3r);
    x1i = ia_add_flt(x1i, x3i);

    tmk = ia_sub_flt(x1r, x3r);
    x3r = ia_sub_flt(tmk, x3r);
    tmk = ia_sub_flt(x1i, x3i);
    x3i = ia_sub_flt(tmk, x3i);

    x0r = ia_add_flt(x0r, x1r);
    x0i = ia_add_flt(x0i, x1i);

    tmk = ia_sub_flt(x0r, x1r);
    x1r = ia_sub_flt(tmk, x1r);
    tmk = ia_sub_flt(x0i, x1i);
    x1i = ia_sub_flt(tmk, x1i);

    x2r = ia_add_flt(x2r, x3i);
    x2i = ia_sub_flt(x2i, x3r);

    tmk = ia_sub_flt(x2r, x3i);
    x3i = ia_sub_flt(tmk, x3i);
    tmk = ia_add_flt(x2i, x3r);
    x3r = ia_add_flt(tmk, x3r);

    *ptr_y++ = x0r;
    *ptr_y++ = x0i;
    *ptr_y++ = x2r;
    *ptr_y++ = x2i;
    *ptr_y++ = x1r;
    *ptr_y++ = x1i;
    *ptr_y++ = x3i;
    *ptr_y++ = x3r;
    
  }
  
  
  ptr_y -= 2 * n_points;
  
  del = 4;
  nodespacing = 64;
  in_loop_cnt = n_points >> 4;
  //printf("%d", n_points >> 4);
  for (i = n_stages - 1; i > 0; i--)
  { 
    const FLOAT32 *twiddles = ptr_w;
    FLOAT32 *data = ptr_y;
    FLOAT32 W1, W2, W3, W4, W5, W6;
    WORD32 sec_loop_cnt;
    
    for (k = in_loop_cnt; k != 0; k--)
    {
      x0r = (*data);
      x0i = (*(data + 1));
      data += (del << 1);
      //printf("%d\n", k);

      x1r = (*data);
      x1i = (*(data + 1));
      data += (del << 1);

      x2r = (*data);
      x2i = (*(data + 1));
      data += (del << 1);

      x3r = (*data);
      x3i = (*(data + 1));
      data -= 3 * (del << 1);

      x0r = ia_add_flt(x0r, x2r);
      x0i = ia_add_flt(x0i, x2i);
      x2r = ia_msu_flt(x0r, x2r, 2);
      x2i = ia_msu_flt(x0i, x2i, 2);
      x1r = ia_add_flt(x1r, x3r);
      x1i = ia_add_flt(x1i, x3i);
      x3r = ia_msu_flt(x1r, x3r, 2);
      x3i = ia_msu_flt(x1i, x3i, 2);

      x0r = ia_add_flt(x0r, x1r);
      x0i = ia_add_flt(x0i, x1i);
      x1r = ia_msu_flt(x0r, x1r, 2);
      x1i = ia_msu_flt(x0i, x1i, 2);
      x2r = ia_add_flt(x2r, x3i);
      x2i = ia_sub_flt(x2i, x3r);
      x3i = ia_msu_flt(x2r, x3i, 2);
      x3r = ia_mac_flt(x2i, x3r, 2);

      *data = x0r;
      *(data + 1) = x0i;
      data += (del << 1);

      *data = x2r;
      *(data + 1) = x2i;
      data += (del << 1);

      *data = x1r;
      *(data + 1) = x1i;
      data += (del << 1);

      *data = x3i;
      *(data + 1) = x3r;
      data += (del << 1);
    }
    data = ptr_y + 2;

    sec_loop_cnt = (nodespacing * del);
    sec_loop_cnt = (sec_loop_cnt / 4) + (sec_loop_cnt / 8) - (sec_loop_cnt / 16) +
                   (sec_loop_cnt / 32) - (sec_loop_cnt / 64) + (sec_loop_cnt / 128) -
                   (sec_loop_cnt / 256);
    
    for (j = nodespacing; j <= sec_loop_cnt; j += nodespacing)
    {
      W1 = *(twiddles + j);
      W4 = *(twiddles + j + 257);
      W2 = *(twiddles + (j << 1));
      W5 = *(twiddles + (j << 1) + 257);
      W3 = *(twiddles + j + (j << 1));
      W6 = *(twiddles + j + (j << 1) + 257);

      // twiddles += nodespacing;
      for (k = in_loop_cnt; k != 0; k--)
      {
        FLOAT32 tmp;
        /*x0 is loaded later to avoid register crunch*/

        data += (del << 1);

        x1r = *data;
        x1i = *(data + 1);
        data += (del << 1);

        x2r = *data;
        x2i = *(data + 1);
        data += (del << 1);

        x3r = *data;
        x3i = *(data + 1);
        data -= 3 * (del << 1);

        tmp = ia_sub_flt(ia_mul_flt(x1r, W1), ia_mul_flt(x1i, W4));
        x1i = ia_mac_flt(ia_mul_flt(x1r, W4), x1i, W1);
        x1r = tmp;

        tmp = ia_sub_flt(ia_mul_flt(x2r, W2), ia_mul_flt(x2i, W5));
        x2i = ia_mac_flt(ia_mul_flt(x2r, W5), x2i, W2);
        x2r = tmp;

        tmp = ia_sub_flt(ia_mul_flt(x3r, W3), ia_mul_flt(x3i, W6));
        x3i = ia_mac_flt(ia_mul_flt(x3r, W6), x3i, W3);
        x3r = tmp;

        x0r = (*data);
        x0i = (*(data + 1));

        x0r = ia_add_flt(x0r, (x2r));
        x0i = ia_add_flt(x0i, (x2i));
        x2r = ia_msu_flt(x0r, x2r, 2);
        x2i = ia_msu_flt(x0i, x2i, 2);
        x1r = ia_add_flt(x1r, x3r);
        x1i = ia_add_flt(x1i, x3i);
        x3r = ia_msu_flt(x1r, x3r, 2);
        x3i = ia_msu_flt(x1i, x3i, 2);

        x0r = ia_add_flt(x0r, (x1r));
        x0i = ia_add_flt(x0i, (x1i));
        x1r = ia_msu_flt(x0r, x1r, 2);
        x1i = ia_msu_flt(x0i, x1i, 2);
        x2r = ia_add_flt(x2r, (x3i));
        x2i = ia_sub_flt(x2i, (x3r));
        x3i = ia_msu_flt(x2r, x3i, 2);
        x3r = ia_mac_flt(x2i, x3r, 2);

        *data = x0r;
        *(data + 1) = x0i;
        data += (del << 1);

        *data = x2r;
        *(data + 1) = x2i;
        data += (del << 1);

        *data = x1r;
        *(data + 1) = x1i;
        data += (del << 1);

        *data = x3i;
        *(data + 1) = x3r;
        data += (del << 1);
      }
      data -= 2 * n_points;
      data += 2;
    }
    //printf("%d\n", j);

    for (; j <= (nodespacing * del) >> 1; j += nodespacing)
    {
      W1 = *(twiddles + j);
      W4 = *(twiddles + j + 257);
      W2 = *(twiddles + (j << 1));
      W5 = *(twiddles + (j << 1) + 257);
      W3 = *(twiddles + j + (j << 1) - 256);
      W6 = *(twiddles + j + (j << 1) + 1);
      // twiddles += nodespacing;

      for (k = in_loop_cnt; k != 0; k--)
      {
        FLOAT32 tmp;
        /*x0 is loaded later to avoid register crunch*/

        data += (del << 1);

        x1r = *data;
        x1i = *(data + 1);
        data += (del << 1);

        x2r = *data;
        x2i = *(data + 1);
        data += (del << 1);

        x3r = *data;
        x3i = *(data + 1);
        data -= 3 * (del << 1);

        tmp = ia_sub_flt(ia_mul_flt(x1r, W1), ia_mul_flt(x1i, W4));
        x1i = ia_mac_flt(ia_mul_flt(x1r, W4), x1i, W1);
        x1r = tmp;

        tmp = ia_sub_flt(ia_mul_flt(x2r, W2), ia_mul_flt(x2i, W5));
        x2i = ia_mac_flt(ia_mul_flt(x2r, W5), x2i, W2);
        x2r = tmp;

        tmp = ia_add_flt(ia_mul_flt(x3r, W6), ia_mul_flt(x3i, W3));
        x3i = ia_add_flt(ia_negate_flt(ia_mul_flt(x3r, W3)), ia_mul_flt(x3i, W6));
        x3r = tmp;

        x0r = (*data /*/ 2 */);
        x0i = (*(data + 1));

        x0r = ia_add_flt(x0r, (x2r));
        x0i = ia_add_flt(x0i, (x2i));
        x2r = ia_msu_flt(x0r, x2r, 2);
        x2i = ia_msu_flt(x0i, x2i, 2);
        x1r = ia_add_flt(x1r, x3r);
        x1i = ia_add_flt(x1i, x3i);
        x3r = ia_msu_flt(x1r, x3r, 2);
        x3i = ia_msu_flt(x1i, x3i, 2);

        x0r = ia_add_flt(x0r, (x1r));
        x0i = ia_add_flt(x0i, (x1i));
        x1r = ia_msu_flt(x0r, x1r, 2);
        x1i = ia_msu_flt(x0i, x1i, 2);
        x2r = ia_add_flt(x2r, (x3i));
        x2i = ia_sub_flt(x2i, (x3r));
        x3i = ia_msu_flt(x2r, x3i, 2);
        x3r = ia_mac_flt(x2i, x3r, 2);

        *data = x0r;
        *(data + 1) = x0i;
        data += (del << 1);

        *data = x2r;
        *(data + 1) = x2i;
        data += (del << 1);

        *data = x1r;
        *(data + 1) = x1i;
        data += (del << 1);

        *data = x3i;
        *(data + 1) = x3r;
        data += (del << 1);
      }
      data -= 2 * n_points;
      data += 2;
    }
    //printf("%d\n", j);
    for (; j <= sec_loop_cnt * 2; j += nodespacing)
    {
      W1 = *(twiddles + j);
      W4 = *(twiddles + j + 257);
      W2 = *(twiddles + (j << 1) - 256);
      W5 = *(twiddles + (j << 1) + 1);
      W3 = *(twiddles + j + (j << 1) - 256);
      W6 = *(twiddles + j + (j << 1) + 1);

      // twiddles += nodespacing;
      for (k = in_loop_cnt; k != 0; k--)
      {
        FLOAT32 tmp;
        /*x0 is loaded later to avoid register crunch*/

        data += (del << 1);

        x1r = *data;
        x1i = *(data + 1);
        data += (del << 1);

        x2r = *data;
        x2i = *(data + 1);
        data += (del << 1);

        x3r = *data;
        x3i = *(data + 1);
        data -= 3 * (del << 1);

        tmp = ia_sub_flt(ia_mul_flt(x1r, W1), ia_mul_flt(x1i, W4));
        x1i = ia_mac_flt(ia_mul_flt(x1r, W4), x1i, W1);
        x1r = tmp;

        tmp = ia_add_flt(ia_mul_flt(x2r, W5), ia_mul_flt(x2i, W2));
        x2i = ia_add_flt(ia_negate_flt(ia_mul_flt(x2r, W2)), ia_mul_flt(x2i, W5));
        x2r = tmp;

        tmp = ia_add_flt(ia_mul_flt(x3r, W6), ia_mul_flt(x3i, W3));
        x3i = ia_add_flt(ia_negate_flt(ia_mul_flt(x3r, W3)), ia_mul_flt(x3i, W6));
        x3r = tmp;

        x0r = (*data);
        x0i = (*(data + 1));

        x0r = ia_add_flt(x0r, (x2r));
        x0i = ia_add_flt(x0i, (x2i));
        x2r = ia_msu_flt(x0r, x2r, 2);
        x2i = ia_msu_flt(x0i, x2i, 2);
        x1r = ia_add_flt(x1r, x3r);
        x1i = ia_add_flt(x1i, x3i);
        x3r = ia_msu_flt(x1r, x3r, 2);
        x3i = ia_msu_flt(x1i, x3i, 2);

        x0r = ia_add_flt(x0r, (x1r));
        x0i = ia_add_flt(x0i, (x1i));
        x1r = ia_msu_flt(x0r, x1r, 2);
        x1i = ia_msu_flt(x0i, x1i, 2);
        x2r = ia_add_flt(x2r, (x3i));
        x2i = ia_sub_flt(x2i, (x3r));
        x3i = ia_msu_flt(x2r, x3i, 2);
        x3r = ia_mac_flt(x2i, x3r, 2);

        *data = x0r;
        *(data + 1) = x0i;
        data += (del << 1);

        *data = x2r;
        *(data + 1) = x2i;
        data += (del << 1);

        *data = x1r;
        *(data + 1) = x1i;
        data += (del << 1);

        *data = x3i;
        *(data + 1) = x3r;
        data += (del << 1);
      }
      data -= 2 * n_points;
      data += 2;
    }
    //printf("%d\n", j);
    for (; j < nodespacing * del; j += nodespacing)
    {
      W1 = *(twiddles + j);
      W4 = *(twiddles + j + 257);
      W2 = *(twiddles + (j << 1) - 256);
      W5 = *(twiddles + (j << 1) + 1);
      W3 = *(twiddles + j + (j << 1) - 512);
      W6 = *(twiddles + j + (j << 1) - 512 + 257);

      // twiddles += nodespacing;
      for (k = in_loop_cnt; k != 0; k--)
      {
        FLOAT32 tmp;
        /*x0 is loaded later to avoid register crunch*/

        data += (del << 1);

        x1r = *data;
        x1i = *(data + 1);
        data += (del << 1);

        x2r = *data;
        x2i = *(data + 1);
        data += (del << 1);

        x3r = *data;
        x3i = *(data + 1);
        data -= 3 * (del << 1);

        tmp = ia_sub_flt(ia_mul_flt(x1r, W1), ia_mul_flt(x1i, W4));
        x1i = ia_mac_flt(ia_mul_flt(x1r, W4), x1i, W1);
        x1r = tmp;

        tmp = ia_add_flt(ia_mul_flt(x2r, W5), ia_mul_flt(x2i, W2));
        x2i = ia_add_flt(ia_negate_flt(ia_mul_flt(x2r, W2)), ia_mul_flt(x2i, W5));
        x2r = tmp;

        tmp = ia_add_flt(ia_negate_flt(ia_mul_flt(x3r, W3)), ia_mul_flt(x3i, W6));
        x3i = ia_mac_flt(ia_mul_flt(x3r, W6), x3i, W3);
        x3r = tmp;

        x0r = (*data);
        x0i = (*(data + 1));

        x0r = ia_add_flt(x0r, (x2r));
        x0i = ia_add_flt(x0i, (x2i));
        x2r = ia_msu_flt(x0r, x2r, 2);
        x2i = ia_msu_flt(x0i, x2i, 2);
        x1r = ia_add_flt(x1r, x3r);
        x1i = ia_sub_flt(x1i, x3i);
        x3r = ia_msu_flt(x1r, x3r, 2);
        x3i = ia_mac_flt(x1i, x3i, 2);

        x0r = ia_add_flt(x0r, (x1r));
        x0i = ia_add_flt(x0i, (x1i));
        x1r = ia_msu_flt(x0r, x1r, 2);
        x1i = ia_msu_flt(x0i, x1i, 2);
        x2r = ia_add_flt(x2r, (x3i));
        x2i = ia_sub_flt(x2i, (x3r));
        x3i = ia_msu_flt(x2r, x3i, 2);
        x3r = ia_mac_flt(x2i, x3r, 2);

        *data = x0r;
        *(data + 1) = x0i;
        data += (del << 1);

        *data = x2r;
        *(data + 1) = x2i;
        data += (del << 1);

        *data = x1r;
        *(data + 1) = x1i;
        data += (del << 1);

        *data = x3i;
        *(data + 1) = x3r;
        data += (del << 1);
      }
      data -= 2 * n_points;
      data += 2;
    }
    nodespacing >>= 2;
    del <<= 2;
    in_loop_cnt >>= 2;
  }
  if (not_power_4)
  {
    const FLOAT32 *twiddles = ptr_w;
    nodespacing <<= 1;

    for (j = del / 2; j != 0; j--)
    {
      FLOAT32 W1 = *twiddles;
      FLOAT32 W4 = *(twiddles + 257);
      FLOAT32 tmp;
      twiddles += nodespacing;

      x0r = *ptr_y;
      x0i = *(ptr_y + 1);
      ptr_y += (del << 1);

      x1r = *ptr_y;
      x1i = *(ptr_y + 1);

      tmp = ia_sub_flt(ia_mul_flt(x1r, W1), ia_mul_flt(x1i, W4));
      x1i = (FLOAT32)ia_mac_flt(ia_mul_flt(x1r, W4), x1i, W1);
      x1r = tmp;

      *ptr_y = ia_sub_flt((x0r), (x1r));
      *(ptr_y + 1) = ia_sub_flt((x0i), (x1i));
      ptr_y -= (del << 1);

      *ptr_y = ia_add_flt((x0r), (x1r));
      *(ptr_y + 1) = ia_add_flt((x0i), (x1i));
      ptr_y += 2;
    }
    twiddles = ptr_w;
    for (j = del / 2; j != 0; j--)
    {
      FLOAT32 W1 = *twiddles;
      FLOAT32 W4 = *(twiddles + 257);
      FLOAT32 tmp;
      twiddles += nodespacing;

      x0r = *ptr_y;
      x0i = *(ptr_y + 1);
      ptr_y += (del << 1);

      x1r = *ptr_y;
      x1i = *(ptr_y + 1);

      tmp = ia_add_flt(ia_mul_flt(x1r, W4), ia_mul_flt(x1i, W1));
      x1i = ia_add_flt(ia_negate_flt(ia_mul_flt(x1r, W1)), ia_mul_flt(x1i, W4));
      x1r = tmp;

      *ptr_y = ia_sub_flt((x0r), (x1r));
      *(ptr_y + 1) = ia_sub_flt((x0i), (x1i));
      ptr_y -= (del << 1);

      *ptr_y = ia_add_flt((x0r), (x1r));
      *(ptr_y + 1) = ia_add_flt((x0i), (x1i));
      ptr_y += 2;
    }
    //printf("%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n", x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i);
  }

  for (i = 0; i < n_points; i++)
  {
    ptr_real[i] = y[2 * i];
    ptr_imag[i] = y[2 * i + 1];
    //printf("%.2f  %.2f\n", y[2 * i], y[2 * i + 1]);
  }
}
