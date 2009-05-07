#ifndef _ECLIPSE_TYPES_H_
#define _ECLIPSE_TYPES_H_

#include "image.h"

typedef float   pixelvalue ;
#define TABSPERPIX      (1000)
#define KERNEL_WIDTH    (2.0)
#define KERNEL_SAMPLES  (1+(int)(TABSPERPIX * KERNEL_WIDTH))
#define PI_NUMB     (3.1415926535897932384626433832795)
#define TANH_STEEPNESS   (5.0)

double * generate_interpolation_kernel(const char * kernel_type);
double * generate_tanh_kernel(double steep);
double sinc(double x);
void reverse_tanh_kernel(double * data, int nn);
void trans (float *res, char ttype, int u, int v, int nx, int ny, float alpha, float beta);

Image * image_warp_generic(const Image *image_in, const char *kernel_type);
void CalcTrans (int inwidth, int inheight, float alpha=0.0, float beta=0.0);
void FreeTrans ();

#endif
