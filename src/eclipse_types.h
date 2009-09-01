#ifndef _ECLIPSE_TYPES_H_
#define _ECLIPSE_TYPES_H_

#include "image.h"

typedef float   pixelvalue ;
#define TABSPERPIX      (1000)
#define KERNEL_WIDTH    (2.0)
#define KERNEL_SAMPLES  (1+(int)(TABSPERPIX * KERNEL_WIDTH))
#define PI_NUMB     (3.1415926535897932384626433832795)
#define TANH_STEEPNESS   (5.0)

float * generate_interpolation_kernel(const char * kernel_type);
float * generate_tanh_kernel(double steep);
double sinc(double x);
void reverse_tanh_kernel(float * data, int nn);
void trans (float *res, char ttype, int u, int v, int nx, int ny, float alpha, float beta);

int function1d_search_value(pixelvalue * x, int len, pixelvalue key, int * found_ptr) ;
int function1d_natural_spline(pixelvalue *, pixelvalue	*, int,	pixelvalue	* 	splx,
    	pixelvalue	* 	sply,
			      int 			spllen);

Image * image_warp_generic(const Image *image_in, Image *image_out=NULL);
void CalcTrans (int inwidth, int inheight, float alpha, float beta, const char *kernel_type);
void FreeTrans ();

#endif
