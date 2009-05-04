#ifndef _ECLIPSE_TYPES_H_
#define _ECLIPSE_TYPES_H_

typedef float   pixelvalue ;
#define TABSPERPIX      (1000)
#define KERNEL_WIDTH    (2.0)
#define KERNEL_SAMPLES  (1+(int)(TABSPERPIX * KERNEL_WIDTH))
#define PI_NUMB     (3.1415926535897932384626433832795)
#define TANH_STEEPNESS   (5.0)
#define MAX_COLUMN_NUMBER  (40000)
#define MAX_LINE_NUMBER    (40000)

double * generate_interpolation_kernel(char * kernel_type);
double * generate_tanh_kernel(double steep);
double sinc(double x);
void reverse_tanh_kernel(double * data, int nn);
void trans (float *res, char ttype, int u, int v, int nx, int ny, float alpha, float beta);
void SaveFits (const float *array, int nx, int ny, const char*name);

float * image_warp_generic(float *image_in, int sizex, int sizey,
			   char *kernel_type, char proj, float alpha, float beta);
void CalcTrans (int inwidth, int inheight, float alpha=0.0, float beta=0.0);
void FreeTrans ();

#endif
