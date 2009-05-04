#include "fitsio.h"
#include "FreeImage.h"

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <cstring>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include <time.h>
#include <sys/time.h>

#include "eclipse_types.h"

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

int main( int argc, char ** argv );
void Help (void);
BYTE CastFloat (float);

float *data[24];
bool verbose = false;
int pad = 3;

// Dimensions of the output image
int outwidth  = 1408;
int outheight = 1408;

// Orientation of the dome 
float alpha = 0.0; // Rotation around z-axis
float beta = 0.0;  // Tilt of the dome

float angle = 90.0; // Not used at the moment!

int main( int argc, char ** argv )
{
  int n;
  int npix_in;
  int pheight;
  int pwidth;
  timeval tv1,tv2;
  long t;

  if (argc == 1) {
    Help ();
    return (EXIT_SUCCESS);
  }

  while ( (n=getopt (argc, argv, "v:f:z:cg:hm:p:a:b:t:o:s:x:y:l:")) != -1) {
    char * temp;
    float tempvalf;
    int tempvali;

    switch (n) {
    case 'l':
      //xlabel = optarg;
      break;
    case 'm':
      //ylabel = optarg;
      break;
    case 'h':
      Help ();
      return (EXIT_SUCCESS);
    }
  }

  if (optind < argc ) {
    int status = 0; 
    fitsfile *fptr;
    int hdutype;

  }

  FreeImage_Initialise ();
  int width,height, bpp;

  gettimeofday (&tv1, NULL);

  for (int image = 0; image<6; image++) {
    FIBITMAP *bitmap1 = FreeImage_Load (FIF_TARGA, argv[image+1]);
    FREE_IMAGE_TYPE imtype1 = FreeImage_GetImageType (bitmap1);
    if (image == 0) {
      bpp = FreeImage_GetBPP (bitmap1);
      width = FreeImage_GetWidth (bitmap1);
      height = FreeImage_GetHeight (bitmap1);
      npix_in = (width+2*pad)*(height+2*pad);
      pwidth = width+2*pad;
      pheight = height+2*pad;
    }
    else {
      int bpp1 = FreeImage_GetBPP (bitmap1);
      int width1 = FreeImage_GetWidth (bitmap1);
      int height1 = FreeImage_GetHeight (bitmap1);
      if (bpp!=bpp1 || width!=width1 || height!=height1) {
	printf ("Incompatible images\n");
	FreeImage_Unload (bitmap1);
	FreeImage_DeInitialise ();
	for (int i=0; i<3*image; i++)
	  delete[] data[i];
      }
    }

    BYTE *bits1 = FreeImage_GetBits (bitmap1);
    
    data[3*image] = new float[npix_in];
    for (int i = 0; i<width; i++) 
      for (int j = 0; j<height; j++) 
	(data[3*image])[i+pad+pwidth*(j+pad)] = bits1[3*(i+j*width)+FI_RGBA_RED];

    data[3*image+1] = new float[npix_in];
    for (int i = 0; i<width; i++) 
      for (int j = 0; j<height; j++) 
	(data[3*image+1])[i+pad+pwidth*(j+pad)] = bits1[3*(i+j*width)+FI_RGBA_GREEN];

    data[3*image+2] = new float[npix_in];
    for (int i = 0; i<width; i++) 
      for (int j = 0; j<height; j++) 
	(data[3*image+2])[i+pad+pwidth*(j+pad)] = bits1[3*(i+j*width)+FI_RGBA_BLUE];
    
    FreeImage_Unload (bitmap1);
  }
  // printf ("Images loaded successfully\n");

  if (pad>0) {
    for (int k=0; k<3; k++)
      for (int i=0; i<pad; i++)
	for (int j=pad; j<pad+height; j++) {
	  (data[6+k])[i+pwidth*j] = (data[9+k])[pwidth-2*pad+i+pwidth*j];
	  (data[0+k])[i+pwidth*j] = (data[12+k])[pwidth-2*pad+i+pwidth*j];
	  (data[6+k])[pwidth-pad+i+pwidth*j] = (data[12+k])[pad+i+pwidth*j];
	  (data[0+k])[pwidth-pad+i+pwidth*j] = (data[9+k])[pad+i+pwidth*j];
	  (data[9+k])[pwidth-pad+i+pwidth*j] = (data[6+k])[pad+i+pwidth*j];
	  (data[9+k])[i+pwidth*j] = (data[0+k])[pwidth-2*pad+i+pwidth*j];
	  (data[12+k])[i+pwidth*j] = (data[6+k])[pwidth-2*pad+i+pwidth*j];
	  (data[12+k])[pwidth-pad+i+pwidth*j] = (data[0+k])[pad+i+pwidth*j];
       }
    for (int k=0; k<3; k++)
      for (int j=0; j<pad; j++)
	for (int i=pad; i<pad+width; i++) {
	  (data[6+k])[i+pwidth*(pheight-pad+j)] = (data[15+k])[i+pwidth*(pheight-2*pad+j)];
	  (data[15+k])[i+pwidth*(pheight-pad+j)] = (data[6+k])[i+pwidth*(pheight-2*pad+j)];
	  (data[0+k])[i+pwidth*(pheight-pad+j)] = (data[15+k])[pwidth-i+pwidth*(pad+j)];
	  (data[15+k])[i+pwidth*j] = (data[0+k])[pwidth-i+pwidth*(pheight-2*pad+j)];
	  (data[9+k])[i+pwidth*(pheight-pad+j)] = (data[15+k])[j+pad+pwidth*(i)];
	  (data[15+k])[j+pwidth*(i)] = (data[9+k])[i+pwidth*(j+pheight-2*pad)];
	  (data[12+k])[i+pwidth*(pheight-pad+j)] = (data[15+k])[pheight-2*pad+j+pwidth*(pheight-i)];
	  (data[15+k])[pheight-pad+j+pwidth*(pheight-i)] = (data[12+k])[i+pwidth*(pheight-2*pad+j)];
	  (data[0+k])[i+pwidth*(j)] = (data[3+k])[i+pwidth*(pheight-2*pad+j)];
	  (data[3+k])[i+pwidth*(pheight-pad+j)] = (data[0+k])[i+pwidth*(j+pad)];
	  (data[6+k])[i+pwidth*(j)] = (data[3+k])[pwidth-i+pwidth*(pad+j)];
	  (data[3+k])[pwidth-i+pwidth*(j)] = (data[6+k])[i+pwidth*(pad+j)];
	  (data[12+k])[i+pwidth*(j)] = (data[3+k])[j+pad+pwidth*(i)];
	  (data[3+k])[j+pwidth*(i)] = (data[12+k])[i+pwidth*(j+pad)];

	  (data[9+k])[i+pwidth*(j)] = (data[3+k])[pheight-2*pad+j+pwidth*(pwidth-i)];
	  (data[3+k])[pheight-pad+j+pwidth*(pwidth-i)] = (data[9+k])[i+pwidth*(j+pad)];
      }

    //HACKY solution to the problem with dark pixels...
    for (int k=0; k<18; k++) {
      for (int i=0; i<pad; i++) {
	(data[k])[i+pwidth*(pad-1)] = (data[0+k])[i+pwidth*(pad)];
	(data[k])[i+pwidth*(pheight-pad)] = (data[0+k])[i+pwidth*(pheight-pad-1)];
	(data[k])[pwidth-i-1+pwidth*(pad-1)] = (data[0+k])[pwidth-i-1+pwidth*(pad)];
	(data[k])[pwidth-1-i+pwidth*(pheight-pad)] = (data[0+k])[pwidth-1-i+pwidth*(pheight-pad-1)];
      }
      for (int j=0; j<pad; j++) {
	(data[k])[pad-1+pwidth*(j)] = (data[k])[pad+pwidth*(j)];
	(data[k])[pwidth-pad+pwidth*(j)] = (data[k])[pwidth-pad-1+pwidth*(j)];
	(data[k])[pad-1+pwidth*(pheight-1-j)] = (data[k])[pad+pwidth*(pheight-1-j)];
	(data[k])[pwidth-pad+pwidth*(pheight-1-j)] = (data[k])[pwidth-pad-1+pwidth*(pheight-1-j)];
      }
    }
  }
  gettimeofday (&tv2, NULL);
  t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
  printf ("Loading images & padding took %ld milliseconds\n", t);

  gettimeofday (&tv1, NULL);
  CalcTrans (width, height, alpha, beta);
  gettimeofday (&tv2, NULL);
  t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
  printf ("Calculating translations took %ld milliseconds\n", t);

  gettimeofday (&tv1, NULL);
  float* outred   =  image_warp_generic(data[0], pwidth, pheight, NULL, 'B', 0.0, 0.0);
  float* outgreen =  image_warp_generic(data[1], pwidth, pheight, NULL, 'B', 0.0, 0.0);
  float* outblue  =  image_warp_generic(data[2], pwidth, pheight, NULL, 'B', 0.0, 0.0);

  float* out1 =  image_warp_generic(data[3], pwidth, pheight, NULL, 'D', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outred[i]+=out1[i];
  delete[] out1;
  out1 =  image_warp_generic(data[4], pwidth, pheight, NULL, 'D', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outgreen[i]+=out1[i];
  delete[] out1;
  out1 =  image_warp_generic(data[5], pwidth, pheight, NULL, 'D', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outblue[i]+=out1[i];
  delete[] out1;

  out1 =  image_warp_generic(data[6], pwidth, pheight, NULL, 'F', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outred[i]+=out1[i];
  delete[] out1;
  out1 =  image_warp_generic(data[7], pwidth, pheight, NULL, 'F', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outgreen[i]+=out1[i];
  delete[] out1;
  out1 =  image_warp_generic(data[8], pwidth, pheight, NULL, 'F', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outblue[i]+=out1[i];
  delete[] out1;

  out1 =  image_warp_generic(data[9], pwidth, pheight, NULL, 'L', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outred[i]+=out1[i];
  delete[] out1;
  out1 =  image_warp_generic(data[10], pwidth, pheight, NULL, 'L', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outgreen[i]+=out1[i];
  delete[] out1;
  out1 =  image_warp_generic(data[11], pwidth, pheight, NULL, 'L', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outblue[i]+=out1[i];
  delete[] out1;

  out1 =  image_warp_generic(data[12], pwidth, pheight, NULL, 'R', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outred[i]+=out1[i];
  delete[] out1;
  out1 =  image_warp_generic(data[13], pwidth, pheight, NULL, 'R', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outgreen[i]+=out1[i];
  delete[] out1;
  out1 =  image_warp_generic(data[14], pwidth, pheight, NULL, 'R', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outblue[i]+=out1[i];
  delete[] out1;

  out1 =  image_warp_generic(data[15], pwidth, pheight, NULL, 'U', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outred[i]+=out1[i];
  delete[] out1;
  out1 =  image_warp_generic(data[16], pwidth, pheight, NULL, 'U', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outgreen[i]+=out1[i];
  delete[] out1;
  out1 =  image_warp_generic(data[17], pwidth, pheight, NULL, 'U', 0.0, 0.0);
  for (int i=0; i<outwidth*outheight; i++) outblue[i]+=out1[i];
  delete[] out1;

  gettimeofday (&tv2, NULL);
  t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
  printf ("Interpolations took %ld milliseconds\n", t);

  // Construct and save output image
  gettimeofday (&tv1, NULL);
  FIBITMAP *outbitmap = FreeImage_Allocate (outwidth, outheight, bpp);
  BYTE *outbits = FreeImage_GetBits (outbitmap);
  for (int i = 0; i<outwidth*outheight; i++) {
    outbits[3*i+FI_RGBA_RED] =   CastFloat (outred[i]);
    outbits[3*i+FI_RGBA_GREEN] = CastFloat (outgreen[i]);
    outbits[3*i+FI_RGBA_BLUE] =  CastFloat (outblue[i]);
  }
  FreeImage_Save (FIF_PNG, outbitmap, "del.png");
  FreeImage_Unload (outbitmap);
  gettimeofday (&tv2, NULL);
  t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
  printf ("Constructing and saving final image took %ld milliseconds\n", t);

  // Cleanup
  delete[] outred;
  delete[] outblue;
  delete[] outgreen;
  for (int i=0; i<18; i++)
    delete[] data[i];
  FreeImage_DeInitialise ();
  FreeTrans ();
  
  return EXIT_SUCCESS;
}

BYTE CastFloat (float f)
{
  if (f<0.0)
    return (0);
  
  if (f>255.0)
    return (255);

  return ((BYTE)roundf(f));
}  

void Help ()
{
  printf ("Try domemaster *ice*\n");
}

void SaveFits (const float *array, int nx, int ny, const char*name)
{
  fitsfile *fptr;
  int status = 0;

  long naxes[2] = {nx, ny};
  long fpixel[2] = {1, 1};
  fits_create_file  (&fptr, name, &status);
  if (status) {
    std::cout << "ERROR: Cannot create file " << name
              << " - it probably exists already" 
              << std::endl;
    return;
  }

  fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status);

  fits_write_pix(fptr, TFLOAT, fpixel, nx*ny,(void *)array, &status);
  fits_close_file (fptr, &status);
}
