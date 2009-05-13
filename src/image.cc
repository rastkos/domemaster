#include "FreeImage.h"
#include "image.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern bool discard_alpha;
extern OUTPUT_TYPE outtype;
extern Output output[4];

Image::Image (int w, int h, int b)
  : width(w), height(h), bpp(b), red(NULL), green(NULL), blue(NULL), 
    alpha(NULL), proj(UNKNOWN), pad(0)
{
  red   = new float[h*w];
  blue  = new float[h*w];
  green = new float[h*w];
  if (b==32) {
    cpi = 4;
    alpha = new float[h*w];
  }
  else
    cpi=3;

  for (int i=0;i<width*height;i++) {
    red[i] = 0.0;
    blue[i] = 0.0;
    green[i] = 0.0;
  }
  if (alpha!=NULL)
    for (int i=0;i<width*height;i++) 
      alpha[i] = 0.0;  
}

Image::Image (const char *name, int npad)
  :width(0), height(0), bpp(0), red(NULL), green(NULL), blue(NULL), 
   alpha(NULL), proj(UNKNOWN), pad(npad)
{
  imname = name;
  if (imname!=NULL) {
    FREE_IMAGE_FORMAT fif  = FreeImage_GetFileType (imname);
    if (fif == FIF_UNKNOWN) {
      printf ("File %s:Unknown filetype\nAborting\n",  imname);
      abort ();
    }
    FIBITMAP *bitmap1 = FreeImage_Load (fif, imname);
    FREE_IMAGE_TYPE imtype1 = FreeImage_GetImageType (bitmap1);
    bool fail = false;
    cpi;
    
    bpp = FreeImage_GetBPP (bitmap1);
    if (bpp==24)
      cpi = 3;
    else if (bpp==32) {
      if (!output[outtype].alpha && !discard_alpha) {
	printf ("The specified file format does not support alpha channel\n"
		"Switching to png instead\n");
	outtype = OUTPUT_PNG;
      }
      cpi = 4;
    }
    else  {
      fail = true;
      printf ("Only 24 or 32 bit images are supported\n");
    }
    width = FreeImage_GetWidth (bitmap1);
    height = FreeImage_GetHeight (bitmap1);
    int pwidth=width+2*npad;
    int pheight=height+2*npad;
    int npix_in=pwidth*pheight;
    
    BYTE *bits1 = FreeImage_GetBits (bitmap1);
    int pitch = FreeImage_GetPitch (bitmap1);
    
    red   = new float[npix_in];
    green = new float[npix_in];
    blue  = new float[npix_in];
    for (int y = 0; y<height; y++) {
      BYTE *pixel = bits1;
      for (int x = 0; x<width; x++) { 
	red[x+npad+pwidth*(y+npad)]   = pixel[FI_RGBA_RED];
	green[x+npad+pwidth*(y+npad)] = pixel[FI_RGBA_GREEN];
	blue[x+npad+pwidth*(y+npad)]  = pixel[FI_RGBA_BLUE];
	pixel += cpi;
      }
      bits1 += pitch;
    }

    if (cpi==4) {
      bits1 = FreeImage_GetBits (bitmap1);
      alpha = new float[npix_in];
      for (int y = 0; y<height; y++) {
	BYTE *pixel = bits1;
	for (int x = 0; x<width; x++) {
 	  alpha[x+npad+pwidth*(y+npad)] = pixel[FI_RGBA_ALPHA];
	  pixel += cpi;
	}
	bits1 += pitch;
      }
    }

    if (!strncmp (imname, "b_", 2))
      proj = BACK;
    else if (!strncmp (imname, "d_", 2))
      proj = DOWN;
    else if (!strncmp (imname, "f_", 2)) 
      proj = FRONT;
    else if (!strncmp (imname, "l_", 2))
      proj = LEFT;
    else if (!strncmp (imname, "r_", 2))
      proj = RIGHT;
    else if (!strncmp (imname, "u_", 2))
      proj = UP;
    
    FreeImage_Unload (bitmap1);
  }
}

Image::~Image ()
{
  if (red!=NULL)
    delete[] red;
  if (green!=NULL)
    delete[] green;
  if (blue!=NULL)
    delete[] blue;
  if (alpha!=NULL)
    delete[] alpha;
}

bool Image::operator == (const Image & d )
{
  // Comparison with empty images is ok
  if (red==NULL || d.red==NULL)
    return (true);

  if (height != d.height || width != d.width || bpp!=d.bpp)
    return (false);

  if ( (alpha!=NULL && d.alpha==NULL) ||
       (alpha==NULL && d.alpha!=NULL) )
    return (false);

  return (true);
}

bool Image::operator += (const Image & d )
{
  if (red !=NULL && d.red != NULL) {
    for (int i=0; i< (height+2*pad)*(width+2*pad); i++) {
      red[i]   += d.red[i];
      green[i] += d.green[i];
      blue[i]  += d.blue[i];
    }
    if (alpha!=NULL && d.alpha != NULL) 
      for (int i=0; i< (height+2*pad)*(width+2*pad); i++) 
	alpha[i] += d.alpha[i];
  }
}
