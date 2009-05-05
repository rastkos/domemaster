#include "FreeImage.h"
#include "image.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern bool discard_alpha;
extern OUTPUT_TYPE outtype;
extern Output output[4];

Image::Image (const char *name, int pad)
  :width(0), height(0), bpp(0), red(NULL), green(NULL), blue(NULL), 
   alpha(NULL), proj(UNKNOWN)
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
    int cpi;
    
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
    int pwidth=width+2*pad;
    int pheight=height+2*pad;
    int npix_in=pwidth*pheight;
    
    BYTE *bits1 = FreeImage_GetBits (bitmap1);
    
    red = new float[npix_in];
    for (int i = 0; i<width; i++) 
      for (int j = 0; j<height; j++) 
	red[i+pad+pwidth*(j+pad)] = bits1[cpi*(i+j*width)+FI_RGBA_RED];
    
    green = new float[npix_in];
    for (int i = 0; i<width; i++) 
      for (int j = 0; j<height; j++) 
	green[i+pad+pwidth*(j+pad)] = bits1[cpi*(i+j*width)+FI_RGBA_GREEN];
    
    blue = new float[npix_in];
    for (int i = 0; i<width; i++) 
      for (int j = 0; j<height; j++) 
	blue[i+pad+pwidth*(j+pad)] = bits1[cpi*(i+j*width)+FI_RGBA_BLUE];
    
    if (cpi==4) {
      alpha = new float[npix_in];
      for (int i = 0; i<width; i++) 
	for (int j = 0; j<height; j++) 
	  alpha[i+pad+pwidth*(j+pad)] = bits1[cpi*(i+j*width)+FI_RGBA_ALPHA];
    }
    else
      alpha = NULL;
  
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
