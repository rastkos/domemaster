#include "FreeImage.h"
#include "image.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cmath>

extern bool discard_alpha;
extern OUTPUT_TYPE outtype;
extern Output output[4];

float GimpContrastBrightness (float value, float databrightness, float datacontrast);

Image::Image (int w, int h, int b)
  : width(w), height(h), bpp(b), red(NULL), green(NULL), blue(NULL),
    alpha(NULL), saturation(NULL), brightness(NULL), hue(NULL), histogram(NULL), 
    histtrans(NULL), proj(UNKNOWN), pad(0)
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

  memset (red, 0, width*height*sizeof(float));
  memset (blue, 0, width*height*sizeof(float));
  memset (green, 0, width*height*sizeof(float));

  if (alpha!=NULL)
    memset (alpha, 0, width*height*sizeof(float));
}

Image::Image (const char *name, int npad)
  :width(0), height(0), bpp(0), red(NULL), green(NULL), blue(NULL), 
   alpha(NULL), saturation(NULL), brightness(NULL), hue(NULL), histogram(NULL), 
   histtrans(NULL), proj(UNKNOWN), pad(npad)
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

  if (saturation!=NULL)
    delete[] saturation;
  if (brightness!=NULL)
    delete[] brightness;
  if (hue!=NULL)
    delete[] hue;

  if (histogram!=NULL)
    delete[] histogram;
  if (histtrans!=NULL)
    delete[] histtrans;
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

/*
The following routine is derived from ImageMagick; adapted by J. Reunanen


Copyright 1999-2009 ImageMagick Studio LLC, a non-profit organization 
dedicated to making software imaging solutions freely available.

   Licensed under the ImageMagick License (the "License"); you may not use
   this file except in compliance with the License.  You may obtain a copy
   of the License at

     http://www.imagemagick.org/script/license.php

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
   WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
   License for the specific language governing permissions and limitations
   under the License.
*/

void Image::ToHSB ()
{
  if (hue==NULL)
    hue = new float[(height+2*pad)*(width+2*pad)];
  if (saturation==NULL)
    saturation = new float[(height+2*pad)*(width+2*pad)];
  if (brightness==NULL)
    brightness = new float[(height+2*pad)*(width+2*pad)];

  for (int i=0; i< (height+2*pad)*(width+2*pad); i++) {
    float min= (red[i] < green[i] ? red[i] : green[i]);
    if ( blue[i] < min)
      min= blue[i];
    float max = (red[i] > green[i] ? red[i] : green[i]);
    if ( blue[i] > max)
      max= blue[i];

    if (max == 0.0) {
      saturation[i] = 0.0;
      brightness[i] = 0.0;
      hue[i] = 0.0;
      continue;
    }
    float range = max - min;

    saturation[i] = range / max;
    brightness[i] = max/255.0; 

    if (range == 0.0) {
      hue[i] = 0.0;
      continue;
    }

    if ( red[i] == max)
      hue[i]=(green[i]- blue[i])/range;
    else
      if ( green[i] == max)
	hue[i]= (2.0+(blue[i]- red[i])/range);
      else
	hue[i]= (4.0+(red[i]- green[i])/range);
    hue[i]/=6.0;
    if (hue[i] < 0.0)
      hue[i]+=1.0;
  }

}

void Image::ToRGB ()
{
  for (int i=0; i< (height+2*pad)*(width+2*pad); i++) {
    if (saturation[i]>1.0)
      saturation[i]=1.0;
    if (saturation[i]<0.0)
      saturation[i]=0.0;
    if (brightness[i]>1.0)
      brightness[i]=1.0;
    if (brightness[i]<0.0)
      brightness[i]=0.0;
    
    if (saturation[i] == 0.0) {
      red[i]=255.0*brightness[i]+0.5;
      green[i]=red[i];
      blue[i]=red[i];
      continue;
    }
    double h=6.0*(hue[i]-floor(hue[i]));
    double f=h-floor((double) h);
    double p=brightness[i]*(1.0-saturation[i]);
    double q=brightness[i]*(1.0-saturation[i]*f);
    double t=brightness[i]*(1.0-(saturation[i]*(1.0-f)));
    
    switch ((int) h)  {
    case 0:
    default: {
      red[i]= 255.0*brightness[i]+0.5;
      green[i]= 255.0*t+0.5;
      blue[i]= 255.0*p+0.5;
      break;
    }
    case 1:    {
      red[i]= 255.0*q+0.5;
      green[i]= 255.0*brightness[i]+0.5;
      blue[i]= 255.0*p+0.5;
      break;
    }
    case 2:    {
      red[i]= 255.0*p+0.5;
      green[i]= 255.0*brightness[i]+0.5;
      blue[i]= 255.0*t+0.5;
      break;
    }
    case 3:    {
      red[i]= 255.0*p+0.5;
      green[i]= 255.0*q+0.5;
      blue[i]= 255.0*brightness[i]+0.5;
      break;
    }
    case 4:   {
      red[i]= 255.0*t+0.5;
      green[i]= 255.0*p+0.5;
      blue[i]= 255.0*brightness[i]+0.5;
      break;
    }
    case 5:   {
      red[i]= 255.0*brightness[i]+0.5;
      green[i]= 255.0*p+0.5;
      blue[i]= 255.0*q+0.5;
      break;
    }}
    // No need to check whether the colour values are legit, as 
    // they are checked later in main.cc

  }
}

void Image::Modulate (float bright, float contrast, float gamma)
{
  float databrightness = bright/255.0;
  float datacontrast = contrast/127.0;
  float contpow;

  if (datacontrast <0.0)
    contpow = 1.0 + datacontrast;
  else
    contpow = (datacontrast == 1.0) ? 127 : 1.0 / (1.0 - datacontrast);

  float ctrans[256];
  for (int i=0; i<256; i++) {
    ctrans[i] = (float)i;
    ctrans[i] = GimpContrastBrightness (ctrans[i]/255.0,  contpow, databrightness);
    if (ctrans[i] >= 0.0)
      ctrans[i] = 255.0*pow(ctrans[i], gamma);
  }

  for (int i=0; i< (height+2*pad)*(width+2*pad); i++) {
    int rid = (int)red[i];
    int gid = (int)green[i];
    int bid = (int)blue[i];

    if (rid>=0 & rid<255) 
      red[i]   = ctrans[rid] + (ctrans[rid+1] - ctrans[rid])*( red[i] - (float)rid);
    if (gid>=0 & gid<255) 
      green[i] = ctrans[gid] + (ctrans[gid+1] - ctrans[gid])*(green[i]- (float)gid);
    if (bid>=0 & bid<255) 
      blue[i]  = ctrans[bid] + (ctrans[bid+1] - ctrans[bid])*( blue[i]- (float)bid);    
  }
}

void Image::ModulateHSB (float sat, float bright, float contrast, float gamma)
{
  float databrightness = bright/255.0;
  float datacontrast = contrast/127.0;
  float contpow;

  if (datacontrast <0.0)
    contpow = 1.0 + datacontrast;
  else
    contpow = (datacontrast == 1.0) ? 127 : 1.0 / (1.0 - datacontrast);

  float ctrans[256];
  for (int i=0; i<256; i++) {
    ctrans[i] = (float)(i);
    ctrans[i] = GimpContrastBrightness (ctrans[i]/255.0,  contpow, databrightness);
    if (ctrans[i] >= 0.0)
      ctrans[i] = 255.0*pow(ctrans[i], gamma);
  }

  if (hue == NULL) {
    ToHSB ();
  }

  for (int i=0; i< (height+2*pad)*(width+2*pad); i++) {    
    saturation[i] = sat*saturation[i];
  }

  if (bright != 0.0 | contrast != 0.0 | gamma!= 1.0) 
    for (int i=0; i< (height+2*pad)*(width+2*pad); i++) {    

    //brightness[i] = (brightness[i]-0.5)*1.3+0.5;
    //brightness[i] = brightness[i]*1.4;
    //brightness[i] = pow(brightness[i], brg2);
    //brightness = brg2*brightness*brightness + brg1*brightness + brg0;
    //brightness = (brightness-brg2)/(brg1-brg2);

    //float lum = 255.0 - brightness;
//     float lum = brightness;
//     if (lum<brg2)
//       lum = lum/128.0*brg2;
//     else
//       lum = (255.0-(2.0*brg2-lum))*128.0/(255.0-brg2);
//     //brightness = 255.0 - lum;
//     brightness = lum;

    //brightness = asinh(brg2*brightness)/asinh(brg2);

    //     brightness[i] = (brightness[i]+0.2/1.2);
    // if (brightness[i]<0.55)
    //   brightness[i] = 0.55+(brightness[i]-0.55)*4;

      int bid = (int)(255.0*brightness[i]);
      if (bid>=0 & bid<255) 
	brightness[i]  = ctrans[bid] + (ctrans[bid+1] - ctrans[bid])*( brightness[i]- (float)bid);
  }
}


void Image::Histogram (float frac)
{
  int cumhist[256];
  if (hue==NULL) {
    printf ("ToHSB\n");
    ToHSB();
  }
  if (histogram == NULL) {
    histogram = new int[256];
    histtrans = new float[256];
  }
  for (int i=0; i<256; i++) {
    histogram[i] = 0;
    histtrans[i] = 0.0;
    cumhist[i] = 0;
  }

  for (int i=pad; i< (width+pad); i++)
    for (int j=pad; j< (height+pad); j++) {
      int n = 255.0*brightness[i+(width+2*pad)*j];
      histogram[n]++;
    }

  cumhist[0] = histogram[0];
  for (int i=1; i<256; i++)
    cumhist[i] = cumhist[i-1]+histogram[i];

  for (int i=0; i<256; i++)
    printf ("%d %d %d\n", i, histogram[i], cumhist[i]);

  int min = 0;
  int max = 255;
  while (histogram[min] == 0 && min<256) min++;
  while (histogram[max] == 0 && min>-1) max--;

  for (int i = min; i<256; i++) {
    histtrans[i] = (float)cumhist[i]/(float)cumhist[255]*255.0;      
  }

  for (int i = 0; i<256; i++) 
    histtrans[i] = frac*histtrans[i] + (1.0-frac)*(float)i;

  for (int i = 0; i< (height+2*pad)*(width+2*pad); i++) {
    float nr=255.0*brightness[i];
    if (nr<0) nr=0;
    if (nr>255) nr=255;
    
    if (nr<254){
      // Does this actually work like I think it should???
      float fr = nr - floor(nr);
      float d = histtrans[(int)nr+1] - histtrans[(int)nr];
      brightness[i] = (histtrans[(int)nr]+d*fr)/255.0; 
    }
    else 
      brightness[i] = histtrans[(int)nr]/255.0;
  }
}


void Image::HistogramEqualize (Image *image, float frac)
{
  if (histogram == NULL)
    Histogram (frac);
  if (image->hue == NULL)
    image->ToHSB ();

  for (int i = 0; i< (image->height+2*pad)*(image->width+2*pad); i++) {
    int nr=255.0*image->brightness[i];
    if (nr<0) nr=0;
    if (nr>255) nr=255;
     
    image->brightness[i] = histtrans[nr]/255.0;
  }
}

float GimpContrastBrightness (float value, float datacontrast, float databrightness)
{
  float nvalue;
  

  /* apply brightness */
  if (databrightness < 0.0)
    value = value * (1.0 + databrightness);
  else
    value = value + ((1.0 - value) * databrightness);
  
  if (value > 0.5) {
    nvalue = 1.0 - value;
    if (nvalue < 0.0)
      nvalue = 0.0;
    nvalue = 0.5 * pow (2.0 * nvalue, datacontrast);
    value = 1.0 - nvalue;
  }
  else {
    if (value < 0.0)
      value = 0.0;
    else 
      value = 0.5 * pow (2.0 * value, datacontrast);
  }
  
  return (value);
}

