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

#include <pthread.h>

#include "main.h"
#include "eclipse_types.h"
#include "../config.h"
#include "image.h"

#ifdef HAVE_GTK
#include <gdk/gdk.h>
#include <pango/pangocairo.h>
#endif

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif


/*
 *  IOThread: Load Images, Save Output. Preference always to save!
 *  ProcessThread: 
 *
 *
 */



int main (int argc, char ** argv );
void Help (void);
void Pad (Image *front, Image *back, Image *left, 
	  Image *right, Image *down, Image *up, 
	  bool padseams);
void * ProcessImages (void * thread_data);

bool verbose = false;
bool discard_alpha = false;
bool padseams = true;
bool interp = true;

int pad = 4;
OUTPUT_TYPE outtype = OUTPUT_PNG;
char buffer[64];

Output output[4] = {
  {".png", FIF_PNG, PNG_DEFAULT, true},
  {".tga", FIF_TARGA, TARGA_DEFAULT, true},
  {".bmp", FIF_BMP, BMP_SAVE_RLE, true},
  {".jpg", FIF_JPEG, 90, false}
};

// Dimensions of the output image
int outwidth  = 1408;
int outheight = 1408;

// Orientation of the dome 
float alpha = 0.0; // Rotation around z-axis
float beta = 0.0;  // Tilt of the dome
float aperture = 180.0;  // Full dome

float angle = 90.0; // Not used at the moment!
float saturation = 1.0; // Not used at the moment!
float gammacorr=1.0;
float brightcorr=0.0;
float contrastcorr=0.0;
const char *kerneltype = "sinc"; // Interpolation kernel


void * ProcessImages (void * thread_data)
{
  Image ** images = (Image **)thread_data;
  timeval tv1,tv2;
  if (verbose)
    gettimeofday (&tv1, NULL);
  
  Image *out = image_warp_generic(images[0]);
  image_warp_generic(images[1], out);
  image_warp_generic(images[2], out);
  image_warp_generic(images[3], out);
  image_warp_generic(images[4], out);
  image_warp_generic(images[5], out);

  if (verbose) {
    gettimeofday (&tv2, NULL);
    long t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
    printf ("Interpolations took %ld milliseconds\n", t);
  }

  // HACK: saturation!=1.0 
  if (saturation != 1.0) {
    if (verbose)
      gettimeofday (&tv1, NULL);
    out->ModulateHSB (saturation, brightcorr, contrastcorr, gammacorr);
    out->ToRGB (); 
    if (verbose) {
      gettimeofday (&tv2, NULL);
      long t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
      printf ("HSB colour modifications took %ld milliseconds\n", t);
    }
  }
  else if (gammacorr!=1.0 | contrastcorr!=0.0 | brightcorr!= 0.0) {
    if (verbose)
      gettimeofday (&tv1, NULL);
    out->Modulate (brightcorr, contrastcorr, gammacorr);
    //out->ToRGB (); // HACK
    if (verbose) {
      gettimeofday (&tv2, NULL);
      long t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
      printf ("RGB colour modifications took %ld milliseconds\n", t);
    }
  }

  pthread_exit((void *)out);
}


int main( int argc, char ** argv )
{
  int cpi; // Channels per image
  int n;
  int npix_in;
  int pheight;
  int pwidth;
  timeval tv1,tv2,tv3;
  long t;
  int nimages = 6;
  char op[6] = {'F', 'B', 'D', 'L', 'R', 'U'};

  if (argc == 1) {
    Help ();
    return (EXIT_SUCCESS);
  }

  while ( (n=getopt (argc, argv, "a:b:c:dg:hik:o:pr:s:t:v")) != -1) {
    switch (n) {
    case 'a':
      aperture=atof (optarg); 
      break;
    case 'b':
      brightcorr=atof (optarg); 
      break;
    case 'c':
      contrastcorr=atof (optarg); 
      break;
    case 'd':
      discard_alpha = true; 
      break;
    case 'g':
      gammacorr=atof (optarg); 
      break;
    case 'h':
      Help ();
      return (EXIT_SUCCESS);
    case 'i': // not fully implemented for now
      interp = false;
      break;
    case 'k':
      // Selection of interpolation kernel
      if (! (!strcmp(optarg, "default") 
	  || !strcmp(optarg, "tanh") 
	  || !strcmp(optarg, "sinc") 
	  || !strcmp(optarg, "sinc2")
	  || !strcmp(optarg, "lanczos")
	  || !strcmp(optarg, "hamming")
	  || !strcmp(optarg, "hann")))
	printf ("Unrecognized interpolation kernel %s\n", optarg);
      else
	kerneltype = optarg;
      break;
    case 'o':
      if (!strcasecmp(optarg, "png"))
	outtype = OUTPUT_PNG;
      else if (!strcasecmp(optarg, "targa") || !strcasecmp(optarg, "tga"))
	outtype = OUTPUT_TARGA;
      else if (!strcasecmp(optarg, "bmp"))
	outtype = OUTPUT_BMP;
      else if (!strcasecmp(optarg, "jpg") || !strcasecmp(optarg, "jpeg"))
	outtype = OUTPUT_JPEG;

      else
	printf ("Unrecognized file type - using PNG\n");
      break;
    case 'p':
      padseams = false; 
      break;
    case 'r':
      alpha=atof (optarg);
      break;
    case 's':
      outwidth=atoi (optarg);
      outheight = outwidth;
      break;
    case 't':
      beta=atof (optarg);
      break;
    case 'v':
      verbose = true;
      break;
    }
  }

  // Fill vectors with file names
  std::vector<const char*> bfiles, dfiles, ffiles, lfiles, rfiles, ufiles;
  if (optind < argc ) {
    for (int i=optind; i<argc; i++) {
      if (!strncmp (argv[i], "b_", 2))
	bfiles.push_back (argv[i]);
      else if (!strncmp (argv[i], "d_", 2))
	dfiles.push_back (argv[i]);
      else if (!strncmp (argv[i], "f_", 2)) 
	ffiles.push_back (argv[i]);
      else if (!strncmp (argv[i], "l_", 2))
	lfiles.push_back (argv[i]);
      else if (!strncmp (argv[i], "r_", 2))
	rfiles.push_back (argv[i]);
      else if (!strncmp (argv[i], "u_", 2))
	ufiles.push_back (argv[i]);
    }
  }

  unsigned int nn = ffiles.size();

  if (nn==0) {
    printf ("No f-files given!\n");
    return (EXIT_FAILURE);
  }

  FreeImage_Initialise ();
  int oldwidth = 0;
  int oldheight = 0;

#ifdef HAVE_GTK
  cairo_surface_t *surface 
    //    = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, outwidth, outheight);
    = cairo_image_surface_create (CAIRO_FORMAT_A8, outwidth, outheight);
  cairo_t *cr = cairo_create (surface);
  cairo_set_source_rgba (cr, 0.0, 0.0, 0.0, 0.0); //BACKGROUND
  cairo_paint (cr);

  char b[32];
  PangoLayout *layout = pango_cairo_create_layout (cr);
  pango_layout_set_alignment (layout, PANGO_ALIGN_LEFT);
  pango_cairo_update_layout (cr, layout);
  
  int fsize = 20*outheight/1408;
  sprintf (&b[0], "Sans Bold %d", fsize);
  PangoFontDescription *desc = pango_font_description_from_string (&b[0]);
  pango_layout_set_font_description (layout, desc);
  pango_font_description_free (desc);
  cairo_set_source_rgb (cr, 1, 0, 0);
  
  sprintf (&b[0], "K:%s", kerneltype);
  pango_layout_set_text (layout, &b[0], -1);
  cairo_move_to (cr, fsize*1.5, outheight-fsize*9);
  pango_cairo_show_layout (cr, layout);
  
  sprintf (&b[0], "A:%.1f", aperture);
  pango_layout_set_text (layout, &b[0], -1);
  cairo_move_to (cr, fsize*1.5, outheight-fsize*7.5);
  pango_cairo_show_layout (cr, layout);
  
  sprintf (&b[0], "T:%.2f R:%.1f", beta, alpha);
  pango_layout_set_text (layout, &b[0], -1);
  cairo_move_to (cr, fsize*1.5, outheight-fsize*6);
  pango_cairo_show_layout (cr, layout);
  
  if (saturation != 1.0)
    sprintf (&b[0], "S:%.2f (HSB)", saturation);
  else
    sprintf (&b[0], "S:%.2f (RGB)", saturation);
  pango_layout_set_text (layout, &b[0], -1);
  cairo_move_to (cr, fsize*1.5, outheight-fsize*4.25);
  pango_cairo_show_layout (cr, layout);
  
  sprintf (&b[0], "G:%.3f B:%.1f C:%.1f", gammacorr, brightcorr, contrastcorr);
  pango_layout_set_text (layout, &b[0], -1);
  cairo_move_to (cr, fsize*1.5, outheight-fsize*2.5);
  pango_cairo_show_layout (cr, layout);

    
  int fsize1 = 30*outheight/1408;
  sprintf (&b[0], "Sans Bold %d", fsize1);
  desc = pango_font_description_from_string (&b[0]);
  pango_layout_set_font_description (layout, desc);
  pango_font_description_free (desc);

  cairo_rectangle (cr, outwidth*0.5, outheight*0.8, outwidth*0.5, outheight*0.2);
  //cairo_clip(cr);
#endif

  for (unsigned index=0; index<ffiles.size(); index++) {
    if (verbose)
      gettimeofday (&tv1, NULL);

    Image *front = new Image(ffiles[index], pad);
    Image *back;
    Image *left;
    Image *right;
    Image *up;
    Image *down;
    
    //front->Histogram(1.0);
    //front->ToRGB();

    if (index<bfiles.size())
      back = new Image (bfiles[index], pad);
    else
      back = new Image (NULL, pad);

    if (index<lfiles.size())
      left = new Image (lfiles[index], pad);
    else
      left = new Image (NULL, pad);

    if (index<rfiles.size())
      right = new Image (rfiles[index], pad);
    else
      right = new Image (NULL, pad);

    if (index<dfiles.size())
      down = new Image (dfiles[index], pad);
    else
      down = new Image (NULL, pad);

    if (index<ufiles.size())
      up = new Image (ufiles[index], pad);
    else
      up = new Image (NULL, pad);

    if ( ! (*front == *back && *front == *left
	    && *front == *right && *front == *up && *front == *down ) )
      printf ("OOH\n");

    Pad (front, back, left, right, down, up, padseams);

    if (verbose) {
      gettimeofday (&tv2, NULL);
      t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
      printf ("Loading images & padding took %ld milliseconds\n", t);
    }

    if (front->width!=oldwidth ||front->height!=oldheight ) {
      oldwidth = front->width;
      oldheight = front->height;
      if (verbose) 
	gettimeofday (&tv1, NULL);
      if (index != 0)
	FreeTrans ();
      CalcTrans (front->width, front->height, alpha, beta, kerneltype);
      if (verbose) {
	gettimeofday (&tv2, NULL);
	t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
	printf ("Calculating translations took %ld milliseconds\n", t);
	
	gettimeofday (&tv1, NULL);
      }
    }

    Image *images[6];
    images[0] = front;
    images[1] = back;
    images[2] = left;
    images[3] = right;
    images[4] = up;
    images[5] = down;

    int rc;
    pthread_t threads;
    rc = pthread_create(&threads, NULL, ProcessImages, (void *)&images[0]);

    void *status;
    rc = pthread_join(threads, &status);
    Image *out = (Image *)status;

    cpi = front->cpi;

    if (verbose)
      gettimeofday (&tv1, NULL);
#ifdef HAVE_GTK
   
    cairo_set_source_rgba (cr, 0.0, 0.0, 0.0, 0.0);
    cairo_set_operator (cr, CAIRO_OPERATOR_CLEAR);
  
    //pango_cairo_update_layout (cr, layout);
    cairo_fill_preserve (cr);
    //cairo_paint (cr);
    cairo_set_operator (cr, CAIRO_OPERATOR_OVER);
    
    cairo_set_source_rgba (cr, 1.0, 0.0, 0.0, 1.0);
    char *dot = strrchr (front->imname, '.');
    int strlength;
    if (dot != NULL)
      strlength = strlen (front->imname) - strlen(dot) - 2;
    else
      strlength = -1;

    pango_layout_set_text (layout, front->imname+2, strlength);
    cairo_move_to (cr, outwidth*0.84, outheight*0.93);
    pango_cairo_show_layout (cr, layout);

    unsigned char*buf = cairo_image_surface_get_data (surface);
    //unsigned long*lbuf = (unsigned long*)buf;
    for (int i=0; i<outwidth*outheight; i++) {
      if (buf[i]!=0) {
	int indx = i % outwidth;
	int indy = i / outwidth;
	int ind = indx + (outheight-indy-1)*outwidth;
	out->red[ind]=buf[i];
	out->green[ind]=buf[i];
	out->blue[ind]=buf[i];
      }
    }
#endif
    
    // Construct and save output image
    // SaveImage (const Image *out, const Image *front, bool)
    FIBITMAP *outbitmap ;
    if (front->bpp==32 && discard_alpha) {
      outbitmap = FreeImage_Allocate (outwidth, outheight, 24);
      BYTE *outbits = FreeImage_GetBits (outbitmap);
      int pitch = FreeImage_GetPitch (outbitmap);
      for (int y = 0; y<outheight; y++) {
	BYTE *pixel = outbits;
	for (int x=0; x<outwidth; x++) {
	  pixel[FI_RGBA_RED] =   CastFloat (out->red[x+y*outwidth]);
	  pixel[FI_RGBA_GREEN] = CastFloat (out->green[x+y*outwidth]);
	  pixel[FI_RGBA_BLUE] =  CastFloat (out->blue[x+y*outwidth]);
	  pixel += 3;
	}
	outbits += pitch;
      }
    } else {
      outbitmap = FreeImage_Allocate (outwidth, outheight, front->bpp);
      BYTE *outbits = FreeImage_GetBits (outbitmap);
      int pitch = FreeImage_GetPitch (outbitmap);
      for (int y = 0; y<outheight; y++) {
	BYTE *pixel = outbits;
	for (int x=0; x<outwidth; x++) {
	  pixel[FI_RGBA_RED] =   CastFloat (out->red[x+y*outwidth]);
	  pixel[FI_RGBA_GREEN] = CastFloat (out->green[x+y*outwidth]);
	  pixel[FI_RGBA_BLUE] =  CastFloat (out->blue[x+y*outwidth]);
	  if (cpi == 4) 
	    pixel[FI_RGBA_ALPHA] = CastFloat (out->alpha[x+y*outwidth]);
	  pixel += cpi;
	}
	outbits += pitch;
      }
    }
    int ss = strlen(front->imname)-1;
   
    //FreeImage_AdjustColors (outbitmap, brightness, contrast, gamma, FALSE);
    // FreeImage_AdjustGamma (outbitmap, 1.3);
    // FreeImage_AdjustBrightness (outbitmap, 0);
    // FreeImage_AdjustContrast (outbitmap, 50);

    while (ss>0) {
      if (front->imname[ss]=='.')
	break;
      ss--;
    }
    for (int i=2; i<ss; i++)
      buffer[i-2] = front->imname[i];
    buffer[ss-2] = 0;
    strcat (&buffer[0], output[outtype].name);

    FreeImage_Save (output[outtype].fif, outbitmap, &buffer[0], output[outtype].flags);
    FreeImage_Unload (outbitmap);
    if (verbose) {
      gettimeofday (&tv2, NULL);
      t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
      printf ("Constructing and saving final image took %ld milliseconds\n", t);
    }

    delete out;
    delete front;
    delete back;
    delete right;
    delete left;
    delete up;
    delete down;
  } 

  FreeImage_DeInitialise ();
  FreeTrans ();
#ifdef HAVE_GTK
  g_object_unref (layout);
  cairo_destroy (cr);
  cairo_surface_destroy (surface);
#endif
  
  return EXIT_SUCCESS;
}

void Pad (Image *front, Image *back, Image *left, Image *right, Image *down, Image *up, 
	  bool padseams)
{
  int pad = front->pad;
  int width = front->width;
  int height = front->height;
  int pwidth = width+2*pad;
  int pheight = height+2*pad;
  float *data[24];
  int cpi = front->cpi;

  data[0] = front->red;
  data[1] = front->green;
  data[2] = front->blue;
  data[3] = front->alpha;

  data[4] = back->red;
  data[5] = back->green;
  data[6] = back->blue;
  data[7] = back->alpha;

  data[8]  = down->red;
  data[9]  = down->green;
  data[10] = down->blue;
  data[11] = down->alpha;

  data[12] = left->red;
  data[13] = left->green;
  data[14] = left->blue;
  data[15] = left->alpha;

  data[16] = right->red;
  data[17] = right->green;
  data[18] = right->blue;
  data[19] = right->alpha;

  data[20] = up->red;
  data[21] = up->green;
  data[22] = up->blue;
  data[23] = up->alpha;
  if (padseams && pad>0) {
    for (int k=0; k<cpi; k++)
      for (int i=0; i<pad; i++)
	for (int j=pad; j<pad+height; j++) {
	  if (data[0]!=NULL && data[12]!=NULL) {
	    (data[0+k])[i+pwidth*j] = (data[12+k])[pwidth-2*pad+i+pwidth*j]; 
	    (data[12+k])[pwidth-pad+i+pwidth*j] = (data[0+k])[pad+i+pwidth*j];
	  }
	  if (data[4]!=NULL && data[16]!=NULL) {
	    (data[4+k])[i+pwidth*j] = (data[16+k])[pwidth-2*pad+i+pwidth*j];
	    (data[16+k])[pwidth-pad+i+pwidth*j] = (data[4+k])[pad+i+pwidth*j];
	  }
	  if (data[0]!=NULL && data[16]!=NULL) {
	    (data[0+k])[pwidth-pad+i+pwidth*j] = (data[16+k])[pad+i+pwidth*j];
	    (data[16+k])[i+pwidth*j] = (data[0+k])[pwidth-2*pad+i+pwidth*j];
	  }
	  if (data[4]!=NULL && data[12]!=NULL) {
	    (data[4+k])[pwidth-pad+i+pwidth*j] = (data[12+k])[pad+i+pwidth*j];
	    (data[12+k])[i+pwidth*j] = (data[4+k])[pwidth-2*pad+i+pwidth*j];
	  }
	}
    for (int k=0; k<cpi; k++)
      for (int j=0; j<pad; j++)
	for (int i=pad; i<pad+width; i++) {
	  if (data[0]!=NULL && data[20]!=NULL) {
	    (data[0+k])[i+pwidth*(pheight-pad+j)] = (data[20+k])[i+pwidth*(pheight-2*pad+j)];
	    (data[20+k])[i+pwidth*(pheight-pad+j)] = (data[0+k])[i+pwidth*(pheight-2*pad+j)];
	  }
	  if (data[4]!=NULL && data[20]!=NULL) {
	    (data[4+k])[i+pwidth*(pheight-pad+j)] = (data[20+k])[pwidth-i+pwidth*(pad+j)];
	    (data[20+k])[i+pwidth*j] = (data[4+k])[pwidth-i+pwidth*(pheight-2*pad+j)];
	  }
	  if (data[12]!=NULL && data[20]!=NULL) {
	    (data[12+k])[i+pwidth*(pheight-pad+j)] = (data[20+k])[j+pad+pwidth*(i)];
	    (data[20+k])[j+pwidth*(i)] = (data[12+k])[i+pwidth*(j+pheight-2*pad)];
	  }
	  if (data[16]!=NULL && data[20]!=NULL) {
	    (data[16+k])[i+pwidth*(pheight-pad+j)] = (data[20+k])[pheight-2*pad+j+pwidth*(pheight-i)];
	    (data[20+k])[pheight-pad+j+pwidth*(pheight-i)] = (data[16+k])[i+pwidth*(pheight-2*pad+j)];
	  }
 	  if (data[4]!=NULL && data[8]!=NULL) {
 	    (data[4+k])[i+pwidth*(j)] = (data[8+k])[pwidth-i+pwidth*(pad+j)];
 	    (data[8+k])[pwidth-i+pwidth*(j)] = (data[4+k])[i+pwidth*(j+pad)];
 	  }
	  if (data[0]!=NULL && data[8]!=NULL) {
	    (data[0+k])[i+pwidth*(j)] = (data[8+k])[i+pwidth*(pheight-2*pad+j)];
	    (data[8+k])[pwidth-i+pwidth*(j)] = (data[0+k])[i+pwidth*(pad+j)];
	  }
	  if (data[16]!=NULL && data[8]!=NULL) {
	    (data[16+k])[i+pwidth*(j)] = (data[8+k])[pheight-2*pad+j+pwidth*(pwidth-i)];
	    (data[8+k])[j+pwidth*(i)] = (data[16+k])[i+pwidth*(j+pad)];	 
	  } 
	  if (data[12]!=NULL && data[8]!=NULL) {
	    (data[12+k])[i+pwidth*(j)] = (data[8+k])[j+pad+pwidth*(i)];
	    (data[8+k])[pheight-pad+j+pwidth*(pwidth-i)] = (data[12+k])[i+pwidth*(j+pad)];
	  }
	}
    
    //HACKY solution to the problem with dark pixels...
    for (int k=0; k<24; k++) {
      if (data[k]!=NULL && !(cpi==3 && k%4==3) ) {
	for (int i=0; i<pad; i++) {
	  (data[k])[i+pwidth*(pad-1)] = (data[k])[i+pwidth*(pad)];
	  (data[k])[i+pwidth*(pheight-pad)] = (data[k])[i+pwidth*(pheight-pad-1)];
	  (data[k])[pwidth-i-1+pwidth*(pad-1)] = (data[k])[pwidth-i-1+pwidth*(pad)];
	  (data[k])[pwidth-1-i+pwidth*(pheight-pad)] = (data[k])[pwidth-1-i+pwidth*(pheight-pad-1)];
	}
	for (int j=0; j<pad; j++) {
	  (data[k])[pad-1+pwidth*(j)] = (data[k])[pad+pwidth*(j)];
	  (data[k])[pwidth-pad+pwidth*(j)] = (data[k])[pwidth-pad-1+pwidth*(j)];
	  (data[k])[pad-1+pwidth*(pheight-1-j)] = (data[k])[pad+pwidth*(pheight-1-j)];
	  (data[k])[pwidth-pad+pwidth*(pheight-1-j)] = (data[k])[pwidth-pad-1+pwidth*(pheight-1-j)];
	}
      }
      if (data[k]!=NULL && pad>=2) {
	if (!(cpi==3 && k%4==3) ) {
	  for (int i=0; i<pad; i++) {
	    (data[k])[i+pwidth*(pad-2)] = (data[k])[i+pwidth*(pad)];
	    (data[k])[i+pwidth*(pheight-pad+1)] = (data[k])[i+pwidth*(pheight-pad-1)];
	    (data[k])[pwidth-i-1+pwidth*(pad-2)] = (data[k])[pwidth-i-1+pwidth*(pad)];
	    (data[k])[pwidth-1-i+pwidth*(pheight-pad+1)] = (data[k])[pwidth-1-i+pwidth*(pheight-pad-1)];
	  }
	  for (int j=0; j<pad; j++) {
	    (data[k])[pad-2+pwidth*(j)] = (data[k])[pad+pwidth*(j)];
	    (data[k])[pwidth-pad+1+pwidth*(j)] = (data[k])[pwidth-pad-1+pwidth*(j)];
	    (data[k])[pad-2+pwidth*(pheight-1-j)] = (data[k])[pad+pwidth*(pheight-1-j)];
	    (data[k])[pwidth-pad+1+pwidth*(pheight-1-j)] = (data[k])[pwidth-pad-1+pwidth*(pheight-1-j)];
	  }
	} 
      }
      if (data[k]!=NULL && pad>=3) {
	if (!(cpi==3 && k%4==3) ) {
	  for (int i=0; i<pad; i++) {
	    (data[k])[i+pwidth*(pad-3)] = (data[k])[i+pwidth*(pad)];
	    (data[k])[i+pwidth*(pheight-pad+2)] = (data[k])[i+pwidth*(pheight-pad-1)];
	    (data[k])[pwidth-i-1+pwidth*(pad-3)] = (data[k])[pwidth-i-1+pwidth*(pad)];
	    (data[k])[pwidth-1-i+pwidth*(pheight-pad+2)] = (data[k])[pwidth-1-i+pwidth*(pheight-pad-1)];
	  }
	  for (int j=0; j<pad; j++) {
	    (data[k])[pad-3+pwidth*(j)] = (data[k])[pad+pwidth*(j)];
	    (data[k])[pwidth-pad+2+pwidth*(j)] = (data[k])[pwidth-pad-1+pwidth*(j)];
	    (data[k])[pad-3+pwidth*(pheight-1-j)] = (data[k])[pad+pwidth*(pheight-1-j)];
	    (data[k])[pwidth-pad+2+pwidth*(pheight-1-j)] = (data[k])[pwidth-pad-1+pwidth*(pheight-1-j)];
	  }
	} 
      }
    }
  }
}

void Help ()
{
  const char *fiver=FreeImage_GetVersion();
  const char *domever = VERSION ;
  printf ("Domemaster %s using FreeImage %s\n\nUSAGE:\n", domever, fiver);
  printf ("  domemaster [options] images\n\nAvailable options:\n");
  printf ("  -a <angle>  Aperture, fulldome is 180 degrees [180]\n");
  printf ("  -b <value>  Modify brightness\n");
  printf ("  -c <value>  Modify contrast\n");
  printf ("  -d          Disable alpha channel\n");
  printf ("  -g <value>  Ably gamma correction to the image\n");
  printf ("  -h          Print this help message\n");
  printf ("  -i          Turn interpolation off\n");
  printf ("  -k <kernel> Select interpolation kernel. Possible values are\n"
	  "              tanh, sinc, sinc2, lanczos, hamming, hann. [sinc]\n");
  printf ("  -o <format> Output file format, possible values png, bmp, targa, jpeg. [png]\n");
  printf ("  -p          Disable padding\n");
  printf ("  -r <angle>  Rotation angle [0.0]\n");
  printf ("  -s <size>   Size of the output image. [1408]\n");
  printf ("  -t <angle>  Tilt angle [0.0]\n");
  printf ("  -v          Verbose messages\n");
}

