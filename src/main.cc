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

#include "main.h"
#include "eclipse_types.h"
#include "../config.h"

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

int main( int argc, char ** argv );
void Help (void);

float *data[24];
bool verbose = false;
bool discard_alpha = false;
bool padseams = true;
int pad = 3;
OUTPUT_TYPE outtype = OUTPUT_PNG;
char buffer[64];
const char *root = "del";

//enum OUTPUT_TYPE  { OUTPUT_PNG=0, OUTPUT_TARGA, OUTPUT_BMP, OUTPUT_JPEG };

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

float angle = 90.0; // Not used at the moment!

int main( int argc, char ** argv )
{
  const char *kernel = NULL; // Interpolation kernel
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

  while ( (n=getopt (argc, argv, "dhk:o:pr:s:t:v")) != -1) {
    switch (n) {
    case 'd':
      discard_alpha = true; 
      break;
    case 'h':
      Help ();
      return (EXIT_SUCCESS);
    case 'i': // not fully implemented for now
      nimages = strlen (optarg);
      if (nimages>6) {
	printf ("Too many images given in -i option\n");
	abort ();
      }
      for (int i=0; i<strlen(optarg); i++) {
	bool fail=false;
	char c=toupper(optarg[i]);
	switch (optarg[i]) {
	case 'B':
	case 'D':
	case 'F':
	case 'L':
	case 'R':
	case 'U':
	  for (int j=0; j<i; j++)
	    if (op[j] == c)
	      fail = true;
	  if (!fail)
	    op[i] = c;
	  break;
	default:
	  fail = true;
	  printf ("Unable to parse imagestring\n");
	  abort ();
	  break;
	}
      }
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
	kernel = optarg;
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
  std::vector<const char*> bfiles, dfiles, ffiles, lfiles, rfiles, ufiles, roots;
  if (optind < argc ) {
    for (int i=optind; i<argc; i++) {
      if (!strncmp (argv[i], "b_", 2))
	bfiles.push_back (argv[i]);
      else if (!strncmp (argv[i], "d_", 2))
	dfiles.push_back (argv[i]);
      else if (!strncmp (argv[i], "f_", 2)) {
	ffiles.push_back (argv[i]);
	roots.push_back( &(argv[i])[2]);
      }
      else if (!strncmp (argv[i], "l_", 2))
	lfiles.push_back (argv[i]);
      else if (!strncmp (argv[i], "r_", 2))
	rfiles.push_back (argv[i]);
      else if (!strncmp (argv[i], "u_", 2))
	ufiles.push_back (argv[i]);
    }
  }

  //
  unsigned int nn = roots.size();

  if (nn==0) {
    printf ("No f-files given!\n");
    return (EXIT_FAILURE);

  if (bfiles.size()==1) {
    printf ("Using %s as b-file for all\n", bfiles[0]);
    for (unsigned int i=1; i<nn; i++)
      bfiles.push_back (bfiles[0]);
  }
  else if (bfiles.size()!=nn && bfiles.size()!=0) 
    printf ("Cannot proceed: incorrect number of b-files\n");


  }
  for (unsigned int i=0; i<roots.size(); i++) 
    printf ("%s\n", roots[i]);


  

  FreeImage_Initialise ();

  
  //  for (unsigned index=0; index<roots.size(); index++) {

    int width,height, bpp;
    gettimeofday (&tv1, NULL);
    for (int image = 0; image<nimages; image++) {
      const char *imname = argv[optind+image];
      printf ("%s\n", imname);
      
      FREE_IMAGE_FORMAT fif  = FreeImage_GetFileType (imname);
      if (fif == FIF_UNKNOWN) {
	printf ("File %s:Unknown filetype\nAborting\n",  imname);
	abort ();
      }
      FIBITMAP *bitmap1 = FreeImage_Load (fif, imname);
      FREE_IMAGE_TYPE imtype1 = FreeImage_GetImageType (bitmap1);
      bool fail = false;
      
      if (image == 0) {
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
	else {
	  fail = true;
	  printf ("Only 24 or 32 bit images are supported\n");
	}
	width = FreeImage_GetWidth (bitmap1);
	height = FreeImage_GetHeight (bitmap1);
	if (width != height) {
	  fail = true;
	  printf ("Only square images are supported - %dx%d\n", width, height);
	}
	npix_in = (width+2*pad)*(height+2*pad);
	pwidth = width+2*pad;
	pheight = height+2*pad;
      }
      else {
	int bpp1 = FreeImage_GetBPP (bitmap1);
	int width1 = FreeImage_GetWidth (bitmap1);
	int height1 = FreeImage_GetHeight (bitmap1);
	if (bpp!=bpp1 || width!=width1 || height!=height1) {
	  printf ("Images do not have identical sizes or bpps\n");
	  fail = true;
	}
      }
      if (fail) {
	FreeImage_Unload (bitmap1);
	FreeImage_DeInitialise ();
	for (int i=0; i<cpi*image; i++)
	  if (data[i]!=NULL)
	    delete[] data[i];
	abort ();
      }
      
      BYTE *bits1 = FreeImage_GetBits (bitmap1);
      
      data[4*image] = new float[npix_in];
      for (int i = 0; i<width; i++) 
	for (int j = 0; j<height; j++) 
	  (data[4*image])[i+pad+pwidth*(j+pad)] = bits1[cpi*(i+j*width)+FI_RGBA_RED];
      
      data[4*image+1] = new float[npix_in];
      for (int i = 0; i<width; i++) 
	for (int j = 0; j<height; j++) 
	  (data[4*image+1])[i+pad+pwidth*(j+pad)] = bits1[cpi*(i+j*width)+FI_RGBA_GREEN];
      
      data[4*image+2] = new float[npix_in];
      for (int i = 0; i<width; i++) 
	for (int j = 0; j<height; j++) 
	  (data[4*image+2])[i+pad+pwidth*(j+pad)] = bits1[cpi*(i+j*width)+FI_RGBA_BLUE];
      
      if (cpi==4) {
	data[4*image+3] = new float[npix_in];
	for (int i = 0; i<width; i++) 
	  for (int j = 0; j<height; j++) 
	    (data[cpi*image+3])[i+pad+pwidth*(j+pad)] = bits1[cpi*(i+j*width)+FI_RGBA_ALPHA];
      }
      else
	data[4*image+3] = NULL;
      
      FreeImage_Unload (bitmap1);
    }
    // printf ("Images loaded successfully\n");
    
    // WAS: B, D, F, L, R, U
    // NOW: F, B, D, L, R, U
    // SO: 0+k -> 4+k
    //     4+k -> 8+k
    //     8+k -> 0+k
    if (pad>0 && padseams) {
      for (int k=0; k<cpi; k++)
	for (int i=0; i<pad; i++)
	  for (int j=pad; j<pad+height; j++) {
	    (data[0+k])[i+pwidth*j] = (data[12+k])[pwidth-2*pad+i+pwidth*j]; 
	    (data[4+k])[i+pwidth*j] = (data[16+k])[pwidth-2*pad+i+pwidth*j];
	    (data[0+k])[pwidth-pad+i+pwidth*j] = (data[16+k])[pad+i+pwidth*j];
	    (data[4+k])[pwidth-pad+i+pwidth*j] = (data[12+k])[pad+i+pwidth*j];
	    (data[12+k])[pwidth-pad+i+pwidth*j] = (data[0+k])[pad+i+pwidth*j];
	    (data[12+k])[i+pwidth*j] = (data[4+k])[pwidth-2*pad+i+pwidth*j];
	    (data[16+k])[i+pwidth*j] = (data[0+k])[pwidth-2*pad+i+pwidth*j];
	    (data[16+k])[pwidth-pad+i+pwidth*j] = (data[4+k])[pad+i+pwidth*j];
	  }
      for (int k=0; k<cpi; k++)
	for (int j=0; j<pad; j++)
	  for (int i=pad; i<pad+width; i++) {
	    (data[0+k])[i+pwidth*(pheight-pad+j)] = (data[20+k])[i+pwidth*(pheight-2*pad+j)];
	    (data[20+k])[i+pwidth*(pheight-pad+j)] = (data[0+k])[i+pwidth*(pheight-2*pad+j)];
	    (data[4+k])[i+pwidth*(pheight-pad+j)] = (data[20+k])[pwidth-i+pwidth*(pad+j)];
	    (data[20+k])[i+pwidth*j] = (data[4+k])[pwidth-i+pwidth*(pheight-2*pad+j)];
	    (data[12+k])[i+pwidth*(pheight-pad+j)] = (data[20+k])[j+pad+pwidth*(i)];
	    (data[20+k])[j+pwidth*(i)] = (data[12+k])[i+pwidth*(j+pheight-2*pad)];
	    (data[16+k])[i+pwidth*(pheight-pad+j)] = (data[20+k])[pheight-2*pad+j+pwidth*(pheight-i)];
	    (data[20+k])[pheight-pad+j+pwidth*(pheight-i)] = (data[16+k])[i+pwidth*(pheight-2*pad+j)];
	    (data[4+k])[i+pwidth*(j)] = (data[8+k])[pwidth-i+pwidth*(pad+j)];
	    (data[8+k])[i+pwidth*(pheight-pad+j)] = (data[4+k])[i+pwidth*(j+pad)];
	    (data[0+k])[i+pwidth*(j)] = (data[8+k])[i+pwidth*(pheight-2*pad+j)];
	    (data[8+k])[pwidth-i+pwidth*(j)] = (data[0+k])[i+pwidth*(pad+j)];
	    (data[16+k])[i+pwidth*(j)] = (data[8+k])[pheight-2*pad+j+pwidth*(pwidth-i)];
	    (data[8+k])[j+pwidth*(i)] = (data[16+k])[i+pwidth*(j+pad)];	  
	    (data[12+k])[i+pwidth*(j)] = (data[8+k])[j+pad+pwidth*(i)];
	    (data[8+k])[pheight-pad+j+pwidth*(pwidth-i)] = (data[12+k])[i+pwidth*(j+pad)];
	  }
      
      //HACKY solution to the problem with dark pixels...
      for (int k=0; k<24; k++) {
	if (!(cpi==3 && k%4==3) ) {
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
    }

    //if (index==0) {
      if (verbose) {
	gettimeofday (&tv2, NULL);
	t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
	printf ("Loading images & padding took %ld milliseconds\n", t);
	
	gettimeofday (&tv1, NULL);
      }
      CalcTrans (width, height, alpha, beta);
      if (verbose) {
	gettimeofday (&tv2, NULL);
	t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
	printf ("Calculating translations took %ld milliseconds\n", t);
	
	gettimeofday (&tv1, NULL);
      }
      //}

    float * outalpha;
    gettimeofday (&tv2, NULL);
    float* outred   =  image_warp_generic(data[0], pwidth, pheight, kernel, op[0], 0.0, 0.0);
    float* outgreen =  image_warp_generic(data[1], pwidth, pheight, kernel, op[0], 0.0, 0.0);
    float* outblue  =  image_warp_generic(data[2], pwidth, pheight, kernel, op[0], 0.0, 0.0);
    if (cpi==4) 
      outalpha = image_warp_generic(data[3], pwidth, pheight, kernel, op[0], 0.0, 0.0);
    
    float *out1;
    for (int j=1; j<nimages; j++) {
      out1 =  image_warp_generic(data[4*j], pwidth, pheight, kernel, op[j], 0.0, 0.0);
      for (int i=0; i<outwidth*outheight; i++) outred[i]+=out1[i];
      delete[] out1;
      out1 =  image_warp_generic(data[4*j+1], pwidth, pheight, kernel, op[j], 0.0, 0.0);
      for (int i=0; i<outwidth*outheight; i++) outgreen[i]+=out1[i];
      delete[] out1;
      out1 =  image_warp_generic(data[4*j+2], pwidth, pheight, kernel, op[j], 0.0, 0.0);
      for (int i=0; i<outwidth*outheight; i++) outblue[i]+=out1[i];
      delete[] out1;
      if (cpi==4) {
	out1 =  image_warp_generic(data[4*j+3], pwidth, pheight, kernel, op[j], 0.0, 0.0);
	for (int i=0; i<outwidth*outheight; i++) outalpha[i]+=out1[i];
	delete[] out1;
      }
    }
    
    if (verbose) {
      gettimeofday (&tv2, NULL);
      t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
      printf ("Interpolations took %ld milliseconds\n", t);
    }
    
    // Construct and save output image
    gettimeofday (&tv1, NULL);
    FIBITMAP *outbitmap ;
    if (bpp==32 && discard_alpha) {
      outbitmap = FreeImage_Allocate (outwidth, outheight, 24);
      BYTE *outbits = FreeImage_GetBits (outbitmap);
      for (int i = 0; i<outwidth*outheight; i++) {
	outbits[3*i+FI_RGBA_RED] =   CastFloat (outred[i]);
	outbits[3*i+FI_RGBA_GREEN] = CastFloat (outgreen[i]);
	outbits[3*i+FI_RGBA_BLUE] =  CastFloat (outblue[i]);
      }
    } else {
      outbitmap = FreeImage_Allocate (outwidth, outheight, bpp);
      BYTE *outbits = FreeImage_GetBits (outbitmap);
      for (int i = 0; i<outwidth*outheight; i++) {
	outbits[cpi*i+FI_RGBA_RED] =   CastFloat (outred[i]);
	outbits[cpi*i+FI_RGBA_GREEN] = CastFloat (outgreen[i]);
	outbits[cpi*i+FI_RGBA_BLUE] =  CastFloat (outblue[i]);
	if (cpi == 4) 
	  outbits[cpi*i+FI_RGBA_ALPHA] =   CastFloat (outalpha[i]);
	
      }
    }
    sprintf (&buffer[0], "%s%s", root, output[outtype].name);
    FreeImage_Save (output[outtype].fif, outbitmap, &buffer[0], output[outtype].flags);
    FreeImage_Unload (outbitmap);
    if (verbose) {
      gettimeofday (&tv2, NULL);
      t = (tv2.tv_sec - tv1.tv_sec)*1000 +(tv2.tv_usec - tv1.tv_usec)/1000;
      printf ("Constructing and saving final image took %ld milliseconds\n", t);
    }
    //}

  // Cleanup
  delete[] outred;
  delete[] outblue;
  delete[] outgreen;
  if (cpi==4)
    delete[] outalpha;
  for (int i=0; i<18; i++)
    if (data[i] != NULL)
      delete[] data[i];
  FreeImage_DeInitialise ();
  FreeTrans ();
  
  return EXIT_SUCCESS;
}

void Help ()
{
  const char *fiver=FreeImage_GetVersion();
  const char *domever = VERSION ;
  printf ("Domemaster %s using FreeImage %s\n\nUSAGE:\n", domever, fiver);
  printf ("  domemaster [options] images\n\nAvailable options:\n");
  printf ("  -d          Disable alpha channel\n");
  printf ("  -h          Print this help message\n");
  printf ("  -k <kernel> Select interpolation kernel. Possible values are\n"
	  "              tanh, sinc, sinc2, lanzcos, hamming, hann. [tanh]\n");
  printf ("  -o <format> Output file format, possible values png, bmp, targa, jpeg. [png]\n");
  printf ("  -p          Disable padding\n");
  printf ("  -r <angle>  Rotation angle (possibly)\n");
  printf ("  -s <size>   Size of the output image. [1408]\n");
  printf ("  -t <angle>  Tilt angle\n");
  printf ("  -v          Verbose messages\n");
}
