#include "eclipse_types.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

extern int outwidth;
extern int outheight;
extern float angle;
extern int pad;

float ap=PI_NUMB;// *180.0/180.0;
float *xcoord[6];
float *ycoord[6];
int nnz[6] = {0,0,0,0,0,0};
std::vector<int> xint;
std::vector<int> yint;

/*-------------------------------------------------------------------------*/
/**
  @brief        Warp an image according to a polynomial transformation.
  @param        image_in                Image to warp.
  @param        kernel_type             Interpolation kernel to use.
  @param        poly_u                  Polynomial transform in U.
  @param        poly_v                  Polynomial transform in V.
  @return       1 newly allocated image.

  Warp an image according to a polynomial transform. Provide two
  polynomials (see poly2d.h for polynomials in this library) Pu and Pv such
  as:

  \begin{verbatim}
  x = poly2d_compute(Pu, u, v)
  y = poly2d_compute(Pv, u, v)
  \end{verbatim}

  Attention! The polynomials define a reverse transform. (u,v) are
  coordinates in the warped image and (x,y) are coordinates in the original
  image. The transform you provide is used to compute from the warped
  image, which pixels contributed in the original image.

  The output image will have strictly the same size as in the input image.
  Beware that for extreme transformations, this might lead to blank images
  as result.

  See the function generate_interpolation_kernel() for possible kernel
  types. If you want to use a default kernel, provide NULL for kernel type.

  The returned image is a newly allocated objet, use image_del() to
  deallocate it.

 */
/*--------------------------------------------------------------------------*/
Image* image_warp_generic (const Image *inimage,
			   const char  *kernel_type)
{
  //image_t    * image_out ;
  Image *outimage;
  int          i, j, k ;
  int          lx_out, ly_out ;
  double       curred, curgreen, curblue, curalpha ;
  double      neighbors_red[16] ;
  double      neighbors_green[16] ;
  double      neighbors_blue[16] ;
  double      neighbors_alpha[16] ;
  double       rsc[8],
    sumrs ;
  double       x, y ;
  int          px, py ;
  int          pos ;
  int          tabx, taby ;
  double     * kernel ;
  int          leaps[16] ;
  
  // if (image_in == NULL) return NULL ;
  if (inimage == NULL) return NULL ;
  
  /* Generate default interpolation kernel */
  kernel = generate_interpolation_kernel(kernel_type) ;
  if (kernel == NULL) {
    fprintf(stderr,"ERROR: image_warp_generic\n"
	    "   cannot generate kernel: aborting resampling\n") ;
    return NULL ;
  }
  
  /* Compute new image size   */
  lx_out = (int)outwidth ;
  ly_out = (int)outheight;
  
  outimage = new Image(lx_out, ly_out, inimage->bpp); 

  float *xs;
  float *ys;
  int n;
  switch (inimage->proj) {
  case BACK:
    xs = xcoord[0]; ys = ycoord[0]; n=nnz[0]; break;
  case DOWN:
    xs = xcoord[1]; ys = ycoord[1]; n=nnz[1]; break;
  case FRONT:
    xs = xcoord[2]; ys = ycoord[2]; n=nnz[2]; break;
  case LEFT:
    xs = xcoord[3]; ys = ycoord[3]; n=nnz[3]; break;
  case RIGHT:
    xs = xcoord[4]; ys = ycoord[4]; n=nnz[4]; break;
  case UP:
    xs = xcoord[5]; ys = ycoord[5]; n=nnz[5]; break;
  default:
    n=0; break;
  }
  
  if (n==0) {
    for (int i=0; i<lx_out*ly_out; i++) {
      outimage->red[i]=0.0;
      outimage->blue[i]=0.0;
      outimage->green[i]=0.0;
    }
    if (inimage->alpha!=NULL)
      for (int i=0; i<lx_out*ly_out; i++) 
	outimage->alpha[i]=0.0;
    return (outimage);
  }

  /* Pre compute leaps for 16 closest neighbors positions */

  int sizex = (inimage->width+2*inimage->pad);
  int sizey = (inimage->height+2*inimage->pad);
  
  leaps[0] = -1 -   sizex; //image_in->lx ;
  leaps[1] =    -   sizex; //image_in->lx ;
  leaps[2] =  1 -   sizex; //image_in->lx ;
  leaps[3] =  2 -   sizex; //image_in->lx ;
  leaps[8] = -1 +   sizex; //image_in->lx ;
  leaps[9] =        sizex; //image_in->lx ;
  leaps[10]=  1 +   sizex; //image_in->lx ;
  leaps[11]=  2 +   sizex; //image_in->lx ;
  leaps[12]= -1 + 2*sizex; //image_in->lx ;
  leaps[13]=      2*sizex; //image_in->lx ;
  leaps[14]=  1 + 2*sizex; //image_in->lx ;
  leaps[15]=  2 + 2*sizex; //image_in->lx ;
  leaps[4] = -1 ;
  leaps[5] =  0 ;
  leaps[6] =  1 ;
  leaps[7] =  2 ;
  
  /* Double loop on the output image  */
  for (i=0 ; i< lx_out ; i++) {
    for (j=0 ; j < ly_out ; j++) {
      /* Compute the original source for this pixel   */
      
      //trans (&res[0], proj, i, j, sizex, sizey, 0.0, 90.0);
      //x = res[0];
      //y = res[1];
      x = xs[i+lx_out*j];
      y = ys[i+lx_out*j];
 
      /* Which is the closest integer positioned neighbor?    */
      px = (int)x ;
      py = (int)y ;
     
      if ((px < 1) ||
	  (px > (sizex-3)) ||
	  (py < 1) ||
	  (py > (sizey-3))) {
	(outimage->red)[i+j*lx_out] = 0.0 ;
	outimage->green[i+j*lx_out] = 0.0 ;
	outimage->blue[i+j*lx_out] = 0.0 ;
	if (outimage->alpha != NULL)
	  outimage->alpha[i+j*lx_out] = 0.0 ;
      }
      else {
	/* Now feed the positions for the closest 16 neighbors  */
	pos = px + py * sizex ;
	for (k=0 ; k<16 ; k++) {
	  neighbors_red[k]   = inimage->red[(int)(pos+leaps[k])] ;
	  neighbors_green[k] = inimage->green[(int)(pos+leaps[k])] ;
	  neighbors_blue[k]  = inimage->blue[(int)(pos+leaps[k])] ;
	}
	if (outimage->alpha != NULL)
	  for (k=0 ; k<16 ; k++) 
	    neighbors_alpha[k] = inimage->alpha[(int)(pos+leaps[k])] ;
	  
	/* Which tabulated value index shall we use?    */
	tabx = (int)((x - (double)px) * (double)(TABSPERPIX)) ;
	taby = (int)((y - (double)py) * (double)(TABSPERPIX)) ;
	
	/* Compute resampling coefficients  */
	/* rsc[0..3] in x, rsc[4..7] in y   */
	
	rsc[0] = kernel[TABSPERPIX + tabx] ;
	rsc[1] = kernel[tabx] ;
	rsc[2] = kernel[TABSPERPIX - tabx] ;
	rsc[3] = kernel[2 * TABSPERPIX - tabx] ;
	rsc[4] = kernel[TABSPERPIX + taby] ;
	rsc[5] = kernel[taby] ;
	rsc[6] = kernel[TABSPERPIX - taby] ;
	rsc[7] = kernel[2 * TABSPERPIX - taby] ;
	
	sumrs = (rsc[0]+rsc[1]+rsc[2]+rsc[3]) *
	  (rsc[4]+rsc[5]+rsc[6]+rsc[7]) ;
	
	/* Compute interpolated pixel now   */
	curred =   rsc[4] * (  rsc[0]*neighbors_red[0] +
			       rsc[1]*neighbors_red[1] +
			       rsc[2]*neighbors_red[2] +
			       rsc[3]*neighbors_red[3] ) +
	  rsc[5] * (  rsc[0]*neighbors_red[4] +
		      rsc[1]*neighbors_red[5] +
		      rsc[2]*neighbors_red[6] +
		      rsc[3]*neighbors_red[7] ) +
	  rsc[6] * (  rsc[0]*neighbors_red[8] +
		      rsc[1]*neighbors_red[9] +
		      rsc[2]*neighbors_red[10] +
		      rsc[3]*neighbors_red[11] ) +
	  rsc[7] * (  rsc[0]*neighbors_red[12] +
		      rsc[1]*neighbors_red[13] +
		      rsc[2]*neighbors_red[14] +
		      rsc[3]*neighbors_red[15] ) ; 
	curblue =   rsc[4] * ( rsc[0]*neighbors_blue[0] +
			       rsc[1]*neighbors_blue[1] +
			       rsc[2]*neighbors_blue[2] +
			       rsc[3]*neighbors_blue[3] ) +
	  rsc[5] * (  rsc[0]*neighbors_blue[4] +
		      rsc[1]*neighbors_blue[5] +
		      rsc[2]*neighbors_blue[6] +
		      rsc[3]*neighbors_blue[7] ) +
	  rsc[6] * (  rsc[0]*neighbors_blue[8] +
		      rsc[1]*neighbors_blue[9] +
		      rsc[2]*neighbors_blue[10] +
		      rsc[3]*neighbors_blue[11] ) +
	  rsc[7] * (  rsc[0]*neighbors_blue[12] +
		      rsc[1]*neighbors_blue[13] +
		      rsc[2]*neighbors_blue[14] +
		      rsc[3]*neighbors_blue[15] ) ; 
	curgreen =   rsc[4] * (rsc[0]*neighbors_green[0] +
			       rsc[1]*neighbors_green[1] +
			       rsc[2]*neighbors_green[2] +
			       rsc[3]*neighbors_green[3] ) +
	  rsc[5] * (  rsc[0]*neighbors_green[4] +
		      rsc[1]*neighbors_green[5] +
		      rsc[2]*neighbors_green[6] +
		      rsc[3]*neighbors_green[7] ) +
	  rsc[6] * (  rsc[0]*neighbors_green[8] +
		      rsc[1]*neighbors_green[9] +
		      rsc[2]*neighbors_green[10] +
		      rsc[3]*neighbors_green[11] ) +
	  rsc[7] * (  rsc[0]*neighbors_green[12] +
		      rsc[1]*neighbors_green[13] +
		      rsc[2]*neighbors_green[14] +
		      rsc[3]*neighbors_green[15] ) ; 
	if (outimage->alpha != NULL) 
	  curalpha =   rsc[4] * (  rsc[0]*neighbors_alpha[0] +
				   rsc[1]*neighbors_alpha[1] +
				   rsc[2]*neighbors_alpha[2] +
				   rsc[3]*neighbors_alpha[3] ) +
	    rsc[5] * (  rsc[0]*neighbors_alpha[4] +
			rsc[1]*neighbors_alpha[5] +
			rsc[2]*neighbors_alpha[6] +
			rsc[3]*neighbors_alpha[7] ) +
	    rsc[6] * (  rsc[0]*neighbors_alpha[8] +
			rsc[1]*neighbors_alpha[9] +
			rsc[2]*neighbors_alpha[10] +
			rsc[3]*neighbors_alpha[11] ) +
	    rsc[7] * (  rsc[0]*neighbors_alpha[12] +
			rsc[1]*neighbors_alpha[13] +
			rsc[2]*neighbors_alpha[14] +
			rsc[3]*neighbors_alpha[15] ) ; 
	
	/* Affect the value to the output image */
	outimage->red[i+j*lx_out] = (pixelvalue)(curred/sumrs) ;
	outimage->green[i+j*lx_out] = (pixelvalue)(curgreen/sumrs) ;
	outimage->blue[i+j*lx_out] = (pixelvalue)(curblue/sumrs) ;
	if (outimage->alpha != NULL) 
	  outimage->alpha[i+j*lx_out] = (pixelvalue)(curalpha/sumrs) ;
	/* done ! */
      }      
    }
  }
  free(kernel) ;
  return outimage ;
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Generate an interpolation kernel to use in this module.
  @param	kernel_type		Type of interpolation kernel.
  @return	1 newly allocated array of doubles.

  Provide the name of the kernel you want to generate. Supported kernel
  types are:

  \begin{tabular}{ll}
  NULL			&	default kernel, currently "tanh" \\
  "default"		&	default kernel, currently "tanh" \\
  "tanh"		&	Hyperbolic tangent \\
  "sinc2"		&	Square sinc \\
  "lanczos"		&	Lanczos2 kernel \\
  "hamming"		&	Hamming kernel \\
  "hann"		&	Hann kernel
  \end{tabular}

  The returned array of doubles is ready of use in the various re-sampling
  functions in this module. It must be deallocated using free().
 */
/*--------------------------------------------------------------------------*/
double * generate_interpolation_kernel(const char * kernel_type)
{
  double  *	tab ;
  int     	i ;
  double  	x ;
  double		alpha ;
  double		inv_norm ;
  int     	samples = KERNEL_SAMPLES ;

  if (kernel_type==NULL) {
    tab = generate_interpolation_kernel("tanh") ;
  } else if (!strcmp(kernel_type, "default")) {
    tab = generate_interpolation_kernel("tanh") ;
  } else if (!strcmp(kernel_type, "sinc")) {
    tab = (double*)malloc(samples * sizeof(double)) ;
    tab[0] = 1.0 ;
    tab[samples-1] = 0.0 ;
    for (i=1 ; i<samples ; i++) {
      x = (double)KERNEL_WIDTH * (double)i/(double)(samples-1) ;
      tab[i] = sinc(x) ;
    }
  } else if (!strcmp(kernel_type, "sinc2")) {
    tab = (double*)malloc(samples * sizeof(double)) ;
    tab[0] = 1.0 ;
    tab[samples-1] = 0.0 ;
    for (i=1 ; i<samples ; i++) {
      x = 2.0 * (double)i/(double)(samples-1) ;
      tab[i] = sinc(x) ;
      tab[i] *= tab[i] ;
    }
  } else if (!strcmp(kernel_type, "lanczos")) {
    tab = (double*)malloc(samples * sizeof(double)) ;
    for (i=0 ; i<samples ; i++) {
      x = (double)KERNEL_WIDTH * (double)i/(double)(samples-1) ;
      if (fabs(x)<2) {
	tab[i] = sinc(x) * sinc(x/2) ;
      } else {
	tab[i] = 0.00 ;
      }
    }
  } else if (!strcmp(kernel_type, "hamming")) {
    tab = (double*)malloc(samples * sizeof(double)) ;
    alpha = 0.54 ;
    inv_norm  = 1.00 / (double)(samples - 1) ;
    for (i=0 ; i<samples ; i++) {
      x = (double)i ;
      if (i<(samples-1)/2) {
	tab[i] = alpha + (1-alpha) * cos(2.0*PI_NUMB*x*inv_norm) ;
      } else {
	tab[i] = 0.0 ;
      }
    }
  } else if (!strcmp(kernel_type, "hann")) {
    tab = (double*)malloc(samples * sizeof(double)) ;
    alpha = 0.50 ;
    inv_norm  = 1.00 / (double)(samples - 1) ;
    for (i=0 ; i<samples ; i++) {
      x = (double)i ;
      if (i<(samples-1)/2) {
	tab[i] = alpha + (1-alpha) * cos(2.0*PI_NUMB*x*inv_norm) ;
      } else {
	tab[i] = 0.0 ;
      }
    }
  } else if (!strcmp(kernel_type, "tanh")) {
    tab = generate_tanh_kernel(TANH_STEEPNESS) ;
  } else {
    fprintf(stderr, "unrecognized kernel type [%s]: aborting generation\n",
	    kernel_type) ;
    return NULL ;
  }
  return tab ;
}


/*-------------------------------------------------------------------------*/
/**
  @brief	Generate a hyperbolic tangent kernel.
  @param	steep	Steepness of the hyperbolic tangent parts.
  @return	1 pointer to a newly allocated array of doubles.

  The following function builds up a good approximation of a box filter. It
  is built from a product of hyperbolic tangents. It has the following
  properties:

  \begin{itemize}
  \item It converges very quickly towards +/- 1.
  \item The converging transition is very sharp.
  \item It is infinitely differentiable everywhere (i.e. smooth).
  \item The transition sharpness is scalable.
  \end{itemize}

  The returned array must be deallocated using free().
 */
/*--------------------------------------------------------------------------*/
#define hk_gen(x,s) (((tanh(s*(x+0.5))+1)/2)*((tanh(s*(-x+0.5))+1)/2))
double * generate_tanh_kernel(double steep)
{
    double  *   kernel ;
    double  *   x ;
    double      width ;
    double      inv_np ;
    double      ind ;
    int         i ;
    int         np ;
    int         samples ;

    width   = (double)TABSPERPIX / 2.0 ; 
    samples = KERNEL_SAMPLES ;
    np      = 32768 ; /* Hardcoded: should never be changed */
    inv_np  = 1.00 / (double)np ;

    /*
     * Generate the kernel expression in Fourier space
     * with a correct frequency ordering to allow standard FT
     */
    x = (double*)malloc((2*np+1)*sizeof(double)) ;
    for (i=0 ; i<np/2 ; i++) {
        ind      = (double)i * 2.0 * width * inv_np ;
        x[2*i]   = hk_gen(ind, steep) ;
        x[2*i+1] = 0.00 ;
    }
    for (i=np/2 ; i<np ; i++) {
        ind      = (double)(i-np) * 2.0 * width * inv_np ;
        x[2*i]   = hk_gen(ind, steep) ;
        x[2*i+1] = 0.00 ;
    }

    /* 
     * Reverse Fourier to come back to image space
     */
    reverse_tanh_kernel(x, np) ;

    /*
     * Allocate and fill in returned array
     */
    kernel = (double*)malloc(samples * sizeof(double)) ;
    for (i=0 ; i<samples ; i++) {
        kernel[i] = 2.0 * width * x[2*i] * inv_np ;
    }
    free(x) ;
    return kernel ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Cardinal sine.
  @param	x	double value.
  @return	1 double.

  Compute the value of the function sinc(x)=sin(pi*x)/(pi*x) at the
  requested x.
 */
/*--------------------------------------------------------------------------*/
double sinc(double x)
{
    if (fabs(x)<1e-4)
        return (double)1.00 ;
    else
        return ((sin(x * (double)PI_NUMB)) / (x * (double)PI_NUMB)) ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Bring a hyperbolic tangent kernel from Fourier to normal space.
  @param	data	Kernel samples in Fourier space.
  @param	nn		Number of samples in the input kernel.
  @return	void

  Bring back a hyperbolic tangent kernel from Fourier to normal space. Do
  not try to understand the implementation and DO NOT MODIFY THIS FUNCTION.
 */
/*--------------------------------------------------------------------------*/
#define KERNEL_SW(a,b) tempr=(a);(a)=(b);(b)=tempr
void reverse_tanh_kernel(double * data, int nn)
{
    unsigned long   n,
					mmax,
					m,
					i, j,
					istep ;
    double  wtemp,
            wr,
            wpr,
            wpi,
            wi,
            theta;
    double  tempr,
            tempi;

    n = (unsigned long)nn << 1;
    j = 1;
    for (i=1 ; i<n ; i+=2) {
        if (j > i) {
            KERNEL_SW(data[j-1],data[i-1]);
            KERNEL_SW(data[j],data[i]);
        }
        m = n >> 1;
        while (m>=2 && j>m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax = 2;
    while (n > mmax) {
        istep = mmax << 1;
        theta = 2 * PI_NUMB / mmax;
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr  = 1.0;
        wi  = 0.0;
        for (m=1 ; m<mmax ; m+=2) {
            for (i=m ; i<=n ; i+=istep) {
                j = i + mmax;
                tempr = wr * data[j-1] - wi * data[j];
                tempi = wr * data[j]   + wi * data[j-1];
                data[j-1] = data[i-1] - tempr;
                data[j]   = data[i]   - tempi;
                data[i-1] += tempr;
                data[i]   += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}
#undef KERNEL_SW

void FreeTrans ()
{
  for (int i=0; i<6; i++) {
    delete[] xcoord[i];
    delete[] ycoord[i];
  }
}

void CalcTrans (int inwidth1, int inheight1, float alpha, float beta)
{
  float np = float(pad); //(angle/180.0-0.5)*(float)inwidth + 
  float inwidth, inheight;
  inwidth = inwidth1;
  inheight=inheight1;

//   np = (angle/180.0-0.5)*(float)inwidth1 + (float)pad;
//   inwidth = inwidth1*(90.0/angle);
//   inheight = inheight1*(90.0/angle);

  // printf ("%f\n", np);
  for (int i=0; i<6; i++) {
    xcoord[i] = new float[outwidth*outheight];
    ycoord[i] = new float[outwidth*outheight];
    for (int j=0; j<outwidth*outheight; j++) {
      (xcoord[i])[j] = -1.0;
      (ycoord[i])[j] = -1.0;
    }
  }

  //printf ("%f %f\n", alpha, beta);

  bool rotate = false;
  float cosa = cos(alpha*0.01745329251);
  float sina = sin(alpha*0.01745329251);
  float cosb = cos(beta*0.01745329251);
  float sinb = sin(beta*0.01745329251);
  if (alpha != 0.0 || beta != 0.0) {
    // printf ("Rotate\n");
    rotate = true;
  }

  for (int u=0; u<outwidth; u++)
    for (int v=0; v<outwidth; v++) {
      // Normalizing the coordinates of the output image between [-1,1]
      float udot = 2.0*(float)u/(float)(outwidth-1.0)-1.0;
      float vdot = 2.0*(float)v/(float)(outheight-1.0)-1.0;
      
      float phi = atan2(vdot,udot);
      float ruv = sqrt(udot*udot+vdot*vdot);
      
      // Azimuth - angle between l.o.s and z-axis
      float theta = ruv*(ap)/2.0;
      
      // Unit vector of l.o.s.
      float x0 = sin(theta)*cos(phi);
      float y0 = sin(theta)*sin(phi);
      float z0 = cos(theta);
      if (rotate) {
	float x1 = x0;
	float y1 = cosb*y0 + sinb*z0;
	float z1 =-sinb*y0 + cosb*z0;
	
	float x2 = cosa*x1 + sina*y1;
	float y2 = -sina*x1+ cosa*y1;
	float z2 = z1;

	x0 = x2; y0=y2; z0= z2;
      }

      float dx,dy,dz,dr;
      
      float xcb=-1.0;
      float xcd=-1.0;
      float xcf=-1.0;
      float xcl=-1.0;
      float xcr=-1.0;
      float xcu=-1.0;
      float ycb=-1.0;
      float ycd=-1.0;
      float ycf=-1.0;
      float ycl=-1.0;
      float ycr=-1.0;
      float ycu=-1.0;
      if (ruv<=1.0) {
	//  B: BACK, i.e. y=-1;
	if (y0<0.0) {
	  dx = x0/y0;
	  dy = -z0/y0;
	  xcb = (float)np+0.5*(dx+1.0)*(float)(inwidth);
	  ycb = (float)np+0.5*(dy+1.0)*(float)(inheight);
	}
	// D: DOWN, i.e. z=-1;
	if (z0<0.0) {
	  dx = -x0/z0;
	  dy = -y0/z0;
	  xcd = (float)np+0.5*(dx+1.0)*(float)(inwidth);
	  ycd = (float)np+0.5*(dy+1.0)*(float)(inheight);
	}
	// F: FRONT, i.e. y=1;
	if (y0>0.0) {
	  dx = x0/y0;
	  dy = z0/y0;
	  xcf = (float)np+0.5*(dx+1.0)*(float)(inwidth);
	  ycf = (float)np+0.5*(dy+1.0)*(float)(inheight);
	}
	// L: LEFT, i.e. x=-1;
	if (x0<0.0) {
	  dx = -y0/x0;
	  dy = -z0/x0;
	  xcl = (float)np+0.5*(dx+1.0)*(float)(inwidth);
	  ycl = (float)np+0.5*(dy+1.0)*(float)(inheight);
	}
	// R: RIGHT, i.e. x=1;
	if (x0>0.0) {
	  dx = -y0/x0;
	  dy = z0/x0;
	  xcr = (float)np+0.5*(dx+1.0)*(float)(inwidth);
	  ycr = (float)np+0.5*(dy+1.0)*(float)(inheight);
	}
	// U: TOP/UP, i.e. z=1
	if (z0>0.0) {
	  dx = x0/z0;
	  dy = y0/z0;
	  xcu = (float)np+0.5*(dx+1.0)*(float)(inwidth);
	  ycu = (float)np+0.5*(dy+1.0)*(float)(inheight);
	}
	if (xcb<1.0 || xcb>=(2.0*np+(float)inwidth-2.0) || ycb<1.0 || ycb>=(2.0*np+(float)inwidth-2.0)) 
	  {xcb=-1.0;ycb=-1.0;} //0
	if (xcd<1.0 || xcd>=(2.0*np+(float)inwidth-2.0) || ycd<1.0 || ycd>=(2.0*np+(float)inwidth-2.0)) 
	  {xcd=-1.0;ycd=-1.0;} //1
	if (xcf<1.0 || xcf>=(2.0*np+(float)inwidth-2.0) || ycf<1.0 || ycf>=(2.0*np+(float)inwidth-2.0)) 
	  {xcf=-1.0;ycf=-1.0;} //2
	if (xcl<1.0 || xcl>=(2.0*np+(float)inwidth-2.0) || ycl<1.0 || ycl>=(2.0*np+(float)inwidth-2.0)) 
	  {xcl=-1.0;ycl=-1.0;} //3
	if (xcr<1.0 || xcr>=(2.0*np+(float)inwidth-2.0) || ycr<1.0 || ycr>=(2.0*np+(float)inwidth-2.0)) 
	  {xcr=-1.0;ycr=-1.0;} //4
	if (xcu<1.0 || xcu>=(2.0*np+(float)inwidth-2.0) || ycu<1.0 || ycu>=(2.0*np+(float)inwidth-2.0)) 
	  {xcu=-1.0;ycu=-1.0;} //5
      } //ruv<1.0
      while (true) {
	(xcoord[2])[u+v*outwidth] = xcf;
	(ycoord[2])[u+v*outwidth] = ycf;
	if (xcf>=0.0 && ycf>=0.0) {
	  nnz[2]++;
	  break;
	}
	(xcoord[3])[u+v*outwidth] = xcl;
	(ycoord[3])[u+v*outwidth] = ycl;
	if (xcl>=0.0 && ycl>=0.0) {
	  nnz[3]++;
	  break;
	}
	(xcoord[4])[u+v*outwidth] = xcr;
	(ycoord[4])[u+v*outwidth] = ycr;
	if (xcr>=0.0 && ycr>=0.0) {
	  nnz[4]++;
	  break;
	}
	(xcoord[5])[u+v*outwidth] = xcu;
	(ycoord[5])[u+v*outwidth] = ycu;
	if (xcu>=0.0 && ycu>=0.0) {
	  nnz[5]++;
	  break;
	}
	(xcoord[0])[u+v*outwidth] = xcb;
	(ycoord[0])[u+v*outwidth] = ycb;
	if (xcb>=0.0 && ycb>=0.0) {
	  nnz[0]++;
	  break;
	}
	(xcoord[1])[u+v*outwidth] = xcd;
	(ycoord[1])[u+v*outwidth] = ycd;
	if (xcd>=0.0 || ycd>=0.0) {
	  nnz[1]++;
	  break;
	}
	if (ruv<=1.0) {
	    xint.push_back (u);
	    yint.push_back (v);
	}
	break;
      }
    }
  if (xint.size()>0) {
    int i = xint.size();
    printf ("Interpolation needed for %d pixels\n", i);
  }
}


/*
 *
 *  Do not remove this function even though it's not used. 
 *  Testing of stiching is easier to do with this than with CalcTrans
 *
 *  At the moment it does not exactly correspond to CalcTrans, though
 */
void trans (float *res, char ttype, int u, int v, int inwidth, int inheight, float alpha, float beta)
{
  float np = 0;
  // Normalizing the coordinates of the output image between [-1,1]
  float udot = 2.0*(float)u/(float)outwidth-1.0;
  float vdot = 2.0*(float)v/(float)outheight-1.0;

  float phi = atan2(vdot,udot);
  float ruv = sqrt(udot*udot+vdot*vdot);

  // Azimuth - angle between l.o.s and z-axis
  float theta = ruv*(ap)/2.0;
 
  bool rotate = false;
  float cosa = cos(alpha*0.01745329251);
  float sina = sin(alpha*0.01745329251);
  float cosb = cos(beta*0.01745329251);
  float sinb = sin(beta*0.01745329251);
  if (alpha != 0.0 || beta != 0.0) {
    // printf ("Rotate\n");
    rotate = true;
  }

  // Unit vector of l.o.s.
  float x0 = sin(theta)*cos(phi);
  float y0 = sin(theta)*sin(phi);
  float z0 = cos(theta);
  if (rotate) {
    float x1 = x0;
    float y1 = cosb*y0 + sinb*z0;
    float z1 =-sinb*y0 + cosb*z0;
    
    float x2 = cosa*x1 + sina*y1;
    float y2 = -sina*x1+ cosa*y1;
    float z2 = z1;

    x0 = x2; y0=y2; z0= z2;
  }

  float dx,dy,dz,dr;

  // Default values 
  res[0] = -1.0; res[1] = -1.0;

  //
  // These needs to be checked!!!
  // 
  if (ruv<=1.0) {
    switch (ttype) {
    case 'U': //  TOP/UP, i.e. z=1
      if (z0>0.0) {
	dx = x0/z0;
	dy = y0/z0;
	res[0] = np+0.5*(dx+1.0)*(float)(inwidth-2.0*np);
	res[1] = np+0.5*(dy+1.0)*(float)(inheight-2.0*np);
      }
      break;
    case 'D': //  BOTTOM/DOWN, i.e. z=-1;
      if (z0<0.0) {
	dx = x0/z0;
	dy = y0/z0;
	res[0] = np+0.5*(dx+1.0)*(float)(inwidth-2.0*np);
	res[1] = np+0.5*(dy+1.0)*(float)(inheight-2.0*np);
      }
      break;
    case 'F': //  FRONT, i.e. y=1;
      if (y0>0.0) {
	dx = x0/y0;
	dy = z0/y0;
	res[0] = np+0.5*(dx+1.0)*(float)(inwidth-2.0*np);
	res[1] = np+0.5*(dy+1.0)*(float)(inheight-2.0*np);
      }
      break;
    case 'B': //  BACK, i.e. y=-1;
      if (y0<0.0) {
	dx = x0/y0;
	dy = -z0/y0;
	res[0] = np+0.5*(dx+1.0)*(float)(inwidth-2.0*np);
	res[1] = np+0.5*(dy+1.0)*(float)(inheight-2.0*np);
      }
      break;
    case 'L': //  LEFT, i.e. x=-1;
      if (x0<0.0) {
	dx = -y0/x0;
	dy = -z0/x0;
	res[0] = np+0.5*(dx+1.0)*(float)(inwidth-2.0*np);
	res[1] = np+0.5*(dy+1.0)*(float)(inheight-2.0*np);
      }
      break;
    case 'R': //  RIGHT, i.e. x=1;
      if (x0>0.0) {
	dx = -y0/x0;
	dy = z0/x0;
	res[0] = np+0.5*(dx+1.0)*(float)(inwidth-2.0*np);
	res[1] = np+0.5*(dy+1.0)*(float)(inheight-2.0*np);
      }
      break;
    default:
      printf ("ERROR\n");
      abort ();
      break;
    }
  }
  if (res[0]>inwidth || res[0]<-1.0)
    res[0] = -1.0;
  if (res[1]>inheight || res[1]<-1.0) 
    res[1] = -1.0;
}

