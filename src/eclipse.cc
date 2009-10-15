/*      Copyright (c) 1995-2001, European Southern Observatory
        All rights reserved.

        Redistribution and use in source and binary forms, with or
        without modification, are permitted provided that the following
        conditions are met:

                Redistributions of source code must retain the above
                copyright notice, this list of conditions and the
                following disclaimer. 

                Redistributions in binary form must reproduce the above
                copyright notice, this list of conditions and the
                following disclaimer in the documentation and/or other
                materials provided with the distribution. 

                Neither the name of the European Southern Observatory
                nor the names of its contributors may be used to endorse
                or promote products derived from this software without
                specific prior written permission. 

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
        CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
        INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
        MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE
        LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
        EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
        TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
        DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
        LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
        IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
        THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
    Adapted by J. Reunanen, 2009
 */

#include "eclipse_types.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

extern int outwidth;
extern int outheight;
extern float angle;
extern float aperture;
extern int pad;
extern bool interp;

float *xcoord[6];
float *ycoord[6];
int umin[6], umax[6], vmin[6], vmax[6];
int nnz[6] = {0,0,0,0,0,0};
std::vector<int> xint;
std::vector<int> yint;

float  * kernel ;

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
Image* image_warp_generic (const Image *inimage, Image *outimage1)
{
  //image_t    * image_out ;
  Image *outimage;
  int          i, j, k ;
  int          lx_out, ly_out ;
  float       curred, curgreen, curblue, curalpha ;
  float      neighbors_red[16] ;
  float      neighbors_green[16] ;
  float      neighbors_blue[16] ;
  float      neighbors_alpha[16] ;
  float       rsc[8], sumrs ;
  float       x, y ;
  int          px, py ;
  int          pos ;
  int          tabx, taby ;
  int          leaps[16] ;
  
  if (inimage == NULL) return NULL ;

  /* Compute new image size   */
  lx_out = (int)outwidth ;
  ly_out = (int)outheight;
  
  if (outimage1==NULL)
    outimage = new Image(lx_out, ly_out, inimage->bpp); 
  else
    outimage = outimage1;

  float *xs;
  float *ys;
  int n, x0, x1, y0, y1;
  switch (inimage->proj) {
  case BACK:
    xs = xcoord[0]; ys = ycoord[0]; n=nnz[0]; 
    x0 = umin[0]; x1=umax[0]; y0=vmin[0]; y1=vmax[0];
    break;
  case DOWN:
    xs = xcoord[1]; ys = ycoord[1]; n=nnz[1];  
    x0 = umin[1]; x1=umax[1]; y0=vmin[1]; y1=vmax[1];
    break;
  case FRONT:
    xs = xcoord[2]; ys = ycoord[2]; n=nnz[2];  
    x0 = umin[2]; x1=umax[2]; y0=vmin[2]; y1=vmax[2];
    break;
  case LEFT:
    xs = xcoord[3]; ys = ycoord[3]; n=nnz[3];  
    x0 = umin[3]; x1=umax[3]; y0=vmin[3]; y1=vmax[3];
    break;
  case RIGHT:
    xs = xcoord[4]; ys = ycoord[4]; n=nnz[4];  
    x0 = umin[4]; x1=umax[4]; y0=vmin[4]; y1=vmax[4];
    break;
  case UP:
    xs = xcoord[5]; ys = ycoord[5]; n=nnz[5];  
    x0 = umin[5]; x1=umax[5]; y0=vmin[5]; y1=vmax[5];
    break;
  default:
    n=0; break;
  }
  
  if (n==0) 
    return (outimage);

  x0--; y0--; x1++; y1++;
  if (x0<0) x0=0;
  if (y0<0) y0=0;
  if (x1>lx_out) x1=lx_out;
  if (y1>ly_out) y1=ly_out;


  int sizex = (inimage->width+2*inimage->pad);
  int sizey = (inimage->height+2*inimage->pad);
  
  /* Pre compute leaps for 16 closest neighbors positions */
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

  if (!interp) {
    for (i=x0 ; i< x1 ; i++) {
      for (j=y0 ; j < y1 ; j++) {
	x = xs[i+lx_out*j];
	y = ys[i+lx_out*j];
	
	px = (int)x ;
	py = (int)y ;
	
	if ((px >= 0) || (py >= 0) ) {
	  outimage->red[i+j*lx_out] = inimage->red[px+py*sizex];
	  outimage->green[i+j*lx_out] = inimage->green[px+py*sizex];
	  outimage->blue[i+j*lx_out] = inimage->blue[px+py*sizex];
	  if (outimage->alpha != NULL)
	    outimage->alpha[i+j*lx_out] = inimage->alpha[px+py*sizex];
	}
      }
    }
    return (outimage);
  }

  /* Double loop on the output image  */
  for (i=x0 ; i< x1 ; i++) {
    for (j=y0 ; j < y1 ; j++) {
      /* Compute the original source for this pixel   */
      
      //trans (&res[0], proj, i, j, sizex, sizey, 0.0, 90.0);
      //x = res[0];
      //y = res[1];
      x = xs[i+lx_out*j];
      y = ys[i+lx_out*j];
 
      /* Which is the closest integer positioned neighbor?    */
      px = (int)x ;
      py = (int)y ;
     
//       if ((px < 1) ||
// 	  (px > (sizex-3)) ||
// 	  (py < 1) ||
// 	  (py > (sizey-3))) {
// 	(outimage->red)[i+j*lx_out] = 0.0 ;
// 	outimage->green[i+j*lx_out] = 0.0 ;
// 	outimage->blue[i+j*lx_out] = 0.0 ;
// 	if (outimage->alpha != NULL)
// 	  outimage->alpha[i+j*lx_out] = 0.0 ;
//       }
//       else {
      if ((px >= 1) && (px <= (sizex-3)) && (py >= 1) && (px <= (sizey-3))) {
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
float * generate_interpolation_kernel(const char * kernel_type)
{
  float  *	tab ;
  int     	i ;
  double  	x ;
  double		alpha ;
  double		inv_norm ;
  int     	samples = KERNEL_SAMPLES ;

  if (kernel_type==NULL) {
    tab = generate_interpolation_kernel("sinc") ;
  } else if (!strcmp(kernel_type, "default")) {
    tab = generate_interpolation_kernel("sinc") ;
  } else if (!strcmp(kernel_type, "sinc")) {
    tab = (float*)malloc(samples * sizeof(float)) ;
    tab[0] = 1.0 ;
    tab[samples-1] = 0.0 ;
    for (i=1 ; i<samples ; i++) {
      x = (double)KERNEL_WIDTH * (double)i/(double)(samples-1) ;
      tab[i] = sinc(x) ;
    }
  } else if (!strcmp(kernel_type, "sinc2")) {
    tab = (float*)malloc(samples * sizeof(float)) ;
    tab[0] = 1.0 ;
    tab[samples-1] = 0.0 ;
    for (i=1 ; i<samples ; i++) {
      x = 2.0 * (double)i/(double)(samples-1) ;
      tab[i] = sinc(x) ;
      tab[i] *= tab[i] ;
    }
  } else if (!strcmp(kernel_type, "lanczos")) {
    tab = (float*)malloc(samples * sizeof(float)) ;
    for (i=0 ; i<samples ; i++) {
      x = (double)KERNEL_WIDTH * (double)i/(double)(samples-1) ;
      if (fabs(x)<2) {
	tab[i] = sinc(x) * sinc(x/2) ;
      } else {
	tab[i] = 0.00 ;
      }
    }
  } else if (!strcmp(kernel_type, "hamming")) {
    tab = (float*)malloc(samples * sizeof(float)) ;
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
    tab = (float*)malloc(samples * sizeof(float)) ;
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
float * generate_tanh_kernel(double steep)
{
    float  *   kernel ;
    float  *   x ;
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
    x = (float*)malloc((2*np+1)*sizeof(float)) ;
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
    kernel = (float*)malloc(samples * sizeof(float)) ;
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
void reverse_tanh_kernel(float * data, int nn)
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

/*-------------------------------------------------------------------------*/
/**
  @brief	Interpolate a vector along new abscissas.
  @param	x		List of x positions.
  @param	y		List of y positions.
  @param	len		Number of samples in x and y.
  @param	splx	Input new list of x positions.
  @param	sply	Output list of interpolated y positions.
  @param	spllen	Number of samples in splx and sply.
  @return	Int 0 if Ok, -1 if error.

  Reference:

  \begin{verbatim}
  	Numerical Analysis, R. Burden, J. Faires and A. Reynolds.
  	Prindle, Weber & Schmidt 1981 pp 112
  \end{verbatim}

  Provide in input a known list of x and y values, and a list where
  you want the signal to be interpolated. The returned signal is
  written into sply.
 */
/*--------------------------------------------------------------------------*/
int function1d_natural_spline( pixelvalue *x, pixelvalue *y,
    	int 			len, pixelvalue	* 	splx, pixelvalue	* 	sply,
    	int 			spllen)
{
  int 			end;
  int 			loc,
    found;
  register int 	i,
    j,
    n;
  double 		*	h;			/* vector of deltas in x */
  double 		*	alpha;
  double 		*	l,
    *	mu,
    *	z,
    *	a,
    *	b,
    *	c,
    *	d,
    v;
  
  end = len - 1;
  
  a = (double *)malloc(sizeof(double) * spllen * 9) ;
  b = a + len;
  c = b + len;
  d = c + len;
  h = d + len;
  l = h + len;
  z = l + len;
  mu = z + len;
  alpha = mu + len;
  
  for (i = 0; i < len; i++) {
    a[i] = (double)y[i];
  }
  
  /* Calculate vector of differences */
  for (i = 0; i < end; i++) {
    h[i] = (double)x[i + 1] - (double)x[i];
    if (h[i] < 0.0) {
      free(a) ;
      return -1;
    }
  }
  
  /* Calculate alpha vector */
  for (n = 0, i = 1; i < end; i++, n++) {
    /* n = i - 1 */
    alpha[i] = 3.0 * ((a[i+1] / h[i]) - (a[i] / h[n]) - (a[i] / h[i]) +
		      (a[n] / h[n]));
  }
  
  /* Vectors to solve the tridiagonal matrix */
  l[0] = l[end] = 1.0;
  mu[0] = mu[end] = 0.0;
  z[0] = z[end] = 0.0;
  c[0] = c[end] = 0.0;
  
  /* Calculate the intermediate results */
  for (n = 0, i = 1; i < end; i++, n++) {
    /* n = i-1 */
    l[i] = 2 * (h[i] + h[n]) - h[n] * mu[n];
    mu[i] = h[i] / l[i];
    z[i] = (alpha[i] - h[n] * z[n]) / l[i];
  }
  for (n = end, j = end - 1; j >= 0; j--, n--) {
    /* n = j + 1 */
    c[j] = z[j] - mu[j] * c[n];
    b[j] = (a[n] - a[j]) / h[j] - h[j] * (c[n] + 2.0 * c[j]) / 3.0;
    d[j] = (c[n] - c[j]) / (3.0 * h[j]);
  }
  
  /* Now calculate the new values */
  for (j = 0; j < spllen; j++) {
    v = (double)splx[j];
    sply[j] = (pixelvalue)0;
    
    /* Is it outside the interval? */
    if ((v < (double)x[0]) || (v > (double)x[end])) {
      continue;
    }
    /* Search for the interval containing v in the x vector */
    loc = function1d_search_value(x, len, (pixelvalue)v, &found);
    if (found) {
      sply[j] = y[loc];
    } else {
      loc--;
      v -= (double)x[loc];
      sply[j] = (pixelvalue)(a[loc]+v*(b[loc]+v*(c[loc]+v*d[loc])));
    }
  }
  free(a) ;
  return 0;
}

/*-------------------------------------------------------------------------*/
/**
  @brief	Conducts a binary search for a value.
  @param	x			Contains the abscissas of interpolation.
  @param	len			Length of the x array.
  @param	key			The value to locate in x.
  @param	found_ptr	Output flag, 1 if value was found, else 0.
  @return	The index of the largest value in x for which x[i]<key.

  This function does a binary search for a value in an array. This
  routine is to be called only if key is in the interval between x[0]
  and x[len-1]. The input x array is supposed sorted.

 */
/*--------------------------------------------------------------------------*/
int function1d_search_value(pixelvalue	*	x,
    	int 			len,
    	pixelvalue 		key,
    	int 		*	found_ptr)
{
  int	high,
    low,
    middle;
  
  low  = 0;
  high = len - 1;
  
  while (high >= low) {
    middle = (high + low) / 2;
    if (key > x[middle]) {
      low = middle + 1;
    } else if (key < x[middle]) {
      high = middle - 1;
    } else {
      *found_ptr = 1;
      return (middle);
    }
  }
  *found_ptr = 0;
  return (low);
}





void FreeTrans ()
{
  for (int i=0; i<6; i++) {
    delete[] xcoord[i];
    delete[] ycoord[i];
  }
  free(kernel) ;
  
}

void CalcTrans (int inwidth1, int inheight1, float alpha, float beta,
		const char  *kernel_type)
{
  float ap = aperture/180.0*PI_NUMB;
  float np = float(pad); //(angle/180.0-0.5)*(float)inwidth + 
  float inwidth, inheight;
  inwidth = inwidth1;
  inheight=inheight1;

  kernel = generate_interpolation_kernel(kernel_type) ;
  if (kernel == NULL) {
    fprintf(stderr,"ERROR: CalcTrans\n"
	    "   cannot generate kernel: aborting resampling\n") ;
  }

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
    umin[i] = outwidth;
    umax[i] = 0;
    vmin[i] = outheight;
    vmax[i] = 0;
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
	int v0 = v; //outwidth-v;
	int u0 = outwidth-u-1; 
	//int v0 = outwidth-v-1;
	(xcoord[2])[u0+v0*outwidth] = xcf;
	(ycoord[2])[u0+v0*outwidth] = ycf;
	if (xcf>=0.0 && ycf>=0.0) {
	  if (umax[2]<u0) umax[2] = u0;
	  if (vmax[2]<v0) vmax[2] = v0;
	  if (umin[2]>u0) umin[2] = u0;
	  if (vmin[2]>v0) vmin[2] = v0;
	  nnz[2]++;
	  break;
	}
	(xcoord[3])[u0+v0*outwidth] = xcl;
	(ycoord[3])[u0+v0*outwidth] = ycl;
	if (xcl>=0.0 && ycl>=0.0) {
	  if (umax[3]<u0) umax[3] = u0;
	  if (vmax[3]<v0) vmax[3] = v0;
	  if (umin[3]>u0) umin[3] = u0;
	  if (vmin[3]>v0) vmin[3] = v0;
	  nnz[3]++;
	  break;
	}
	(xcoord[4])[u0+v0*outwidth] = xcr;
	(ycoord[4])[u0+v0*outwidth] = ycr;
	if (xcr>=0.0 && ycr>=0.0) {
	  if (umax[4]<u0) umax[4] = u0;
	  if (vmax[4]<v0) vmax[4] = v0;
	  if (umin[4]>u0) umin[4] = u0;
	  if (vmin[4]>v0) vmin[4] = v0;
	  nnz[4]++;
	  break;
	}
	(xcoord[5])[u0+v0*outwidth] = xcu;
	(ycoord[5])[u0+v0*outwidth] = ycu;
	if (xcu>=0.0 && ycu>=0.0) {
	  if (umax[5]<u0) umax[5] = u0;
	  if (vmax[5]<v0) vmax[5] = v0;
	  if (umin[5]>u0) umin[5] = u0;
	  if (vmin[5]>v0) vmin[5] = v0;
	  nnz[5]++;
	  break;
	}
	(xcoord[0])[u0+v0*outwidth] = xcb;
	(ycoord[0])[u0+v0*outwidth] = ycb;
	if (xcb>=0.0 && ycb>=0.0) {
	  if (umax[0]<u0) umax[0] = u0;
	  if (vmax[0]<v0) vmax[0] = v0;
	  if (umin[0]>u0) umin[0] = u0;
	  if (vmin[0]>v0) vmin[0] = v0;
	  nnz[0]++;
	  break;
	}
	(xcoord[1])[u0+v0*outwidth] = xcd;
	(ycoord[1])[u0+v0*outwidth] = ycd;
	if (xcd>=0.0 || ycd>=0.0) {
	  if (umax[1]<u0) umax[1] = u0;
	  if (vmax[1]<v0) vmax[1] = v0;
	  if (umin[1]>u0) umin[1] = u0;
	  if (vmin[1]>v0) vmin[1] = v0;
	  nnz[1]++;
	  break;
	}
	if (ruv<=1.0) {
	    xint.push_back (u0);
	    yint.push_back (v0);
	}
	break;
      }
    }
  if (xint.size()>0) {
    int i = xint.size();
    printf ("***BUG***\n  Interpolation needed for %d pixels\n", i);
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
  float ap = aperture/180.0*PI_NUMB;

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


