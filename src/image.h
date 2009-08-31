#ifndef IMAGE_H
#define IMAGE_H

#include "FreeImage.h"

enum PROJTYPE {UNKNOWN=0,FRONT, BACK, LEFT, RIGHT, UP, DOWN};
enum OUTPUT_TYPE  { OUTPUT_PNG=0, OUTPUT_TARGA, OUTPUT_BMP, OUTPUT_JPEG };

class Output {
 public:
  const char *name;
  FREE_IMAGE_FORMAT fif;
  int flags;
  bool alpha;
};

class Image {
 public:
  Image (const char *, int);
  Image (int, int, int);
  ~Image ();
  void Modulate (float, float, float, float);
  void Histogram (float);
  void ToHSB ();
  void ToRGB ();
  void HistogramEqualize (Image *, float);

  int width;
  int height;
  int bpp;

  float *red;
  float *green;
  float *blue;  
  float *alpha;
  float *saturation;
  float *brightness;
  float *hue;

  int *histogram;
  float *histtrans;
  PROJTYPE proj;
  const char* imname;
  int cpi;
  int pad;

  bool operator == (const Image & d );
  bool operator += (const Image & d );
 private:
  Image ();
};

#endif
