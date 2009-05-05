#ifndef IMAGE_H
#define IMAGE_H

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
  ~Image ();

  int width;
  int height;
  int bpp;

  float *red;
  float *green;
  float *blue;
  float *alpha;
  PROJTYPE proj;
  const char* imname;

 private:
  Image ();
};

#endif
