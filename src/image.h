#ifndef IMAGE_H
#define IMAGE_H

class Image {
 public:
  Image (const char *);
  ~Image ();

  int width;
  int height;
  int bpp;

  float *red;
  float *green;
  float *blue;
  float *alpha;

 private:
  Image ();
};

#endif
