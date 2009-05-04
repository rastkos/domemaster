#ifndef MAIN_H
#define MAIN_H

#include "FreeImage.h"

inline BYTE CastFloat (float);

inline BYTE CastFloat (float f)
{
  if (f<0.0)
    return (0);
  
  if (f>255.0)
    return (255);

  return ((BYTE)roundf(f));
}  

enum OUTPUT_TYPE  { OUTPUT_PNG=0, OUTPUT_TARGA, OUTPUT_BMP, OUTPUT_JPEG };

class Output {
 public:
  const char *name;
  FREE_IMAGE_FORMAT fif;
  int flags;
  bool alpha;
};
  

#endif
