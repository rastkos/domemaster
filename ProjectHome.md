Rendering images for planetarium shows requires rendering software—e.g. Maya, 3dsmax—that supports fisheye projection. With the popular open source project Blender such option does not exists (at least as yet). One option is to use spherical reflective surface, but this causes problems with some of the visual effects. Alternatively one can use a six camera setup ("hemicube"), and stitch them together to create a 180 degree hemisphere. Domemaster is a free open source tool to do the latter.

Domemaster:
  * uses [FreeImage](http://freeimage.sourceforge.net/) for loading and saving images
  * developed on Linux, but possibly compiles on Mac without modifications
  * handles alpha channel
  * is reasonably fast even though not multithreaded
  * allows modifying the saturation of colours

## Usage ##
The program recognizes camera positions by file prefix:
  * f`_`  front, b`_` back, u`_` up, d`_` down, l`_` left, r`_` right
  * only the front image is required. The file name after the prefix determines the name of the output image: e.g if the input frame is f\_test0005.jpeg, then the output image is called test0005.png (png is the default output format)
  * all input images should be square and cover exactly 90 degrees vertically and horizontally

Examples:
  * domemaster `*`.png
    * The command shell will probably choke on large number of input files... eventually domemaster will need to parse directories itself.
  * domemaster -h
```
USAGE:
  domemaster [options] images

Available options:
  -a <angle>  Aperture, fulldome is 180 degrees [180]
  -b <value>  Modify brightness
  -c <value>  Modify contrast
  -d          Disable alpha channel
  -g <value>  Apply gamma correction to the image
  -h          Print this help message
  -i          Turn interpolation off
  -k <kernel> Select interpolation kernel. Possible values are
              tanh, sinc, sinc2, lanczos, hamming, hann. [sinc]
  -o <format> Output file format, possible values png, bmp, targa, jpeg. [png]
  -p          Disable padding
  -r <angle>  Rotation angle [0.0]
  -s <size>   Size of the output image. [1408]
  -t <angle>  Tilt angle [0.0]
  -v          Verbose messages
```