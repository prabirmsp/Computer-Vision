#include <stdlib.h>
class MonoImage {
public:
  MonoImage(R2Image& image);
  MonoImage(MonoImage& image);
  ~MonoImage(void);

  double* operator[](int x);
  double get(int x, int y);
  void set(int x, int y, double val);

private:
  int width;
  int height;
  double * pixels;
};

MonoImage::
MonoImage(R2Image& image) {
  width = image.Width();
  height = image.Height();
  const int size = width * height;
  pixels = (double *) malloc(sizeof(double) * size);
  assert(pixels);

  double* c ;
  double sum;
  // compute averge greyscale
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      c = image[i][j].Components();
      sum = 0;
      for (int k = 0; k < 3; k++) {
        sum += c[k];
      }
      pixels[i*height + j] = sum / 3;
    }
  }
}

MonoImage::
MonoImage(MonoImage& image) {
  width = image.width;
  height = image.height;
  const int size = width * height;
  pixels = new double[size];
  assert(pixels);
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      pixels[i*height + j] = image[i][j];
    }
  }
}

MonoImage::
~MonoImage(void) {
  if (pixels) free(pixels);
}

double* MonoImage::
operator[](int x) {
  return &pixels[x*height];
}

double MonoImage::
get(int x, int y) {
  return pixels[x*height + y];
}

void MonoImage::
set(int x, int y, double val) {
  pixels[x* height + y] = val;
}
