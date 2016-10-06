
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
  double* pixels;
};

MonoImage::
MonoImage(R2Image& image) {
  width = image.Width();
  height = image.Height();
  pixels = (double *)malloc(sizeof(double) * width * height);
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
      pixels[i*width + j] = sum / 3;
    }
  }
}

MonoImage::
MonoImage(MonoImage& image) {
  width = image.width;
  height = image.height;
  pixels = (double *)malloc(sizeof(double) * width * height);
  assert(pixels);
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      pixels[i*width + j] = image[i][j];
    }
  }
}

MonoImage::
~MonoImage(void) {
  if (pixels) delete pixels;
}

double* MonoImage::
operator[](int x) {
  return &pixels[x*width];
}

double MonoImage::
get(int x, int y) {
  return pixels[x*width + y];
}

void MonoImage::
set(int x, int y, double val) {
  pixels[x* width + y] = val;
}
