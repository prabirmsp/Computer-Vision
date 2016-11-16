// Source file for image class



// Include files

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include <iostream> // Include for I/O streams
//#include <tgmath.h>
#include <math.h>
#include "MonoImage.cpp"
#include <cmath>
#include <limits>

#include <queue>
#include <vector>
#include <unordered_map>
#include <cstdlib>
#include <ctime>


////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0),
    height(0)
{
}



R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0),
    height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width),
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width),
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width),
    height(image.height)

{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}


void R2Image::
svdTest(void)
{
	// fit a 2D conic to five points
	R2Point p1(1.2,3.5);
	R2Point p2(2.1,2.2);
	R2Point p3(0.2,1.6);
	R2Point p4(0.0,0.5);
	R2Point p5(-0.2,4.2);

	// build the 5x6 matrix of equations
	double** linEquations = dmatrix(1,5,1,6);

	linEquations[1][1] = p1[0]*p1[0];
	linEquations[1][2] = p1[0]*p1[1];
	linEquations[1][3] = p1[1]*p1[1];
	linEquations[1][4] = p1[0];
	linEquations[1][5] = p1[1];
	linEquations[1][6] = 1.0;

	linEquations[2][1] = p2[0]*p2[0];
	linEquations[2][2] = p2[0]*p2[1];
	linEquations[2][3] = p2[1]*p2[1];
	linEquations[2][4] = p2[0];
	linEquations[2][5] = p2[1];
	linEquations[2][6] = 1.0;

	linEquations[3][1] = p3[0]*p3[0];
	linEquations[3][2] = p3[0]*p3[1];
	linEquations[3][3] = p3[1]*p3[1];
	linEquations[3][4] = p3[0];
	linEquations[3][5] = p3[1];
	linEquations[3][6] = 1.0;

	linEquations[4][1] = p4[0]*p4[0];
	linEquations[4][2] = p4[0]*p4[1];
	linEquations[4][3] = p4[1]*p4[1];
	linEquations[4][4] = p4[0];
	linEquations[4][5] = p4[1];
	linEquations[4][6] = 1.0;

	linEquations[5][1] = p5[0]*p5[0];
	linEquations[5][2] = p5[0]*p5[1];
	linEquations[5][3] = p5[1]*p5[1];
	linEquations[5][4] = p5[0];
	linEquations[5][5] = p5[1];
	linEquations[5][6] = 1.0;

	printf("\n Fitting a conic to five points:\n");
	printf("Point #1: %f,%f\n",p1[0],p1[1]);
	printf("Point #2: %f,%f\n",p2[0],p2[1]);
	printf("Point #3: %f,%f\n",p3[0],p3[1]);
	printf("Point #4: %f,%f\n",p4[0],p4[1]);
	printf("Point #5: %f,%f\n",p5[0],p5[1]);

	// compute the SVD
	double singularValues[7]; // 1..6
	double** nullspaceMatrix = dmatrix(1,6,1,6);
	svdcmp(linEquations, 5, 6, singularValues, nullspaceMatrix);

	// get the result
	printf("\n Singular values: %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6]);

	// find the smallest singular value:
	int smallestIndex = 1;
	for(int i=2;i<7;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

	// solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
	printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex]);

	// make sure the solution is correct:
	printf("Equation #1 result: %f\n",	p1[0]*p1[0]*nullspaceMatrix[1][smallestIndex] +
										p1[0]*p1[1]*nullspaceMatrix[2][smallestIndex] +
										p1[1]*p1[1]*nullspaceMatrix[3][smallestIndex] +
										p1[0]*nullspaceMatrix[4][smallestIndex] +
										p1[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #2 result: %f\n",	p2[0]*p2[0]*nullspaceMatrix[1][smallestIndex] +
										p2[0]*p2[1]*nullspaceMatrix[2][smallestIndex] +
										p2[1]*p2[1]*nullspaceMatrix[3][smallestIndex] +
										p2[0]*nullspaceMatrix[4][smallestIndex] +
										p2[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #3 result: %f\n",	p3[0]*p3[0]*nullspaceMatrix[1][smallestIndex] +
										p3[0]*p3[1]*nullspaceMatrix[2][smallestIndex] +
										p3[1]*p3[1]*nullspaceMatrix[3][smallestIndex] +
										p3[0]*nullspaceMatrix[4][smallestIndex] +
										p3[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #4 result: %f\n",	p4[0]*p4[0]*nullspaceMatrix[1][smallestIndex] +
										p4[0]*p4[1]*nullspaceMatrix[2][smallestIndex] +
										p4[1]*p4[1]*nullspaceMatrix[3][smallestIndex] +
										p4[0]*nullspaceMatrix[4][smallestIndex] +
										p4[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #5 result: %f\n",	p5[0]*p5[0]*nullspaceMatrix[1][smallestIndex] +
										p5[0]*p5[1]*nullspaceMatrix[2][smallestIndex] +
										p5[1]*p5[1]*nullspaceMatrix[3][smallestIndex] +
										p5[0]*nullspaceMatrix[4][smallestIndex] +
										p5[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	R2Point test_point(0.34,-2.8);

	printf("A point off the conic: %f\n",	test_point[0]*test_point[0]*nullspaceMatrix[1][smallestIndex] +
											test_point[0]*test_point[1]*nullspaceMatrix[2][smallestIndex] +
											test_point[1]*test_point[1]*nullspaceMatrix[3][smallestIndex] +
											test_point[0]*nullspaceMatrix[4][smallestIndex] +
											test_point[1]*nullspaceMatrix[5][smallestIndex] +
											nullspaceMatrix[6][smallestIndex]);

	return;
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
  // Brighten the image by multiplying each pixel component by the factor.
  // This is implemented for you as an example of how to access and set pixels
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();
    }
  }
}

// Weights for the 3x3 kernel
int SobelWeight(int x, int y) {
  return x * ((y == 0) ? 2 : 1);
}

void R2Image::
SobelX(void)
{
  // Apply the Sobel oprator to the image in X direction
  R2Image oldR2Image(*this);
  R2Pixel halfPixel(0.5,0.5,0.5,1);
  for (int i = 1; i < width - 1; i++) {
    for (int j = 1; j < height - 1; j++) {
      Pixel(i,j).Reset(0,0,0,1);
      for (int lx = -1; lx < 2; lx++) {
        for (int ly = -1; ly < 2; ly++ ) {
          Pixel(i,j) += oldR2Image.Pixel(i + lx, j + ly) * SobelWeight(lx, ly);
        }
      }
      //Pixel(i, j) += halfPixel; // add half a pixel for viewing purpose
      //Pixel(i, j).Clamp(); // clamp for viewing purpose
    }
  }
  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  //fprintf(stderr, "SobelX() not implemented\n");
}

void R2Image::
SobelY(void)
{
	// Apply the Sobel oprator to the image in Y direction
  R2Image oldR2Image(width, height, Pixels());
  R2Pixel halfPixel(0.5,0.5,0.5,1);
  for (int i = 1; i < width - 1; i++) {
    for (int j = 1; j < height - 1; j++) {
      Pixel(i,j).Reset(0,0,0,1);
      for (int lx = -1; lx < 2; lx++) {
        for (int ly = -1; ly < 2; ly++ ) {
          Pixel(i,j) += oldR2Image.Pixel(i + lx, j + ly) * SobelWeight(ly, lx);
        }
      }
      //Pixel(i, j) += halfPixel; // add half a pixel for viewing purpose
      //Pixel(i, j).Clamp();
    }
  }
  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  //fprintf(stderr, "SobelY() not implemented\n");
}

void R2Image::
LoG(void)
{
  // Apply the LoG oprator to the image

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
  // Changes the saturation of an image
  // Find a formula that changes the saturation without affecting the image brightness

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "ChangeSaturation(%g) not implemented\n", factor);
}


// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{
  // Gaussian blur of the image. Separable solution is preferred

  // Initialize weights
  int sigmaInt = (int) (sigma + 0.5);
  int kernelSize = (6 * sigmaInt) + 1;
  double weights[kernelSize];
  double sum = 0;
  for (int i = 0; i < kernelSize; i++) {
    int x = i - (3 * sigmaInt);
    double num = exp(-1 * x * x / (2 * sigma * sigma));
    double denom = pow(2.0 * M_PI, 0.5) * sigma;
    double weight = num / denom;
    weights[i] = weight;
    sum += weight;
  }

  // Normalize kernel
  for (int i = 0; i < kernelSize; i++) {
    weights[i] /= sum;
  }

  // Temporary image
  R2Image firstPass(width, height);

  // First pass - x
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      R2Pixel cur;
      for(int lx = -3 * sigmaInt; lx <= 3 * sigmaInt; lx++) {
        int x = std::min(std::max(i + lx, 0), width - 1);
        cur += Pixel(x, j) * weights[lx + (3 * sigmaInt)];
      }
      firstPass.SetPixel(i, j, cur);
    }
  }

  // Second pass - y
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      R2Pixel cur;
      for(int ly = -3 * sigmaInt; ly <= 3 * sigmaInt; ly++) {
        int y = std::min(std::max(j + ly, 0), height - 1);
        cur += firstPass.Pixel(i, y) * weights[ly + (3 * sigmaInt)];
      }
      SetPixel(i, j, cur);
    }
  }
}

R2Image ComputeHarrisImage(R2Image* orig, double sigma) {
  int width = orig->Width();
  int height = orig->Height();
  R2Image img(width,height);

  // Initialize the temporary images
  R2Image Ix_sq(*orig);
  R2Image Iy_sq(*orig);
  R2Image Ix_Iy(width, height);

  Ix_sq.SobelX();
  Iy_sq.SobelY();

  R2Pixel ix, iy;
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      ix = Ix_sq[i][j];
      iy = Iy_sq[i][j];
      Ix_Iy[i][j] = ix * iy;
      Ix_sq[i][j] = ix * ix;
      Iy_sq[i][j] = iy * iy;
    }
  }

  // Blur (weighted average)
  Ix_sq.Blur(sigma);
  Iy_sq.Blur(sigma);
  Ix_Iy.Blur(sigma);
  R2Pixel grey(0.5,0.5,0.5,1);
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      R2Pixel temp = (Ix_sq[i][j] * Iy_sq[i][j])
                    - (Ix_Iy[i][j] * Ix_Iy[i][j])
                    - (0.04 * (Ix_sq[i][j] + Iy_sq[i][j])
                            * (Ix_sq[i][j] + Iy_sq[i][j]));
      //              + grey;
      //temp.Clamp();
      img.SetPixel(i, j, temp);
    }
  }
  return img;
}

typedef struct Point {
  int x;
  int y;
  float val;
} Point;


struct ComparePoint {
    bool operator()(Point const & p1, Point const & p2) {
        // return "true" if "p1" is ordered before "p2", for example:
        return p1.val < p2.val;
    }
};


void mark (R2Image &image, Point p, R2Pixel color)  {
  const int radius = 5;
  R2Pixel white (1, 1, 1, 1);
  for (int i = -1* radius; i <= radius; i++) {
    for (int j = -1* radius; j <= radius; j++) {
      int x = p.x + i;
      int y = p.y + j;
      if (x >= 0 && x < image.Width() && y >= 0 && y < image.Height()) {
        if(std::abs(i) < radius - 1 && std::abs(j) < radius - 1 )
          image.SetPixel(x, y, color);
        else
          image.SetPixel(x, y, white);
      }
    }
  }
}

void line (R2Image &image, Point p1, Point p2, R2Pixel color)  {
  if (p1.x == p2.x) {
    for (int i = std::min(p1.y, p2.y); i <= std::max(p1.y, p2.y); i++) {
      image.SetPixel(p1.x, i, color);
    }
  } else {
    double m = ((double) (p2.y - p1.y)) / ((double) (p2.x - p1.x));
    for (int i = std::min(p1.x, p2.x); i <= std::max(p1.x, p2.x); i++) {
      for (int j = std::min(p1.y, p2.y); j <= std::max(p1.y, p2.y); j++) {
        double dist = abs((double)(j - p1.y) - (m * (i - p1.x)));
        if (dist < 2.0) {
          image.SetPixel(i, j, color);
        }
      }
    }
  }
}


void getFeaturePoints(R2Image* harris, std::vector<Point> &points, int numFeaturePoints, int borderx, int bordery) {
  int width = harris->Width();
  int height = harris->Height();
  std::priority_queue<Point, std::vector<Point>, ComparePoint> q;
  bool valid[width][height];
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      valid[i][j] = true;
    }
  }

  for (int i = borderx; i < width - borderx; i++) {
    for (int j = bordery; j < height - bordery; j++) {
      R2Pixel cur = harris->Pixel(i,j);
      double sum = cur.Red() + cur.Green() + cur.Blue();
      Point p;
      p.x = i;
      p.y = j;
      p.val = sum;
      q.push(p);
    }
  }

  const int invalidRadius = 10;
  int pointCount = 0;
  while(pointCount < numFeaturePoints) {
    Point p = q.top();
    if (valid[p.x][p.y]) {
      points[pointCount] = p;
      //printf("-%d,%d\n", p.x,p.y);
      for (int i = -1 * invalidRadius; i <= invalidRadius; i++) {
        for (int j = -1 * invalidRadius; j <= invalidRadius; j++) {
          int x = std::min(std::max(0, p.x + i), width - 1);
          int y = std::min(std::max(0, p.y + j), height - 1);
          valid[x][y] = false;
        }
      }
      pointCount++;
    }
    q.pop();
  }
}

void R2Image::
Harris(double sigma)
{
    // Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
	// Output should be 50% grey at flat regions, white at corners and black/dark near edges
  R2Image harris = ComputeHarrisImage(this, sigma);
  printf("Harris Image Computed\n");

  const int numFeaturePoints = 150;
  std::vector<Point> points(numFeaturePoints);
  getFeaturePoints(&harris, points, numFeaturePoints, 0, 0);
  for(unsigned int i = 0; i < numFeaturePoints; i++) {
    mark(*this, points[i], R2Pixel(1, 0, 0, 1));
  }
  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  //fprintf(stderr, "Harris(%g) not implemented\n", sigma);
}

double computeSSD(MonoImage& I_0, MonoImage& I_1, Point p, int searchx, int searchy, const int ssdWindowRadius) {
  double sum = 0;
  for (int i = -1 * ssdWindowRadius; i <= ssdWindowRadius; i++) {
    for (int j = -1 * ssdWindowRadius; j <= ssdWindowRadius; j++) {
      double diff = I_0[p.x + i][p.y + j]
                   - I_1[searchx + i][searchy + j];
      sum += diff * diff;
    }
  }
  return sum;
}


void  track(R2Image * thisImage, R2Image * otherImage, int numFeaturePoints, std::vector<Point> &points, std::unordered_map<int,Point> &matchedPoints) {
  int width = thisImage->Width();
  int height = thisImage->Height();

  // Options to change
  const int sigma = 2;
  const int ssdWindowRadius = 6 * sigma + 1;
  const int windowx = (int) (0.2f * width);
  const int windowy = (int) (0.2f * height);

  printf("Computing Harris Image... ");
  R2Image harris = ComputeHarrisImage(thisImage, 2);
  printf("Completed\n");

  printf("Finding Feature Points... ");
  getFeaturePoints(&harris, points, numFeaturePoints, 0, 0);
  //getFeaturePoints(&harris, points, numFeaturePoints, windowx / 2, windowy / 2);
  printf("Completed\n");

  MonoImage curMono(*thisImage);
  MonoImage otherMono(*otherImage);

  const double infinity = std::numeric_limits<double>::infinity();

  for(int ip = 0; ip < numFeaturePoints; ip++) {
    printf("Tracking points... %d%%", ip * 100 / numFeaturePoints);
    Point p = points[ip];
    int startx = std::max(ssdWindowRadius, p.x - (windowx / 2));
    int starty = std::max(ssdWindowRadius, p.y - (windowy / 2));
    int endx = std::min(otherImage->Width() - ssdWindowRadius, p.x + (windowx / 2));
    int endy = std::min(otherImage->Height() - ssdWindowRadius, p.y + (windowy / 2));
    Point best;
    best.val = infinity;
    for (int i = startx; i < endx; i++) {
      for (int j = starty; j < endy; j++) {
        double ssd = computeSSD(curMono, otherMono, p, i, j, ssdWindowRadius);
        if (ssd < best.val) {
          best.val = ssd;
          best.x = i;
          best.y = j;
        }
      }
    }
    if(best.val != infinity)
      matchedPoints[ip] = best;
    std::cout<<'\r';
    std::cout.flush();
  }
  printf("Tracking Points... Completed\n");
  return;
}


void R2Image::
trackMovement(R2Image * otherImage) {

  const int numFeaturePoints = 100;

  std::vector<Point> points(numFeaturePoints);
  std::unordered_map<int, Point> matchedPoints(numFeaturePoints);
  track(this, otherImage, numFeaturePoints, points, matchedPoints);
  *this = *otherImage;

  R2Pixel red (1, 0, 0, 1);
  for (int i = 0; i < numFeaturePoints; i++) {
    if (matchedPoints.count(i) > 0) {
      line(*this, points[i], matchedPoints[i], red);
      mark(*this, matchedPoints[i], red);
    }
  }
}

void ransac() {

}

void R2Image::
trackMovementRansac(R2Image * otherImage) {

  const int numFeaturePoints = 100;

  std::vector<Point> points(numFeaturePoints);
  std::unordered_map<int, Point> matchedPoints(numFeaturePoints);
  track(this, otherImage, numFeaturePoints, points, matchedPoints);
  *this = *otherImage;

  const int numIterations = 100;
  const int numSubsetPoints = 4;
  const double acceptThreshold = 5;

  int bestNumMatches = -1;
  double bestAverageDx = 0;
  double bestAverageDy = 0;

  for (int iteration = 0; iteration < numIterations; iteration++) {
    int subset [numSubsetPoints];
    double averageDx = 0;
    double averageDy = 0;
    for (int i = 0; i < numSubsetPoints; i++) {
      do {
        subset[i] = std::rand() % numFeaturePoints;
      } while (matchedPoints.count(subset[i]) == 0);
      averageDx += matchedPoints[subset[i]].x - points[subset[i]].x;
      averageDy += matchedPoints[subset[i]].y - points[subset[i]].y;
    }
    averageDx /= numSubsetPoints;
    averageDy /= numSubsetPoints;

    int numMatches = 0;
    for (int i = 0; i < numFeaturePoints; i++) {
      if (matchedPoints.count(subset[i]) > 0) {
        double deltax = matchedPoints[subset[i]].x - points[subset[i]].x;
        double deltay = matchedPoints[subset[i]].y - points[subset[i]].y;
        if (abs(deltax - averageDx) <= acceptThreshold
         && abs(deltay - averageDy) <= acceptThreshold) {
           numMatches ++;
         }
      }
    }
    if (numMatches > bestNumMatches) {
      bestNumMatches = numMatches;
      bestAverageDx = averageDx;
      bestAverageDy = averageDy;
    }
  }

  const R2Pixel red(1, 0, 0, 1);
  const R2Pixel green(0, 1, 0, 1);
  for (int i = 0; i < numFeaturePoints; i++) {
    if (matchedPoints.count(i) > 0) {
      double deltax = matchedPoints[i].x - points[i].x;
      double deltay = matchedPoints[i].y - points[i].y;
      if (abs(deltax - bestAverageDx) <= acceptThreshold
       && abs(deltay - bestAverageDy) <= acceptThreshold) {
         line(*this, points[i], matchedPoints[i], green);
         mark(*this, matchedPoints[i], green);
       } else{
         line(*this, points[i], matchedPoints[i], red);
         mark(*this, matchedPoints[i], red);
       }
    }
  }
}

void dlt(int numSubsetPoints, int * subset, std::vector<Point> &points, std::unordered_map<int,Point> &matchedPoints, double * H) {
  double** leq = dmatrix(1, 2 * numSubsetPoints, 1, 9);
  for (int i = 0; i < numSubsetPoints; i++) {
    // first line
    leq[2*(i+1)-1][1] = 0;
    leq[2*(i+1)-1][2] = 0;
    leq[2*(i+1)-1][3] = 0;
    leq[2*(i+1)-1][4] = -1 * points[subset[i]].x;
    leq[2*(i+1)-1][5] = -1 * points[subset[i]].y;
    leq[2*(i+1)-1][6] = -1;
    leq[2*(i+1)-1][7] = matchedPoints[subset[i]].y * points[subset[i]].x;
    leq[2*(i+1)-1][8] = matchedPoints[subset[i]].y * points[subset[i]].y;
    leq[2*(i+1)-1][9] = matchedPoints[subset[i]].y;
    // second line
    leq[2*(i+1)][1] = points[subset[i]].x;
    leq[2*(i+1)][2] = points[subset[i]].y;
    leq[2*(i+1)][3] = 1;
    leq[2*(i+1)][4] = 0;
    leq[2*(i+1)][5] = 0;
    leq[2*(i+1)][6] = 0;
    leq[2*(i+1)][7] = -1 * points[subset[i]].x * matchedPoints[subset[i]].x;
    leq[2*(i+1)][8] = -1 * matchedPoints[subset[i]].x * points[subset[i]].y;
    leq[2*(i+1)][9] = -1 * matchedPoints[subset[i]].x;
  }

  double singularValues[10];
	double** nullspaceMatrix = dmatrix(1,9,1,9);
	svdcmp(leq, 2 * numSubsetPoints, 9, singularValues, nullspaceMatrix);

  // find the smallest singular value:
  int smallestIndex = 1;
  for(int i=2;i<10;i++)
    if(singularValues[i]<singularValues[smallestIndex])
      smallestIndex=i;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      H[i*3 + j] = nullspaceMatrix[i*3+j+1][smallestIndex];
    }
  }
}

void R2Image::
dltTest() {

      std::vector<Point> ps(4);
      Point p1;
      p1.x = 0;
      p1.y = 0;
      ps[0] = p1;
      Point p2;
      p2.x = 1;
      p2.y = 0;
      ps[1] = p2;
      Point p3;
      p3.x = 1;
      p3.y = 1;
      ps[2] = p3;
      Point p4;
      p4.x = 0;
      p4.y = 1;
      ps[3] = p4;
      std::unordered_map<int, Point> mps(4);
      Point p5;
      p5.x = 1;
      p5.y = 2;
      mps[0] = p5;
      Point p6;
      p6.x = 1;
      p6.y = 1;
      mps[1] = p6;
      Point p7;
      p7.x = 3;
      p7.y = 1;
      mps[2] = p7;
      Point p8;
      p8.x = 3;
      p8.y = 2;
      mps[3] = p8;
      double HH[9];
      int ss [4];
      for (int i = 0; i < 4; i ++)
      ss[i] = i;
      dlt(4, ss, ps, mps, HH);


        printf("Matrix H:\n");
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            printf("%2.2lf\t", HH[i*3+j]);
          }
          printf("\n");
        }


          Point p = p1;
          double w = HH[6] * p.x + HH[7] * p.y + HH[8];
          double y = (HH[3] * p.x + HH[4] * p.y + HH[5]) / w;
          double x = (HH[0] * p.x + HH[1] * p.y + HH[2]) / w;

          printf("%lf, %lf, %lf\n", w,x,y);
          printf("%d, %d\n", (int)(x + 0.5),(int) (y + 0.5));
        return;




  R2Point points[4][2] =
    {{ R2Point(0,0), R2Point(1, 2)},
     { R2Point(1,0), R2Point(1,1)},
     { R2Point(1,1), R2Point(3,1)},
     { R2Point(0,1), R2Point(3,2)}};
  double** leq = dmatrix(1, 8, 1, 9);
  for (int i = 1; i <= 4; i++) {
    // first line
    leq[2*i-1][1] = 0;
    leq[2*i-1][2] = 0;
    leq[2*i-1][3] = 0;
    leq[2*i-1][4] = -1 * points[i-1][0][0];
    leq[2*i-1][5] = -1 * points[i-1][0][1];
    leq[2*i-1][6] = -1;
    leq[2*i-1][7] = points[i-1][1][1] * points[i-1][0][0];
    leq[2*i-1][8] = points[i-1][1][1] * points[i-1][0][1];
    leq[2*i-1][9] = points[i-1][1][1];
    // second line
    leq[2*i][1] = points[i-1][0][0];
    leq[2*i][2] = points[i-1][0][1];
    leq[2*i][3] = 1;
    leq[2*i][4] = 0;
    leq[2*i][5] = 0;
    leq[2*i][6] = 0;
    leq[2*i][7] = -1 * points[i-1][0][0] * points[i-1][1][0];
    leq[2*i][8] = -1 * points[i-1][1][0] * points[i-1][0][1];
    leq[2*i][9] = -1 * points[i-1][1][0];
  }

  double singularValues[10];
	double** nullspaceMatrix = dmatrix(1,9,1,9);
	svdcmp(leq, 8, 9, singularValues, nullspaceMatrix);

  // find the smallest singular value:
  int smallestIndex = 1;
  for(int i=2;i<10;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;
  printf("Singular Values:\n");
  for(int i = 1; i < 10; i++) {
    printf("%2.2lf\t", singularValues[i]);
  }
  printf("Smallest singular value: %lf i:%d\n", singularValues[smallestIndex], smallestIndex);

  printf("Matrix H:\n");
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      printf("%2.2lf\t", nullspaceMatrix[i*3+j+1][smallestIndex]);
    }
    printf("\n");
  }
}

int rounddtoi (double x) {
  return (int) (x + 0.5);
}

void R2Image::
trackMovementDltRansac(R2Image * otherImage) {

  const int numFeaturePoints = 100;

  std::vector<Point> points(numFeaturePoints);
  std::unordered_map<int, Point> matchedPoints(numFeaturePoints);
  track(this, otherImage, numFeaturePoints, points, matchedPoints);
  *this = *otherImage;

  const int numIterations = 100;
  const int numSubsetPoints = 4;
  const double acceptThreshold = 5; // pixels

  int bestNumMatches = -1;
  double bestH [9];

  for (int iteration = 0; iteration < numIterations; iteration++) {
    int subset [numSubsetPoints];
    for (int i = 0; i < numSubsetPoints; i++) {
      do {
        subset[i] = std::rand() % numFeaturePoints;
      } while (matchedPoints.count(subset[i]) == 0);
    }

    double H [9];
    dlt(numSubsetPoints, subset, points, matchedPoints, H);
    int numMatches = 0;
    for (int i = 0; i < numFeaturePoints; i++) {
      if (matchedPoints.count(subset[i]) > 0) {
        Point p = points[subset[i]];
        double w = H[6] * p.x + H[7] * p.y + H[8];
        int x = rounddtoi((H[0] * p.x + H[1] * p.y + H[2]) / w);
        int y = rounddtoi((H[3] * p.x + H[4] * p.y + H[5]) / w);
        if (abs(x - matchedPoints[subset[i]].x) <= acceptThreshold
         && abs(y - matchedPoints[subset[i]].y) <= acceptThreshold) {
           numMatches ++;
         }
      }
    }
    if (numMatches > bestNumMatches) {
      bestNumMatches = numMatches;
      for (int i = 0; i < 9; i++)
        bestH[i] = H[i];
    }
  }

  const R2Pixel red(1, 0, 0, 1);
  const R2Pixel green(0, 1, 0, 1);
  for (int i = 0; i < numFeaturePoints; i++) {
    if (matchedPoints.count(i) > 0) {
      Point p = points[i];
      double w = bestH[6] * p.x + bestH[7] * p.y + bestH[8];
      int x = rounddtoi((bestH[0] * p.x + bestH[1] * p.y + bestH[2]) / w);
      int y = rounddtoi((bestH[3] * p.x + bestH[4] * p.y + bestH[5]) / w);
      if (abs(x - matchedPoints[i].x) <= acceptThreshold
       && abs(y - matchedPoints[i].y) <= acceptThreshold) {
         line(*this, points[i], matchedPoints[i], green);
         mark(*this, matchedPoints[i], green);
       } else{
         line(*this, points[i], matchedPoints[i], red);
         mark(*this, matchedPoints[i], red);
       }
    }
  }
}


/*
Kernel:

 0 -1  0
-1  5 -1
 0 -1  0  */
int SharpenKernel(int x, int y) {
  if (x == 0 && y == 0)
    return 5;
  else if (x == 0 || y == 0)
    return -1;
  return 0;
}

void R2Image::
Sharpen()
{
  // Sharpen an image using a linear filter. Use a kernel of your choosing.
  R2Image oldR2Image(width, height, Pixels());
  R2Pixel halfPixel(0.5,0.5,0.5,1);
  for (int i = 1; i < width - 1; i++) {
    for (int j = 1; j < height - 1; j++) {
      Pixel(i,j).Reset(0,0,0,1);
      for (int lx = -1; lx < 2; lx++) {
        for (int ly = -1; ly < 2; ly++ ) {
          Pixel(i,j) += oldR2Image.Pixel(i + lx, j + ly) * SharpenKernel(lx, ly);
        }
      }
      Pixel(i, j).Clamp();
    }
  }
  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  //fprintf(stderr, "Sharpen() not implemented\n");
}


void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching translation (pixel precision is OK), and blend the translated "otherImage"
	// into this image with a 50% opacity.
	fprintf(stderr, "fit other image using translation and blend imageB over imageA\n");
	return;
}

void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching homography, and blend the transformed "otherImage" into this image with a 50% opacity.
	fprintf(stderr, "fit other image using a homography and blend imageB over imageA\n");
	return;
}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp);
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);

  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);

  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);

  // Check info header
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }

  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }

  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
	// Read pixel values
	int red, green, blue;
	if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
	  fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
	  fclose(fp);
	  return 0;
	}

	// Assign pixel values
	double r = (double) red / max_value;
	double g = (double) green / max_value;
	double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }

    // Close file
    fclose(fp);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" {
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}




int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width; 	/* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 95, TRUE);
  jpeg_start_compress(&cinfo, TRUE);

  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}
