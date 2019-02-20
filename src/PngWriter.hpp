#ifndef PNG_WRITER_HPP
#define PNG_WRITER_HPP

// ======================================================================
// PngWriter.hpp
// 
// This class can write a png image using libpng routines (i.e. png.h).
//
// Sample usage:
//
//    PngWriter png(nx,ny);
//
//    // then one or more calls to...
//    
//    png.set(i, j, red, green, blue); // 0 <= red,green,blue <= 255
// 
//    png.write("myfile.png");
//
//    // at this point you can change the image and write again...
//
//    png.set(i, j, red, green, blue);
//
//    png.write("myfile2.png");
// 
//
// History:
//  authors: Frank Ham & Phuc Quang - July 2013
// ======================================================================

// ======================================================================
// Copyright (c) 2013 Frank Ham and Phuc Quang
//
// License: MIT (http://opensource.org/licenses/MIT)
//
// use for any purpose, commercial or otherwise.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
// ======================================================================

#include <png.h>
#include "CurvilinearInterpolator.hpp"
#include "AbstractCSolver.hpp"

class PngWriter {
  
  public:

    int mpiRank;  

    unsigned char (*buffer)[3]; // 0 <= r,g,b < 255 
    int nx,ny;
    int dumpInterval;

    double *fieldPtr;
    string varName;
    int planeInd;
    double fraction;

    double x_min[3], x_max[3];
    bool bboxFlag;

    bool interpolatorFlag;
    CurvilinearInterpolator *ci;  

    bool valFlag;
    double valMax, valMin;

    string timeStepString;

    enum Colormap {GREYSCALE,
		   RAINBOW, 
		   BWR};

    Colormap cm;

  //General constructor w/ floating value bounds
  PngWriter(int dumpInterval, const int width, const int height, double *fieldPtr, string varName, int planeInd, double fraction, Colormap cm){
    nx = width;
    ny = height;
    buffer = new unsigned char[nx*ny][3];

    // fill buffer with a "nice" cyan [73,175,205] -- eventually a designer should choose this ;)
    for (int i = 0; i < nx*ny; ++i) {
      buffer[i][0] = 73;
      buffer[i][1] = 175;
      buffer[i][2] = 205;
    }

    this->cm = cm;

    this->fieldPtr = fieldPtr;
    this->varName  = varName;
    this->planeInd = planeInd;
    this->fraction = fraction;

    interpolatorFlag = false;
    ci = NULL;

    valFlag = false;
    valMax  = 0.0;
    valMin  = 0.0;

    bboxFlag = false;

    this->dumpInterval = dumpInterval;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  }

  //Constructor for fixed value bounds
  PngWriter(int dumpInterval, const int width, const int height, double *fieldPtr, string varName, int planeInd, double fraction, double valMin, double valMax, Colormap cm){
    nx = width;
    ny = height;
    buffer = new unsigned char[nx*ny][3];

    // fill buffer with a "nice" cyan [73,175,205] -- eventually a designer should choose this ;)
    for (int i = 0; i < nx*ny; ++i) {
      buffer[i][0] = 73;
      buffer[i][1] = 175;
      buffer[i][2] = 205;
    }

    this->fieldPtr = fieldPtr;
    this->varName  = varName;
    this->planeInd = planeInd;
    this->fraction = fraction;

    this->cm = cm;

    interpolatorFlag = false;
    ci = NULL;

    valFlag = true;
    this->valMax  = valMax;
    this->valMin  = valMin;

    bboxFlag = false;

    this->dumpInterval = dumpInterval;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  }

  //Constructor for fixed value bounds and fixed bounds
  PngWriter(int dumpInterval, const int width, const int height, double *fieldPtr, string varName, int planeInd, double fraction, double valMin, double valMax, double x_min[3], double x_max[3], Colormap cm){
    nx = width;
    ny = height;
    buffer = new unsigned char[nx*ny][3];

    // fill buffer with a "nice" cyan [73,175,205] -- eventually a designer should choose this ;)
    for (int i = 0; i < nx*ny; ++i) {
      buffer[i][0] = 73;
      buffer[i][1] = 175;
      buffer[i][2] = 205;
    }

    this->fieldPtr = fieldPtr;
    this->varName  = varName;
    this->planeInd = planeInd;
    this->fraction = fraction;

    this->cm = cm;

    interpolatorFlag = false;
    ci = NULL;

    valFlag = true;
    this->valMax  = valMax;
    this->valMin  = valMin;

    this->dumpInterval = dumpInterval;

   
    bboxFlag = true;
    for(int ip = 0; ip < 3; ip++){
        this->x_min[ip] = x_min[ip];
        this->x_max[ip] = x_max[ip];
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  }



  ~PngWriter() {

    delete[] buffer;

  }
  
  void set(const int i,const int j,const unsigned char r,const unsigned char g,const unsigned char b) {
    // recall that for png files, the pixels are ordered from the top left, so modify 
    // this set routine's passed j so that zero is at the bottom left...
    buffer[(ny-j-1)*nx+i][0] = r;
    buffer[(ny-j-1)*nx+i][1] = g;
    buffer[(ny-j-1)*nx+i][2] = b;
  }
  
  void write(const char * filename) {

    // note: we completely skip any error handling treatment here for simplicity.
    
    FILE * fp = fopen(filename,"wb");
    if (!fp) {
      std::cout << "Warning: could not open png file: " << filename << std::endl;
      return;
    }
    
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
						  (png_voidp)NULL,NULL,NULL);
    if (!png_ptr) {
      fclose(fp);
      std::cout << "Warning: could not create png_ptr" << std::endl;
      return;
    }
    
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
      png_destroy_write_struct(&png_ptr,(png_infopp)NULL);
      std::cout << "Warning: could not create info_ptr" << std::endl;
      return;
    }
    
    png_init_io(png_ptr, fp);

    png_set_IHDR(png_ptr, info_ptr,
		 nx, ny,             // width, height
		 8,                  // bits per pixel -- 16 does not work with blockbuster
		 PNG_COLOR_TYPE_RGB, // non-alpha options are PNG_COLOR_TYPE_RGB,PNG_COLOR_TYPE_GRAY,
		 PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_DEFAULT,
		 PNG_FILTER_TYPE_DEFAULT);

    // Some bits per pixel notes: 16 does not work with blockbuster, and there are also 
    // issues with PNG_COLOR_TYPE_GRAY interpretation, so stick to 8 and PNG_COLOR_TYPE_RGB
    // for now. Note that if you do use 16, pay attention to MSB/LSB order. Endian is 
    // flipped on my linux workstation...
    
    png_write_info(png_ptr, info_ptr);
    
    // set up row pointers to point into the raw image data in buffer...
    
    png_byte * row_pointers[ny];
    for (int i = 0; i < ny; ++i)
      row_pointers[i] = (png_byte*)(buffer + i*nx);
    
    png_write_image(png_ptr, row_pointers);

    png_write_end(png_ptr, NULL);

    png_destroy_write_struct(&png_ptr, &info_ptr);

    fclose(fp);

    // leave the image data in buffer to allow the calling process to change it and/or
    // write another image. It is deleted in the destructor.

  }
  
  void getRainbowColormap(double f, int &r, int &g, int &b){

	int r_map[5] = {0,   0,   0,   255, 255};
	int g_map[5] = {0,   255, 255, 255, 0};
	int b_map[5] = {255, 255, 0,   0,   0};

	int chunk = (int) floor(f/0.2500001);

	double interp = (f-0.2500001*(double)chunk)/0.2500001;

	r = (int) ( (1.0-interp)*r_map[chunk] + interp*r_map[chunk+1]);
	g = (int) ( (1.0-interp)*g_map[chunk] + interp*g_map[chunk+1]);
	b = (int) ( (1.0-interp)*b_map[chunk] + interp*b_map[chunk+1]);
	
  }

  void getPARAVIEW_BWR(double f, int &r, int &g, int &b){

	double r_mapf[33] = {0.07514311,  0.247872569, 0.339526309, 0.409505078, 0.468487184, 0.520796675, 0.568724526, 0.613686735, \
		             0.656658579, 0.698372844, 0.739424025, 0.780330104, 0.821573924, 0.863634967, 0.907017747, 0.936129275, \
			     0.943467973, 0.990146732, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.961891484, 0.916482116};
;
	double g_mapf[33] = {0.468049805, 0.498782363, 0.528909511, 0.558608486, 0.588057293, 0.617435078, 0.646924167, 0.676713218, \
			     0.707001303, 0.738002964, 0.769954435, 0.803121429, 0.837809045, 0.874374691, 0.913245283, 0.938743558, \
			     0.943498599, 0.928791426, 0.88332677,  0.833985467, 0.788626485, 0.746206642, 0.70590052,  0.667019783, \
			     0.6289553,   0.591130233, 0.552955184, 0.513776083, 0.472800903, 0.428977855, 0.380759558, 0.313155629, 0.236630659};

	double b_mapf[33] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.983038586, 0.943398095, 0.917447482, 0.861943246, 0.803839606, \
			     0.750707739, 0.701389973, 0.654994046, 0.610806959, 0.568237474, 0.526775617, 0.485962266, 0.445364274, \
			     0.404551679, 0.363073592, 0.320428137, 0.265499262, 0.209939162};
;
	int r_map[33], g_map[33], b_map[33];

	for(int i = 0; i < 33; i++){
	    r_map[i] = (int)(255.0*r_mapf[i]);
	    g_map[i] = (int)(255.0*g_mapf[i]);
	    b_map[i] = (int)(255.0*b_mapf[i]);
	}

	double dd = 0.031250001;

	int chunk = (int) floor(f/dd);

	double interp = (f-dd*(double)chunk)/dd;

	r = (int) ( (1.0-interp)*r_map[chunk] + interp*r_map[chunk+1]);
	g = (int) ( (1.0-interp)*g_map[chunk] + interp*g_map[chunk+1]);
	b = (int) ( (1.0-interp)*b_map[chunk] + interp*b_map[chunk+1]);
	
  }



};

#endif
