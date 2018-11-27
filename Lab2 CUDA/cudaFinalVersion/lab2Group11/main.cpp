//===================================================================================================================
// File:        main.cpp
// Created on:  27-10-11
// Authors:     Dirk Vos, Mark Schrauwen, Michiel van der Vlag
//
// Description: The main file for calling the CUDA code. This code processes al the images in folder 'images' in the
// 				upper folder. It calculates the average processing time based on the value of AVERAGE_FILTER_SIZE.
//				This code first calls the sequential code and afterwards the CUDA code. The processed information is
//				put in the 'run' folder in the stdout.log file.
//===================================================================================================================

#include <CImg.h>
#include <Timer.hpp>
#include <iostream>
#include <iomanip>
#include <cstring>

//-------------------------------------------------------------------------------------------------------------------
// Local definitions, macros
//-------------------------------------------------------------------------------------------------------------------
// Use commenting to disable/enable functions
#define AVERAGE_FILTER_SIZE		5					// This number represent the size of the average calculation
#define NUMBER_IMAGES			16					// The total number of images that need to be processed
//#define IMAGE_SAVING
#define PATH					"../images/"		// Path with all the original images
#define PATHOUT					"./imagesOutput/"	// Output path for processed images

using cimg_library::CImg;
using LOFAR::NSTimer;
using std::cout;
using std::cerr;
using std::endl;
using std::fixed;
using std::setprecision;

// Constants
const unsigned int HISTOGRAM_SIZE     = 256;
const unsigned int CONTRAST_THRESHOLD = 80;

extern void getCudaDeviceInformation();
extern void rgb2gray(			unsigned char *inputImage, unsigned char *grayImage, const int width, const int height, double *totaltime);
extern void histogram1D(		unsigned char *grayImage, const int width, const int height,
                                unsigned int *histogram, const unsigned int HISTOGRAM_SIZE, double *totaltime);
extern void contrast1D(			unsigned char *grayImage, const int width, const int height, unsigned int *histogram,
                                const 	unsigned int HISTOGRAM_SIZE, const unsigned int CONTRAST_THRESHOLD, double *totaltime);
extern void triangularSmooth(	unsigned char *grayImage, unsigned char *smoothImage, const int width, const int height, double *totaltime);
extern void imageProcess(		unsigned char *inputImage, unsigned char *smoothImage, const int width, const int height,
                                const unsigned int HISTOGRAM_SIZE, const unsigned int CONTRAST_THRESHOLD, double *totaltime);

int main(int argc, char *argv[])
{
	// variables
	double totaltime[AVERAGE_FILTER_SIZE][6];
	double speedup;
	unsigned int i, j;

	const char* imageNames[NUMBER_IMAGES] = {	
												"image00.jpg",
	                                            "image01.jpg",
	                                            "image02.jpg",
	                                            "image03.jpg",
	                                            "image04.bmp",
	                                            "image05.jpg",
	                                            "image06.jpg",
	                                            "image07.jpg",
	                                            "image08.jpg",
	                                            "image09.bmp",
	                                            "image10.jpg",
	                                            "image11.jpg",
	                                            "image12.jpg",
	                                            "image13.jpg",
	                                            "image14.jpg",
	                                            "image15.jpg"
	                                        };

	CImg< unsigned char > smoothImage;
	CImg< unsigned char > grayImage;
	CImg< unsigned char > inputImage;
	unsigned int* histogram;

	// check arguments
	if ( argc != 1 ) {
		cerr << "Program does not take parameters " << endl;
		return 1;
	}

	// get device information and print it to output
	getCudaDeviceInformation();

	// Proces every image
	for (i = 0; i < NUMBER_IMAGES; i++) {
		// Concatenate correct imagename
		char imageName[100];
		strcpy(imageName, PATH);
		strcat(imageName, imageNames[i]);

		// Load the input image
		inputImage.assign(imageName);

		// check image
		if ( inputImage.spectrum() != 3 ) {
			cerr << "The input must be a color image." << endl;
			return 1;
		}

		// output image size
		cout << fixed << setprecision(2);
		cout << "image size " << ((double)inputImage.width()*inputImage.height()) / 1000000 << " Mpixels, ";

		//-------------------------------------------------------------------------------------------------------------------
		// Execution only on CPU
		//-------------------------------------------------------------------------------------------------------------------
		// display text
		cout << "cpu times are ";
		cout << fixed << setprecision(6);

		// average measurements
		for (j = 0; j < AVERAGE_FILTER_SIZE; j++) {
			// create gray scale image space
			grayImage.assign(inputImage.width(), inputImage.height(), 1, 1);

			// create space for histogram
			histogram = new unsigned int [HISTOGRAM_SIZE];

			// create space for smooth image
			smoothImage.assign(grayImage.width(), grayImage.height(), 1, 1);

			// initialize timing
			totaltime[j][0] = totaltime[j][1] = totaltime[j][2] = totaltime[j][3] = totaltime[j][4] = totaltime[j][5] = 0.0;

			// convert to gray scale
			rgb2gray(inputImage.data(), grayImage.data(), inputImage.width(), inputImage.height(), &totaltime[j][1]);

			// Compute histogram
			histogram1D(grayImage.data(), grayImage.width(), grayImage.height(), histogram, HISTOGRAM_SIZE, &totaltime[j][2]);

			// Compute Contrast
			contrast1D(grayImage.data(), grayImage.width(), grayImage.height(), histogram, HISTOGRAM_SIZE, CONTRAST_THRESHOLD, &totaltime[j][3]);

			// Compute smoothing
			triangularSmooth(grayImage.data(), smoothImage.data(), grayImage.width(), grayImage.height(), &totaltime[j][4]);

			// clean up
			delete [] histogram;

			// Give compare image a unique name for comparison
			char str[100];
			strcpy(str, PATHOUT);
			strcat(str, "cpu");
			strcat(str, imageNames[i]);
#ifdef IMAGE_SAVING
			smoothImage.save(str);
#endif

			grayImage.clear();
			smoothImage.clear();
		}
		for (j = 1; j < AVERAGE_FILTER_SIZE; j++) {
			totaltime[0][0] += totaltime[j][0];
			totaltime[0][1] += totaltime[j][1];
			totaltime[0][2] += totaltime[j][2];
			totaltime[0][3] += totaltime[j][3];
			totaltime[0][4] += totaltime[j][4];
			totaltime[0][5] += totaltime[j][5];
		}

		// display timings
		cout << totaltime[0][0] << ", ";
		cout << totaltime[0][1] << ", ";
		cout << totaltime[0][2] << ", ";
		cout << totaltime[0][3] << ", ";
		cout << totaltime[0][4] << ", ";
		cout << totaltime[0][5] << ", ";
		cout << "tot " << totaltime[0][0] + totaltime[0][1] + totaltime[0][2] + totaltime[0][3] + totaltime[0][4] + totaltime[0][5] << ", ";
		speedup = totaltime[0][0] + totaltime[0][1] + totaltime[0][2] + totaltime[0][3] + totaltime[0][4] + totaltime[0][5];

		//-------------------------------------------------------------------------------------------------------------------
		// Execution also on GPU
		//-------------------------------------------------------------------------------------------------------------------
		// display text
		cout << "GPU times are ";

		// average measurements
		for (j = 0; j < AVERAGE_FILTER_SIZE; j++) {
			// create space for smooth image
			smoothImage.assign(inputImage.width(), inputImage.height(), 1, 1);

			// initialize timing
			totaltime[j][0] = totaltime[j][1] = totaltime[j][2] = totaltime[j][3] = totaltime[j][4] = totaltime[j][5] = 0.0;

			// Compute all steps on GPU, Calling the CUDA code
			imageProcess(inputImage.data(), smoothImage.data(), inputImage.width(), inputImage.height(),
			             HISTOGRAM_SIZE, CONTRAST_THRESHOLD, totaltime[j]);

			// Give compare image a unique name for comparison
			char str[100];
			strcpy(str, PATHOUT);
			strcat(str, "GPU");
			strcat(str, imageNames[i]);
#ifdef IMAGE_SAVING
			smoothImage.save(str);
#endif
			// clean up
			smoothImage.clear();
		}
		for (j = 1; j < AVERAGE_FILTER_SIZE; j++) {
			totaltime[0][0] += totaltime[j][0];
			totaltime[0][1] += totaltime[j][1];
			totaltime[0][2] += totaltime[j][2];
			totaltime[0][3] += totaltime[j][3];
			totaltime[0][4] += totaltime[j][4];
			totaltime[0][5] += totaltime[j][5];
		}

		// display timings
		cout << totaltime[0][0] << ", ";
		cout << totaltime[0][1] << ", ";
		cout << totaltime[0][2] << ", ";
		cout << totaltime[0][3] << ", ";
		cout << totaltime[0][4] << ", ";
		cout << totaltime[0][5] << ", ";
		cout << "tot " << totaltime[0][0] + totaltime[0][1] + totaltime[0][2] + totaltime[0][3] + totaltime[0][4] + totaltime[0][5] << ", ";
		speedup = speedup / (totaltime[0][0] + totaltime[0][1] + totaltime[0][2] + totaltime[0][3] + totaltime[0][4] + totaltime[0][5]);
		cout << "speedup is " << speedup << endl;

		// clean up
		inputImage.clear();
	}
	return 0;
}