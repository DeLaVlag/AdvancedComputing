//===================================================================================================================
// File:        main.cpp
// Created on:  27-10-11
// Authors:     Dirk Vos, Mark Schrauwen, Michiel van der Vlag
//
// Description: The filter file called from main.cpp. This file contains the actual CUDA code. Different functions
//				of the code can be enabled\disabled by (un)commenting the macros as defined below.
//===================================================================================================================

#include <Timer.hpp>
#include <iostream>
#include <iomanip>

//-------------------------------------------------------------------------------------------------------------------
// Global definitions, macros
//-------------------------------------------------------------------------------------------------------------------
// Use commenting to disable/enable functions
#define MAX_BLOCKSIZE	32			// The size of one square block. A block will have size MAX_BLOCKSIZE * MAX_BLOCKSIZE
#define DATA_SIZE		7			// The number of variables of Memory Mapping part.
#define SHARED_MEM					// If uncommented, shared memory will be used.
#define TEXTURE_MEM					// If uncommented, Texture memory will be used.

using LOFAR::NSTimer;
using std::cout;
using std::cerr;
using std::endl;
using std::fixed;
using std::setprecision;

// global information
static unsigned int numMultiProcessors;
static unsigned int numThreadsPerBlock[2];

// define 2D textures
texture<unsigned char, cudaTextureType2D, cudaReadModeElementType> texture2DRed;
texture<unsigned char, cudaTextureType2D, cudaReadModeElementType> texture2DGreen;
texture<unsigned char, cudaTextureType2D, cudaReadModeElementType> texture2DBlue;


//-------------------------------------------------------------------------------------------------------------------
// Memory Mapping
//-------------------------------------------------------------------------------------------------------------------
/**
 *   unsigned int *dev_data is filled in the following way:
 *	  index		content
 *		0      	width
 *		1		height
 *		2		HISTOGRAM_SIZE
 *		3		CONTRAST_THRESHOLD
 *		4		pitch red image
 *		5		pitch green image
 *		6		pitch blue image
 *		7		min			(filled on device)
 *		8		max			(filled on device)
 *		9		histogram	(filled on device)
 */


//-------------------------------------------------------------------------------------------------------------------
// Device discovery
//-------------------------------------------------------------------------------------------------------------------
__host__ void getCudaDeviceInformation(void) {
	// get GPU device
	int device;
	if (cudaGetDevice(&device) != cudaSuccess) {
		cout << "main - cuda get device failed" << endl;
		exit(1);
	}

	// get GPU properties
	cudaDeviceProp prop;
	if (cudaGetDeviceProperties (&prop, device) != cudaSuccess) {
		cout << "main - cuda get device properties failed" << endl;
		exit(1);
	}

	// save to global information
	numMultiProcessors    = prop.multiProcessorCount;
	numThreadsPerBlock[0] = MAX_BLOCKSIZE;
	numThreadsPerBlock[1] = MAX_BLOCKSIZE;
}


//-------------------------------------------------------------------------------------------------------------------
// helper functions
//-------------------------------------------------------------------------------------------------------------------
__device__ unsigned char getElementRed(unsigned int x, unsigned int y) {
	return tex2D(texture2DRed, x, y);
}

__device__ unsigned char getElementGreen(unsigned int x, unsigned int y) {
	return tex2D(texture2DGreen, x, y);
}

__device__ unsigned char getElementBlue(unsigned int x, unsigned int y) {
	return tex2D(texture2DBlue, x, y);
}

__device__ unsigned int getWidth(unsigned int *data) {
	return data[0];
}

__device__ unsigned int getHeight(unsigned int *data) {
	return data[1];
}

__device__ unsigned int getHistogramSize(unsigned int *data) {
	return data[2];
}

__device__ unsigned int getContrastThreshold(unsigned int *data) {
	return data[3];
}

__device__ unsigned int getPitchRedImage(unsigned int *data) {
	return data[4];
}

__device__ unsigned int getPitchGreenImage(unsigned int *data) {
	return data[5];
}

__device__ unsigned int getPitchBlueImage(unsigned int *data) {
	return data[6];
}

__device__ unsigned int getMin(unsigned int *data) {
	return data[7];
}

__device__ unsigned int getMax(unsigned int *data) {
	return data[8];
}

__device__ void setMin(unsigned int *data, unsigned int value) {
	data[7] = value;
}

__device__ void setMax(unsigned int *data, unsigned int value) {
	data[8] = value;
}

__device__ unsigned int * getHistogram(unsigned int *data) {
	return &data[9];
}

//-------------------------------------------------------------------------------------------------------------------
// rgb2gray
//-------------------------------------------------------------------------------------------------------------------
#ifdef TEXTURE_MEM
__global__ void rgb2grayCudaKernel(unsigned char *image, unsigned int *data) {
	// load width and height
	unsigned int width  = getWidth(data);
	unsigned int height = getHeight(data);

	// load block ID's. Thread ids
	int y = threadIdx.y + blockIdx.y * blockDim.y;
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned int pitch = getPitchRedImage(data);

#ifdef SHARED_MEM
	// only run threads that are in the image
	if (x < width && y < height) {
		// allocate shared mem for this block
		__shared__ unsigned char redImage[MAX_BLOCKSIZE][MAX_BLOCKSIZE];
		__shared__ unsigned char greenImage[MAX_BLOCKSIZE][MAX_BLOCKSIZE];
		__shared__ unsigned char blueImage[MAX_BLOCKSIZE][MAX_BLOCKSIZE];

		// copy from global mem to shared mem
		redImage[threadIdx.x][threadIdx.y]   = getElementRed(x, y);
		greenImage[threadIdx.x][threadIdx.y] = getElementGreen(x, y);
		blueImage[threadIdx.x][threadIdx.y]  = getElementBlue(x, y);

		__syncthreads();

		// execute grey scaling code
		float grayPix = 0.0f;
		float r       = static_cast< float >(redImage[threadIdx.x][threadIdx.y]);
		float g       = static_cast< float >(greenImage[threadIdx.x][threadIdx.y]);
		float b       = static_cast< float >(blueImage[threadIdx.x][threadIdx.y]);
		grayPix       = (0.3f * r) + (0.59f * g) + (0.11f * b);

		__syncthreads();

		// write back to global mem
		image[(y * pitch) + x] = static_cast< unsigned char >(grayPix);
	}
#else
	// only run threads that are in the image
	if (x < width && y < height) {
		// execute grey scaling code without shared mem
		float grayPix = 0.0f;
		float r       = static_cast< float >(getElementRed(x, y));
		float g       = static_cast< float >(getElementGreen(x, y));
		float b       = static_cast< float >(getElementBlue(x, y));
		grayPix       = (0.3f * r) + (0.59f * g) + (0.11f * b);

		// write back to global mem
		image[(y * pitch) + x] = static_cast< unsigned char >(grayPix);
	}
#endif
}
#else
__global__ void rgb2grayCudaKernel(unsigned char *red, unsigned char *green, unsigned char *blue, unsigned int *data) {
	// load width and height
	unsigned int width  = getWidth(data);
	unsigned int height = getHeight(data);

	// load block ID's. Thread ids
	int y = threadIdx.y + blockIdx.y * blockDim.y;
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned int redpitch   = getPitchRedImage(data);
	unsigned int greenpitch = getPitchGreenImage(data);
	unsigned int bluepitch  = getPitchBlueImage(data);

#ifdef SHARED_MEM
	// only run threads that are in the image
	if (x < width && y < height) {
		// allocate shared mem for this block
		__shared__ unsigned char redImage[MAX_BLOCKSIZE][MAX_BLOCKSIZE];
		__shared__ unsigned char greenImage[MAX_BLOCKSIZE][MAX_BLOCKSIZE];
		__shared__ unsigned char blueImage[MAX_BLOCKSIZE][MAX_BLOCKSIZE];

		// copy from global mem to shared mem
		redImage[threadIdx.x][threadIdx.y]   = red[(y * redpitch) + x];
		greenImage[threadIdx.x][threadIdx.y] = green[(y * greenpitch) + x];
		blueImage[threadIdx.x][threadIdx.y]  = blue[(y * bluepitch) + x];

		__syncthreads();

		// execute grey scaling code
		float grayPix = 0.0f;
		float r = static_cast< float >(redImage[threadIdx.x][threadIdx.y]);
		float g = static_cast< float >(greenImage[threadIdx.x][threadIdx.y]);
		float b = static_cast< float >(blueImage[threadIdx.x][threadIdx.y]);
		grayPix = (0.3f * r) + (0.59f * g) + (0.11f * b);

		__syncthreads();

		// write back to global mem
		red[(y * redpitch) + x] = static_cast< unsigned char >(grayPix);
	}
#else
	// only run threads that are in the image
	if (x < width && y < height) {
		// execute grey scaling code without shared mem
		float grayPix = 0.0f;
		float r = static_cast< float >(red[(y * redpitch) + x]);
		float g = static_cast< float >(green[(y * greenpitch) + x]);
		float b = static_cast< float >(blue[(y * bluepitch) + x]);
		grayPix = (0.3f * r) + (0.59f * g) + (0.11f * b);

		// write back to global mem
		red[(y * redpitch) + x] = static_cast< unsigned char >(grayPix);
	}
#endif
}
#endif

__host__ void rgb2gray(unsigned char *inputImage, unsigned char *grayImage, const int width, const int height, double *totaltime) {
	NSTimer kernelTime = NSTimer("kernelTime", false, false);

	kernelTime.start();
	for ( int y = 0; y < height; y++ )
	{
		for ( int x = 0; x < width; x++ )
		{
			float grayPix = 0.0f;
			float r = static_cast< float >(inputImage[(y * width) + x]);
			float g = static_cast< float >(inputImage[(width * height) + (y * width) + x]);
			float b = static_cast< float >(inputImage[(2 * width * height) + (y * width) + x]);

			grayPix = (0.3f * r) + (0.59f * g) + (0.11f * b);

			grayImage[(y * width) + x] = static_cast< unsigned char >(grayPix);
		}
	}
	kernelTime.stop();
	*totaltime = kernelTime.getElapsed();
}

//-------------------------------------------------------------------------------------------------------------------
// histogram1D
//-------------------------------------------------------------------------------------------------------------------
#ifdef TEXTURE_MEM
__global__ void histogram1DCudaKernel(unsigned int *data) {
	// load width and height
	unsigned int width  = getWidth(data);
	unsigned int height = getHeight(data);

	// load block ID's. Thread ids
	int y = threadIdx.y + blockIdx.y * blockDim.y;
	int x = threadIdx.x + blockIdx.x * blockDim.x;

	// get histogram pointer
	unsigned int *histogram = getHistogram(data);
	unsigned int index;

#ifdef SHARED_MEM
	// only run threads that are in the image
	if (x < width && y < height) {
		// allocate shared mem for this block
		__shared__ unsigned char grayImage[MAX_BLOCKSIZE][MAX_BLOCKSIZE];

		// copy from global mem to shared mem
		grayImage[threadIdx.x][threadIdx.y] = getElementRed(x, y);

		__syncthreads();

		// get histogram index
		index = static_cast< unsigned int >(grayImage[threadIdx.x][threadIdx.y]);

		__syncthreads();

		// add pixel to histogram in one threadsafe operation
		atomicAdd(&histogram[index], 1);
	}
#else
	// only run threads that are in the image
	if (x < width && y < height) {
		// get histogram index
		index = static_cast< unsigned int >(getElementRed(x, y));

		// add pixel to histogram in one threadsafe operation
		atomicAdd(&histogram[index], 1);
	}
#endif
}
#else
__global__ void histogram1DCudaKernel(unsigned char* image, unsigned int *data) {
	// load width and height
	unsigned int width  = getWidth(data);
	unsigned int height = getHeight(data);

	// load block ID's. Thread ids
	int y = threadIdx.y + blockIdx.y * blockDim.y;
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned int pitch = getPitchRedImage(data);

	// get histogram pointer
	unsigned int *histogram = getHistogram(data);
	unsigned int index;

#ifdef SHARED_MEM
	// only run threads that are in the image
	if (x < width && y < height) {
		// allocate shared mem for this block
		__shared__ unsigned char grayImage[MAX_BLOCKSIZE][MAX_BLOCKSIZE];

		// copy from global mem to shared mem
		grayImage[threadIdx.x][threadIdx.y] = image[(y * pitch) + x];

		__syncthreads();

		// get histogram index
		index = static_cast< unsigned int >(grayImage[threadIdx.x][threadIdx.y]);

		__syncthreads();

		// add pixel to histogram in one threadsafe operation
		atomicAdd(&histogram[index], 1);
	}
#else
	// only run threads that are in the image
	if (x < width && y < height) {
		// get histogram index
		index = static_cast< unsigned int >(image[(y * pitch) + x]);

		// add pixel to histogram in one threadsafe operation
		atomicAdd(&histogram[index], 1);
	}
#endif
}
#endif

__host__ void histogram1D(unsigned char *grayImage, const int width, const int height,
                          unsigned int *histogram, const unsigned int HISTOGRAM_SIZE, double *totaltime)
{
	NSTimer kernelTime = NSTimer("kernelTime", false, false);

	memset(reinterpret_cast< void * >(histogram), 0, HISTOGRAM_SIZE * sizeof(int));

	kernelTime.start();
	// Kernel
	for ( int y = 0; y < height; y++ )
	{
		for ( int x = 0; x < width; x++ )
		{
			histogram[static_cast< unsigned int >(grayImage[(y * width) + x])] += 1;
		}
	}
	// /Kernel
	kernelTime.stop();
	*totaltime = kernelTime.getElapsed();
}

//-------------------------------------------------------------------------------------------------------------------
// contrast1D
//-------------------------------------------------------------------------------------------------------------------
__global__ void contrastMinKernel(unsigned int *data) {
	// load width and height, histogramsize and contrast threshold
	unsigned int histogramSize     = getHistogramSize(data);
	unsigned int contrastThreshold = getContrastThreshold(data);

	// get histogram pointer
	unsigned int *histogram = getHistogram(data);

	// find minimum
	unsigned int i = 0;
	while ( (i < histogramSize) && (histogram[i] < contrastThreshold) ) {
		i++;
	}
	setMin(data, i);
}

__global__ void contrastMaxKernel(unsigned int *data) {
	// load width and height, histogramsize and contrast threshold
	unsigned int histogramSize     = getHistogramSize(data);
	unsigned int contrastThreshold = getContrastThreshold(data);
	unsigned int min               = getMin(data);

	// get histogram pointer
	unsigned int *histogram = getHistogram(data);

	// find maximum
	unsigned int i = histogramSize - 1;
	while ( (i > min) && (histogram[i] < contrastThreshold) ) {
		i--;
	}
	setMax(data, i);
}

#ifdef TEXTURE_MEM
__global__ void contrast1DKernel(unsigned char *image, unsigned int *data) {
	// load width and height
	unsigned int width  = getWidth(data);
	unsigned int height = getHeight(data);

	// load block ID's. Thread ids
	int y = threadIdx.y + blockIdx.y * blockDim.y;
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned int pitch = getPitchRedImage(data);

	// calculate difference
	float diff = getMax(data) - getMin(data);

	// get pixel
	unsigned int min = getMin(data);
	unsigned int max = getMax(data);

#ifdef SHARED_MEM
	// only run threads that are in the image
	if (x < width && y < height) {
		// allocate shared mem for this block
		__shared__ unsigned char grayImage[MAX_BLOCKSIZE][MAX_BLOCKSIZE];

		// copy from global mem to shared mem
		grayImage[threadIdx.x][threadIdx.y] = getElementRed(x, y);

		__syncthreads();

		// apply contrast enhancement
		unsigned char pixel = grayImage[threadIdx.x][threadIdx.y];
		if ( pixel < min )	{
			pixel = 0;
		}
		else if ( pixel > max )	{
			pixel = 255;
		}
		else	{
			pixel = static_cast< unsigned char >(255.0f * (pixel - min) / diff);
		}

		__syncthreads();

		// write back pixel
		image[(y * pitch) + x] = pixel;
	}
#else
	// only run threads that are in the image
	if (x < width && y < height) {
		// apply contrast enhancement
		unsigned char pixel = getElementRed(x, y);
		if ( pixel < min )	{
			pixel = 0;
		}
		else if ( pixel > max )	{
			pixel = 255;
		}
		else	{
			pixel = static_cast< unsigned char >(255.0f * (pixel - min) / diff);
		}

		// write back pixel
		image[(y * pitch) + x] = pixel;
	}
#endif
}
#else
__global__ void contrast1DKernel(unsigned char *image, unsigned int *data) {
	// load width and height
	unsigned int width  = getWidth(data);
	unsigned int height = getHeight(data);

	// load block ID's. Thread ids
	int y = threadIdx.y + blockIdx.y * blockDim.y;
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned int pitch = getPitchRedImage(data);

	// calculate difference
	float diff = getMax(data) - getMin(data);

	// get pixel
	unsigned int min = getMin(data);
	unsigned int max = getMax(data);

#ifdef SHARED_MEM
	// only run threads that are in the image
	if (x < width && y < height) {
		// allocate shared mem for this block
		__shared__ unsigned char grayImage[MAX_BLOCKSIZE][MAX_BLOCKSIZE];

		// copy from global mem to shared mem
		grayImage[threadIdx.x][threadIdx.y] = image[(y * pitch) + x];

		__syncthreads();

		// apply contrast enhancement
		unsigned char pixel = grayImage[threadIdx.x][threadIdx.y];
		if ( pixel < min )	{
			pixel = 0;
		}
		else if ( pixel > max )	{
			pixel = 255;
		}
		else	{
			pixel = static_cast< unsigned char >(255.0f * (pixel - min) / diff);
		}

		__syncthreads();

		// write back pixel
		image[(y * pitch) + x] = pixel;
	}
#else
	// only run threads that are in the image
	if (x < width && y < height) {
		// apply contrast enhancement
		unsigned char pixel = image[(y * pitch) + x];
		if ( pixel < min )	{
			pixel = 0;
		}
		else if ( pixel > max )	{
			pixel = 255;
		}
		else	{
			pixel = static_cast< unsigned char >(255.0f * (pixel - min) / diff);
		}

		// write back pixel
		image[(y * pitch) + x] = pixel;
	}
#endif
}
#endif

__host__ void contrast1D(unsigned char *grayImage, const int width, const int height,
                         unsigned int *histogram, const unsigned int HISTOGRAM_SIZE,
                         const unsigned int CONTRAST_THRESHOLD, double *totaltime)
{
	unsigned int i = 0;
	NSTimer kernelTime = NSTimer("kernelTime", false, false);

	while ( (i < HISTOGRAM_SIZE) && (histogram[i] < CONTRAST_THRESHOLD) )
	{
		i++;
	}
	unsigned int min = i;

	i = HISTOGRAM_SIZE - 1;
	while ( (i > min) && (histogram[i] < CONTRAST_THRESHOLD) )
	{
		i--;
	}
	unsigned int max = i;
	float diff = max - min;

	kernelTime.start();
	// Kernel
	for ( int y = 0; y < height; y++ )
	{
		for (int x = 0; x < width; x++ )
		{
			unsigned char pixel = grayImage[(y * width) + x];

			if ( pixel < min )
			{
				pixel = 0;
			}
			else if ( pixel > max )
			{
				pixel = 255;
			}
			else
			{
				pixel = static_cast< unsigned char >(255.0f * (pixel - min) / diff);
			}

			grayImage[(y * width) + x] = pixel;
		}
	}
	// /Kernel
	kernelTime.stop();
	*totaltime = kernelTime.getElapsed();
}

//-------------------------------------------------------------------------------------------------------------------
// triangularSmooth
//-------------------------------------------------------------------------------------------------------------------
#ifdef TEXTURE_MEM
__global__ void triangularSmoothKernel(unsigned char *image, unsigned int *data) {
	// load width and height
	unsigned int width  = getWidth(data);
	unsigned int height = getHeight(data);

	// load block ID's. Thread ids
	int y = threadIdx.y + blockIdx.y * blockDim.y;
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned int pitch = getPitchGreenImage(data);

	// only run threads that are in the image
	if (x < width && y < height) {
		// declare variables
		unsigned int filterItem = 0;
		float filterSum = 0.0f;
		float smoothPix = 0.0f;
		unsigned char value;
		int fy, fx;
		const float filter[] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
		                        1.0f, 2.0f, 2.0f, 2.0f, 1.0f,
		                        1.0f, 2.0f, 3.0f, 2.0f, 1.0f,
		                        1.0f, 2.0f, 2.0f, 2.0f, 1.0f,
		                        1.0f, 1.0f, 1.0f, 1.0f, 1.0f
		                       };

		// do the smoothing
		for ( fy = y - 2; fy < y + 3; fy++ )
		{
			for ( fx = x - 2; fx < x + 3; fx++ )
			{
				if ( ((fy < 0) || (fy >= height)) || ((fx < 0) || (fx >= width)) )
				{
					filterItem++;
					continue;
				}

				smoothPix += getElementRed(fx, fy) * filter[filterItem];
				filterSum += filter[filterItem];
				filterItem++;
			}
		}
		smoothPix /= filterSum;
		value = static_cast< unsigned char >(smoothPix);

		// write back to global mem
		image[(y * pitch) + x] = value;
	}
}
#else
__global__ void triangularSmoothKernel(unsigned char *red, unsigned char *green, unsigned int *data) {
	// load width and height
	unsigned int width  = getWidth(data);
	unsigned int height = getHeight(data);

	// load block ID's. Thread ids
	int y = threadIdx.y + blockIdx.y * blockDim.y;
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned int redpitch   = getPitchRedImage(data);
	unsigned int greenpitch = getPitchGreenImage(data);

	// only run threads that are in the image
	if (x < width && y < height) {
		// declare variables
		unsigned int filterItem = 0;
		float filterSum = 0.0f;
		float smoothPix = 0.0f;
		unsigned char value;
		int fy, fx;
		const float filter[] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
		                        1.0f, 2.0f, 2.0f, 2.0f, 1.0f,
		                        1.0f, 2.0f, 3.0f, 2.0f, 1.0f,
		                        1.0f, 2.0f, 2.0f, 2.0f, 1.0f,
		                        1.0f, 1.0f, 1.0f, 1.0f, 1.0f
		                       };

		// do the smoothing
		for ( fy = y - 2; fy < y + 3; fy++ )
		{
			for ( fx = x - 2; fx < x + 3; fx++ )
			{
				if ( ((fy < 0) || (fy >= height)) || ((fx < 0) || (fx >= width)) )
				{
					filterItem++;
					continue;
				}

				smoothPix += red[(fy * redpitch) + fx] * filter[filterItem];
				filterSum += filter[filterItem];
				filterItem++;
			}
		}
		smoothPix /= filterSum;
		value = static_cast< unsigned char >(smoothPix);

		// write back to global mem
		green[(y * greenpitch) + x] = value;
	}
}
#endif

__host__ void triangularSmooth(unsigned char *grayImage, unsigned char *smoothImage, const int width, const int height, double *totaltime)
{
	NSTimer kernelTime = NSTimer("kernelTime", false, false);

	kernelTime.start();
	// Kernel
	for ( int y = 0; y < height; y++ )
	{
		for ( int x = 0; x < width; x++ )
		{
			unsigned int filterItem = 0;
			float filterSum = 0.0f;
			float smoothPix = 0.0f;
			const float filter[] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
			                        1.0f, 2.0f, 2.0f, 2.0f, 1.0f,
			                        1.0f, 2.0f, 3.0f, 2.0f, 1.0f,
			                        1.0f, 2.0f, 2.0f, 2.0f, 1.0f,
			                        1.0f, 1.0f, 1.0f, 1.0f, 1.0f
			                       };

			for ( int fy = y - 2; fy < y + 3; fy++ )
			{
				for ( int fx = x - 2; fx < x + 3; fx++ )
				{
					if ( ((fy < 0) || (fy >= height)) || ((fx < 0) || (fx >= width)) )
					{
						filterItem++;
						continue;
					}

					smoothPix += grayImage[(fy * width) + fx] * filter[filterItem];
					filterSum += filter[filterItem];
					filterItem++;
				}
			}

			smoothPix /= filterSum;
			smoothImage[(y * width) + x] = static_cast< unsigned char >(smoothPix);
		}
	}
	// /Kernel
	kernelTime.stop();
	*totaltime = kernelTime.getElapsed();
}

__host__ void imageProcess(unsigned char *inputImage, unsigned char *smoothImage, const int width, const int height,
                           const unsigned int HISTOGRAM_SIZE, const unsigned int CONTRAST_THRESHOLD, double *totaltime) {

//-------------------------------------------------------------------------------------------------------------------
// calculation occupancy GPU
//-------------------------------------------------------------------------------------------------------------------
	// variables
	unsigned int numBlocks[2];
	unsigned int numThreads[2];
	unsigned int widthCounter, heightCounter;

	// calculate number of blocks in both directions
	widthCounter = 1;
	while (width > numMultiProcessors * widthCounter * numThreadsPerBlock[0]) {
		widthCounter++;
	}
	numBlocks[0] = numMultiProcessors * widthCounter;
	heightCounter = 1;
	while (height > numMultiProcessors * heightCounter * numThreadsPerBlock[1]) {
		heightCounter++;
	}
	numBlocks[1] = numMultiProcessors * heightCounter;
	dim3 blockGrid(numBlocks[0], numBlocks[1]);

	// calculate number of threads per block in both directions
	numThreads[0] = width / numBlocks[0];
	if (width % numBlocks[0] > 0)
		numThreads[0]++;
	numThreads[1] = height / numBlocks[1];
	if (height % numBlocks[1] > 0)
		numThreads[1]++;
	dim3 threadGrid(numThreads[0], numThreads[1]);

//-------------------------------------------------------------------------------------------------------------------
// Write Memory to GPU device
//-------------------------------------------------------------------------------------------------------------------
	// device variables
	unsigned char *dev_redimage, *dev_greenimage, *dev_blueimage;
	unsigned int  *dev_data;
	// host variables
	NSTimer kernelTime = NSTimer("kernelTime", false, false);

	// initialize pitch values (will be filled in by the malloc function later)
	size_t redimage_pitch   = 0;
	size_t greenimage_pitch = 0;
	size_t blueimage_pitch  = 0;

	// allocate 2D memory on GPU
	if (cudaMallocPitch<unsigned char>(&dev_redimage, &redimage_pitch, width, height) != cudaSuccess) {
		cout << "imageProcess - cuda pitch malloc red failed" << endl;
		exit(1);
	}
	if (cudaMallocPitch<unsigned char>(&dev_greenimage, &greenimage_pitch, width, height) != cudaSuccess) {
		cout << "imageProcess - cuda pitch malloc green failed" << endl;
		exit(1);
	}
	if (cudaMallocPitch<unsigned char>(&dev_blueimage, &blueimage_pitch, width, height) != cudaSuccess) {
		cout << "imageProcess - cuda pitch malloc blue failed" << endl;
		exit(1);
	}
	if (cudaMalloc(&dev_data, HISTOGRAM_SIZE * sizeof(int) + 9 * sizeof(int)) != cudaSuccess) {
		cout << "imageProcess - cuda malloc data failed" << endl;
		exit(1);
	}
	unsigned int data[DATA_SIZE] = {width, height, HISTOGRAM_SIZE, CONTRAST_THRESHOLD, redimage_pitch, greenimage_pitch, blueimage_pitch};

	// start timing
	kernelTime.start();

	// copy data from host to device.
	if (cudaMemcpy2D(dev_redimage, redimage_pitch, &inputImage[0], width, width, height, cudaMemcpyHostToDevice) != cudaSuccess) {
		cout << "imageProcess - cuda 2D mem cpy red failed" << endl;
		exit(1);
	}
	if (cudaMemcpy2D(dev_greenimage, greenimage_pitch, &inputImage[width * height], width, width, height, cudaMemcpyHostToDevice) != cudaSuccess) {
		cout << "imageProcess - cuda 2D mem cpy green failed" << endl;
		exit(1);
	}
	if (cudaMemcpy2D(dev_blueimage, blueimage_pitch, &inputImage[2 * width * height], width, width, height, cudaMemcpyHostToDevice) != cudaSuccess) {
		cout << "imageProcess - cuda 2D mem cpy blue failed" << endl;
		exit(1);
	}
	if (cudaMemcpy(dev_data, data, DATA_SIZE * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) {
		cout << "imageProcess - cuda 2D mem cpy data failed" << endl;
		exit(1);
	}
	if (cudaMemset(&dev_data[9], 0, HISTOGRAM_SIZE * sizeof(int)) != cudaSuccess) {
		cout << "imageProcess - cuda mem set histogram failed" << endl;
		exit(1);
	}

	// stop timing
	kernelTime.stop();
	totaltime[0] = kernelTime.getElapsed();

#ifdef TEXTURE_MEM
	//Bind the image to the texture. Now the kernel will read the input image through the texture cache.
	if (cudaBindTexture2D(NULL, texture2DRed, dev_redimage, width, height, redimage_pitch) != cudaSuccess) {
		cout << "imageProcess - cuda bind 2D texture red failed" << endl;
		exit(1);
	}
	if (cudaBindTexture2D(NULL, texture2DGreen, dev_greenimage, width, height, greenimage_pitch) != cudaSuccess) {
		cout << "imageProcess - cuda bind 2D texture green failed" << endl;
		exit(1);
	}
	if (cudaBindTexture2D(NULL, texture2DBlue, dev_blueimage, width, height, blueimage_pitch) != cudaSuccess) {
		cout << "imageProcess - bind 2D texture blue failed" << endl;
		exit(1);
	}

	// set border access to zero
	texture2DRed.addressMode[0]   = texture2DRed.addressMode[1] 	= cudaAddressModeBorder;
	texture2DGreen.addressMode[0] = texture2DGreen.addressMode[1] 	= cudaAddressModeBorder;
	texture2DBlue.addressMode[0]  = texture2DBlue.addressMode[1] 	= cudaAddressModeBorder;
#endif

//-------------------------------------------------------------------------------------------------------------------
// RGB to gray scale conversion
//-------------------------------------------------------------------------------------------------------------------
	// start timing
	kernelTime.reset();
	kernelTime.start();

	// create grid for kernel functions, execute the kernel on the GPU
#ifdef TEXTURE_MEM
	rgb2grayCudaKernel <<< blockGrid, threadGrid>>>(dev_redimage, dev_data);
#else
	rgb2grayCudaKernel <<< blockGrid, threadGrid>>>(dev_redimage, dev_greenimage, dev_blueimage, dev_data);
#endif
	if (cudaGetLastError() != cudaSuccess) {
		cout << "imageProcess - cuda start rgb2gray kernel on device failed" << endl;
		exit(1);
	}

	// stop timing
	kernelTime.stop();
	totaltime[1] = kernelTime.getElapsed();

//-------------------------------------------------------------------------------------------------------------------
// Creating histogram
//-------------------------------------------------------------------------------------------------------------------
	// start timing
	kernelTime.reset();
	kernelTime.start();

	// create grid for kernel functions, execute the kernel on the GPU
#ifdef TEXTURE_MEM
	histogram1DCudaKernel <<< blockGrid, threadGrid>>>(dev_data);
#else
	histogram1DCudaKernel <<< blockGrid, threadGrid>>>(dev_redimage, dev_data);
#endif
	if (cudaGetLastError() != cudaSuccess) {
		cout << "imageProcess - cuda start histogram1D kernel on device failed" << endl;
		exit(1);
	}

	// stop timing
	kernelTime.stop();
	totaltime[2] = kernelTime.getElapsed();

//-------------------------------------------------------------------------------------------------------------------
// Computing contrast
//-------------------------------------------------------------------------------------------------------------------
	// start timing
	kernelTime.reset();
	kernelTime.start();

	// create grid for kernel functions, execute the kernel on the GPU
	contrastMinKernel <<< 1, 1>>>(dev_data);
	if (cudaGetLastError() != cudaSuccess) {
		cout << "imageProcess - cuda start contrastMin kernel on device failed" << endl;
		exit(1);
	}

	// create grid for kernel functions, execute the kernel on the GPU
	contrastMaxKernel <<< 1, 1>>>(dev_data);
	if (cudaGetLastError() != cudaSuccess) {
		cout << "imageProcess - cuda start contrastMax kernel on device failed" << endl;
		exit(1);
	}

	// create grid for kernel functions, execute the kernel on the GPU
#ifdef TEXTURE_MEM
	contrast1DKernel <<< blockGrid, threadGrid>>>(dev_redimage, dev_data);
#else
	contrast1DKernel <<< blockGrid, threadGrid>>>(dev_redimage, dev_data);
#endif
	if (cudaGetLastError() != cudaSuccess) {
		cout << "imageProcess - cuda start contrast1D kernel on device failed" << endl;
		exit(1);
	}

	// stop timing
	kernelTime.stop();
	totaltime[3] = kernelTime.getElapsed();

//-------------------------------------------------------------------------------------------------------------------
// Computing Smoothing
//-------------------------------------------------------------------------------------------------------------------
	// start timing
	kernelTime.reset();
	kernelTime.start();

	// create grid for kernel functions, execute the kernel on the GPU
#ifdef TEXTURE_MEM
	triangularSmoothKernel <<< blockGrid, threadGrid>>>(dev_greenimage, dev_data);
#else
	triangularSmoothKernel <<< blockGrid, threadGrid>>>(dev_redimage, dev_greenimage, dev_data);
#endif
	if (cudaGetLastError() != cudaSuccess) {
		cout << "triangularSmoothCuda - cuda start kernels on device failed" << endl;
		exit(1);
	}

	// stop timing
	kernelTime.stop();
	totaltime[4] = kernelTime.getElapsed();

	// start timing
	kernelTime.reset();
	kernelTime.start();

	// read back result from GPU
	if (cudaMemcpy2D(smoothImage, width, dev_greenimage, greenimage_pitch, width, height, cudaMemcpyDeviceToHost) != cudaSuccess) {
		cout << "triangularSmoothCuda - cuda mem copy smooth image to host failed" << endl;
		exit(1);
	}

	// stop timing
	kernelTime.stop();
	totaltime[5] = kernelTime.getElapsed();

#ifdef TEXTURE_MEM
	//Release the texture
	cudaUnbindTexture(texture2DRed);
	cudaUnbindTexture(texture2DGreen);
	cudaUnbindTexture(texture2DBlue);
#endif

	// free memory on GPU
	cudaFree(dev_redimage);
	cudaFree(dev_greenimage);
	cudaFree(dev_blueimage);
	cudaFree(dev_data);
}