#include "matrix_cl.h"

// defines used
#define device_id				0		// always select device 0
#define average_filter_size		3		// execute a run multiple times to test consistency

/****************************************************
the following function calculate the below equation
   vector_out = vector_in x matrix_in
 ***************************************************/
void matrix_mult_sq(int size, cl_double *vector_in,
               cl_double *matrix_in, cl_double *vector_out){
  int rows, cols;
  int j;
  for(cols=0; cols<size; cols++){
    vector_out[cols] = 0.0;
    for(j=0,rows=0; rows<size; j++,rows++)
      vector_out[cols] += vector_in[j] * matrix_in[rows*size + cols];
  }
}

/*****************************************************
the following function generates a "size"-element vector
and a "size x size" matrix
 ****************************************************/
void matrix_vector_gen(cl_int size, cl_double *matrix, cl_double *vector){
  int i;
  for(i=0; i<size; i++)
    vector[i] = ((double)rand())/65535.0;
  for(i=0; i<size*size; i++)
    matrix[i] = ((double)rand())/5307.0;
}

int main(int argc, char *argv[]){

	// check argument count
    if(argc < 4){
        printf("Usage: %s (matrix/vector_size) (local group size start) (local group size end)\n", argv[0]);
        return 0;
    }

	// load arguments into values
    cl_int size = atoi(argv[1]);
	cl_int localSizeStart = atoi(argv[2]);
    cl_int localSizeEnd = atoi(argv[3]);
	cl_int localSize;

	// check values
    if((size == 0) || (localSizeStart == 0) || (localSizeEnd == 0)){
        printf("incorrect arguments, make sure the arguments are integers greater than zero\n");
        exit(-1);
	}
	
	// variables for time measurement
	double time_sq, average_time_sq;
	double time_opencl, average_time_opencl;
	int i, j;
	int cols;
	double average_filter[2][average_filter_size];

	// parameters used by opencl
	cl_event mulDone;
	cl_int status;
	cl_uint numPlatforms;
	cl_platform_id *platforms;
	cl_uint numDevices;
	cl_device_id *devices;
	cl_context context;
	cl_command_queue cmdQueue;
	cl_mem bufferVectorIn, bufferMatrixIn, bufferVectorOut;
	char *mulBuffer;
	FILE *mulFile;
	size_t mulSize;
	char *mulFileName = "vectorMatrixMul.cl";
	cl_program program;
	const char options[] = "-cl-std=CL1.2";
	cl_kernel mulKernel;
	size_t localWorkSize[1];
	size_t globalWorkSize[1];  
	
	// allocate space for input matrices and output matrices
	cl_double *vector = (double *)malloc(sizeof(cl_double)*size);
	cl_double *matrix = (double *)malloc(sizeof(cl_double)*size*size);
	cl_double *result_sq = (double *)malloc(sizeof(cl_double)*size);
	cl_double *result_pl = (double *)malloc(sizeof(cl_double)*size);
	
	// acquire opencl program pointer
	mulFile = fopen(mulFileName, "r");
	if(mulFile == NULL){
		printf("cannot open .cl file\n");
		printf("current path: %s\n", mulFileName);
		exit(-1);
	}
	fseek(mulFile, 0, SEEK_END);
	mulSize = ftell(mulFile);
	rewind(mulFile);

	// read kernel source into buffer
	mulBuffer = (char*) malloc(mulSize + 1);
	mulBuffer[mulSize] = '\0';
	fread(mulBuffer, sizeof(char), mulSize, mulFile);
	fclose(mulFile);
	
	//-----------------------------------------------------
	// STEP 1: Discover and initialize the platforms
	//-----------------------------------------------------
	// Use clGetPlatformIDs() to retrieve the number of platforms
	numPlatforms = 0;
	if(clGetPlatformIDs(0, NULL, &numPlatforms) != CL_SUCCESS){
		printf("error getting number of platform IDs\n");
		exit(-1);
	}
	
	// Allocate enough space for each platform
	platforms = (cl_platform_id*)malloc(numPlatforms*sizeof(cl_platform_id));
 
	// Fill in platforms with clGetPlatformIDs()
	if(clGetPlatformIDs(numPlatforms, platforms, NULL) != CL_SUCCESS){
		printf("error getting platform IDs\n");
		exit(-1);
	}

	//-----------------------------------------------------
	// STEP 2: Discover and initialize the devices
	//-----------------------------------------------------
	// Use clGetDeviceIDs() to retrieve the number of devices present
	numDevices = 0;
	if(clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices) != CL_SUCCESS){
		printf("error getting number of device IDs\n");
		exit(-1);
	}
	
	// Allocate enough space for each device
	devices = (cl_device_id*)malloc(numDevices*sizeof(cl_device_id));

	// Fill in devices with clGetDeviceIDs()
	if(clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_ALL, numDevices, devices, NULL) != CL_SUCCESS){
		printf("error getting device IDs\n");
		exit(-1);
	}
	
	// check if device is present
	if(device_id >= numDevices){
		printf("error device not found\n");
		exit(-1);
	}

	//-----------------------------------------------------
	// STEP 3: Create a context
	//-----------------------------------------------------
	// Create a context using clCreateContext() and associate it with the devices
	context = NULL;
	context = clCreateContext(NULL, 1, &devices[device_id], NULL, NULL, &status);
	if(status != CL_SUCCESS){
		printf("error creating context\n");
		exit(-1);
	}

	//-----------------------------------------------------
	// STEP 4: Create a command queue
	//-----------------------------------------------------
	// Create a command queue using clCreateCommandQueue(),
	// and associate it with the device you want to execute 
	// on
	cmdQueue = clCreateCommandQueue(context, devices[device_id], 0, &status);
	if(status != CL_SUCCESS){
		printf("error creating command queue\n");
		exit(-1);
	}
	
	//-----------------------------------------------------
	// STEP 5: Create and compile the program
	//----------------------------------------------------- 
	program = clCreateProgramWithSource(context, 1, (const char**) &mulBuffer, &mulSize, &status);
	if(status != CL_SUCCESS){
		printf("error creating program\n");
		exit(-1);
	}

	// Build (compile) the program for the devices with clBuildProgram()
	if(clBuildProgram(program, 1, &devices[device_id], options, NULL, NULL) != CL_SUCCESS){
		printf("error building program\n");
		exit(-1);
	}

	//-----------------------------------------------------
	// STEP 6: Create the kernel
	//----------------------------------------------------- 
	// Use clCreateKernel() to create a kernel from the program
	mulKernel = NULL;
	mulKernel = clCreateKernel(program, "mul_kernel", &status);
	if(status != CL_SUCCESS){
		printf("error in step 7\n");
		exit(-1);
	}
		
	// execute for every matrix size between start and end value being a multiple of localsize
	for(localSize = localSizeStart; localSize <= localSizeEnd; localSize++){
		// check if size is a multiple of local size
		if(size % localSize == 0){
			// do every variation 5 times for consistency
			for(i=0; i<average_filter_size; i++){
			
				// create random vector and matrix
				matrix_vector_gen(size, matrix, vector);

				// execute sequential part
				time_sq = omp_get_wtime();
				matrix_mult_sq(size, vector, matrix, result_sq);
				time_sq = omp_get_wtime() - time_sq;
				
				// add to average filter
				average_filter[0][i] = time_sq;

				// start measuring parallel part executed with opencl
				time_opencl = omp_get_wtime();

				//-----------------------------------------------------
				// STEP 7: Create device buffers, images and copy data to buffers
				//-----------------------------------------------------
				// create and fill buffers
				bufferVectorIn = clCreateBuffer(context, CL_MEM_READ_ONLY, size*sizeof(cl_double), NULL, &status);
				if(status != CL_SUCCESS){
					printf("error creating buffer for vectorIn\n");
					exit(-1);
				}
				
				bufferMatrixIn = clCreateBuffer(context, CL_MEM_READ_ONLY, size*size*sizeof(cl_double), NULL, &status);
				if(status != CL_SUCCESS){
					printf("error creating buffer for matrixIn\n");
					exit(-1);
				}

				bufferVectorOut = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size*sizeof(cl_double), NULL, &status);
				if(status != CL_SUCCESS){
					printf("error creating buffer for vector out\n");
					exit(-1);
				}

				if(clEnqueueWriteBuffer (cmdQueue, bufferVectorIn, CL_FALSE, 0, size*sizeof(cl_double), vector, 0, NULL, NULL) != CL_SUCCESS){
					printf("error enqueue-ing write buffer vectorIn\n");
					exit(-1);
				}

				if(clEnqueueWriteBuffer (cmdQueue, bufferMatrixIn, CL_FALSE, 0,	size*size*sizeof(cl_double), matrix, 0,	NULL, NULL) != CL_SUCCESS){
					printf("error enqueue-ing write buffer matrixIn\n");
					exit(-1);
				}

				//-----------------------------------------------------
				// STEP 8: Set the kernel arguments
				//----------------------------------------------------- 
				// Associate the input and output buffers with the kernel using clSetKernelArg()
				if(clSetKernelArg(mulKernel, 0, sizeof(cl_mem), &bufferVectorIn) != CL_SUCCESS){
					printf("error settings 1nd kernel argument\n");
					exit(-1);
				}
					
				if(clSetKernelArg(mulKernel, 1, sizeof(cl_mem), &bufferMatrixIn) != CL_SUCCESS){
					printf("error settings 2nd kernel argument\n");
					exit(-1);
				}
					
				if(clSetKernelArg(mulKernel, 2, sizeof(cl_mem), &bufferVectorOut) != CL_SUCCESS){
					printf("error settings 3nd kernel argument\n");
					exit(-1);
				}
					
				if(clSetKernelArg(mulKernel, 3, sizeof(cl_int), &size) != CL_SUCCESS){
					printf("error settings 4nd kernel argument\n");
					exit(-1);
				}

				//-----------------------------------------------------
				// STEP 9: Configure the work-item structure
				//----------------------------------------------------- 
				// Define an index space (global work size) of work items for execution. 
				// A workgroup size (local work size) is not required, but can be used.
				globalWorkSize[0] = size;		// total number of work items
				localWorkSize[0] = localSize;	// number of work items done multithreaded per computation unit 

				if(clEnqueueNDRangeKernel(cmdQueue, mulKernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &mulDone) != CL_SUCCESS){
					printf("error enqueue-ing kernel\n");
					exit(-1);
				}
				
				if(clEnqueueReadBuffer(cmdQueue, bufferVectorOut, CL_TRUE, 0, size*sizeof(cl_double), result_pl, 1, &mulDone, NULL) != CL_SUCCESS){
					printf("error enqueue-ing read buffer\n");
					exit(-1);
				}

				// end measuring time parallel execution with opencl
				time_opencl = omp_get_wtime() - time_opencl;

				// add to average filter
				average_filter[1][i] = time_opencl;
				
				//check vector validity
				for(cols=0; cols<size; cols++){
					if((int) result_sq[cols] != (int) result_pl[cols]){
						printf("wrong at col %d\n", cols);
						return 0;
					}
				}
				
				//-----------------------------------------------------
				// STEP 10: Release OpenCL resources
				//----------------------------------------------------- 
				// Free OpenCL resources
				clReleaseMemObject(bufferVectorIn);
				clReleaseMemObject(bufferMatrixIn);
				clReleaseMemObject(bufferVectorOut);
			} // end for(i=0; i<average_filter_size; i++)
			
			// calculate averaged values
			average_time_sq = 0;
			average_time_opencl = 0;
			for(j=0; j<average_filter_size; j++){
				average_time_sq += average_filter[0][j];
				average_time_opencl += average_filter[1][j];
			}
			average_time_sq = average_time_sq / average_filter_size;
			average_time_opencl = average_time_opencl / average_filter_size;
			
			// show averaged value
			printf("matrix size %5d, workgroup size %3d, seq exec %1.6f, parallel exec %1.6f, seq/parallel %1.6f\n", size, localSize, average_time_sq, average_time_opencl, (average_time_sq/average_time_opencl));
			
		} // end if(size % localSize == 0)
	} // end for(localSize = localSizeStart; localSize <= localSizeEnd; localSize++)
	
	// free opencl resources
	clReleaseContext(context);
	clReleaseCommandQueue(cmdQueue);
	clReleaseProgram(program);
	clReleaseKernel(mulKernel);
	
	//Free up memory
	free(platforms);
	free(devices);
	free(mulBuffer);
	free(vector);
	free(matrix);
	free(result_sq);
	free(result_pl);

    return EXIT_SUCCESS;
}

