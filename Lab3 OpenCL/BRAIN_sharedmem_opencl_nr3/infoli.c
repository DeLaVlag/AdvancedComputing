#include "infoli.h"
#include "init.h"
#include "ioFile.h"

// defines used
#define device_id				0		// always select device 0

// timing
#define averaging_number		3		// execute a run multiple times to test consistency
#define timing_variables		6		// number of timing variables to print

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp (){
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

int main(int argc, char *argv[]){
	// variables
    cl_uint i, j, k, n, p, q;
	cl_int return_value;
    int simSteps = 0;
    int inputFromFile = 0;
    int initSteps;
    int previApp;
    cl_mod_prec *cellStatePtr;
	unsigned int sim_step;
	cl_bool* state_ptr;
    cl_mod_prec* iApp; 
	cl_mem dev_cellStatePtr;
    cl_mem dev_iApp, dev_i;

    timestamp_t t0, t1, usecs, tComputeStart, tComputeEnd, tUpdateStart, tUpdateEnd, tReadStart, tReadEnd, tWriteEnd, tWriteStart, tInitStart, tInitEnd, tLoopStart, tLoopEnd, tWriteFileStart, tWriteFileEnd;
    timestamp_t tCompute, tRead, tWrite, tUpdate, tWriteFile, tInit, tLoop;
    tCompute = tRead = tWrite = tUpdate, tWriteFile = tInit = tLoop = 0;

	// parameters used by opencl
	cl_int status;
	cl_uint numPlatforms;
	cl_platform_id *platforms;
	cl_event writeiAppDone, writeiDone, writeCellStatesDone;
	cl_event computeDone;
	cl_event readiAppDone, readCellStatesDone;
	cl_uint numDevices;
	cl_device_id *devices;
	cl_context context;
	cl_command_queue cmdQueue;
	//const char options[] = "-cl-std=CL1.2 -cl-opt-disable";
	const char options[] = "-cl-std=CL1.2";
	size_t localWorkSize[2];
	size_t globalWorkSize[2];
	FILE *computeFile;
    char* computeFileName = "compute_kernel.cl";
	char *computeBuffer;
	size_t computeSize;
	cl_program computeProgram;
	cl_kernel computeKernel;
	
	// timing variables
	timestamp_t avg_Values[timing_variables][averaging_number];
	unsigned int sim_time, networksize_min, networksize_max, blocksize;
	unsigned int exec_count, network_count;
	unsigned int sum;

	// check number of arguments
	if (argc != 5) {
		printf("Usage: %s \n simtime(ms) BLOCKSIZE NETWORKSIZE_MIN NETWORKSIZE_MAX \n", argv[0]);
		return 0;
	}
	
	// get arguments
	sim_time = atoi(argv[1]);
	//printf("sim time (ms):\t %d \n", sim_time);
	blocksize = atoi(argv[2]);
	//printf("block size:\t %d \n", blocksize);
	networksize_min = atoi(argv[3]);
	//printf("network size min:\t %d x %d \n", networksize_min, networksize_min);
	networksize_max = atoi(argv[4]);
	//printf("network size max:\t %d x %d \n", networksize_max, networksize_max);
	if(networksize_min > networksize_max){
		printf("NETWORKSIZE_MAX should be bigger then NETWORKSIZE_MIN\n");
		return 0;
	}
	
	// calculate sim steps
	simSteps = ceil(sim_time/DELTA);
	
	// execute an execution multiple times to average out result
	for(network_count=networksize_min; network_count<=networksize_max; network_count++){
	
		// check if networksize is a multiple of blocksize, if not then skip execution
		if((network_count%blocksize) != 0)
			continue;
	
		// execute an execution multiple times to average out result
		for(exec_count=0; exec_count<averaging_number; exec_count++){
			// get initial timing
			t0 = get_timestamp();
			if(EXTRA_TIMING)
				tInitStart = get_timestamp();

			DEBUG_PRINT(("Inferior Olive Model (%d x %d cell mesh)\n", network_count, network_count));

			//Allocate space for the cellStates, the state pointer and the iapp parameter
			cellStatePtr = (cl_mod_prec*)malloc(2*network_count*network_count*CELL_STATE_SIZE*sizeof(cl_mod_prec));
			if(cellStatePtr==NULL){
				printf("Error: Couldn't malloc for cellStatePtr\n");
				exit(EXIT_FAILURE);
			}
			state_ptr = (cl_bool*)malloc(sizeof(cl_bool));
			if(state_ptr==NULL){
				printf("Error: Couldn't malloc for state_ptr\n");
				exit(EXIT_FAILURE);
			}
			iApp = (cl_mod_prec*)malloc(sizeof(cl_mod_prec));
			if(iApp==NULL){
				printf("Error: Couldn't malloc for iApp\n");
				exit(EXIT_FAILURE);
			}
			
			// initialize iApp
			*state_ptr = CL_FALSE;
			*iApp = previApp = 0;

			//Write initial state values
			InitState(cellStatePtr,network_count);

			//Initialize g_CaL
			init_g_CaL(cellStatePtr,network_count);
			
			// initialize output files
			if(WRITE_OUTPUT){
				if(initializeOutputFiles() != 0)
					exit(-1);
			}

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
			if(platforms==NULL){
				printf("Error: Couldn't malloc for opencl platforms\n");
				exit(EXIT_FAILURE);
			}
		 
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
			
			DEBUG_PRINT(("number of gpu devices found: %d\n", numDevices));
			
			// Allocate enough space for each device
			devices = (cl_device_id*)malloc(numDevices*sizeof(cl_device_id));
			if(devices==NULL){
				printf("Error: Couldn't malloc for opencl devices\n");
				exit(EXIT_FAILURE);
			}

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
			// STEP 5: Create device buffers
			//-----------------------------------------------------
			// create device buffers
			dev_cellStatePtr = clCreateBuffer(context, CL_MEM_READ_WRITE , 2*network_count*network_count*CELL_STATE_SIZE*sizeof(cl_mod_prec), NULL, &status);
			if(status != CL_SUCCESS){
				printf("error creating buffer for cell states\n");
				exit(-1);
			}
			dev_iApp = clCreateBuffer(context, CL_MEM_READ_ONLY , sizeof(cl_mod_prec), NULL, &status);
			if(status != CL_SUCCESS){
				printf("error creating buffer for iApp\n");
				exit(-1);
			}
			dev_i = clCreateBuffer(context, CL_MEM_READ_ONLY , sizeof(cl_bool), NULL, &status);
			if(status != CL_SUCCESS){
				printf("error creating buffer for simulation steps\n");
				exit(-1);
			}
			
			//-----------------------------------------------------
			// STEP 6: Write host data to device buffers
			//-----------------------------------------------------
			if(EXTRA_TIMING)
				tWriteStart = get_timestamp();
			
			// write host data to device
			if(clEnqueueWriteBuffer (cmdQueue, dev_cellStatePtr, CL_FALSE, 0, 2*network_count*network_count*CELL_STATE_SIZE*sizeof(cl_mod_prec), cellStatePtr, 0, NULL, &writeCellStatesDone) != CL_SUCCESS){
				printf("error writing cell states to device memory\n");
				exit(-1);
			}
			if(clEnqueueWriteBuffer (cmdQueue, dev_i, CL_FALSE, 0, sizeof(cl_bool), state_ptr, 0, NULL, &writeiDone) != CL_SUCCESS){
				printf("error writing state pointer to device memory\n");
				exit(-1);
			}
			if(clEnqueueWriteBuffer (cmdQueue, dev_iApp, CL_FALSE, 0, sizeof(cl_mod_prec), iApp, 0, NULL, &writeiAppDone) != CL_SUCCESS){
				printf("error writing iApp to device memory\n");
				exit(-1);
			}
			
			if(EXTRA_TIMING){
				tWriteEnd = get_timestamp();
				tWrite = (tWriteEnd - tWriteStart);
			}

			//-----------------------------------------------------
			// STEP 7: Create and compile the programs
			//-----------------------------------------------------
			// open compute kernel code
			computeFile = fopen(computeFileName, "r");
			if(computeFile == NULL){
				printf("cannot open compute file\n");
				printf("current path: %s\n", computeFileName);
				exit(EXIT_FAILURE);
			}
			fseek(computeFile, 0, SEEK_END);
			computeSize = ftell(computeFile);
			rewind(computeFile);
			computeBuffer = (char*) malloc(computeSize + 1);
			computeBuffer[computeSize] = '\0';
			fread(computeBuffer, sizeof(char), computeSize, computeFile);
			fclose(computeFile);
			// create program from kernel code
			computeProgram = clCreateProgramWithSource(context, 1, (const char**) &computeBuffer, &computeSize, &status);
			if(status != CL_SUCCESS){
				printf("error creating compute kernel program\n");
				exit(-1);
			}
			// build program from kernel code
			if(clBuildProgram(computeProgram, 1, &devices[device_id], options, NULL, NULL) != CL_SUCCESS){
				// Build error log
				size_t len = 0;
				cl_int ret = CL_SUCCESS;
				ret = clGetProgramBuildInfo(computeProgram, *devices, CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
				char *buildlog = calloc(len, sizeof(char));
				ret = clGetProgramBuildInfo(computeProgram, *devices, CL_PROGRAM_BUILD_LOG, len, buildlog, NULL);
				printf("\n\n compute kernel buildlog:   %s\n\n",buildlog);
				switch(ret){
					case CL_INVALID_DEVICE:
						printf("CL_INVALID_DEVICE\n");
						break;
					case CL_INVALID_VALUE:
						printf("CL_INVALID_VALUE\n");
						break;
					case CL_INVALID_PROGRAM:
						printf("CL_INVALID_PROGRAM\n");
						break;
					case CL_OUT_OF_RESOURCES:
						printf("CL_OUT_OF_RESOURCES\n");
						break;
					case CL_OUT_OF_HOST_MEMORY:
						printf("CL_OUT_OF_HOST_MEMORY\n");
						break;
					default:
						printf("default\n");
						break;
				}
				printf("error building compute kernel program\n");
				exit(-1);
			}

			//-----------------------------------------------------
			// STEP 8: Create the kernel
			//-----------------------------------------------------
			// create a kernel from the builded program
			computeKernel = NULL;
			computeKernel = clCreateKernel(computeProgram, "compute_kernel", &status);
			if(status != CL_SUCCESS){
				printf("error creating compute kernel\n");
				exit(-1);
			}

			//-----------------------------------------------------
			// STEP 9: Set the kernel arguments
			//-----------------------------------------------------
			// for compute kernel	
			if(clSetKernelArg(computeKernel, 0, sizeof(cl_mem), &dev_cellStatePtr) != CL_SUCCESS){
				printf("error settings 1st compute kernel argument\n");
				exit(-1);
			}
			if(clSetKernelArg(computeKernel, 1, sizeof(cl_mem), &dev_i) != CL_SUCCESS){
				printf("error settings 3th compute kernel argument\n");
				exit(-1);
			}
			if(clSetKernelArg(computeKernel, 2, sizeof(cl_mem), &dev_iApp) != CL_SUCCESS){
				printf("error settings 4th compute kernel argument\n");
				exit(-1);
			}

			//-----------------------------------------------------
			// STEP 10: Configure the work-item structure
			//-----------------------------------------------------
			globalWorkSize[0] = network_count;		// total number of work items in x direction
			globalWorkSize[1] = network_count;		// total number of work items in y direction
			localWorkSize[0] = blocksize;			// number of hardware threads in x direction
			localWorkSize[1] = blocksize;			// number of hardware threads in y direction
			
			// make sure writes are finished
			if(clWaitForEvents(1, &writeCellStatesDone) != CL_SUCCESS){
				printf("error waiting for cell state write to device memory\n");
				exit(-1);
			}
			if(clWaitForEvents(1, &writeiDone) != CL_SUCCESS){
				printf("error waiting for sim step write to device memory\n");
				exit(-1);
			}
			if(clWaitForEvents(1, &writeiAppDone) != CL_SUCCESS){
				printf("error waiting for iApp write to device memory\n");
				exit(-1);
			}
			
			if(EXTRA_TIMING){
				tInitEnd = get_timestamp();
				tInit = (tInitEnd - tInitStart);
			}
			
			// write output to file
			if(WRITE_OUTPUT){
				writeCellStates(cellStatePtr, network_count, -1);
				writeCellStatesOneCell(cellStatePtr, network_count, -1, 1, 1);
			}
			
			if(EXTRA_TIMING)
				tLoopStart = get_timestamp();
			
			/*-----------------------------------------------------------------------------------------------------------
															simulation loop
			-------------------------------------------------------------------------------------------------------------*/
			for(sim_step=0;sim_step<simSteps;sim_step++){
				//Compute one sim step for all cells
				if(sim_step>20000-1 && sim_step<20500-1){ (*iApp) = 6;} // start @ 1 because skipping initial values
				else{ (*iApp) = 0;}
				
				//-----------------------------------------------------
				// STEP 11.0: Write host data to device buffers
				//-----------------------------------------------------
				if(EXTRA_TIMING)
					tUpdateStart = get_timestamp();
						
				// check if iApp changed and thus needs to be written to the GPU device
				if (previApp != (*iApp)){	
					previApp = (*iApp);
					if(clEnqueueWriteBuffer (cmdQueue, dev_iApp, CL_FALSE, 0, sizeof(cl_mod_prec), iApp, 0, NULL, &writeiAppDone) != CL_SUCCESS){
						printf("error writing cell states to device memory\n");
						exit(-1);
					}
					// wait for write to finish
					if(clWaitForEvents(1, &writeiAppDone) != CL_SUCCESS){
						printf("error waiting for iApp write to device memory\n");
						exit(-1);
					}
				}
				
				if(EXTRA_TIMING){
					tUpdateEnd = get_timestamp();
					tUpdate = (tUpdateEnd - tUpdateStart);
				}
				
				//-----------------------------------------------------
				// STEP 11: Run compute kernel
				//-----------------------------------------------------
				if(EXTRA_TIMING)
					tComputeStart = get_timestamp();
				
				// start neighbour kernel
				if(clEnqueueNDRangeKernel(cmdQueue, computeKernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &computeDone) != CL_SUCCESS){
					printf("error starting compute kernel, probably too much work items requested\n");
					exit(-1);
				}
				
				// wait for kernel to finish
				if(clWaitForEvents(1, &computeDone) != CL_SUCCESS){
					printf("error waiting for compute kernel to finish\n");
					exit(-1);
				}
				if(EXTRA_TIMING){
					tComputeEnd = get_timestamp();
					tCompute = (tComputeEnd - tComputeStart);
				}

				//-----------------------------------------------------
				// STEP 11.3: Read output data from device
				//-----------------------------------------------------
				// transfer data from device to CPU
				if(WRITE_OUTPUT){
					if(EXTRA_TIMING)
						tReadStart = get_timestamp();
						
					// read data from device memory to host memory
					if(clEnqueueReadBuffer(cmdQueue, dev_cellStatePtr, CL_FALSE, 0, 2*network_count*network_count*CELL_STATE_SIZE*sizeof(cl_mod_prec), cellStatePtr, 0, NULL, &readCellStatesDone) != CL_SUCCESS){
						printf("error enqueue-ing read cell states buffer\n");
						exit(-1);
					}

					// wait for read to finish from device to host
					if(clWaitForEvents(1, &readCellStatesDone) != CL_SUCCESS){
						printf("error waiting for reading cell states from device memory\n");
						exit(-1);
					}
					
					if(EXTRA_TIMING){
						tReadEnd = get_timestamp();
						tRead = (tReadEnd - tReadStart);
					}

					// write output to file
					if(EXTRA_TIMING)
						tWriteFileStart = get_timestamp();
					
					// log data
					writeCellStates(cellStatePtr, network_count, sim_step);
					writeCellStatesOneCell(cellStatePtr, network_count, sim_step, 1, 1);
					
					// end timing
					if(EXTRA_TIMING){
						tWriteFileEnd = get_timestamp();
						tWriteFile = (tWriteFileEnd - tWriteFileStart);
					}
				}
				
				if(EXTRA_TIMING)
					tUpdateStart = get_timestamp();
					
				//update state pointer
				*state_ptr = (*state_ptr)^1;
				if(clEnqueueWriteBuffer (cmdQueue, dev_i, CL_FALSE, 0, sizeof(cl_bool), state_ptr, 0, NULL, &writeiDone) != CL_SUCCESS){
					printf("error writing state pointer to device memory\n");
					exit(-1);
				}
				if(clWaitForEvents(1, &writeiDone) != CL_SUCCESS){
					printf("error waiting for sim step write to device memory\n");
					exit(-1);
				}
				
				if(EXTRA_TIMING){
					tUpdateEnd = get_timestamp();
					tUpdate += (tUpdateEnd - tUpdateStart);
				}
			}
			/*-----------------------------------------------------------------------------------------------------------
															End simulation loop
			-------------------------------------------------------------------------------------------------------------*/
			
			// transfer data from device to CPU
			if(EXTRA_TIMING)
				tReadStart = get_timestamp();
				
			// read data from device memory to host memory
			if(clEnqueueReadBuffer(cmdQueue, dev_cellStatePtr, CL_FALSE, 0, 2*network_count*network_count*CELL_STATE_SIZE*sizeof(cl_mod_prec), cellStatePtr, 0, NULL, &readCellStatesDone) != CL_SUCCESS){
				printf("error enqueue-ing read cell states buffer\n");
				exit(-1);
			}

			// wait for read to finish from device to host
			if(clWaitForEvents(1, &readCellStatesDone) != CL_SUCCESS){
				printf("error waiting for reading cell states from device memory\n");
				exit(-1);
			}
			
			if(EXTRA_TIMING){
				tReadEnd = get_timestamp();
				tRead = (tReadEnd - tReadStart);
			}
			
			if(EXTRA_TIMING){
				tLoopEnd = get_timestamp();
				tLoop = (tLoopEnd - tLoopStart);
			}
			
			t1 = get_timestamp();
			usecs = (t1 - t0);// / 1000000;
			
			// get timing values
			avg_Values[0][exec_count] = tInit;
			avg_Values[1][exec_count] = tLoop;
			avg_Values[2][exec_count] = tWrite;
			avg_Values[3][exec_count] = tCompute;
			avg_Values[4][exec_count] = tRead;
			avg_Values[5][exec_count] = tUpdate + tCompute;
			
			DEBUG_PRINT(("%d ms of brain time in %d simulation steps\n", sim_time, simSteps));
			DEBUG_PRINT((" %lld usecs real time \n", usecs));

			if(EXTRA_TIMING){
				// print timing info
				DEBUG_PRINT(("\n"));
				DEBUG_PRINT(("----------------------------------\n"));
				DEBUG_PRINT(("tInit: \t\t %lld \n", tInit));
				DEBUG_PRINT(("tLoop: \t\t %lld \n", tLoop));
				DEBUG_PRINT(("\ttWrite: \t %lld \n", tWrite));
				DEBUG_PRINT(("\ttCompute: \t %lld \n", tCompute));
				DEBUG_PRINT(("\ttRead: \t\t %lld \n", tRead));
				DEBUG_PRINT(("\ttWriteFile:\t %lld \n", tWriteFile));
				DEBUG_PRINT(("\t----------- + \n"));
				DEBUG_PRINT(("\ttSumLoop: \t %lld \n", (tWriteFile + tCompute + tRead)));
				DEBUG_PRINT(("----------------------------------\n"));
				DEBUG_PRINT(("tSum: \t %lld \n", (tInit + tLoop)));
			}
			
			//-----------------------------------------------------
			// STEP 12: Release OpenCL resources
			//----------------------------------------------------- 
			// free opencl resources
			clReleaseMemObject(dev_cellStatePtr);
			clReleaseMemObject(dev_iApp);
			clReleaseMemObject(dev_i);
			clReleaseContext(context);
			clReleaseCommandQueue(cmdQueue);
			clReleaseProgram(computeProgram);
			clReleaseKernel(computeKernel);
			
			//Free up memory
			free(platforms);
			free(devices);
			free(computeBuffer);

			//Free up memory and close files
			free(cellStatePtr);
			free(iApp);
			free(state_ptr);
			if(WRITE_OUTPUT)
				deleteFiles();
			
		} // end of averaging for loop
		
		// compute averages
		for(i=0;i<timing_variables; i++){
			sum = 0;
			for(j=0; j<averaging_number; j++)
				sum += avg_Values[i][j];
			avg_Values[i][0] = sum/averaging_number;
		}
		
		// display timing values
		#ifndef DEBUG
			printf("network size %dx%d", network_count, network_count);
			printf(",block size %dx%d", blocksize, blocksize);
			printf(",sim time %lld ms", sim_time);
			printf(",tInit %lld", avg_Values[0][0]);
			printf(",tLoop %lld", avg_Values[1][0]);
			printf(",tWrite %lld", avg_Values[2][0]);
			printf(",tCompute %lld", avg_Values[3][0]);
			printf(",tRead %lld", avg_Values[4][0]);
			printf(",tLoopCycle %lld", avg_Values[5][0]);
			printf(",tTotal %lld", avg_Values[0][0]+avg_Values[1][0]+avg_Values[4][0]);	// loop + init + read
			printf("\n");
		#endif
		
	} // end of network dimension for loop

    return EXIT_SUCCESS;
}