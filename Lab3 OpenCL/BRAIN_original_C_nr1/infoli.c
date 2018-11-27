/*
 *
 * Copyright (c) 2012, Neurasmus B.V., The Netherlands,
 * web: www.neurasmus.com email: info@neurasmus.com
 *
 * Any use or reproduction in whole or in parts is prohibited
 * without the written consent of the copyright owner.
 *
 * All Rights Reserved.
 *
 *
 * Author: Sebastian Isaza
 * Created: 19-01-2012
 * Modified: 07-08-2012
 *
 * Description: Top source file of the Inferior Olive model, originally written
 * in Matlab by Jornt De Gruijl. It contains the implementation of all functions.
 * The main function allocates the necessary memory, initializes the system
 * state and runs the model calculations.
 *
 */
#include "infoli.h"

// timing
#define averaging_number		3		// execute a run multiple times to test consistency
#define timing_variables		6		// number of timing variables to print

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp ()
{
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

int main(int argc, char *argv[]){

    char *outFileName = "InferiorOlive_Output.txt";
    FILE *pOutFile;
    int i, j, k, p, q;
    int simSteps = 0;
    int simTime = 0;
    int inputFromFile = 0;
    int initSteps;
    cellState ***cellStatePtr;
    cellCompParams **cellCompParamsPtr;
    int seedvar;
    char temp[100];//warning: this buffer may overflow
    mod_prec iApp;
    timestamp_t t0, t1, usecs, tNeighbourStart, tNeighbourEnd, tComputeStart, tComputeEnd, tInitStart, tInitEnd, tLoopStart, tLoopEnd, tWriteFileStart, tWriteFileEnd;
    timestamp_t tNeighbour, tCompute, tWrite, tUpdate, tRead, tWriteFile, tInit, tLoop;
    tNeighbour = tCompute = tWriteFile = tUpdate = tWrite = tRead = tInit = tLoop = 0;
	
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
	for(network_count=networksize_min; network_count<=networksize_max; network_count+=blocksize){
	
		// execute an execution multiple times to average out result
		for(exec_count=0; exec_count<averaging_number; exec_count++){
		
			t0 = get_timestamp();
			if(EXTRA_TIMING){
				tInitStart = get_timestamp();    
			}
			DEBUG_PRINT(("Inferior Olive Model (%d x %d cell mesh)\n", network_count, network_count));

			//Open output file
			if(WRITE_OUTPUT){
				pOutFile = fopen(outFileName,"w");
				if(pOutFile==NULL){
					printf("Error: Couldn't create %s\n", outFileName);
					exit(EXIT_FAILURE);
				}
				writeOutput(temp, ("#simSteps Time(ms) Input(Iapp) Output(V_axon)\n"), pOutFile);
			}


			//Malloc for the array of cellStates and cellCompParams
			mallocCells(&cellCompParamsPtr, &cellStatePtr, network_count);
			//Write initial state values
			InitState(cellStatePtr[0], network_count);
			//Initialize g_CaL
			init_g_CaL(cellStatePtr, network_count);
			
			//Random initialization: put every cell in a different oscillation state
			if(RAND_INIT){
				random_init(cellCompParamsPtr, cellStatePtr, network_count);
			}

			if(EXTRA_TIMING){
				tInitEnd = get_timestamp();
				tLoopStart = get_timestamp();  
			}

			for(i=0;i<simSteps;i++){
				//Compute one sim step for all cells
				if(i>20000-1 && i<20500-1){ iApp = 6;} // start @ 1 because skipping initial values
				else{ iApp = 0;}
				if(WRITE_OUTPUT){
					sprintf(temp, "%d %.2f %.1f ", i+1, i*0.05, iApp); // start @ 1 because skipping initial values
					fputs(temp, pOutFile);
				}
				for(j=0;j<network_count;j++){
					for(k=0;k<network_count;k++){
						if(EXTRA_TIMING){
							tNeighbourStart = get_timestamp();
						}
						neighbors(cellCompParamsPtr, cellStatePtr, i, j, k, network_count);
						if(EXTRA_TIMING){
							tNeighbourEnd = get_timestamp();
							tNeighbour = (tNeighbourEnd - tNeighbourStart);
							tComputeStart = get_timestamp(); 
						}
						compute(cellCompParamsPtr, cellStatePtr, iApp, i, j, k);
						if(EXTRA_TIMING){
							tComputeEnd = get_timestamp();
							tCompute = (tComputeEnd - tComputeStart);
							tWriteFileStart = get_timestamp();
						}

						if(WRITE_OUTPUT)
							writeOutputDouble(temp, cellStatePtr[(i%2)^1][j][k].axon.V_axon, pOutFile);
							
						if(EXTRA_TIMING){
							tWriteFileEnd = get_timestamp();
							tWriteFile = (tWriteFileEnd - tWriteFileStart);
						}
					}
				}
				if(EXTRA_TIMING){
					tWriteFileStart = get_timestamp();
				}
				
				if(WRITE_OUTPUT)
					writeOutput(temp, ("\n"), pOutFile);
					
				if(EXTRA_TIMING){
					tWriteFileEnd = get_timestamp();
					tWriteFile = (tWriteFileEnd - tWriteFileStart);
				}
			}
			if(EXTRA_TIMING){
				tLoopEnd = get_timestamp();
			}

			t1 = get_timestamp();
			usecs = (t1 - t0);// / 1000000;
			
			// get timing values
			avg_Values[0][exec_count] = tInit;
			avg_Values[1][exec_count] = tLoop;
			avg_Values[2][exec_count] = tWrite;
			avg_Values[3][exec_count] = tCompute + tNeighbour;
			avg_Values[4][exec_count] = tRead;
			avg_Values[5][exec_count] = tUpdate + tCompute + tNeighbour;
			
			DEBUG_PRINT(("%d ms of brain time in %d simulation steps\n", simTime, simSteps));
			DEBUG_PRINT((" %lld usecs real time \n", usecs));

			if(EXTRA_TIMING){
				tInit = (tInitEnd - tInitStart);
				tLoop = (tLoopEnd - tLoopStart);
				
				DEBUG_PRINT(("\n"));
				DEBUG_PRINT(("----------------------------------\n"));
				DEBUG_PRINT(("tInit: \t\t %lld \n", tInit));
				DEBUG_PRINT(("tLoop: \t\t %lld \n", tLoop));
				DEBUG_PRINT(("\ttNeighbour: \t %lld \n", tNeighbour));
				DEBUG_PRINT(("\ttCompute: \t %lld \n", tCompute));
				DEBUG_PRINT(("\ttWriteFile: \t %lld \n", tWriteFile));
				DEBUG_PRINT(("\t----------- + \n"));
				DEBUG_PRINT(("\ttSumLoop: \t %lld \n", (tWriteFile + tCompute + tNeighbour)));
				DEBUG_PRINT(("----------------------------------\n"));
				DEBUG_PRINT(("tSum: \t %lld \n", (tInit + tLoop)));
			}


			//Free up memory and close files
			free(cellStatePtr[0]);
			free(cellStatePtr[1]);
			free(cellStatePtr);
			free(cellCompParamsPtr);
			if(WRITE_OUTPUT)
				fclose (pOutFile);
	
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
			printf(",block size %dx%d", 1, 1);
			printf(",sim time %lld ms", sim_time);
			printf(",tInit %lld", avg_Values[0][0]);
			printf(",tLoop %lld", avg_Values[1][0]);
			printf(",tWrite %lld", avg_Values[2][0]);
			printf(",tCompute %lld", avg_Values[3][0]);
			printf(",tRead %lld", avg_Values[4][0]);
			printf(",tLoopCycle %lld", avg_Values[5][0]);
			printf(",tTotal %lld", avg_Values[0][0]+avg_Values[1][0]);	// loop + init
			printf("\n");
		#endif
		
	} // end of network dimension for loop

    return EXIT_SUCCESS;
}

/**
Input: cellCompParamsPtr, cellStatePtr, i, j, k
cellCompParamsPtr: Array of struct which stores values of neighbours for each cell.
cellStatePtr: Array with values for each cell.
i: current simulation step
j: current position in dimension 1 of the IO network
k: current position in dimension 2 of the IO network

Retreive the voltage of the dendrite (V_dend) from each neighbour
**/
void neighbors(cellCompParams **cellCompParamsPtr, cellState ***cellStatePtr, int i, int j, int k, unsigned int network_dimension){
    int n, p, q ;
    n = 0;
    for(p=j-1;p<=j+1;p++){
        for(q=k-1;q<=k+1;q++){
            if(((p!=j)||(q!=k)) && ((p>=0)&&(q>=0)) && ((p<network_dimension)&&(q<network_dimension))){
                cellCompParamsPtr[j][k].neighVdend[n++] = cellStatePtr[i%2][p][q].dend.V_dend;
            }else if(p==j && q==k){
                ;   // do nothing, this is the cell itself
            }else{
                //store same V_dend so that Ic becomes zero by the subtraction
                cellCompParamsPtr[j][k].neighVdend[n++] = cellStatePtr[i%2][j][k].dend.V_dend;
            }
        }
    }
}

/**
Input: cellCompParamsPtr, cellStatePtr, iApp ,i, j, k
cellCompParamsPtr: Array of struct which stores values of neighbours for each cell.
cellStatePtr: Array with values for each cell.
iApp: Extenal input of the dendrite
i: Current simulation step
j: Current position in dimension 1 of the IO network
k: Current position in dimension 2 of the IO network

Retreive the external input of the dedrite 
and update the previous and new state of the current cell.
Then Compute the new variables of the current cell with ComputeOneCell.
**/
void compute(cellCompParams **cellCompParamsPtr, cellState ***cellStatePtr, int iApp, int i, int j, int k){
    cellCompParamsPtr[j][k].iAppIn = iApp;
    cellCompParamsPtr[j][k].prevCellState = &cellStatePtr[i%2][j][k];
    cellCompParamsPtr[j][k].newCellState = &cellStatePtr[(i%2)^1][j][k];
    //Compute one Cell...
    ComputeOneCell(&cellCompParamsPtr[j][k]);
}

