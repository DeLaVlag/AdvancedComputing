#include "ioFile.h"
#include "variables.h"

// private variables
static char temp[100];		//warning: this buffer may overflow
static FILE *pOutFileCellstates, *pOutFileCellstatesOne;
static FILE *pOutFileCellcomp, *pOutFileCellcompOne;

int initializeOutputFiles(){
	// variables
	char *outFileNameCellStates = "cellstates_opencl.txt";
	char *outFileNameCellStatesOne = "cellstates_singlecell_opencl.txt";
	char *outFileNameCellComp = "cellcomp_opencl.txt";
	char *outFileNameCellCompOne = "cellcomp_singlecell_opencl.txt";
	int result = 0;
	
	// open cell states file
    pOutFileCellstates = fopen(outFileNameCellStates,"w");
    if(pOutFileCellstates==NULL){
        printf("Error: Couldn't create %s\n", outFileNameCellStates);
        result = -1;
    }
    writeOutput(temp, ("#simSteps statesConfig(0/1) AxonVoltages(V)\n"), pOutFileCellstates);
	
	// open cell state single cell file
	pOutFileCellstatesOne = fopen(outFileNameCellStatesOne,"w");
    if(pOutFileCellstatesOne==NULL){
        printf("Error: Couldn't create %s\n", outFileNameCellStatesOne);
        result = -1;
    }
    writeOutput(temp, ("#simSteps statesConfig(0/1) States(-)\n"), pOutFileCellstatesOne);
	
	// open cell comp parameters file
	pOutFileCellcomp = fopen(outFileNameCellComp,"w");
    if(pOutFileCellcomp==NULL){
        printf("Error: Couldn't create %s\n", outFileNameCellComp);
        result = -1;
    }
    writeOutput(temp, ("#simSteps iApp(-) statesConfig(0/1) AxonVoltages(V)\n"), pOutFileCellcomp);
	
	// open cell comp parameters single cell file
	pOutFileCellcompOne = fopen(outFileNameCellCompOne,"w");
    if(pOutFileCellcompOne==NULL){
        printf("Error: Couldn't create %s\n", outFileNameCellCompOne);
        result = -1;
    }
    writeOutput(temp, ("#simSteps iApp(-) Neighbours statesConfig(0/1) States(-)\n"), pOutFileCellcompOne);
	
	// return result
	return result;
}

void deleteFiles(){
	fclose(pOutFileCellstates);
    fclose(pOutFileCellstatesOne);
	fclose(pOutFileCellcomp);
	fclose(pOutFileCellcompOne);
}

void writeOutput(char *temp, char* s, FILE *pOutFile){
    if(WRITE_OUTPUT){
        sprintf(temp, "%s", s);
        fputs(temp, pOutFile);
    }

}

void writeOutputDouble(char *temp, cl_mod_prec d, FILE *pOutFile){
    if(WRITE_OUTPUT){
    sprintf(temp, "%.8f ", d);
    fputs(temp, pOutFile);
    }
}

/**
 *	cellStatePtr the pointer to the cell states
 *	simStep the simulation step
 * 			columns are #simSteps, statesConfig(0/1), AxonVoltages(V)
 */
void writeCellStates(cl_mod_prec* cellStatePtr, unsigned int network_dimension, unsigned int simStep){
	// variables
	unsigned int j,k;
	
    if(WRITE_OUTPUT){
		sprintf(temp, "%d %d ", simStep, 0);
		fputs(temp, pOutFileCellstates);
		for(j = 0; j < network_dimension; j++){
			for(k = 0; k < network_dimension; k++){
				writeOutputDouble(temp, cellStatePtr[CELLSTATE_BASEADDR(0,j,k,network_dimension) + AXON_V], pOutFileCellstates);
			}                
		}
		writeOutput(temp, ("\n"), pOutFileCellstates);
		sprintf(temp, "%d %d ", simStep, 1);
		fputs(temp, pOutFileCellstates);
		for(j = 0; j < network_dimension; j++){
			for(k = 0; k < network_dimension; k++){
				writeOutputDouble(temp, cellStatePtr[CELLSTATE_BASEADDR(1,j,k,network_dimension) + AXON_V], pOutFileCellstates);
			}                
		}
		writeOutput(temp, ("\n"), pOutFileCellstates);
    }
}

/**
 *	cellStatePtr the pointer to the cell states
 *	simStep the simulation step
 *	x location of cell in dim 1
 *	y location of cell in dim 2
 * 			columns are #simSteps, statesConfig(0/1), States(-)
 */
void writeCellStatesOneCell(cl_mod_prec* cellStatePtr, unsigned int network_dimension, unsigned int simStep, unsigned int x, unsigned int y){
	// variables
	unsigned int n;
	
    if(WRITE_OUTPUT){
		sprintf(temp, "%d %d ", simStep, 0);
		fputs(temp, pOutFileCellstatesOne);
		for(n=0; n<CELL_STATE_SIZE; n++){
			writeOutputDouble(temp, cellStatePtr[CELLSTATE_BASEADDR(0,x,y,network_dimension) + n], pOutFileCellstatesOne);         
		}
		writeOutput(temp, ("\n"), pOutFileCellstatesOne);
		sprintf(temp, "%d %d ", simStep, 1);
		fputs(temp, pOutFileCellstatesOne);
		for(n=0; n<CELL_STATE_SIZE; n++){
			writeOutputDouble(temp, cellStatePtr[CELLSTATE_BASEADDR(1,x,y,network_dimension) + n], pOutFileCellstatesOne);         
		}
		writeOutput(temp, ("\n"), pOutFileCellstatesOne);
    }
}

/**
 *	cellCompParamsPtr the pointer to the cell compute parameters
 *	simStep the simulation step
 * 	iApp the weird variable that nobody knows what it does
 * 			columns are #simSteps, iApp(-), statesConfig(0/1), AxonVoltages(V)
 */
void writeCellCompParameters(cl_mod_prec* cellCompParamsPtr, unsigned int network_dimension, unsigned int simStep, cl_mod_prec iApp){
	// variables
	unsigned int j,k;
	
    if(WRITE_OUTPUT){
		sprintf(temp, "%d %d %.2f", simStep, 0, iApp);
		fputs(temp, pOutFileCellcomp);
		for(j = 0; j < network_dimension; j++){
			for(k = 0; k < network_dimension; k++){
				writeOutputDouble(temp, cellCompParamsPtr[CELLCOMP_BASEADDR(j,k,network_dimension) + PREVSTATESTARTADD + AXON_V], pOutFileCellcomp);
			}                
		}
		writeOutput(temp, ("\n"), pOutFileCellcomp);
		sprintf(temp, "%d %d %.2f", simStep, 1, iApp);
		fputs(temp, pOutFileCellcomp);
		for(j = 0; j < network_dimension; j++){
			for(k = 0; k < network_dimension; k++){
				writeOutputDouble(temp, cellCompParamsPtr[CELLCOMP_BASEADDR(j,k,network_dimension) + NEXTSTATESTARTADD + AXON_V], pOutFileCellcomp);
			}                
		}
		writeOutput(temp, ("\n"), pOutFileCellcomp);
    }
}

/**
 *	cellCompParamsPtr the pointer to the cell compute parameters
 *	simStep the simulation step
 * 	iApp the weird variable that nobody knows what it does
 *	x location of cell in dim 1
 *	y location of cell in dim 2
 * 			columns are #simSteps, iApp(-), Neighbours, statesConfig(0/1), States(-)
 */
void writeCellCompParametersOneCell(cl_mod_prec* cellCompParamsPtr, unsigned int network_dimension, unsigned int simStep, cl_mod_prec iApp, unsigned int x, unsigned int y){
	// variables
	unsigned int n;
	
    if(WRITE_OUTPUT){
		sprintf(temp, "%d %.2f ", simStep, iApp);
		fputs(temp, pOutFileCellcompOne);
		for(n=0; n<15; n++){
			writeOutputDouble(temp, cellCompParamsPtr[CELLCOMP_BASEADDR(x,y,network_dimension) + VNEIGHSTARTADD + n], pOutFileCellcompOne);         
		}
		writeOutput(temp, (" \t0 "), pOutFileCellcompOne);
		for(n=0; n<CELL_STATE_SIZE; n++){
			writeOutputDouble(temp, cellCompParamsPtr[CELLCOMP_BASEADDR(x,y,network_dimension) + PREVSTATESTARTADD + n], pOutFileCellcompOne);         
		}
		writeOutput(temp, ("\n"), pOutFileCellcompOne);
		sprintf(temp, "%d %.2f ", simStep, iApp);
		fputs(temp, pOutFileCellcompOne);
		for(n=0; n<15; n++){
			writeOutputDouble(temp, cellCompParamsPtr[CELLCOMP_BASEADDR(x,y,network_dimension) + VNEIGHSTARTADD + n], pOutFileCellcompOne);         
		}
		writeOutput(temp, (" \t1 "), pOutFileCellcompOne);
		for(n=0; n<CELL_STATE_SIZE; n++){
			writeOutputDouble(temp, cellCompParamsPtr[CELLCOMP_BASEADDR(x,y,network_dimension) + NEXTSTATESTARTADD + n], pOutFileCellcompOne);         
		}
		writeOutput(temp, ("\n"), pOutFileCellcompOne);
		sprintf(temp, "%.2f \n", cellCompParamsPtr[CELLCOMP_BASEADDR(x,y,network_dimension) + IAPP_IN]);
		fputs(temp, pOutFileCellcompOne);
    }
}
