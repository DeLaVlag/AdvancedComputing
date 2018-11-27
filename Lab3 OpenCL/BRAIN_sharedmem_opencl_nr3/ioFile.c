#include "ioFile.h"
#include "variables.h"

// private variables
static char temp[100];		//warning: this buffer may overflow
static FILE *pOutFileCellstates, *pOutFileCellstatesOne;

int initializeOutputFiles(){
	// variables
	char *outFileNameCellStates = "cellstates_opencl.txt";
	char *outFileNameCellStatesOne = "cellstates_singlecell_opencl.txt";
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
	
	// return result
	return result;
}

void deleteFiles(){
	fclose(pOutFileCellstates);
    fclose(pOutFileCellstatesOne);
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
