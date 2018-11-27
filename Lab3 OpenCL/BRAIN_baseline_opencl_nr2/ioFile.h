#include "infoli.h"

int initializeOutputFiles();
void writeOutput(char *,char *,FILE *pOutFile);
void writeOutputDouble(char *, cl_mod_prec, FILE *pOutFile);
void writeCellStates(cl_mod_prec* cellStatePtr, unsigned int network_dimension, unsigned int simStep);
void writeCellStatesOneCell(cl_mod_prec* cellStatePtr, unsigned int network_dimension, unsigned int simStep, unsigned int x, unsigned int y);
void writeCellCompParameters(cl_mod_prec* cellCompParamsPtr, unsigned int network_dimension, unsigned int simStep, cl_mod_prec iApp);
void writeCellCompParametersOneCell(cl_mod_prec* cellCompParamsPtr, unsigned int network_dimension, unsigned int simStep, cl_mod_prec iApp, unsigned int x, unsigned int y);
void deleteFiles();