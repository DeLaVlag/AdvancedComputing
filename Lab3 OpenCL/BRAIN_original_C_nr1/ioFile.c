#include "infoli.h"

void writeOutput(char *temp, char* s, FILE *pOutFile){
    if(WRITE_OUTPUT){
        sprintf(temp, "%s", s);
        fputs(temp, pOutFile);
    }
}

void writeOutputDouble(char *temp, double d, FILE *pOutFile){
        if(WRITE_OUTPUT){
        sprintf(temp, "%.8f ", d);
        fputs(temp, pOutFile);
    }
}