#include "init.h" 

void InitState(cl_mod_prec *cellStatePtr, unsigned int network_dimension){
	// variables
    int j,k,n;
    cl_mod_prec  cellStateInit[CELL_STATE_SIZE];
	
    //Initial dendritic parameters
    cellStateInit[DEND_V]   = -60;  
    cellStateInit[DEND_H]   = 0.0337836;
    cellStateInit[DEND_CAL] = 0.0112788;
    cellStateInit[DEND_P]   = 0.0049291;
    cellStateInit[DEND_I]   = 0.5;
    cellStateInit[DEND_CA2] = 3.7152;

    cellStateInit[SOMA_G]   = 0.68;
    cellStateInit[SOMA_V]   = -60;
    cellStateInit[SOMA_SM]  = 1.0127807;
    cellStateInit[SOMA_SH]  = 0.3596066;
    cellStateInit[SOMA_CK]  = 0.7423159;
    cellStateInit[SOMA_CL]  = 0.0321349;
    cellStateInit[SOMA_PN]  = 0.2369847;
    cellStateInit[SOMA_PP]  = 0.2369847;
    cellStateInit[SOMA_PXS] = 0.1;

    cellStateInit[AXON_V]   = -60;
    cellStateInit[AXON_SM]  = 0.003596066;
    cellStateInit[AXON_SH]  = 0.9;
    cellStateInit[AXON_P]   = 0.2369847;

    //Copy init sate to all cell states
	for(j=0;j<network_dimension;j++){
        for(k=0;k<network_dimension;k++){
			for(n=0; n<CELL_STATE_SIZE; n++){
				cellStatePtr[CELLSTATE_BASEADDR(0,j,k,network_dimension)+n] = cellStateInit[n];
			}
        }
    }

    return;
}

void init_g_CaL(cl_mod_prec *cellStatePtr, unsigned int network_dimension){
	// variables
    int seedvar,j,k;
    seedvar = 1;
	
	//Copy init sate to all cell states
	for(j=0;j<network_dimension;j++){
        for(k=0;k<network_dimension;k++){
			srand(seedvar++);   // use this for debugging, now there is difference
			cellStatePtr[CELLSTATE_BASEADDR(0,j,k,network_dimension)+SOMA_G] = cellStatePtr[CELLSTATE_BASEADDR(1,j,k,network_dimension)+SOMA_G] = 0.68;
        }
    }
}
