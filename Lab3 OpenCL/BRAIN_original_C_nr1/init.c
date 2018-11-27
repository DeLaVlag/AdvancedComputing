#include "infoli.h" 

void mallocCells(cellCompParams ***cellCompParamsPtr, cellState ****cellStatePtr, unsigned int network_dimension){
	int k;
    DEBUG_PRINT(("cellStatePtr: %luB\n", 2*network_dimension*network_dimension*sizeof(cellState)));
    //Two cell state structs are needed so as to avoid having to synchronize all consumers before they start rewriting the cell state.
    (*cellStatePtr) = malloc(2*sizeof(cellState *));//current and next state
    if((*cellStatePtr)==NULL){
        printf("Error: Couldn't malloc for cellStatePtr\n");
        exit(EXIT_FAILURE);
    }
    (*cellStatePtr)[0] = malloc(network_dimension*sizeof(cellState *));
    if(cellStatePtr[0]==NULL){
        printf("Error: Couldn't malloc for cellStatePtr[0]\n");
        exit(EXIT_FAILURE);
    }
    for(k=0;k<network_dimension;k++){
        (*cellStatePtr)[0][k] = malloc(network_dimension*sizeof(cellState));
        if((*cellStatePtr)[0][k]==NULL){
            printf("Error: Couldn't malloc for cellStatePtr[0][k]\n");
            exit(EXIT_FAILURE);
        }
    }
    (*cellStatePtr)[1] = malloc(network_dimension*sizeof(cellState));
    if((*cellStatePtr)[1]==NULL){
        printf("Error: Couldn't malloc for cellStatePt[1]r\n");
        exit(EXIT_FAILURE);
    }
    for(k=0;k<network_dimension;k++){
        (*cellStatePtr)[1][k] = malloc(network_dimension*sizeof(cellState));
        if((*cellStatePtr)[1][k]==NULL){
            printf("Error: Couldn't malloc for cellStatePtr[1][k]\n");
            exit(EXIT_FAILURE);
        }
    }

    DEBUG_PRINT(("cellCompParamsPtr: %luB\n", network_dimension*network_dimension*sizeof(cellCompParams)));
    (*cellCompParamsPtr) = malloc(network_dimension*sizeof(cellCompParams *));
    if((*cellCompParamsPtr) ==NULL){
        printf("Error: Couldn't malloc for cellCompParamsPtr\n");
        exit(EXIT_FAILURE);
    }
    for(k=0;k<network_dimension;k++){
        (*cellCompParamsPtr) [k] = malloc(network_dimension*sizeof(cellCompParams));
        if((*cellCompParamsPtr) [k]==NULL){
            printf("Error: Couldn't malloc for cellCompParamsPtr[k]\n");
            exit(EXIT_FAILURE);
        }
    }
}

void InitState(cellState **cellStatePtr, unsigned int network_dimension){
    int j, k;
    cellState initState;
    //Initial dendritic parameters
    initState.dend.V_dend = -60;
    initState.dend.Calcium_r = 0.0112788;// High-threshold calcium
    initState.dend.Potassium_s = 0.0049291;// Calcium-dependent potassium
    initState.dend.Hcurrent_q = 0.0337836;// H current
    initState.dend.Ca2Plus = 3.7152;// Calcium concentration
    initState.dend.I_CaH   = 0.5;// High-threshold calcium current
    //Initial somatic parameters
    initState.soma.g_CaL = 0.68; //default arbitrary value but it should be randomized per cell
    initState.soma.V_soma = -60;
    initState.soma.Sodium_m = 1.0127807;// Sodium (artificial)
    initState.soma.Sodium_h = 0.3596066;
    initState.soma.Potassium_n = 0.2369847;// Potassium (delayed rectifier)
    initState.soma.Potassium_p = 0.2369847;
    initState.soma.Potassium_x_s = 0.1;// Potassium (voltage-dependent)
    initState.soma.Calcium_k = 0.7423159;// Low-threshold calcium
    initState.soma.Calcium_l = 0.0321349;
    // Initial axonal parameters
    initState.axon.V_axon = -60;
    //sisaza: Sodium_m_a doesn't have a state, therefore this assignment doesn'thave any effect
    initState.axon.Sodium_m_a = 0.003596066;// Sodium (thalamocortical)
    initState.axon.Sodium_h_a = 0.9;
    initState.axon.Potassium_x_a = 0.2369847;// Potassium (transient)

    //Copy init sate to all cell states
    for(j=0;j<network_dimension;j++){
        for(k=0;k<network_dimension;k++){
            memcpy(&cellStatePtr[j][k], &initState, sizeof(cellState));
        }
    }

    return;
}

void init_g_CaL(cellState ***cellStatePtr, unsigned int network_dimension){
	int seedvar, j, k;
	seedvar = 1;
    for(j=0;j<network_dimension;j++){
        for(k=0;k<network_dimension;k++){
            srand(seedvar++);   // use this for debugging, now there is difference
            cellStatePtr[1][j][k].soma.g_CaL = cellStatePtr[0][j][k].soma.g_CaL = 0.68;
            // Uncomment the next two lines to assign different soma conductances to each cell.
            //cellStatePtr[0][j][k].soma.g_CaL = 0.6+(0.2*(rand()%100)/100);
            //cellStatePtr[1][j][k].soma.g_CaL = cellStatePtr[0][j][k].soma.g_CaL;

        }
    }
}

void random_init(cellCompParams **cellCompParamsPtr, cellState ***cellStatePtr, unsigned int network_dimension){
	int seedvar, i, j, k, initSteps;
	seedvar = 1;
    for(j=0;j<network_dimension;j++){
        for(k=0;k<network_dimension;k++){
            //Put each cell at a different random state
            //srand(time(NULL));//Initialize random seed - Too fast when called in a loop.
            srand(seedvar++);   // use this for debugging, now there is difference
            initSteps = rand()%(int)ceil(100/DELTA);
            initSteps = initSteps | 0x00000001;//make it odd, so that the final state is in prevCellState
            printf("%d iterations - ",initSteps);
            for(i=0;i<initSteps;i++){
                //Arrange inputs
                cellCompParamsPtr[j][k].iAppIn = 0;//No stimulus
                cellCompParamsPtr[j][k].prevCellState = &cellStatePtr[i%2][j][k];
                cellCompParamsPtr[j][k].newCellState = &cellStatePtr[(i%2)^1][j][k];
                ComputeOneCell(&cellCompParamsPtr[j][k]);
            }
            printf("Random initialization of the cell states finished.\n");
        }
    }
}