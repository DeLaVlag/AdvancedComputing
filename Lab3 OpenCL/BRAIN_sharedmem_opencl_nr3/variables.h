/*** MACROS ***/
#define RAND_INIT 0 		// make it zero to facilitate debugging
//#define SIMTIME   1 		// in ms, for when no input file is provided
#define DELTA 0.05			// simulation step size

//IO network size is IO_NETWORK_DIM1*IO_NETWORK_DIM2
//#define IO_NETWORK_DIM1 	128
//#define IO_NETWORK_DIM2 	128
//#define IO_NETWORK_SIZE 	IO_NETWORK_DIM1*IO_NETWORK_DIM2
// SM block dimensions
#define MAX_BLOCKSIZEX 			8
#define MAX_BLOCKSIZEY 			8
//#define BLOCKSIZE 			BLOCKSIZEX*BLOCKSIZEY

// model precision
//#define DOUBLE_PRECISION	// comment to use single precision

#ifdef DOUBLE_PRECISION
	typedef double mod_prec;
#else
	typedef float mod_prec;
#endif

// printing define
#define IAPP_MAX_CHARS 6 	//2 integer, the dot, 2 decimals and the delimiter

// Cell properties
#define CONDUCTANCE 0.04	//Conductance for neighbors' coupling
// Capacitance
#define C_M 1
// Somatic conductances (mS/cm2)
#define G_NA_S 150      	// Na gate conductance (=90 in Schweighofer code, 70 in paper) 120 too little
#define G_KDR_S 9.0    		// K delayed rectifier gate conductance (alternative value: 18)
#define G_K_S 5      		// Voltage-dependent (fast) potassium
#define G_LS 0.016  		// Leak conductance (0.015)
// Dendritic conductances (mS/cm2)
#define G_K_CA 35       	// Potassium gate conductance (35)
#define G_CAH 4.5     		// High-threshold Ca gate conductance (4.5)
#define G_LD 0.016   		// Dendrite leak conductance (0.015)
#define G_H 0.125    		// H current gate conductance (1.5) (0.15 in SCHWEIGHOFER 2004)
// Axon hillock conductances (mS/cm2)
#define G_NA_A 240      	// Na gate conductance (according to literature: 100 to 200 times as big as somatic conductance)
#define G_NA_R 0      		// Na (resurgent) gate conductance
#define G_K_A 20      		// K voltage-dependent
#define G_LA 0.016  		// Leak conductance
// Cell morphology
#define P1 0.25        		// Cell surface ratio soma/dendrite (0.2)
#define P2 0.15      		// Cell surface ratio axon(hillock)/soma (0.1)
#define G_INT 0.13       	// Cell internal conductance (0.13)
// Reversal potentials
#define V_NA 55       		// Na reversal potential (55)
#define V_K -75       		// K reversal potential
#define V_CA 120       		// Ca reversal potential (120)
#define V_H -43       		// H current reversal potential
#define V_L 10       		// leak current

// program parameters
//#define DEBUG
#define EXTRA_TIMING 1
#define WRITE_OUTPUT 0

// defines cellCompParams mapping
#define CELL_COMP_PARAMS_SIZE	54
#define IAPP_IN					0		// 1 double
#define VNEIGHSTARTADD 			1		// 15 double
#define PREVSTATESTARTADD 		16		// 19 double
#define NEXTSTATESTARTADD 		35		// 19 double

// defines cellState mapping
#define CELL_STATE_SIZE 19
#define DEND_V 			0		// 1 double
#define DEND_H 			1		// 1 double
#define DEND_CAL 		2		// 1 double
#define DEND_P 			3		// 1 double
#define DEND_I 			4		// 1 double
#define DEND_CA2 		5		// 1 double
#define SOMA_G 			6		// 1 double
#define SOMA_V 			7		// 1 double
#define SOMA_SM 		8		// 1 double
#define SOMA_SH 		9		// 1 double
#define SOMA_CK 		10		// 1 double
#define SOMA_CL 		11		// 1 double
#define SOMA_PN 		12		// 1 double
#define SOMA_PP 		13		// 1 double
#define SOMA_PXS 		14		// 1 double
#define AXON_V 			15		// 1 double
#define AXON_SM 		16		// 1 double
#define AXON_SH 		17		// 1 double
#define AXON_P 			18		// 1 double

// addressing macro
#define CELLSTATE_BASEADDR(n,x,y,stride)	((n * stride * stride * CELL_STATE_SIZE) + (x * stride * CELL_STATE_SIZE + y * CELL_STATE_SIZE))
#define CELLSTATE_BASEADDR_OPENCL(n,x,y)	((n * get_global_size(0) * get_global_size(1) * CELL_STATE_SIZE) + (x * get_global_size(1) * CELL_STATE_SIZE + y * CELL_STATE_SIZE))