#include "kernel.h"

/**
	Input: cellCompParamsPtr, cellStatePtr, i
	cellCompParamsPtr: Array of struct which stores values of neighbours for each cell.
	cellStatePtr: Array with values for each cell.
	i: current simulation step


	Retreive the voltage of the dendrite (V_dend) from each neighbour
**/
__kernel void neighbor_kernel(global mod_prec *cellStatePtr, global mod_prec *cellCompParamsPtr,  global uint *i) {
	// initialize variables
	int j = get_global_id(0); //current position in dimension 1 of the IO network
	int k = get_global_id(1); //current position in dimension 2 of the IO network
	int n = 0, p, q;

	// get neighbouring cells
	for (p = j - 1; p <= j + 1; p++) {
		for (q = k - 1; q <= k + 1; q++) {
			if (((p != j) || (q != k)) && ((p >= 0) && (q >= 0)) && ((p < get_global_size(0)) && (q < get_global_size(1)))) {
				// cellCompParamsPtr[j][k].neighVdend[n++] = cellStatePtr[i % 2][p][q].dend.V_dend;
				cellCompParamsPtr[CELLCOMP_BASEADDR_OPENCL(j,k) + VNEIGHSTARTADD + (n++)] = cellStatePtr[CELLSTATE_BASEADDR_OPENCL(((*i) % 2),p,q) + DEND_V];
			}
			else if (p == j && q == k) {
			}
			else {
				// cellCompParamsPtr[j][k].neighVdend[n++] = cellStatePtr[i % 2][j][k].dend.V_dend;
				cellCompParamsPtr[CELLCOMP_BASEADDR_OPENCL(j,k) + VNEIGHSTARTADD + (n++)] = cellStatePtr[CELLSTATE_BASEADDR_OPENCL(((*i) % 2),j,k) + DEND_V];
			}
		}
	}

}