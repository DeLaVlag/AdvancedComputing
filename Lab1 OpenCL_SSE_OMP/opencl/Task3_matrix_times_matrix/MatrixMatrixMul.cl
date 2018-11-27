__kernel void mul_kernel(global double *matrix1_in, global double *matrix2_in, global double *matrix_out, int size){
	int rows1 = get_global_id(0); 
	int cols = get_global_id(1); 
	double value = 0.0;
	int rows2;

	for(rows2=0; rows2<size; rows2++)
	  value += matrix1_in[rows1*size + rows2] * matrix2_in[rows2*size + cols];
	  
	matrix_out[rows1*size + cols] = value;
}
