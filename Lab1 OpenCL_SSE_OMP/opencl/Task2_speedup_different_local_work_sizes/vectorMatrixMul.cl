__kernel void mul_kernel(global double *vector_in, global double *matrix_in, global double *vector_out, int size){
  int id = get_global_id(0); 
  double value = 0;
  int k;

  for(k = 0; k < size; k++){
    value += matrix_in[(k* size) + id] * vector_in[k];
  }

  vector_out[id] = value;

}
