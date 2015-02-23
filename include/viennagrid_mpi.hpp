


//MPI_Datatype MPI_VIENNAWD_PARTICLE;

//int           type_size  = 3;
//int           blocklen[] = {1, 7, 5};
//MPI_Datatype  types[]    = {MPI_UNSIGNED_CHAR, MPI_SHORT, MPI_FLOAT};
//MPI_Aint      disp[type_size];

//MPI_Address(&(particle->active),   &disp[0]);
//MPI_Address(&(particle->sign),     &disp[1]);
//MPI_Address(&(particle->position), &disp[2]);

//for(int i = type_size-1; i >= 0; i--)
//disp[i] -= disp[0];

//MPI_Type_struct(type_size, blocklen, disp, types, MPI_VIENNAWD_PARTICLE);
//MPI_Type_commit(MPI_VIENNAWD_PARTICLE);



//MPI_Type_free(&MPI_VIENNAWD_PARTICLE);
