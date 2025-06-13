#define ROUND_TO_INT(x) ( (int)(x + 0.5) )

/*
THINGS TO CHECK:
    -- 2D and 3D
    -- Serial and Parallel MPI
    -- ASCII and BINARY
    -- float or double

    COMBINATIONS
    1) 2D + Serial + ASCII
    2) [ok] 3D + Serial + ASCII
    3) [ok] 2D + Parallel + ASCII
    4) [ok] 3D + Parallel + ASCII

    5) 2D + Serial + BINARY + float
    6) [ok] 3D + Serial + BINARY + float
    7) 2D + Parallel + BINARY + float
    8) [ok] 3D + Parallel + BINARY + float

    9) 2D + Serial + BINARY + double
    10) [ok] 3D + Serial + BINARY + double
    11) 2D + Parallel + BINARY + double
    12) [ok] 3D + Parallel + BINARY + double

ORDER OF THINGS TO DO:
    -- Try to run a 3D version on cluster in parallel without getting that divergence things
    -- Check that the ASCII mesh function prints correctly in sequential and parallel
    -- Check that the BINARY mesh function prints correctly in sequenttial and parallel
*/

char *folder_name = NULL;
FILE *log_file = NULL;
FILE *error_file = NULL;


/** 
Checks a condition and, if false, prints an error message and ends execution if requested
(Similar to the standard "assert" macro). */
void ErrorMessage(int Condition, const char *Message, const char *FunctionName, char EndExecution)
{   
  // If the condition is satisfied, we do nothing
  if( Condition )
    return;

  printf("  ===== ERROR in function %s ...\n", (FunctionName) ? FunctionName : "NOT-SPECIFIED");
  printf("  ===== %s\n", Message ? Message : "");
  if( EndExecution ) {
    printf("  ===== Ending execution...\n\n");
    exit(0);
  }

  return;
}

/**
Returns a command line argument given by the user when calling the program
Input Parameter: the index (starting from zero) of the argument to be returned
Output Parameter: string location to store the argument. */
void GetCommandLineArgument(int ArgumentIndex, char *ReturnArgument, int argc, char *argv[])
{   
  /// Checking if the user has provided the argument with the requested index
  if( argc<=ArgumentIndex+1 ) {
    printf("  ===== ERROR in function GetCommandLineArgument ...\n");
    printf("  ===== The requested parameter %d was not provided when calling the program.\n", ArgumentIndex);
    printf("  ===== Ending execution...\n\n");
    exit(0);
  }

  /// Copying the argument into the return string
  strcpy(ReturnArgument, argv[ArgumentIndex + 1]);
  return;
}

/* Parameter: 
-> attempt_restart: tries to restart the simulation from a checkpoint, if available in the folder
-> the name of the folder to be created. You can use the same syntax as "printf" 
to format this folder name */
void OpenSimulationFolder(bool attempt_restart, const char* format, ...)
{
  char *folder_temp = (char *)malloc(1000*sizeof(char));
  folder_name = (char *)malloc(1100*sizeof(char));

  // === Formatting the string with the function parameters
  va_list argptr;
  va_start(argptr, format);
  vsprintf(folder_temp, format, argptr);
  va_end(argptr);

  // === Putting the folder name in the global string for future use
  sprintf(folder_name, "outputs/%s/", folder_temp);

  // === Only root process continues below here
  if( pid() )
    return;

  // === Creating the base "outputs" folder
  system("mkdir outputs/");

  // === Creating the simulation folder
  char command[1000];
  sprintf(command, "mkdir %s", folder_name);
  system(command);

  FILE *file_checkpoint = NULL;
  if( attempt_restart ) {
    char checkpoint_file_name[900];
    sprintf(checkpoint_file_name, "%s/checkpoint", folder_name);
    file_checkpoint = fopen(checkpoint_file_name, "r");
  }

  /// === Opening the new log file
  if( log_file ) 
    fclose(log_file);
  sprintf(command, "%s/log_file.txt", folder_name);
  log_file = fopen(command, file_checkpoint ? "at" : "wt");
  if( !log_file ) {
    printf("  ===== ERROR in function OpenSimulationFolder ...\n");
    printf("  ===== Unable to create the folder and log_file. Do you have permission for creating this folder?\n");
    printf("  ===== Ending execution...\n\n");
    exit(1);
  }

  /// === Opening the new error file
  if( error_file ) 
    fclose(error_file);
  sprintf(command, "%s/error_file.txt", folder_name);
  error_file = fopen(command, file_checkpoint ? "at" : "wt");
  if( !error_file ) {
    printf("  ===== ERROR in function OpenSimulationFolder ...\n");
    printf("  ===== Unable to create the folder and error_file. Do you have permission for creating this folder?\n");
    printf("  ===== Ending execution...\n\n");
    exit(1);
  }

  // Changing the basilisk standard error output "ferr" to be our custom file
  ferr = error_file;
  
  free(folder_temp);
  if( file_checkpoint ) {
    fclose(file_checkpoint);
  }
  return;
}

/**
Closes the current simulation folder and output files currently open. */
void CloseSimulationFolder() {

  // === Only root process does this
  if( pid() )
    return;

  if( log_file )
    fclose(log_file);
  log_file = NULL;

  return;
}


void CreateBasiliskCheckpoint() {
  char file_name[900];
  sprintf(file_name, "%s/checkpoint", folder_name);
  dump(file_name);
}

bool RestoreBasiliskCheckpoint() {
  char file_name[900];
  sprintf(file_name, "%s/checkpoint", folder_name);
  return restore (file = file_name);
}

/**
Read parameters from a file
Parameters: 
1. The name of the file containing the parameters.
2. A list of parameter names to look for in the file. Separated by semi-colon.
3. Pointers to variables where each of the parameters will be stored
NOTE: currently I'm reading all values as "double" variables. 
I can change this in future if there is interest for it */
// void ReadParametersFromFile(const char *file_name, const char *parameters_names, ...)
// {
//   char names_copy[1000], *token;
//   int number_of_parameters = 0;
//   FILE *file;

//   // == Counting how many parameters the user wants to read
//   strcpy(names_copy, parameters_names);
//   token = strtok(names_copy, ";");
//   while( token ) {
//     number_of_parameters++;
//     token = strtok(NULL, ";");
//   }

//   // == Looking for the parameters in the text file
//   if( !(file=fopen(file_name, "rt")) ) {
//     printf("  ===== ERROR in function ReadParametersFromFile ...\n");
//     printf("  ===== Unable to open the input text file \"%s\".\n", file_name);
//     printf("  ===== Ending execution...\n\n");
//     exit(0);
//   }
//   strcpy(names_copy, parameters_names);
//   token = strtok(names_copy, ";");
//   int i_token = 0;
//   while( token ) {
//     char name_read[1000];
//     double value_read;
    
//     // Scanning the entire file from the beginning loking for this token
//     rewind(file);
//     int found_parameter = 0;
//     while( 2==fscanf(file, "%s %lf\n", name_read, &value_read) ) {
//       if( !strcmp(name_read, token) ) {
//         found_parameter = 1;
//         break;
//       }
//     }
//     if( !found_parameter ) {
//       printf("  ===== ERROR in function ReadParametrsFromFile ...");
//       printf("  ===== Could not find the parameter [%s] in the input file.\n", token);
//       printf("  ===== Ending execution...\n\n");
//       exit(0);
//     }

//     // Finding the pointer in the list of pointers given by the user in the function parameters
//     va_list ap;
//     va_start(ap, parameters_names);
//     int i;
//     double *pointer;
//     for( i=0; i<=i_token; i++ )
//       pointer = va_arg(ap, double*);
//     va_end(ap);

//     // Setting the value read in the pointer
//     *pointer = value_read;

//     token = strtok(NULL, ";");
//     i_token++;
//   }
//   fclose( file );

//   return;
// }

/**
Prints something to the log file of this simulation */
void PrintLog(bool print_stdout, const char* format, ...)
{
  // === Only root process will print
  if( pid() )
    return;
  
  if( log_file==NULL )
    return;

  // Formatted printing
  va_list argptr;
  va_start(argptr, format);
  vfprintf(log_file, format, argptr);
  va_end(argptr);
  fflush(log_file);

  if( print_stdout ) {
    va_start(argptr, format);
    vfprintf(stdout, format, argptr);
    va_end(argptr);
  }
  
  return;
}



/**
This is basically a wrapper for the MPI_Gatherv function.
It automatically counts how much data will be sent by each processor */
#if _MPI
void MPI_Gather_Uneven(const void *sendbuf, int sendcount, MPI_Datatype sendreceivetype, void *recvbuf, int root)
{
  int array_counts[npe()];
  int displs[npe()];

  if( pid()==root ) {
    int i;
    for( i=0; i<npe(); i++ ) {
      if( i==root )
        array_counts[i] = sendcount;
      else
        MPI_Recv(array_counts+i, 1, MPI_INT, i, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      displs[i] = (i==0) ? 0 : displs[i-1] + array_counts[i-1];
    }
  }
  else
    MPI_Send(&sendcount, 1, MPI_INT, 0, 1000, MPI_COMM_WORLD);

  MPI_Gatherv(sendbuf, sendcount, sendreceivetype, recvbuf, array_counts, displs, sendreceivetype, root, MPI_COMM_WORLD);
}
#endif

// This function transforms the little-endian entries of an array into big-endian (or vice-versa). This will be used in the BINARY vtk printing functions
// IMPORTANT: currently, this function will also swap the order of the array entries as well as the byte ordering of each individual entry
void SwapArrayBytes(void *Array, size_t number_of_bytes)
{
    int i;
    size_t half_of_bytes = number_of_bytes/2;

    unsigned char *byte_array = (unsigned char *)Array;
    for( i=0; i<half_of_bytes; i++ ) {
        unsigned char temporary = byte_array[i];
        byte_array[i] = byte_array[number_of_bytes - i - 1];
        byte_array[number_of_bytes - i - 1] = temporary;
    }

    return;
}

/**
Prints mesh and (optionally) scalars to a single-precision BINARY VTK file. */
void PrintMeshVTK_Binary_Float(int n, float time, bool print_only_droplet,
                              scalar *list_scalar_data, const char **list_scalar_names,
                              vector *list_vector_data = NULL, const char **list_vector_names = NULL)
{
  FILE *arq;
  char nomeArq[900];

  // ===  Cell is either a square (2D) or cube (3D), meaning either 4 or 8 vertices per cell
  int vertices_per_cell = (dimension==2) ? 4 : 8;

  // === VTK cell code that represents voxels
  // === (check Figure 2 here to understand: https://kitware.github.io/vtk-examples/site/VTKFileFormats/)
  int vtk_cell_type = (dimension==2) ? 9 : 11;

  // === Counting how many local cells we have in the mesh
  // === NOTE: whenever I say "local", i mean things that are locally in this processor (in case of parallel simulation)
  int local_num_cells = 0;
  foreach(serial) {
    if( (!print_only_droplet) || ((f1[] + f2[])>0.01) )
      local_num_cells++;
  }

  // === Allocatting memory for the local vertex arrays
  // === Note: the x,y,z coordinates will all go into the same array, 
  // === so the array will be [x1, y1, z1, x2, y2, z2, ...]
  float *local_vertices = (float *)malloc( 3*vertices_per_cell*local_num_cells*sizeof(float) );

  // === Allocating memory for ALL the local scalar data arrays and then all the local vector data arrays after
  int number_of_scalar_fields = list_len(list_scalar_data);
  int number_of_vector_fields = vectors_len(list_vector_data);
  typedef float** floatpp; // qcc gets uncomfortable if i dont hide the double pointer (?????)
  floatpp local_data = number_of_scalar_fields || number_of_vector_fields ? (float **)malloc( (number_of_scalar_fields + 2*number_of_vector_fields)*sizeof(float *) ) : NULL;
  for( int k=0; k<number_of_scalar_fields + 2*number_of_vector_fields; k++ )
    local_data[k] = (float *)malloc( local_num_cells*sizeof(float) );

  // === Storing all the vertices coordinates in the arrays
  int i = 0, i_cell = 0;
  foreach(serial) {
    if( (!print_only_droplet) || ((f1[] + f2[])>0.01) ) {

      /// Using a macro conditional to avoid checking dimension==2 every loop during execution time...
      #if dimension==2
        local_vertices[i++] = 0.0; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = 0.0; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        local_vertices[i++] = 0.0; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        local_vertices[i++] = 0.0; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x-0.5*Delta;
      #else // dimension
        local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        
        local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x+0.5*Delta;
      #endif

      int list_index = 0;
      for(scalar s in list_scalar_data) {
        local_data[list_index][i_cell] = (float) s[];
        SwapArrayBytes(&local_data[list_index][i_cell], sizeof(float));
        list_index++;
      }

      for(vector v in list_vector_data) {
        local_data[list_index][i_cell] = (float) v.x[];
        local_data[list_index + 1][i_cell] = (float) v.y[];
        SwapArrayBytes(&local_data[list_index][i_cell], sizeof(float));
        SwapArrayBytes(&local_data[list_index + 1][i_cell], sizeof(float));
        list_index += 2;
      }

      i_cell++;
    }
  }

  // === Getting how many cells we have in total between all processes
  // === And then gathering all the local vertices into the root process
  int total_num_cells;
  float *vertices = NULL;
  floatpp data = NULL;
  #if _MPI
    // === Total number of cells
    MPI_Allreduce(&local_num_cells, &total_num_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // === Gathering all the local vertices and the scalar fields data into the root process
    if( pid()==0 ) {
      vertices = (float *)malloc( 3*total_num_cells*vertices_per_cell*sizeof(float) );
      data = (float **)malloc( (number_of_scalar_fields + 2*number_of_vector_fields)*sizeof(float*) );
      for(int k=0; k<(number_of_scalar_fields + 2*number_of_vector_fields); k++)
        data[k] = (float *)malloc( total_num_cells*sizeof(float) );
    }
    MPI_Gather_Uneven(local_vertices, 3*vertices_per_cell*local_num_cells, MPI_FLOAT, vertices, 0);
    int list_index = 0;
    for( list_index=0; list_index<(number_of_scalar_fields + 2*number_of_vector_fields); list_index++ )
      MPI_Gather_Uneven(local_data[list_index], local_num_cells, MPI_FLOAT, (pid()==0) ? data[list_index] : NULL, 0);

    // === Releasing local memory
    free(local_vertices);
    for(int k=0; k<(number_of_scalar_fields + 2*number_of_vector_fields); k++)
      free(local_data[k]);
    if( local_data )
      free(local_data);
  #else
    total_num_cells = local_num_cells;
    vertices = local_vertices;
    data = local_data;
  #endif

  
  
  // === From now on, we only do file-printing stuff. Only the root process will do it
  if( pid()==0 ) {
    // === Opening the file to print the mesh
    sprintf(nomeArq, "%s/Mesh_%04d.vtk", folder_name, n);
    arq = fopen(nomeArq, "wt");

    // === Printing the VTK header information (printed as ASCII text)
    fprintf(arq, "# vtk DataFile Version 2.0\n");
    fprintf(arq, "MESH. step %d time %f\n", n, time);
    fprintf(arq, "BINARY\n");
    fprintf(arq, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(arq, "FIELD FieldData 1\n");
    fprintf(arq, "time 1 1 float\n");
    SwapArrayBytes(&time, sizeof(float));
    fwrite(&time, sizeof(float), 1, arq);
    fprintf(arq, "\n");

    // === Printing all the vertices coordinates (as BINARY)
    fprintf(arq, "POINTS %d float\n", total_num_cells*vertices_per_cell);
    SwapArrayBytes(vertices, 3*vertices_per_cell*total_num_cells*sizeof(float));
    fwrite(vertices, sizeof(float), 3*vertices_per_cell*total_num_cells, arq);
    fprintf(arq, "\n");

    // === Printing all the cells (each cell contains 4 (or 8) indices referring to the vertices printed above)
    fprintf(arq, "CELLS %d %d\n", total_num_cells, (vertices_per_cell + 1)*total_num_cells);
    int *array_cell_indices = malloc( (vertices_per_cell + 1)*total_num_cells*sizeof(int) );
    int offset = 0, vertex_index = 0;
    for( i=0; i<total_num_cells; i++ ) {
      array_cell_indices[offset] = vertex_index;
      array_cell_indices[offset + 1] = vertex_index + 1;
      array_cell_indices[offset + 2] = vertex_index + 2;
      array_cell_indices[offset + 3] = vertex_index + 3;
      #if dimension==3
        array_cell_indices[offset + 4] = vertex_index + 4;
        array_cell_indices[offset + 5] = vertex_index + 5;
        array_cell_indices[offset + 6] = vertex_index + 6;
        array_cell_indices[offset + 7] = vertex_index + 7;
      #endif
      array_cell_indices[offset + vertices_per_cell] = vertices_per_cell;
      offset += vertices_per_cell + 1;
      vertex_index += vertices_per_cell;
    }
    SwapArrayBytes(array_cell_indices, (vertices_per_cell + 1)*total_num_cells*sizeof(int));
    fwrite(array_cell_indices, sizeof(int), (vertices_per_cell + 1)*total_num_cells, arq);
    fprintf(arq, "\n");
        

    fprintf(arq, "CELL_TYPES %d\n", total_num_cells);
    SwapArrayBytes(&vtk_cell_type, sizeof(int));
    for( i=0; i<total_num_cells; i++ )
      array_cell_indices[i] = vtk_cell_type;
    fwrite(array_cell_indices, sizeof(int), total_num_cells, arq);
    fprintf(arq, "\n");


    // === Printing the actual simulation data that is stored in the cells
    int list_index = 0;
    fprintf(arq, "CELL_DATA %d\n", total_num_cells);
    for( scalar s in list_scalar_data ) {
      fprintf(arq, "SCALARS %s float 1\n", list_scalar_names[list_index]);
      fprintf(arq, "LOOKUP_TABLE default\n");
      fwrite(data[list_index], sizeof(float), total_num_cells, arq);
      fprintf(arq, "\n");
      list_index++;
    }

    int index_vector = 0;
    for( vector v in list_vector_data ) {
      fprintf(arq, "VECTORS %s float\n", list_vector_names[index_vector]);
      // fprintf(arq, "LOOKUP_TABLE default\n");

      float *array_vector_data = (float *)malloc( 3*total_num_cells*sizeof(float) );
      for(int i=0; i<total_num_cells; i++) {
        array_vector_data[3*i] = data[list_index][i];
        array_vector_data[3*i + 1] = data[list_index+1][i];
        array_vector_data[3*i + 2] = 0.0;
      }
      fwrite(array_vector_data, sizeof(float), 3*total_num_cells, arq);
      free(array_vector_data);
      fprintf(arq, "\n");

      list_index += 2;
      index_vector++;
    }

    free(array_cell_indices);
    free(vertices);
    for( list_index=0; list_index<(number_of_scalar_fields + 2*number_of_vector_fields); list_index++ )
      free(data[list_index]);
    free(data);
    fclose(arq);
  }

  return;
}


/** 
  I blatantly copied this function from here: [draw.h](http://basilisk.fr/src/draw.h) */
static bool cfilter (Point point, scalar c, double cmin)
{
  double cmin1 = 4.*cmin;
  if (c[] <= cmin) {
    foreach_dimension()
      if (c[1] >= 1. - cmin1 || c[-1] >= 1. - cmin1)
	return true;
    return false;
  }
  if (c[] >= 1. - cmin) {
    foreach_dimension()
      if (c[1] <= cmin1 || c[-1] <= cmin1)
	return true;
    return false;
  }
  int n = 0;
  double min = HUGE, max = - HUGE;
  foreach_neighbor(1) {
    if (c[] > cmin && c[] < 1. - cmin && ++n >= (1 << dimension))
      return true;
    if (c[] > max) max = c[];
    if (c[] < min) min = c[];
  }
  return max - min > 0.5;
}

/** 
  Prints a VTK file with the current interface.
  NOTE: Currently, I am printing the interface in ASCII format.
  I will make the Binary version soon. */
void PrintInterfaceVTK(int n, float time)
{
  double fmin = 1e-3; // do not reconstruct fragments smaller than this

  // === Counting how many line segments we will draw and how many vertices
  int local_count_polys = 0;
  int local_count_vertices = 0;
  for( scalar s in interfaces ) {
    foreach(serial) {
      if( cfilter (point, s, fmin) ) {
        local_count_polys++;
        coord n = interface_normal (point, s);
        double alpha = plane_alpha (s[], n);
        #if dimension==2
          coord v[2];
          int m = facets (n, alpha, v);
        #else
          coord v[12];
          int m = facets (n, alpha, v, 1.1);
        #endif
        local_count_vertices += m;
      }
    }
  }

  // === Alocating memory for the local vertices
  int *local_poly_count_vertices = (int *)malloc( local_count_polys*sizeof(int) );
  double *local_vertices_x = (double *)malloc( local_count_vertices*sizeof(double) );
  double *local_vertices_y = (double *)malloc( local_count_vertices*sizeof(double) );
  double *local_vertices_z = (double *)malloc( local_count_vertices*sizeof(double) );
  
  double *local_normals_x = (double *)malloc( local_count_polys*sizeof(double) );
  double *local_normals_y = (double *)malloc( local_count_polys*sizeof(double) );
  double *local_normals_z = (double *)malloc( local_count_polys*sizeof(double) );

  // === Calculating all polygons again (bad) and saving the vertices to the local arrays
  int index_polygon = 0, index_vertex = 0, index_normal = 0;
  for( scalar s in interfaces ) {
    foreach(serial) {
      if( cfilter (point, s, fmin) ) {
        coord n = interface_normal (point, s);
        double alpha = plane_alpha (s[], n);

        double n_magnitude = sqrt(sq(n.x) + sq(n.y));
        local_normals_x[index_normal] = n.x/n_magnitude;
        local_normals_y[index_normal] = n.y/n_magnitude;
        local_normals_z[index_normal++] = 0.0;

        #if dimension==2
          coord v[2];
          int m = facets (n, alpha, v);
        #else
          coord v[12];
          int m = facets (n, alpha, v, 1.1);
        #endif
        local_poly_count_vertices[index_polygon++] = m;
        for( int i=0; i<m; i++ ) {
          local_vertices_x[index_vertex] = x + v[i].x*Delta;
          local_vertices_y[index_vertex] = y + v[i].y*Delta;
          #if dimension==2
            local_vertices_z[index_vertex] = 0.0;
          #else
            local_vertices_z[index_vertex] = z + v[i].z*Delta;
          #endif
          index_vertex++;
        }
      }
    }
  }

  // === Counting how many polys and vertices we have in total between all processes
  int total_count_vertices = 0, total_count_polys = 0;
  #if _MPI
    MPI_Allreduce(&local_count_vertices, &total_count_vertices, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_count_polys, &total_count_polys, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  #else
    total_count_vertices = local_count_vertices;
    total_count_polys = local_count_polys;
  #endif

  // === Sending all the vertices to the global vertex array in process id=0
  int *poly_count_vertices = NULL;
  double *vertices_x = NULL, *vertices_y = NULL, *vertices_z = NULL;
  double *normals_x = NULL, *normals_y = NULL, *normals_z = NULL;
  #if _MPI
    if( pid()==0 ) {
      poly_count_vertices = (int *)malloc( total_count_polys*sizeof(int) );
      vertices_x = (double *)malloc( total_count_vertices*sizeof(double) );
      vertices_y = (double *)malloc( total_count_vertices*sizeof(double) );
      vertices_z = (double *)malloc( total_count_vertices*sizeof(double) );

      normals_x = (double *)malloc( total_count_polys*sizeof(double) );
      normals_y = (double *)malloc( total_count_polys*sizeof(double) );
      normals_z = (double *)malloc( total_count_polys*sizeof(double) );
    }

    MPI_Gather_Uneven(local_poly_count_vertices, local_count_polys, MPI_INT, poly_count_vertices, 0);
    MPI_Gather_Uneven(local_vertices_x, local_count_vertices, MPI_DOUBLE, vertices_x, 0);
    MPI_Gather_Uneven(local_vertices_y, local_count_vertices, MPI_DOUBLE, vertices_y, 0);
    MPI_Gather_Uneven(local_vertices_z, local_count_vertices, MPI_DOUBLE, vertices_z, 0);

    MPI_Gather_Uneven(local_normals_x, local_count_polys, MPI_DOUBLE, normals_x, 0);
    MPI_Gather_Uneven(local_normals_y, local_count_polys, MPI_DOUBLE, normals_y, 0);
    MPI_Gather_Uneven(local_normals_z, local_count_polys, MPI_DOUBLE, normals_z, 0);

    // === Releasing local memory
    free(local_poly_count_vertices);
    free(local_vertices_x);
    free(local_vertices_y);
    free(local_vertices_z);
    free(local_normals_x);
    free(local_normals_y);
    free(local_normals_z);
  #else
    poly_count_vertices = local_poly_count_vertices;
    vertices_x = local_vertices_x;
    vertices_y = local_vertices_y;
    vertices_z = local_vertices_z;
    normals_x = local_normals_x;
    normals_y = local_normals_y;
    normals_z = local_normals_z;
  #endif

  // === From this point we will only do file-writing stuff. Only process ZERO will do it
  if( pid()!=0 )
    return;

  // === Opening the vtk file
  char nomeArq[900];
  sprintf(nomeArq, "%s/Interface_%04d.vtk", folder_name, n);
  FILE *arq = fopen(nomeArq, "wt");
  if( !arq ) {
    printf("\n\n PrintInterfaceVTK_2D: Problem opening file... \n\n");
    exit(0);
  }

  // === Writing the VTK header
  fprintf(arq, "# vtk DataFile Version 2.0\n");
  fprintf(arq, "INTERFACE. step %d time %lf\n", n, time);
  fprintf(arq, "ASCII\n");
  fprintf(arq, "DATASET POLYDATA\n");

  // fprintf(arq, "FIELD FieldData 1\n");
  // fprintf(arq, "time 1 1 float\n");
  // SwapArrayBytes(&time, sizeof(float));
  // fwrite(&time, sizeof(float), 1, arq);
  // fprintf(arq, "\n");

  fprintf(arq, "POINTS %d float\n", total_count_vertices);

  
  // === Writing all the surface vertices
  for( index_vertex=0; index_vertex<total_count_vertices; index_vertex++ )
    fprintf(arq, "%e %e %e\n", vertices_x[index_vertex], vertices_y[index_vertex], vertices_z[index_vertex]);
  

  // === Writing the polygons conectivity
  int global_index_vertex = 0;
  #if dimension==2
    fprintf(arq, "LINES %d %d\n", total_count_polys, total_count_polys + total_count_vertices);
  #else
    fprintf(arq, "POLYGONS %d %d\n", total_count_polys, total_count_polys + total_count_vertices);
  #endif
  for( int index_polygon=0; index_polygon<total_count_polys; index_polygon++ ) {
    fprintf(arq, "%d", poly_count_vertices[index_polygon]);
    int index_vertex;
    for( index_vertex=0; index_vertex<poly_count_vertices[index_polygon]; index_vertex++ )
      fprintf(arq, " %d", global_index_vertex++);
    fprintf(arq, "\n");
  }

  fprintf(arq, "CELL_DATA %d\n", total_count_polys);
  fprintf(arq, "NORMALS normals float\n");
  for( int index_polygon=0; index_polygon<total_count_polys; index_polygon++ ) {
    fprintf(arq, "%lf %lf %lf\n", normals_x[index_polygon], normals_y[index_polygon], normals_z[index_polygon]);
  }
  
  // === Releasing memory
  free(poly_count_vertices);
  free(vertices_x);
  free(vertices_y);
  free(vertices_z);
  free(normals_x);
  free(normals_y);
  free(normals_z);
  fclose(arq);
  return;
}

#ifdef INCLUDED_VOF_TRACERS
void PrintTracerParticlesVTK(int n, float time, Particles Pin)
{
  // === Counting the number of local particles
  int local_number_particles = 0;
  foreach_particle_in(Pin)
    local_number_particles++;

  // === Allocating memory for the local particles
  double *local_particles_x = (double *)malloc( local_number_particles*sizeof(double) );
  double *local_particles_y = (double *)malloc( local_number_particles*sizeof(double) );

  int i_particle = 0;
  foreach_particle_in(Pin) {
    local_particles_x[i_particle] = x;
    local_particles_y[i_particle] = y;
    i_particle++;
  }
  
  // === Counting how many particles we have in total between all processes
  int total_number_particles = 0;
  #if _MPI
    MPI_Allreduce(&local_number_particles, &total_number_particles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  #else
    total_number_particles = local_number_particles;
  #endif

  // === Sending all the particle positions to the global particle array in process id=0
  double *particles_x = NULL, *particles_y = NULL;
  #if _MPI
    if( pid()==0 ) {
      particles_x = (double *)malloc( total_number_particles*sizeof(double) );
      particles_y = (double *)malloc( total_number_particles*sizeof(double) );
    }

    MPI_Gather_Uneven(local_particles_x, local_number_particles, MPI_DOUBLE, particles_x, 0);
    MPI_Gather_Uneven(local_particles_y, local_number_particles, MPI_DOUBLE, particles_y, 0);

    // === Releasing local memory
    free(local_particles_x);
    free(local_particles_y);
  #else
    particles_x = local_particles_x;
    particles_y = local_particles_y;
  #endif

  // === From this point we will only do file-writing stuff. Only process ZERO will do it
  if( pid()!=0 )
    return;

  // === Opening the vtk file
  char nomeArq[900];
  sprintf(nomeArq, "%s/TracerParticles_%04d.vtk", folder_name, n);
  FILE *arq = fopen(nomeArq, "wt");
  if( !arq ) {
    printf("\n\n PrintInterfaceVTK_2D: Problem opening file... \n\n");
    exit(0);
  }

  // === Writing the VTK header
  fprintf(arq, "# vtk DataFile Version 2.0\n");
  fprintf(arq, "INTERFACE. step %d time %lf\n", n, time);
  fprintf(arq, "ASCII\n");
  fprintf(arq, "DATASET POLYDATA\n");

  fprintf(arq, "POINTS %d float\n", total_number_particles);

  
  // === Writing all the surface vertices
  for( int i=0; i<total_number_particles; i++ )
    fprintf(arq, "%e %e %e\n", particles_x[i], particles_y[i], 0.0);
  
  fprintf(arq, "POLYGONS 1 %d\n", total_number_particles + 1);
  fprintf(arq, "%d", total_number_particles);
  for( int i=0; i<total_number_particles; i++ )
    fprintf(arq, " %d", i);

  
  // === Releasing memory
  free(particles_x);
  free(particles_y);
  fclose(arq);
  return;
}
#endif

#define PRINT_DATA_DUMP
#ifdef PRINT_DATA_DUMP
void PrintUniformMeshDataDump(int n, double time, int nx, int ny, 
                              scalar *list_scalar_data, const char **list_scalar_names,
                              double clip_xMin = nodata, double clip_xMax = nodata, 
                              double clip_yMin = nodata, double clip_yMax = nodata,
                              bool inside_only = false)
{
  int i, j;
  static const float bad_data = 1e+10;

  /// === Sometimes the interpolate function that we use further below gets stuck without updatting boundaries
  boundary(list_scalar_data);
  
  /// === Find the limits of the printable domain box
  double xMin = 1e+10, xMax = -1e+10, yMin = 1e+10, yMax = -1e+10;
  foreach(reduction(min:xMin) reduction(min:yMin) reduction(max:xMax) reduction(max:yMax)) {
    if( (!inside_only) || ((f1[] + f2[])>0.01) ) {
      if( (x - 0.5*Delta) < xMin )
        xMin = x - 0.5*Delta;
      if( (y - 0.5*Delta) < yMin )
        yMin = y - 0.5*Delta;
      if( (x + 0.5*Delta) > xMax )
        xMax = x + 0.5*Delta;
      if( (y + 0.5*Delta) > yMax )
        yMax = y + 0.5*Delta;
    }
  }

  /// === Changing the printable box limits if the user requested it
  xMin = (clip_xMin!=nodata && xMin<clip_xMin) ? clip_xMin : xMin;
  xMax = (clip_xMax!=nodata && xMax>clip_xMax) ? clip_xMax : xMax;
  yMin = (clip_yMin!=nodata && yMin<clip_yMin) ? clip_yMin : yMin;
  yMax = (clip_yMax!=nodata && yMax>clip_yMax) ? clip_yMax : yMax;

  int num_cells = nx*ny;
  double dx = (xMax - xMin)/nx;
  double dy = (yMax - yMin)/ny;

  // Counting how many scalars we want to print
  int number_of_scalar_fields = list_len(list_scalar_data);

  // Allocating memory for all the scalar fields (plus x and y arrays)
  float **values = (float **)malloc((number_of_scalar_fields+2)*sizeof(float *));
  for(i=0; i<number_of_scalar_fields+2; i++)
    values[i] = (float *)malloc(num_cells*sizeof(float));

  // printf("bb %d\n", pid());

  int cell_index = 0;
  for(i=0; i<nx; i++) {
    for(j=0; j<ny; j++) {
      float x_grid = clamp(xMin + (i+0.5)*dx, xMin, xMax);
      float y_grid = clamp(yMin + (j+0.5)*dy, yMin, yMax);

      values[0][cell_index] = x_grid;
      values[1][cell_index] = y_grid;

      float vol_fraction = inside_only ? clamp(interpolate(f1, x_grid, y_grid) + interpolate(f2, x_grid, y_grid), 0.0, 1.0) : 0.5;

      /// === Value of the velocity at the center of the channel
      int scalar_index = 2;
      for(scalar s in list_scalar_data)
        // values[scalar_index++][cell_index] = vol_fraction;
        values[scalar_index++][cell_index] = vol_fraction>0.01 && vol_fraction<1.5  ? interpolate(s, x_grid, y_grid) : bad_data;
      
      cell_index++;
    }
  }

  // PrintLog("bc kkk %d\n", pid());

  /// Allocating memory for the arrays coming from the remote processes into the root
  #if _MPI
    float *remote_values[npe()];
    if( pid()==0 ) {
      for(int id_proc=1; id_proc<npe(); id_proc++)
        remote_values[id_proc] = (float *)malloc(num_cells*sizeof(float));
    }
  
    int scalar_index = 2;
    for(scalar s in list_scalar_data) {
      // printf("a %d kkk %d\n", scalar_index, pid());
      if( pid() )
        MPI_Send(values[scalar_index], num_cells, MPI_FLOAT, 0, 1000 + pid(), MPI_COMM_WORLD);
      else {
        for(int id_proc=1; id_proc<npe(); id_proc++)
          MPI_Recv(remote_values[id_proc], num_cells, MPI_FLOAT, id_proc, 1000 + id_proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          
        cell_index = 0;
        for(i=0; i<nx; i++) {
          for(j=0; j<ny; j++) {
            for(int id_proc=1; id_proc<npe() && values[scalar_index][cell_index]==bad_data; id_proc++ )
              values[scalar_index][cell_index] = remote_values[id_proc][cell_index];
            cell_index++;
          }
        }
      }

      // printf("b %d kkk %d\n", scalar_index, pid());
      scalar_index++;
    }
  
    if( pid()==0 ) {
      for(int id_proc=1; id_proc<npe(); id_proc++)
        free(remote_values[id_proc]);
    }
  #endif

  
  if( pid()==0 ) {
    char fileName[900];
    sprintf(fileName, "%s/MeshDump_%04d.bin", folder_name, n);
    FILE *file_out = (FILE *)fopen(fileName, "w");

    fprintf(file_out, "time %lf\n", time);

    // Printing a header with name of scalars and the dimensiona of the mesh (nx, ny)
    int scalar_index = 0;
    fprintf(file_out, "%d %d\n", nx, ny);
    fprintf(file_out, "x;y");
    for(scalar s in list_scalar_data)
      fprintf(file_out, ";%s", list_scalar_names[scalar_index++]);
    fprintf(file_out, "\n");

    fwrite(values[0], sizeof(float), num_cells, file_out);
    fwrite(values[1], sizeof(float), num_cells, file_out);
    scalar_index = 2;
    for(scalar s in list_scalar_data)
      fwrite(values[scalar_index++], sizeof(float), num_cells, file_out);
    fclose(file_out);
  }

  // Releasing allocated memory
  for(i=0; i<number_of_scalar_fields+2; i++)
    free(values[i]);
  free(values);

  // exit(0);
  return;
}
#endif

// void CalculateDistanceFunction(scalar d, bool signed_distance = false)
// {
//   // === Alocating memory for the local vertices
//   double *local_vertices_x = (double *)malloc( 5000*sizeof(double) );
//   double *local_vertices_y = (double *)malloc( 5000*sizeof(double) );

//   // === Calculating all polygons again (bad) and saving the vertices to the local arrays
//   int index_vertex = 0;
//   foreach(serial) {
//     if( cfilter (point, f, 1e-05) ) {
//       coord n = interface_normal (point, f);
//       double alpha = plane_alpha (f[], n);
//       coord v[2];
//       int m = facets (n, alpha, v);

//       local_vertices_x[index_vertex] = local_vertices_y[index_vertex] = 0.0;
//       for( int i=0; i<m; i++ ) {
//         local_vertices_x[index_vertex] += x + v[i].x*Delta;
//         local_vertices_y[index_vertex] += y + v[i].y*Delta;
//       }
//       local_vertices_x[index_vertex] /= (double)m;
//       local_vertices_y[index_vertex] /= (double)m;

//       // Ignore vertices too close to the solid surface
//       // if( local_vertices_x[index_vertex]>0.001 )
//         index_vertex++;
//     }
//   }

//   // === Counting how many polys and vertices we have in total between all processes
//   int local_count_vertices = index_vertex;
//   int total_count_vertices = 0;
//   #if _MPI
//     MPI_Allreduce(&local_count_vertices, &total_count_vertices, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//   #else
//     total_count_vertices = local_count_vertices;
//   #endif

//   // === Sending all the vertices to the global vertex array in process id=0
//   double *vertices_x = NULL, *vertices_y = NULL;
//   #if _MPI
//     vertices_x = (double *)malloc( total_count_vertices*sizeof(double) );
//     vertices_y = (double *)malloc( total_count_vertices*sizeof(double) );

//     MPI_Gather_Uneven(local_vertices_x, local_count_vertices, MPI_DOUBLE, vertices_x, 0);
//     MPI_Gather_Uneven(local_vertices_y, local_count_vertices, MPI_DOUBLE, vertices_y, 0);

//     // Sending the total array to other processors as well
//     if( pid()==0 ) {
//       for(int i=1; i<npe(); i++) {
//         MPI_Send(vertices_x, total_count_vertices, MPI_DOUBLE, i, 3000 + i, MPI_COMM_WORLD);
//         MPI_Send(vertices_y, total_count_vertices, MPI_DOUBLE, i, 4000 + i, MPI_COMM_WORLD);
//       }
//     }
//     else {
//       MPI_Recv(vertices_x, total_count_vertices, MPI_DOUBLE, 0, 3000 + pid(), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//       MPI_Recv(vertices_y, total_count_vertices, MPI_DOUBLE, 0, 4000 + pid(), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//     }

//     // === Releasing local memory
//     free(local_vertices_x);
//     free(local_vertices_y);
//   #else
//     vertices_x = local_vertices_x;
//     vertices_y = local_vertices_y;
//   #endif

//   foreach() {
//     d[] = 1e+10;
//     for( int i=0; i<total_count_vertices; i++ ) {
//       double local_distance = (x - vertices_x[i])*(x - vertices_x[i]) + (y - vertices_y[i])*(y - vertices_y[i]);
//       local_distance = sqrt(local_distance);
//       if( local_distance<d[] )
//         d[] = local_distance;
//     }
//     d[] *= (f[]<1e-05) && (signed_distance) ? -1.0 : 1.0;
//   }
// }

// double calculate_skewness(bool automatic_center = true)
// {
//   // We begin by calculating the center of mass (xcm, ycm)
//   double wt = 0., xcm = 0., ycm = 0.;
//   foreach( reduction(+:wt), reduction(+:xcm), reduction(+:ycm) ) {
//     wt += f[]*sq(Delta);
//     xcm += f[]*x*sq(Delta);
//     ycm += f[]*y*sq(Delta);
//   }
//   xcm /= wt;
//   ycm /= wt;

//   // Use the origin as center
//   if( automatic_center==false )
//     xcm = ycm = 0.0;

//   // Calculating the skewness properties
//   double x_squared = 0.0, x_cubed = 0.0;
//   foreach( reduction(+:x_squared), reduction(+:x_cubed) ) {
//     x_squared += f[]*sq(Delta)*sq(x - xcm)/wt;
//     x_cubed += f[]*sq(Delta)*cube(x - xcm)/wt;
//   }

//   double coskew_xxx = x_cubed/pow(x_squared, 1.5);
//   return coskew_xxx;
// }

// void calculate_energy_budget(double *kinetic, double *surface, double *dissipated,
//                             double factor_kinetic, double factor_surface, double factor_dissipated)
// {
//   double kinetic_temp = 0.0;
//   double surface_temp = 0.0;
//   double dissipated_temp = 0.0;

//   foreach(reduction(+:kinetic_temp) reduction(+:surface_temp) reduction(+:dissipated_temp)) {

//     // ==== Calculating kinetic energy
//     double local_density = clamp(f[], 0, 1)*(rho1 - rho2) + rho2;
//     kinetic_temp += local_density*0.5*(sq(u.x[]) + sq(u.y[]))*sq(Delta);

//     // ==== Calculating surface energy
//     // face vector s; // need s for normal construction (see below)
//     // s.x.i = -1; // just book keeping -- ignore
//     double fmin = 1e-03;
//     if( cfilter (point, f, fmin) ) {
//       coord n = interface_normal (point, f);
//       double alpha = plane_alpha (f[], n);
//       coord v[2];
//       int m = facets (n, alpha, v);
      
//       if( m==2 ) {
//         double p1x = x + v[0].x*Delta;
//         double p1y = y + v[0].y*Delta;
//         double p2x = x + v[1].x*Delta;
//         double p2y = y + v[1].y*Delta;

//         surface_temp += sqrt( sq(p1x - p2x) + sq(p1y - p2y) );
//       }
//     }

//     // ==== Calculating viscously dissipated energy
//     // Note: this is only the dissipated energy at the current timestep
//     // Note: If you want the total dissipation over the simulation, you have to integrate this over time
//     double dudx = (u.x[1, 0] - u.x[-1, 0])/(2.0*Delta);
//     double dudy = (u.x[0, 1] - u.x[0, -1])/(2.0*Delta);
//     double dvdx = (u.y[1, 0] - u.y[-1, 0])/(2.0*Delta);
//     double dvdy = (u.y[0, 1] - u.y[0, -1])/(2.0*Delta);

//     double viscosity_ratio = 1e-02;
//     double local_viscosity = clamp(f[], 0, 1)*(1.0 - viscosity_ratio) + viscosity_ratio;
//     dissipated_temp += local_viscosity*sq(Delta)*(2.0*dudx*dudx + 2.0*dvdy*dvdy + (dvdx  + dudy)*(dvdx  + dudy));
//   }

//   *kinetic = factor_kinetic*kinetic_temp;
//   *surface = factor_surface*surface_temp;
//   *dissipated += dt*factor_dissipated*dissipated_temp;

//   return;
// }
