#ifndef VTK_OUTPUT_TYPE
  #define VTK_OUTPUT_TYPE float
  #define MPI_VTK_OUTPUT_TYPE MPI_FLOAT
#endif
#define VTK_STRING_TYPE_HELPER(x) #x
#define VTK_STRING_TYPE(x) VTK_STRING_TYPE_HELPER(x)

#define ROUND_TO_INT(x) ( (int)(x + 0.5) )

#include "tag.h"

char *folder_name = NULL;
FILE *log_file = NULL;
FILE *error_file = NULL;

int restart_index = 0; // How many times this simulation has been restarted using dump/restore

double *list_printed_vtk_times; 

// Hiding the double pointer with a typedef because qcc gets really annoyed if I dont do this (why??)
typedef VTK_OUTPUT_TYPE** vtk_output_type_pp; 

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

  /// === Opening the new log file
  if( log_file ) 
    fclose(log_file);
  sprintf(command, "%s/log_file_0.txt", folder_name);
  log_file = fopen(command, attempt_restart ? "at" : "wt");
  if( !log_file ) {
    printf("  ===== ERROR in function OpenSimulationFolder ...\n");
    printf("  ===== Unable to create the folder and log_file. Do you have permission for creating this folder?\n");
    printf("  ===== Ending execution...\n\n");
    exit(1);
  }

  /// === Opening the new error file
  if( error_file ) 
    fclose(error_file);
  sprintf(command, "%s/error_file_0.txt", folder_name);
  error_file = fopen(command, attempt_restart ? "at" : "wt");
  if( !error_file ) {
    printf("  ===== ERROR in function OpenSimulationFolder ...\n");
    printf("  ===== Unable to create the folder and error_file. Do you have permission for creating this folder?\n");
    printf("  ===== Ending execution...\n\n");
    exit(1);
  }

  // Changing the basilisk standard error output "ferr" to be our custom file
  ferr = error_file;

  list_printed_vtk_times = (double *)malloc( 20000*sizeof(double) );
  list_printed_vtk_times[0] = - 1.0;
  
  free(folder_temp);
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

void CreateBasiliskDumpSnapshot(const char *output_name = NULL) {
  char file_name[900];
  if( output_name )
    sprintf(file_name, "%s/%s", folder_name, output_name);  
  else
    sprintf(file_name, "%s/dump_snapshot_%.5lf.bsk", folder_name, t);
    
  dump(file_name);
  return;
}

double RestoreBasiliskCheckpoint(char *input_name = NULL, double snapshot_time = -1.0) {
  char file_name[900];
  if( input_name )
    sprintf(file_name, "%s/%s", folder_name, input_name);
  else if( snapshot_time>=0 )
    sprintf(file_name, "%s/dump_snapshot_%.5lf.bsk", folder_name, snapshot_time);
  else {
    // TO DO: find the most recent snapshot automatically and load it
  }
  
  bool success_restore = restore (file = file_name);


  // Loading the time of the restore
  double restore_time = - 1.0;
  if( success_restore ) {
    FILE *fp = NULL;
    if (!fp && (fp = fopen (file_name, "r")) == NULL)
      return false;
    assert (fp);
    struct DumpHeader header = {0};
    if (fread (&header, sizeof(header), 1, fp) < 1) {
      fprintf (ferr, "restore(): error: expecting header\n");
      exit (1);
    }
    fclose(fp);
    restore_time = header.t;
  }


  if( success_restore && (pid()==0) ) { 
    
    fclose(log_file);
    fclose(error_file);

    restart_index = -1;
    do {
      restart_index++;
      sprintf(file_name, "%s/log_file_%d.txt", folder_name, restart_index);
      log_file = fopen(file_name, "rt");

      if( log_file!=NULL )
        fclose(log_file);
    } while( log_file!=NULL );

    sprintf(file_name, "%s/log_file_%d.txt", folder_name, restart_index);
    log_file = fopen(file_name, "wt");
    sprintf(file_name, "%s/error_file_%d.txt", folder_name, restart_index);
    error_file = fopen(file_name, "wt");

    ferr = error_file;
  }

  return restore_time;
}

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
    fflush(stdout);
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
Prints mesh and (optionally) scalars to a BINARY VTK file. */
void PrintMeshVTK(float time, bool print_only_droplet, scalar *list_scalar_data, const char **list_scalar_names)
{
  char name_file[900];
  double f_threshold = 1e-05;

  // ===  Cell is either a square (2D) or cube (3D), meaning either 4 or 8 vertices per cell
  int vertices_per_cell = (dimension==2) ? 4 : 8;

  // === VTK cell code that represents quads or voxels
  // === check Figure 2 here to understand: https://docs.vtk.org/en/latest/vtk_file_formats/vtk_legacy_file_format.html
  int vtk_cell_type = (dimension==2) ? 9 : 11;

  // === Counting how many local cells we have in the mesh
  // === NOTE: whenever I say "local", i mean things that are locally in this processor (in case of parallel simulation)
  int local_num_cells = 0;
  foreach(serial) {
    if( (!print_only_droplet) || (f[]>f_threshold) )
      local_num_cells++;
  }

  // === Allocatting memory for the local vertex arrays
  // === Note: the x,y,z coordinates will all go into the same array, 
  // === so the array will be [x1, y1, z1, x2, y2, z2, ...]
  VTK_OUTPUT_TYPE *local_vertices = (VTK_OUTPUT_TYPE *)malloc( 3*vertices_per_cell*local_num_cells*sizeof(VTK_OUTPUT_TYPE) );

  // === Allocating memory for ALL the local scalar data arrays
  int number_of_scalar_fields = list_len(list_scalar_data);
  vtk_output_type_pp local_data = number_of_scalar_fields ? (VTK_OUTPUT_TYPE **)malloc( number_of_scalar_fields*sizeof(VTK_OUTPUT_TYPE *) ) : NULL;
  for( int k=0; k<number_of_scalar_fields; k++ )
    local_data[k] = (VTK_OUTPUT_TYPE *)malloc( local_num_cells*sizeof(VTK_OUTPUT_TYPE) );

  // === Storing all the vertices coordinates in the arrays
  int i = 0, i_cell = 0;
  foreach(serial) {
    if( (!print_only_droplet) || (f[]>f_threshold) ) {

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
        local_data[list_index][i_cell] = (VTK_OUTPUT_TYPE) s[];
        SwapArrayBytes(&local_data[list_index][i_cell], sizeof(VTK_OUTPUT_TYPE));
        list_index++;
      }

      i_cell++;
    }
  }

  // === Getting how many cells we have in total between all processes
  // === And then gathering all the local vertices into the root process
  int total_num_cells;
  VTK_OUTPUT_TYPE *vertices = NULL;
  vtk_output_type_pp data = NULL;
  #if _MPI
    // === Total number of cells
    MPI_Allreduce(&local_num_cells, &total_num_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // === Gathering all the local vertices and the scalar fields data into the root process
    if( pid()==0 ) {
      vertices = (VTK_OUTPUT_TYPE *)malloc( 3*total_num_cells*vertices_per_cell*sizeof(VTK_OUTPUT_TYPE) );
      data = (VTK_OUTPUT_TYPE **)malloc( number_of_scalar_fields*sizeof(VTK_OUTPUT_TYPE*) );
      for(int k=0; k<number_of_scalar_fields; k++)
        data[k] = (VTK_OUTPUT_TYPE *)malloc( total_num_cells*sizeof(VTK_OUTPUT_TYPE) );
    }
    MPI_Gather_Uneven(local_vertices, 3*vertices_per_cell*local_num_cells, MPI_VTK_OUTPUT_TYPE, vertices, 0);
    int list_index = 0;
    for( list_index=0; list_index<number_of_scalar_fields; list_index++ )
      MPI_Gather_Uneven(local_data[list_index], local_num_cells, MPI_VTK_OUTPUT_TYPE, (pid()==0) ? data[list_index] : NULL, 0);

    // === Releasing local memory
    free(local_vertices);
    for(int k=0; k<number_of_scalar_fields; k++)
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
    // sprintf(nomeArq, "%s/Mesh_%04d.vtk", folder_name, n);
    sprintf(name_file, "%s/Mesh_%.5lf.vtk", folder_name, time);
    FILE *file_pointer = fopen(name_file, "wt");

    // === Printing the VTK header information (printed as ASCII text)
    fprintf(file_pointer, "# vtk DataFile Version 2.0\n");
    fprintf(file_pointer, "MESH. time %f\n", time);
    fprintf(file_pointer, "BINARY\n");
    fprintf(file_pointer, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(file_pointer, "FIELD FieldData 1\n");
    fprintf(file_pointer, "time 1 1 %s\n", VTK_STRING_TYPE(VTK_OUTPUT_TYPE));
    SwapArrayBytes(&time, sizeof(VTK_OUTPUT_TYPE));
    fwrite(&time, sizeof(VTK_OUTPUT_TYPE), 1, file_pointer);
    fprintf(file_pointer, "\n");

    // === Printing all the vertices coordinates (as BINARY)
    fprintf(file_pointer, "POINTS %d %s\n", total_num_cells*vertices_per_cell, VTK_STRING_TYPE(VTK_OUTPUT_TYPE));
    SwapArrayBytes(vertices, 3*vertices_per_cell*total_num_cells*sizeof(VTK_OUTPUT_TYPE));
    fwrite(vertices, sizeof(VTK_OUTPUT_TYPE), 3*vertices_per_cell*total_num_cells, file_pointer);
    fprintf(file_pointer, "\n");

    // === Printing all the cells (each cell contains 4 (or 8) indices referring to the vertices printed above)
    fprintf(file_pointer, "CELLS %d %d\n", total_num_cells, (vertices_per_cell + 1)*total_num_cells);
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
    fwrite(array_cell_indices, sizeof(int), (vertices_per_cell + 1)*total_num_cells, file_pointer);
    fprintf(file_pointer, "\n");
        

    fprintf(file_pointer, "CELL_TYPES %d\n", total_num_cells);
    SwapArrayBytes(&vtk_cell_type, sizeof(int));
    for( i=0; i<total_num_cells; i++ )
      array_cell_indices[i] = vtk_cell_type;
    fwrite(array_cell_indices, sizeof(int), total_num_cells, file_pointer);
    fprintf(file_pointer, "\n");


    // === Printing the actual simulation data that is stored in the cells
    int list_index = 0;
    fprintf(file_pointer, "CELL_DATA %d\n", total_num_cells);
    for( scalar s in list_scalar_data ) {
      fprintf(file_pointer, "SCALARS %s %s 1\n", list_scalar_names[list_index], VTK_STRING_TYPE(VTK_OUTPUT_TYPE));
      fprintf(file_pointer, "LOOKUP_TABLE default\n");
      fwrite(data[list_index], sizeof(VTK_OUTPUT_TYPE), total_num_cells, file_pointer);
      fprintf(file_pointer, "\n");
      list_index++;
    }

    free(array_cell_indices);
    free(vertices);
    for( list_index=0; list_index<number_of_scalar_fields; list_index++ )
      free(data[list_index]);
    free(data);
    fclose(file_pointer);
  }

  return;
}


void PrintMesh2DCrossSection(float time, int chosen_axis, bool print_only_droplet,
                            scalar *list_scalar_data = NULL, const char **list_scalar_names = NULL)
{
  FILE *arq;
  char nomeArq[900];

  coord box[2] = {{-10.0*(chosen_axis!=0), -10.0*(chosen_axis!=1), -10.0*(chosen_axis!=2)}, 
                  {10.0*(chosen_axis!=0), 10.0*(chosen_axis!=1), 10.0*(chosen_axis!=2)}};
  coord region_sampling = {(chosen_axis==0) ? 1 : 10000, (chosen_axis==1) ? 1 : 10000, (chosen_axis==2) ? 1 : 10000};
  coord p;
  

  // ===  Cell is either a square (2D) or cube (3D), meaning either 4 or 8 vertices per cell
  int vertices_per_cell = 4;

  // === VTK cell code that represents voxels
  // === (check Figure 2 here to understand: https://kitware.github.io/vtk-examples/site/VTKFileFormats/)
  int vtk_cell_type = 9;

  // === Counting how many local cells we have in the mesh
  // === NOTE: whenever I say "local", i mean things that are locally in this processor (in case of parallel simulation)
  int local_num_cells = 0;
  coord previous_cell = {1e+10, 1e+10, 1e+10};
  foreach_region(p, box, region_sampling, serial) { 
    // Avoiding repeated cells (I guess they're ordered...?)
    if( previous_cell.x!=x || previous_cell.y!=y || previous_cell.z!=z ) {
      if( (!print_only_droplet) || (f[]>0.01) ) 
        local_num_cells++;
      previous_cell.x = x;
      previous_cell.y = y;
      previous_cell.z = z;
    }
  }


  // === Allocatting memory for the local vertex arrays
  // === Note: the x,y,z coordinates will all go into the same array, 
  // === so the array will be [x1, y1, z1, x2, y2, z2, ...]
  float *local_vertices = (float *)malloc( 3*vertices_per_cell*local_num_cells*sizeof(float) );

  // === Allocating memory for ALL the local scalar data arrays
  int number_of_scalar_fields = list_len(list_scalar_data);
  typedef float** floatpp; // qcc gets uncomfortable if i dont hide the double pointer (?????)
  floatpp local_data = number_of_scalar_fields ? (float **)malloc( number_of_scalar_fields*sizeof(float *) ) : NULL;
  for( int k=0; k<number_of_scalar_fields; k++ )
    local_data[k] = (float *)malloc( local_num_cells*sizeof(float) );

  // === Storing all the vertices coordinates in the arrays
  int i = 0, i_cell = 0;
  previous_cell.x = 1e+10;
  previous_cell.y = 1e+10;
  previous_cell.z = 1e+10;
  foreach_region(p, box, region_sampling, serial) { 
    // Avoiding repeated cells (I guess they're ordered...?)
    if( previous_cell.x!=x || previous_cell.y!=y || previous_cell.z!=z ) {
      previous_cell.x = x;
      previous_cell.y = y;
      previous_cell.z = z;
      
      if( (!print_only_droplet) || (f[]>0.01) ) {
        /// Using a macro conditional to avoid checking dimension==2 every loop during execution time...
        // local_vertices[i++] = 0.0; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        // local_vertices[i++] = 0.0; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        // local_vertices[i++] = 0.0; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x+0.5*Delta;
        // local_vertices[i++] = 0.0; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        
        if( chosen_axis==0 ) {
          local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = 0.0;
          local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = 0.0;
          local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = 0.0;
          local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = 0.0;
        }
        else if( chosen_axis==1 ) {
          local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = 0.0; local_vertices[i++] = x-0.5*Delta;
          local_vertices[i++] = z-0.5*Delta; local_vertices[i++] = 0.0; local_vertices[i++] = x+0.5*Delta;
          local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = 0.0; local_vertices[i++] = x+0.5*Delta;
          local_vertices[i++] = z+0.5*Delta; local_vertices[i++] = 0.0; local_vertices[i++] = x-0.5*Delta;
        }
        else if( chosen_axis==2 ) {
          local_vertices[i++] = 0.0; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x-0.5*Delta;
          local_vertices[i++] = 0.0; local_vertices[i++] = y-0.5*Delta; local_vertices[i++] = x+0.5*Delta;
          local_vertices[i++] = 0.0; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x+0.5*Delta;
          local_vertices[i++] = 0.0; local_vertices[i++] = y+0.5*Delta; local_vertices[i++] = x-0.5*Delta;
        }

        int list_index = 0;
        for(scalar s in list_scalar_data) {
          local_data[list_index][i_cell] = (float) s[];
          SwapArrayBytes(&local_data[list_index][i_cell], sizeof(float));
          list_index++;
        }

        i_cell++;
      }
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
      data = (float **)malloc( number_of_scalar_fields*sizeof(float*) );
      for(int k=0; k<number_of_scalar_fields; k++)
        data[k] = (float *)malloc( total_num_cells*sizeof(float) );
    }
    MPI_Gather_Uneven(local_vertices, 3*vertices_per_cell*local_num_cells, MPI_FLOAT, vertices, 0);
    int list_index = 0;
    for( list_index=0; list_index<number_of_scalar_fields; list_index++ )
      MPI_Gather_Uneven(local_data[list_index], local_num_cells, MPI_FLOAT, (pid()==0) ? data[list_index] : NULL, 0);

    // === Releasing local memory
    free(local_vertices);
    for(int k=0; k<number_of_scalar_fields; k++)
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
    // sprintf(nomeArq, "%s/Mesh_%04d.vtk", folder_name, n);
    sprintf(nomeArq, "%s/MeshSection_%.5lf.vtk", folder_name, time);
    arq = fopen(nomeArq, "wt");

    // === Printing the VTK header information (printed as ASCII text)
    fprintf(arq, "# vtk DataFile Version 2.0\n");
    fprintf(arq, "MESH. time %f\n", time);
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

    free(array_cell_indices);
    free(vertices);
    for( list_index=0; list_index<number_of_scalar_fields; list_index++ )
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

void PrintInterfaceVTK( VTK_OUTPUT_TYPE time, 
                        scalar *list_scalar_data = NULL, const char **list_scalar_names = NULL,
                        bool tag_droplets = false)
{
  double fmin = 1e-3; // do not reconstruct fragments smaller than this
  int number_of_scalar_fields = list_len(list_scalar_data);
  int list_index;

  scalar scalar_tags[];

  // === Counting how many line segments we will draw and how many vertices
  int local_count_polys = 0;
  int local_count_vertices = 0;
  foreach(serial) {

    // Preparing the threshold scalar for the tag function later
    scalar_tags[] = tag_droplets ? f[] > 1e-04 : 1.0;

    if( cfilter (point, f, fmin) ) {
      local_count_polys++;
      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);
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

  int number_tags = tag_droplets ? tag(scalar_tags) : 1;  

  // === Alocating memory for the local vertices
  int *local_poly_count_vertices = (int *)malloc( local_count_polys*sizeof(int) );
  int *local_poly_tags = (int *)malloc( local_count_polys*sizeof(int) );
  VTK_OUTPUT_TYPE *local_vertices = (VTK_OUTPUT_TYPE *)malloc( 3*local_count_vertices*sizeof(VTK_OUTPUT_TYPE) );
  VTK_OUTPUT_TYPE *local_normals = (VTK_OUTPUT_TYPE *)malloc( 3*local_count_polys*sizeof(VTK_OUTPUT_TYPE) );
  
  // === Allocating memory for the scalar data
  vtk_output_type_pp local_data = number_of_scalar_fields ? (vtk_output_type_pp)malloc( number_of_scalar_fields*sizeof(VTK_OUTPUT_TYPE *) ) : NULL;
  for( int k=0; k<number_of_scalar_fields; k++ )
    local_data[k] = (VTK_OUTPUT_TYPE *)malloc( local_count_polys*sizeof(VTK_OUTPUT_TYPE) );

  // === Calculating all polygons again (bad) and saving the vertices to the local arrays
  int index_polygon = 0, index_vertex_array = 0;//, index_normal = 0;
  int index_interface_cell = 0;
  foreach(serial) {
    if( cfilter (point, f, fmin) ) {
      local_poly_tags[index_polygon] = ROUND_TO_INT( scalar_tags[] );

      coord n = interface_normal (point, f);
      double alpha = plane_alpha (f[], n);
      

      // double n_magnitude = sqrt(sq(n.x) + sq(n.y));
      // local_normals_x[index_normal] = n.x/n_magnitude;
      // local_normals_y[index_normal] = n.y/n_magnitude;
      // local_normals_z[index_normal++] = 0.0;

      #if dimension==2
        coord v[2];
        int m = facets (n, alpha, v);
      #else
        coord v[12];
        int m = facets (n, alpha, v, 1.1);
      #endif

      int list_index = 0;
      for(scalar s in list_scalar_data)
        local_data[list_index++][index_interface_cell] = (VTK_OUTPUT_TYPE) s[];
      index_interface_cell++;

      local_poly_count_vertices[index_polygon++] = m;
      for( int i=0; i<m; i++ ) {
        #if dimension==3
          local_vertices[index_vertex_array++] = z + v[i].z*Delta;
        #else
          local_vertices[index_vertex_array++] = 0.0;
        #endif

        local_vertices[index_vertex_array++] = y + v[i].y*Delta;
        local_vertices[index_vertex_array++] = x + v[i].x*Delta;
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
  int *poly_tags = NULL;
  VTK_OUTPUT_TYPE *vertices = NULL;
  VTK_OUTPUT_TYPE *normals = NULL;
  vtk_output_type_pp data = NULL;
  #if _MPI
    if( pid()==0 ) {
      poly_count_vertices = (int *)malloc( total_count_polys*sizeof(int) );
      poly_tags = (int *)malloc( total_count_polys*sizeof(int) );
      vertices = (VTK_OUTPUT_TYPE *)malloc( 3*total_count_vertices*sizeof(VTK_OUTPUT_TYPE) );
      normals = (VTK_OUTPUT_TYPE *)malloc( 3*total_count_polys*sizeof(VTK_OUTPUT_TYPE) );
      data = (VTK_OUTPUT_TYPE **)malloc( number_of_scalar_fields*sizeof(VTK_OUTPUT_TYPE*) );
      for(int k=0; k<number_of_scalar_fields; k++)
        data[k] = (VTK_OUTPUT_TYPE *)malloc( total_count_polys*sizeof(VTK_OUTPUT_TYPE) );
    }

    MPI_Gather_Uneven(local_poly_count_vertices, local_count_polys, MPI_INT, poly_count_vertices, 0);
    MPI_Gather_Uneven(local_poly_tags, local_count_polys, MPI_INT, poly_tags, 0);
    MPI_Gather_Uneven(local_vertices, 3*local_count_vertices, MPI_VTK_OUTPUT_TYPE, vertices, 0);

    // MPI_Gather_Uneven(local_normals_x, local_count_polys, MPI_DOUBLE, normals_x, 0);
    // MPI_Gather_Uneven(local_normals_y, local_count_polys, MPI_DOUBLE, normals_y, 0);
    // MPI_Gather_Uneven(local_normals_z, local_count_polys, MPI_DOUBLE, normals_z, 0);

    list_index = 0;
    for( list_index=0; list_index<number_of_scalar_fields; list_index++ )
      MPI_Gather_Uneven(local_data[list_index], local_count_polys, MPI_VTK_OUTPUT_TYPE, (pid()==0) ? data[list_index] : NULL, 0);

    


    // === Releasing local memory
    free(local_poly_count_vertices);
    free(local_poly_tags);
    free(local_vertices);
    free(local_normals);
    for(int k=0; k<number_of_scalar_fields; k++)
      free(local_data[k]);
    if( local_data )
      free(local_data);
  #else
    poly_count_vertices = local_poly_count_vertices;
    poly_tags = local_poly_tags;
    vertices = local_vertices;
    normals = local_normals;
    data = local_data;
  #endif

  // === From this point we will only do file-writing stuff. Only process ZERO will do it
  if( pid()!=0 )
    return;

  // === Opening the vtk file
  char nomeArq[900];
  // sprintf(nomeArq, "%s/Interface_%04d.vtk", folder_name, n);
  sprintf(nomeArq, "%s/Interface_%.5lf.vtk", folder_name, time);
  FILE *arq = fopen(nomeArq, "wt");
  if( !arq ) {
    printf("\n\n PrintInterfaceVTK_2D: Problem opening file... \n\n");
    exit(0);
  }

  // === Writing the VTK header
  fprintf(arq, "# vtk DataFile Version 2.0\n");
  fprintf(arq, "INTERFACE. time %lf\n", time);
  fprintf(arq, "BINARY\n");
  fprintf(arq, "DATASET POLYDATA\n");

  // === Flipping the vertices and data arrays because of VTK conventions...
  for(int i=0; i<total_count_polys/2; i++) {
    int temp_swap = poly_count_vertices[i];
    poly_count_vertices[i] = poly_count_vertices[total_count_polys-i-1];
    poly_count_vertices[total_count_polys-i-1] = temp_swap;

    temp_swap = poly_tags[i];
    poly_tags[i] = poly_tags[total_count_polys-i-1];
    poly_tags[total_count_polys-i-1] = temp_swap;

    list_index = 0;
    for( scalar s in list_scalar_data ) {
      double temp_swap_double = data[list_index][i];
      data[list_index][i] = data[list_index][total_count_polys-i-1];
      data[list_index][total_count_polys-i-1] = temp_swap_double;
      list_index++;
    }
  }

  // === Writing all the surface vertices
  fprintf(arq, "POINTS %d %s\n", total_count_vertices, VTK_STRING_TYPE(VTK_OUTPUT_TYPE));
  SwapArrayBytes(vertices, 3*total_count_vertices*sizeof(VTK_OUTPUT_TYPE));
  fwrite(vertices, sizeof(VTK_OUTPUT_TYPE), 3*total_count_vertices, arq);
  fprintf(arq, "\n");
  

  // === Writing the polygons conectivity
  #if dimension==2
    fprintf(arq, "LINES %d %d\n", total_count_polys, total_count_polys + total_count_vertices);
  #else
    fprintf(arq, "POLYGONS %d %d\n", total_count_polys, total_count_polys + total_count_vertices);
  #endif

  // === Printing all the cells (each cell contains 4 (or 8) indices referring to the vertices printed above)
  int *array_cell_indices = malloc( (total_count_polys + total_count_vertices)*sizeof(int) );
  int offset = 0, global_index_vertex = 0;
  for( int index_polygon=0; index_polygon<total_count_polys; index_polygon++ ) {
    int count_vertices_this_poly = poly_count_vertices[index_polygon];

    for( int index_vertex=0; index_vertex<count_vertices_this_poly; index_vertex++ )
      array_cell_indices[offset++] = global_index_vertex++;

    array_cell_indices[offset++] = count_vertices_this_poly;
  }
  SwapArrayBytes(array_cell_indices, (total_count_polys + total_count_vertices)*sizeof(int));
  fwrite(array_cell_indices, sizeof(int), total_count_polys + total_count_vertices, arq);
  fprintf(arq, "\n");
  
  // fprintf(arq, "CELL_DATA %d\n", total_count_polys);
  // fprintf(arq, "NORMALS normals float\n");
  // for( int index_polygon=0; index_polygon<total_count_polys; index_polygon++ ) {
  //   fprintf(arq, "%lf %lf %lf\n", normals_x[index_polygon], normals_y[index_polygon], normals_z[index_polygon]);
  // }


  // === Printing the actual simulation data that is stored in the cells
  list_index = 0;
  fprintf(arq, "CELL_DATA %d\n", total_count_polys);
  for( scalar s in list_scalar_data ) {
    fprintf(arq, "SCALARS %s %s 1\n", list_scalar_names[list_index], VTK_STRING_TYPE(VTK_OUTPUT_TYPE));
    fprintf(arq, "LOOKUP_TABLE default\n");
    SwapArrayBytes(data[list_index], total_count_polys*sizeof(VTK_OUTPUT_TYPE));
    fwrite(data[list_index], sizeof(VTK_OUTPUT_TYPE), total_count_polys, arq);
    fprintf(arq, "\n");
    list_index++;
  }

  if( tag_droplets ) {
    fprintf(arq, "SCALARS tags int 1\n");
    fprintf(arq, "LOOKUP_TABLE default\n");
    SwapArrayBytes(poly_tags, total_count_polys*sizeof(int));
    fwrite(poly_tags, sizeof(int), total_count_polys, arq);
    fprintf(arq, "\n");
  }
  
  // === Releasing memory
  free(poly_count_vertices);
  free(poly_tags);
  free(vertices);
  free(normals);
  free(array_cell_indices);
  for( int list_index=0; list_index<number_of_scalar_fields; list_index++ )
    free(data[list_index]);
  free(data);
  fclose(arq);
  return;
}

#ifdef INCLUDED_VOF_TRACERS
void PrintTracerParticlesVTK(VTK_OUTPUT_TYPE time, Particles Pin)
{
  // === Counting the number of local particles
  int local_number_particles = 0;
  foreach_particle_in(Pin)
    local_number_particles++;

  // === Allocating memory for the local particles
  VTK_OUTPUT_TYPE *local_particles_x = (VTK_OUTPUT_TYPE *)malloc( local_number_particles*sizeof(VTK_OUTPUT_TYPE) );
  VTK_OUTPUT_TYPE *local_particles_y = (VTK_OUTPUT_TYPE *)malloc( local_number_particles*sizeof(VTK_OUTPUT_TYPE) );

  int i_particle = 0;
  foreach_particle_in(Pin) {
    local_particles_x[i_particle] = (VTK_OUTPUT_TYPE)x;
    local_particles_y[i_particle] = (VTK_OUTPUT_TYPE)y;
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
  VTK_OUTPUT_TYPE *particles_x = NULL, *particles_y = NULL;
  #if _MPI
    if( pid()==0 ) {
      particles_x = (VTK_OUTPUT_TYPE *)malloc( total_number_particles*sizeof(VTK_OUTPUT_TYPE) );
      particles_y = (VTK_OUTPUT_TYPE *)malloc( total_number_particles*sizeof(VTK_OUTPUT_TYPE) );
    }

    MPI_Gather_Uneven(local_particles_x, local_number_particles, MPI_VTK_OUTPUT_TYPE, particles_x, 0);
    MPI_Gather_Uneven(local_particles_y, local_number_particles, MPI_VTK_OUTPUT_TYPE, particles_y, 0);

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
  sprintf(nomeArq, "%s/TracerParticles_%.5lf.vtk", folder_name, time);
  FILE *arq = fopen(nomeArq, "wt");
  if( !arq ) {
    printf("\n\n TracerParticles: Problem opening file... \n\n");
    exit(0);
  }

  // === Writing the VTK header
  fprintf(arq, "# vtk DataFile Version 2.0\n");
  fprintf(arq, "TRACERS. time %f\n", (float)time);
  fprintf(arq, "ASCII\n");
  fprintf(arq, "DATASET POLYDATA\n");

  fprintf(arq, "POINTS %d %s\n", total_number_particles, VTK_STRING_TYPE(VTK_OUTPUT_TYPE));

  
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

void UpdateMeshOutputSeriesFile(float time)
{
  // Only process zero will perform file writing stuff
  if( pid() )
    return;

  char file_name[900];

  // Opening the interface file
  sprintf(file_name, "%s/Interface.vtk.series", folder_name);
  FILE *file_interface_series = fopen(file_name, "wt");
  fprintf(file_interface_series, "{\n  \"file-series-version\": \"1.0\",\n  \"files\": [\n");

  // Opening the mesh file
  sprintf(file_name, "%s/Mesh.vtk.series", folder_name);
  FILE *file_mesh_series = fopen(file_name, "wt");
  fprintf(file_mesh_series, "{\n  \"file-series-version\": \"1.0\",\n  \"files\": [\n");

  // Opening the TRACERS file
  sprintf(file_name, "%s/Tracers.vtk.series", folder_name);
  FILE *file_tracers_series = fopen(file_name, "wt");
  fprintf(file_tracers_series, "{\n  \"file-series-version\": \"1.0\",\n  \"files\": [\n");

  int i;
  for( i=0; (list_printed_vtk_times[i]!=-1.0) && (i<20000); i++ ) {
    fprintf(file_interface_series, "    {\"name\": \"Interface_%.5f.vtk\", \"time\": %.5f},\n", list_printed_vtk_times[i], list_printed_vtk_times[i]);
    fprintf(file_mesh_series, "    {\"name\": \"Mesh_%.5f.vtk\", \"time\": %.5f},\n", list_printed_vtk_times[i], list_printed_vtk_times[i]);
    fprintf(file_tracers_series, "    {\"name\": \"TracerParticles_%.5f.vtk\", \"time\": %.5f},\n", list_printed_vtk_times[i], list_printed_vtk_times[i]);
  }
  fprintf(file_interface_series, "    {\"name\": \"Interface_%.5f.vtk\", \"time\": %.5f},\n", time, time);
  fprintf(file_mesh_series, "    {\"name\": \"Mesh_%.5f.vtk\", \"time\": %.5f},\n", time, time);
  fprintf(file_tracers_series, "    {\"name\": \"TracerParticles_%.5f.vtk\", \"time\": %.5f},\n", time, time);
  
  list_printed_vtk_times[i] = time;
  list_printed_vtk_times[i + 1] = -1.0;

  // Closing the interface file
  fprintf(file_interface_series, "  ]\n}\n");
  fclose(file_interface_series);

  // Closing the mesh file
  fprintf(file_mesh_series, "  ]\n}\n");
  fclose(file_mesh_series);

  // Closing the Tracers file
  fprintf(file_tracers_series, "  ]\n}\n");
  fclose(file_tracers_series);
  
  return;
}

#define PRINT_DATA_DUMP
#ifdef PRINT_DATA_DUMP
void PrintUniformMeshData(int n, double time, int mesh_level, 
                              scalar *list_scalar_data, const char **list_scalar_names,
                              double clip_xMin = nodata, double clip_xMax = nodata, 
                              double clip_yMin = nodata, double clip_yMax = nodata,
                              bool inside_only = false)
{
  int i, j;
  static const float bad_data = 1e+10;

  /// === Sometimes the interpolate function that we use further below gets stuck without updatting boundaries
  boundary(list_scalar_data);

  // Finding the minimum dx and limits of the domain
  int base_nx = 1 << mesh_level;
  double dx = L0/base_nx;
  
  /// === Find the limits of the rectangular domain we want to print
  double xMin = 1e+10, xMax = -1e+10, yMin = 1e+10, yMax = -1e+10;
  foreach(reduction(min:xMin) reduction(min:yMin) reduction(max:xMax) reduction(max:yMax)) {
    if( (!inside_only) || f[]>0.01 ) {
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

  xMin = (clip_xMin!=nodata && xMin<clip_xMin) ? clip_xMin : xMin;
  xMax = (clip_xMax!=nodata && xMax>clip_xMax) ? clip_xMax : xMax;
  yMin = (clip_yMin!=nodata && yMin<clip_yMin) ? clip_yMin : yMin;
  yMax = (clip_yMax!=nodata && yMax>clip_yMax) ? clip_yMax : yMax;

  int nx = round( (xMax - xMin)/dx );
  int ny = round( (yMax - yMin)/dx );
  int num_cells = nx*ny;
  dx = (xMax - xMin)/nx;
  double dy = (yMax - yMin)/ny;

  // Counting how many scalars we want to print
  int number_of_scalar_fields = list_len(list_scalar_data);

  // printf("ba kkk %d\n", pid());

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

      float vol_fraction = inside_only ? interpolate(f, x_grid, y_grid) : 0.5;

      /// === Value of the velocity at the center of the channel
      int scalar_index = 2;
      for(scalar s in list_scalar_data)
        // values[scalar_index++][cell_index] = vol_fraction;
        values[scalar_index++][cell_index] = vol_fraction>0.01 && vol_fraction<1.5  ? interpolate(s, x_grid, y_grid) : bad_data;
      
      cell_index++;
    }
  }


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


void PrintUniformMeshDataDump(double time, int nx, int ny, 
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
    if( (!inside_only) || f[]>0.01 ) {
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

      float vol_fraction = inside_only ? interpolate(f, x_grid, y_grid) : 0.5;

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
    sprintf(fileName, "%s/MeshDump_%.5lf.bin", folder_name, t);
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

