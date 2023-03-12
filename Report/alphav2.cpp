// mpic++ alphav2.cpp -o run1
// mpirun -np 4 ./run1

#include <mpi.h>
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG_MODE 0
#define MASTER_PROCESSOR_RANK 0
#define DELTA_T 10
#define INPUT_FILE_DATA "inputfile.txt"




using namespace std;


double GRAVITATIONAL_CONSTANT = 6.6674 * pow(10,-11);
// FILE *fp4;



int rank_of_processor,num_of_processors;

MPI_Datatype aggregat_Type;

double * velocities = NULL;
long long int NUMBER_OF_PARTICLES = 2000;
long long int TOTAL_TIME_RUN = 3000;


// Function prototypes for calculating forces, updating position & velocity
void compute_force(double masses[], double positions[], double forces_each_proc[], int rank_of_processor, int BODIES_per_proc);
void update_positions_velocities(double positions[], double forces_each_proc[], double velocities_per_proc[], int rank_of_processor, int BODIES_per_proc);

// Main function.
int main(int argc, char * argv[]) {

  // Array for all positions and velocities & forces for particles belonging to current processors.
  double * masses;
  double * positions;
  double * velocities_per_proc;
  double * forces_each_proc;

  MPI_Init( & argc, & argv);
  MPI_Comm_size(MPI_COMM_WORLD, & num_of_processors);
  MPI_Comm_rank(MPI_COMM_WORLD, & rank_of_processor);

  // Get number of particles per processor (fancy calculation is to get the ceiling).
  int BODIES_per_proc = (NUMBER_OF_PARTICLES + num_of_processors - 1) / num_of_processors;

  if (rank_of_processor == MASTER_PROCESSOR_RANK) {
    velocities = (double *)malloc(2 * num_of_processors * BODIES_per_proc * sizeof(double));
  }

  masses = (double *)malloc(NUMBER_OF_PARTICLES * sizeof(double));
  positions = (double *)calloc(2 * num_of_processors * BODIES_per_proc, sizeof(double));
  velocities_per_proc = (double *)malloc(2 * BODIES_per_proc * sizeof(double));
  forces_each_proc = (double *)malloc(2 * BODIES_per_proc * sizeof(double));



  // This is to communicate positions and velocity of each chunk (so 4 total doubles for each particle in chunk).
  MPI_Type_contiguous(2 * BODIES_per_proc, MPI_DOUBLE, & aggregat_Type);
  MPI_Type_commit( & aggregat_Type);

  if (rank_of_processor == MASTER_PROCESSOR_RANK) {
    FILE * fp = fopen(INPUT_FILE_DATA, "r");
    if (!fp) {
      printf("Error opening input file.\n");
      exit(1);
    }
    long long int number_offile_list;
    int counter_masses = 0;
    int counter_positions = 0;
    int counter_velocities = 0;
    int line = 0;
      fscanf(fp, "%lld", &number_offile_list);
      cout<<"MAXIMUM NUMBER of objects"<<number_offile_list<<endl;
        FILE *fp4 = fopen("justcoordinates.txt","w+");
        
        fprintf(fp4,"%lld",NUMBER_OF_PARTICLES);

        fprintf(fp4,"\n");
        fprintf(fp4,"%lld",TOTAL_TIME_RUN);
        fprintf(fp4,"\n");

        fclose(fp4);

    for (line = 0; line < NUMBER_OF_PARTICLES; line++) {
      fscanf(fp, "%lf %lf %lf", & positions[counter_positions], & positions[counter_positions + 1], & masses[counter_masses] );
    //   printf("%lf %lf %lf\n ",  positions[counter_positions], positions[counter_positions + 1], masses[counter_masses]);

      velocities[counter_velocities] =0;
      velocities[counter_velocities + 1] =0;




      counter_masses += 1;
      counter_positions += 2;
      counter_velocities += 2;
    }




    if (DEBUG_MODE>=2)
            {
                    counter_masses = 0;
                    counter_positions = 0;
                    counter_velocities = 0;

                        for (line = 0; line < NUMBER_OF_PARTICLES; line++) {
            
                                // fscanf(fp, "%lf %lf %lf ",  & positions[counter_positions], & positions[counter_positions + 1], & masses[counter_masses]);
                                printf("%lf %lf %lf\n ",  positions[counter_positions], positions[counter_positions + 1], masses[counter_masses]);


                            counter_masses += 1;
                            counter_positions += 2;
                            counter_velocities += 2;


                            }
            }





    fclose(fp);
  }

  MPI_Bcast(masses, NUMBER_OF_PARTICLES, MPI_DOUBLE, MASTER_PROCESSOR_RANK, MPI_COMM_WORLD);//send same information
  MPI_Bcast(positions, 2 * num_of_processors * BODIES_per_proc, MPI_DOUBLE, MASTER_PROCESSOR_RANK, MPI_COMM_WORLD);
  MPI_Scatter(velocities, 2 * BODIES_per_proc, MPI_DOUBLE, velocities_per_proc, 2 * BODIES_per_proc, MPI_DOUBLE, MASTER_PROCESSOR_RANK, MPI_COMM_WORLD);
    //send chink of information
  double start_time = MPI_Wtime();

  // We simulate for the specified number of steps.
  int steps = 1;
  for (steps = 1; steps <= TOTAL_TIME_RUN; steps++) {
    compute_force(masses, positions, forces_each_proc, rank_of_processor, BODIES_per_proc);
    update_positions_velocities(positions, forces_each_proc, velocities_per_proc, rank_of_processor, BODIES_per_proc);
    MPI_Allgather(MPI_IN_PLACE, 1, aggregat_Type, positions, 1, aggregat_Type, MPI_COMM_WORLD);






    //just coordinates written in file for simualtion purpose
        if(rank_of_processor == MASTER_PROCESSOR_RANK)
        {

                   int     counter_masses1 = 0;
                    int counter_positions1 = 0;
                     int counter_velocities1 = 0;



        FILE *fp4 = fopen("justcoordinates.txt","a");

                         for (long long int line = 0; line < NUMBER_OF_PARTICLES; line++) 
                            {
            

                            // fscanf(fp, "%lf %lf %lf ",  & positions[counter_positions], & positions[counter_positions + 1], & masses[counter_masses]);
                            // printf("%lf %lf %lf\n ",  positions[counter_positions], positions[counter_positions + 1], masses[counter_masses]);
                            fprintf(fp4 ,"%lf ",  positions[counter_positions1]);

                            counter_masses1 += 1;
                            counter_positions1 += 2;
                            counter_velocities1 += 2;

                            }

                            fprintf(fp4,"\n");
                    counter_masses1 = 0;
                    counter_positions1 = 0;
                    counter_velocities1 = 0;

                         for (long long int line = 0; line < NUMBER_OF_PARTICLES; line++) 
                            {
            

                            // fscanf(fp, "%lf %lf %lf ",  & positions[counter_positions], & positions[counter_positions + 1], & masses[counter_masses]);
                            // printf("%lf %lf %lf\n ",  positions[counter_positions], positions[counter_positions + 1], masses[counter_masses]);
                            fprintf(fp4 ,"%lf ",  positions[counter_positions1+1]);

                            counter_masses1 += 1;
                            counter_positions1 += 2;
                            counter_velocities1 += 2;

                            }

                            fprintf(fp4,"\n");
                    counter_masses1 = 0;
                    counter_positions1 = 0;
                    counter_velocities1 = 0;





        fclose(fp4);
    }
    
  }
//inverse of MPI_GAther at upper line
  MPI_Gather(velocities_per_proc, 1, aggregat_Type, velocities, 1, aggregat_Type, MASTER_PROCESSOR_RANK, MPI_COMM_WORLD);

  if (DEBUG_MODE >= 1 && rank_of_processor == MASTER_PROCESSOR_RANK  ) 
  {
      FILE * final_state = fopen("final_state.txt", "w+");
      if (!final_state) {
        printf("Error creating output file.\n");
        exit(1);
      }

      int particle = 0;
      for (particle = 0; particle < NUMBER_OF_PARTICLES; particle++) {
        fprintf(final_state, "%lf %lf %lf %lf %lf\n", masses[particle], positions[2 * particle], positions[2 * particle + 1], velocities[2 * particle], velocities[2 * particle + 1]);
      }

      fclose(final_state);
  }

  double end_time = MPI_Wtime();
  if (rank_of_processor == MASTER_PROCESSOR_RANK) 
  {
    printf("Time take = %lf s.\n", end_time - start_time);
  }

  MPI_Type_free( & aggregat_Type);
  free(masses);
  free(positions);
  free(velocities_per_proc);
  free(forces_each_proc);
  if (rank_of_processor == MASTER_PROCESSOR_RANK) {
    free(velocities);
  }

  MPI_Finalize();
  return 0;
}









void compute_force(double masses[], double positions[], double forces_each_proc[], int rank_of_processor, int BODIES_per_proc) 
{

  // Starting and ending particle for the current processor.
  int starting_index = rank_of_processor * BODIES_per_proc;
  int ending_index = starting_index + BODIES_per_proc - 1;

  if (starting_index >= NUMBER_OF_PARTICLES) 
  {
    return;
  }
  else if (ending_index >= NUMBER_OF_PARTICLES) 
  {
    ending_index = NUMBER_OF_PARTICLES - 1;
  }


  int particle = starting_index;


  for (particle = starting_index; particle <= ending_index; particle++) 
  {
      double force_x = 0;
      double force_y = 0;
      int i = 0;
      for (i = 0; i < NUMBER_OF_PARTICLES; i++) 
      {
        if (particle == i) 
        {
           continue;
        }
        double x_diff = positions[2 * i] - positions[2 * particle];
        double y_diff = positions[2 * i + 1] - positions[2 * particle + 1];
        double distance = sqrt(x_diff * x_diff + y_diff * y_diff);
        double distance_cubed = distance * distance * distance;

        double force_total = GRAVITATIONAL_CONSTANT * masses[i] / distance;
        force_x += GRAVITATIONAL_CONSTANT * masses[i] * x_diff / distance_cubed;
        force_y += GRAVITATIONAL_CONSTANT * masses[i] * y_diff / distance_cubed;
      }

      forces_each_proc[2 * (particle - starting_index)] = force_x;
      forces_each_proc[2 * (particle - starting_index) + 1] = force_y;
      
      if (DEBUG_MODE >= 2) 
      {
          printf("Force on particle %i = %.3f  %.3f\n", particle, force_x, force_y);
      }

  }



}

void update_positions_velocities(double positions[], double forces_each_proc[], double velocities_per_proc[], int rank_of_processor, int BODIES_per_proc) 
{

  // Starting and ending particle for the current processor.
  int starting_index = rank_of_processor * BODIES_per_proc;
  int ending_index = starting_index + BODIES_per_proc - 1;

  if (starting_index >= NUMBER_OF_PARTICLES) {
    return;
  } else if (ending_index >= NUMBER_OF_PARTICLES) {
    ending_index = NUMBER_OF_PARTICLES - 1;
  }

  int particle = starting_index;
  for (particle = starting_index; particle <= ending_index; particle++) 
  {
      //   s = ut+1at^2
      positions[2 * particle] += velocities_per_proc[2 * (particle - starting_index)] * DELTA_T + (forces_each_proc[2 * (particle - starting_index)] * DELTA_T * DELTA_T / 2);
      positions[2 * particle + 1] += velocities_per_proc[2 * (particle - starting_index) + 1] * DELTA_T + (forces_each_proc[2 * (particle - starting_index) + 1] * DELTA_T * DELTA_T / 2);

      // v = u+at 
      velocities_per_proc[2 * (particle - starting_index)] += forces_each_proc[2 * (particle - starting_index)] * DELTA_T;
      velocities_per_proc[2 * (particle - starting_index) + 1] += forces_each_proc[2 * (particle - starting_index) + 1] * DELTA_T;
      
      
      if (DEBUG_MODE >= 2) 
      {
          printf("Position of particle %i = %.3f  %.3f\n", particle, positions[2*particle], positions[2*particle + 1]);
      }
  }


}
