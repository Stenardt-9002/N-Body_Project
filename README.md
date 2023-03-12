# N Body Simulation using C++/ C-Cuda /C-OpenMP /C-MPI 
The simulation is done using freeglutdev OpenGL (legacy)(See requirements).

---

The following solution was part of activity in HPC course offered in fall of 2020 @ IIITDM Kancheepuram

---

<br>

**Content**
1. [Introduction](#introduction)
1. [OpenMP](#openmp)
1. [MPI](#mpi)
    1. [Approach](#approach)
    1. [Params](#parametres)
1. [CUDA](#cuda)
    1. [Strategy](#strategy)
    1. [Specification](#specification)
1. [Simulation (Only Cuda)](#simulation)
1. [Test](#)
<br>
<br>


## Introduction
<br>

The N-body problem is one of famous problems in classical physics for predicting the motion of
n celestial body that interacts gravitationally in free space.It is a problem for predicting individual
motion of bodies starting from a quasi state.

The problem has been motivation to understand motions of sun ,planets and other celestial
bodies in global clusters . Consider general relativity the problem is difficult to solve and still is
an open problem .See ​two-body problem​ and restricted ​three-body problem​ which have been
solved .

The below problem solution is based on problem of simulation of random 10000 masses
ranging from 34000,25000000 kg on ​2d​ plane on (-1000,1000) on both x and y coordinates
where initial states of bodies are in quasi state (initial velocity and acceleration are 0 in both x
axis and y axis) and initial position are ((-500,500)|(-400,600)) in x and y axis respectively .

>Assumption:
For simulation purpose and ease of calculation classical newton laws are used for computing
velocity and acceleration of individual bodies .

$$ F(i,j)(t) =  {{Gm_i m_j} {|r_i(t)-r_j(t)|}\over{{|r_i(t)-r_j(t)|}^3}}$$


Where i and j the body are applying f(i,j) force on each other. (Newton 3rd law) where G =6.67 ×
10​-11​ Newtons kg​-2​ m​2

NOTE: While calculating position and velocity of actual planetary motion Newton Laws are no
longer valid ,due to fixation of barycente

<br>
<br>
<br>


## OpenMP
<br>

The algorithm uses block distribution as symmetry of Force , F(i,j) = F(j,i) and cyclic distribution
for force calculations ,since all bodies force value computation is calculated and then updated at
end simultaneously .

Since OpenMP uses shared memory for parallelisation thus cyclic distribution computation
parallelisation is easy.

Thus #pragma omp parallel for directive is added for force computation from all n-1 bodies and
workload is distributed on static scheduling where chunks of data are equally distributed to the
number of threads that are available .

To avoid race conditions each data is saved in individual rvector position ,velocity data in vector
data type in c++ so that each thread in parallel constructs ensures access of one element at
index .

Parallelization is done on force of n-1 body computation is independent of each other thus no
memory race condition is introduced .
Fi = ​Σ​Fj (j=1:n j!=i)



```c
For each particle i do:
foreach particlej > i do
ompsetlock(& locks[i])
Fi(t)​+=​fi,j(t).
ompunsetlock(& locks[i])
ompsetlock(& locks[j])
Fj(t)​-=​fi,j(t).
ompunsetlock(& locks[j])
end for
end for
```


```c++
for​ (​lld​ ​timestep​ = ​0​; ​timestep​*​deltat​ < ToalTimeStep ; ​timestep​++)
fprintf​(​fp2​,​"%lf"​,​timestep​*​deltat​);
fprintf​(​fp2​,​"​\n​"​);
lld​ ​j1​;
double​ ​tbegin​ = ​omp_get_wtime​();
// #pragma omp for
#pragma​ ​omp​ ​parallel​ ​for​ ​schedule​(​static​)
private​(​j1​,​r1vector​,​r2vector​,​m1​,​m2​,​Forcex1​,​Forcey2​,​rdiffector​,​acceleration
x​,​accerlerationy​,​velex​,​veley​)
for​ (​lld​ ​i1​ = ​0​; ​i1​ < ​n​; ​i1​++)
{
//each body
// current status then update each body wrt static position of
m1​ = ​mass​[​i1​];
newposxi1​ = ​posx​[​i1​];​newposyi1​ = ​posy​[​i1​];
Forcex1​ = ​0​;
Forcey2​ = ​0​;
//consider all bodies except i1 th to effect n(i1)
for​ (​j1​ = ​0​; ​j1​ < ​n​ ; ​j1​++)
{
if​ (​j1​!=​i1​)

{
m2​ = ​mass​[​j1​];
r1vector​ = ​posx​[​j1​]-​posx​[​i1​] ;​ //wrt to current body
r2vector​ = ​posy​[​j1​]-​posy​[​i1​] ;​ //wrt to current body
rdiffector​ = ​r1vector​*​r1vector​+​r2vector​*​r2vector​;
rdiffector​ = ​pow​(​rdiffector​,​1.5​);
Forcex1​ = ​Forcex1​+

(((​Gconstant​*​m1​*​m2​)/​rdiffector​)*​r1vector​) ;
Forcey2​ = ​Forcey2​+

(((​Gconstant​*​m1​*​m2​)/​rdiffector​)*​r2vector​) ;
}
}
// time difference delatat delta
accelerationx​ = ​Forcex1​/​m1​;
tempvx1​ = ​velex​[​i1​] + ​accelerationx​*​deltat​;
posx​[​i1​] += (​velex​[​i1​]*​deltat​+(​0.5​*​accelerationx​*​deltat​*​deltat​)
);
velex​[​i1​] = ​tempvx1​;
//moement in y direction
accerlerationy​ = ​Forcey2​/​m1​;
tempvy2​ = ​veley​[​i1​] + ​accerlerationy​*​deltat​;
posy​[​i1​] +=
(​veley​[​i1​]*​deltat​+(​0.5​*​accerlerationy​*​deltat​*​deltat​) );
veley​[​i1​] = ​tempvy2​;
}

```





<br>
<br>
<br>



## MPI

<br>

### Approach

Parallelisation with MPI:

Although using similar methodology used in openmp parallelisation for the same problem might
work but will incur much communication overhead,thus block partitioning methodology is used.

Block partitioning methodology is used for MPI thus utilising less memory and improving load
balancing. Cache misses are still prominent but load balance improvement makes up for it.

Strategy :

- Parallelisation of computing of force is done by distributing the array of mass vector to p
nodes reducing the workload to n/p.


- Initial vector position and mass of n bodies information is sent to each worker (because
computation of i th requires all rest n-1 bodies force).

- The velocities of the bodies are distributed among p workers each given their velocity for
usage of updation of individual rvector and velocity value.

- To avoid race condition rvector ,velocity ,acceleration(force) is stored in individual array
data type and updation is done once force is computed for rest n-1 bodies.

MPI_Bcast(masses, NUMBER_OF_PARTICLES, MPI_DOUBLE,
MASTER_PROCESSOR_RANK, MPI_COMM_WORLD);//send same information

MPI_Bcast(positions, 2 * num_of_processors * particles_per_processor,
MPI_DOUBLE, MASTER_PROCESSOR_RANK, MPI_COMM_WORLD);

>Note velocities are scatter (equally distributed)

MPI_Scatter(velocities, 2 * particles_per_processor, MPI_DOUBLE,
curr_proc_velocities, 2 * particles_per_processor, MPI_DOUBLE,
MASTER_PROCESSOR_RANK, MPI_COMM_WORLD);

- For each worker force is calculated and once force(acceleration )is known those specific
bodies velocities are updated.Specific distributed velocities(current proc_velcoties ) are collected again MPI_GATHER.

- Then all arrays of each quantities are freed.
```
//READ DATA
Initialize MPI and call n workers
//Share rvector and masses to all workers .
//Distribute velocity
Parallelise p workers
For each time step:
For each body data of worker :
Compute force from rest bodies under that specific
worker and update it
For each body in worker:
Compute and update Velocity
Compute and update rvector
//Wait until all positions are updated
//MPI_Allgather
Gather velocities of each worker after the individual worker is
done with the task .
//Free Mass
//Free(Velocity)
//Free(Curr_proc_velocity)
//Free(Forces)
End MPI_routine

```

###  Parametres
- Time step incremented by deltaT
> NOTE:Varying deltaT and making it small(1-5) will make simulation smooth but computation
time taken will be extremely large for number of bodies

- Number_of_particles: Increasing number of bodies (maximum 111002) will increase
computation but make simulation uniform .
- Total_Time_Run : The number of timestep until which the simulation will run .
> NOTE:Increasing the Total_Time_Run will increase the runtime of simulation.

- counter_velocities :A single dimension array of length 2n containing Velocity-x and Velocity-y of
each n body at contiguous memory location.

>OpenGL is used for simulating the bodies (look into simulation.cpp)

>NOTE: The simulation is attempted only when all the computation is done for all specific
bodies.It is not done simultaneously when calculating individual velocity and position of bodies
at each timestamp.

- DEBUG_MODE : variable to set 1 for light information display and 2 for all information displayed
at each time interval .



<br>
<br>
<br>




## CUDA
<br>
Parallelisation is achieved by making each body force computation parallel by calling kernel.

Although calling kernel multiple times causes overhead ,dividing all bodies force computation
requires less kernel call but overloads memory limitations of kernel.

Initial vector position and mass of n bodies information is sent to the kernel (because
computation of i th requires all rest n-1 bodies force).

To avoid race condition rvector ,velocity ,acceleration(force) is stored in individual array data
type and updation is done once force is computed from rest n-1 bodies.

### Strategy
```c++
cudaMemcpy(dposx,posx,array_bytes,cudaMemcpyHostToDevice);
cudaMemcpy(dposy,posy,array_bytes,cudaMemcpyHostToDevice);
cudaMemcpy(dmass,mass,array_bytes,cudaMemcpyHostToDevice);
cudaMemcpy(dvelex,velex,array_bytes,cudaMemcpyHostToDevice);
cudaMemcpy(dveley,veley,array_bytes,cudaMemcpyHostToDevice);
kernelcomputenewvel<<<1,threads_in_block>>>(dmass,dvelex,dveley,dposx,dposy,forece
arrayx,forecearrayy,i1,threads_in_block,nobj);
// __global__ void kernelcomputenewvel(double idx ,double masslist[], double
velocitylistx[],double velocitylisty[], double postionx[],double postiony[],lld threadcunt,double
Forcex1[],double Forcey2[] )
cudaMemcpy(hoforecearrayx,forecearrayx,threadsize,cudaMemcpyDeviceToHost);
```



### Specification
<br>

> GPU usage (Tesla P100-PCIE-16GB,memory available 16 G.B.(distributed on server))
//devicequery output
```bash
Device 0: "Tesla P100-PCIE-16GB"
CUDA Driver Version / Runtime Version
10.1 / 10.1
CUDA Capability Major/Minor version number: 6.0
Total amount of global memory:
16281 MBytes (17071734784 bytes)
(56) Multiprocessors, ( 64) CUDA Cores/MP: 3584 CUDA Cores
GPU Max Clock rate:
1329 MHz (1.33 GHz)
Memory Clock rate:
715 Mhz
Memory Bus Width:
4096-bit
L2 Cache Size:
4194304 bytes
Maximum Texture Dimension Size (x,y,z)
1D=(131072), 2D=(131072, 65536), 3D=(16384,
16384, 16384)
Maximum Layered 1D Texture Size, (num) layers 1D=(32768), 2048 layers
Maximum Layered 2D Texture Size, (num) layers 2D=(32768, 32768), 2048 layers
Total amount of constant memory:
65536 bytes
Total amount of shared memory per block:
49152 bytes
Total number of registers available per block: 65536
Warp size:
32
Maximum number of threads per multiprocessor: 2048
Maximum number of threads per block:
1024
Max dimension size of a thread block (x,y,z): (1024, 1024, 64)
Max dimension size of a grid size (x,y,z): (2147483647, 65535, 65535)
```









<br>

## Simulation

<br>

run ```./run1.sh``` to generate simulation file 
Specify Number of threads and number of bodies respectively.

maincuda.cu contains n body parallelised solution .

run ```./simulate1.sh``` to run simulation (requires c++ legacy opengl)

Read entire report  ```Report\main.html```



