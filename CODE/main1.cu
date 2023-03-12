


#include <stdio.h>
#include <stdlib.h>
// #include <cutil.h>
#include<time.h>
#include <iostream>
#include <string>
#include <vector>



#define DEBUG_MODE 0
#define MASTER_PROCESSOR_RANK 0
#define DELTA_T 10
#define INPUT_FILE_DATA "inputfile.txt"

#define ToalTimeStep 10000
const float DEG2RAD = 3.14159/180;

using namespace std;
double  Gconstant = 6.6674 * pow(10,-11);
// #define Gconstant 6.6674 * pow(10,-11)


typedef long long int lld;

// lld n;
lld numeroflines ;


// double *arrx1 *arry2;


// __global__ void kernelcomputenewvel(double* idx ,double *masslist, double *velocitylistx,double *velocitylisty, double *newvelocitylistx,double *newvelocitylisty,double *postionx,double *postiony,double *newpostionx,double *newpostiony ,lld threadcunt )
__global__ void kernelcomputenewvel(double masslist[], double velocitylistx[],double velocitylisty[], double postionx[],double postiony[],double Forcex1[],double Forcey2[] ,lld idx ,lld threadcunt,lld n )
{
    // lld threadid = blockDim.x * blockIdx.x + threadIdx.x;
    lld threadid =  threadIdx.x;

    double m2 ,r1vector,r2vector,rdiffector;
    // double m1 = masslist[idx];
    // double curentpos_1x = postionx[idx];
    // double curentpos_1y = postiony[idx];

    // double curentvel_1x = velocitylistx[idx];
    // double curentvel_1y = velocitylisty[idx];

    Forcex1[threadid] = 0;
    Forcey2[threadid] = 0 ;


    for(lld i1 = threadid;i1<n;i1+=threadcunt)
    {
        // c1temp[i1] = a1[i1]+b1[i1];
        if(i1!=idx)
        {
            m2 = masslist[i1] ;
            r1vector = postionx[i1] - postionx[idx];
            r2vector = postiony[i1] - postiony[idx];

            rdiffector = r1vector*r1vector+r2vector*r2vector;
            rdiffector = pow(rdiffector,1.5);

            // Forcey2[threadid] = Forcey2[threadid] + (((Gconstant*masslist[idx]*m2)/rdiffector)*r2vector) ;
            // Forcex1[threadid] = Forcex1[threadid] + (((Gconstant*masslist[idx]*m2)/rdiffector)*r1vector) ;
             
            // Forcey2[threadid] = Forcey2[threadid] + (((masslist[idx]*m2)/rdiffector)*r2vector) ;
            // Forcex1[threadid] = Forcex1[threadid] + (((masslist[idx]*m2)/rdiffector)*r1vector) ;


            Forcey2[threadid] = Forcey2[threadid] + ((m2/rdiffector)*r2vector) ;
            Forcex1[threadid] = Forcex1[threadid] + ((m2/rdiffector)*r1vector) ;

        }
    }

    // __syncthreads(); // barrier point 


    //update all 
    // if(threadid ==0 )
    // {

    // }

}









int main(int argc, char **argv)
{
    int x1,x2;
    int deltavartobeused;
    // cout<<argc;
    double deltat = DELTA_T;

    if (argc-1==2)
    {
        //argv[1] number of threads ,number of bodies
         x1 = atoi(argv[1]);
         x2 = atoi(argv[2]);
        // cout<<x1;
        // if(argc-1==3)
        // {
        //     deltat = atoi
        // }

    }




    FILE * fp = fopen(INPUT_FILE_DATA, "r");
    if (!fp) {
      printf("Error opening input file.\n");
      exit(1);
    }

    FILE *fp2,*fp3,*fp4;

    // fp4 = fopen("justcoordinates.txt","w+");
    fp2 = fopen("outputtxt.txt","w");
    fp3 = fopen("outputtonlytimedone.txt","w");
    long long int nobj;
     double m1,newposxi1,newposyi1,m2,r1vector,r2vector,rdiffector;
    double Forcex1,Forcey2;
    
    double accelerationx,accerlerationy;
    double tempvx1,tempvy2;
    double average_kernel_time = 0;



    //host particle storage 
    double posx[120000];
    double posy[120000];
    double mass[120000];
    double velex[120000];
    double veley[120000];
    double *forecearrayx,*forecearrayy;
    double *hoforecearrayx,*hoforecearrayy;



    //device particle storage
    double *dposx;
    double *dposy;
    double *dmass;
    double *dvelex;
    double *dveley;

    lld maxnumberofobjs ;
    //temporary store
    // double *dtempposx;
    // double *dtempposy;
    // double *dtempvelex;
    // double *dtempveley;
    
    // lld threads_in_block = 1024;
    lld threads_in_block = x1;




    // hoforecearrayx = new double[]  ;

    long long int array_siez,numer_of_lines;
    
    numer_of_lines = ToalTimeStep/deltat;

      fscanf(fp, "%lld", &maxnumberofobjs);
      cout<<"MAXIMUM NUMBER of objects "<<maxnumberofobjs<<endl;
    //   nobj = 2048;
    //   nobj = 5040;
    nobj = x2 ; 

    //   n = nobj ;
      cout<<"Selected number of objects "<<nobj<<endl;
      cout<<"Number of Threads "<<threads_in_block<<endl;
      cout<<"Deltat "<<deltat<<endl;


        fp4 = fopen("justcoordinates.txt","w+");
        
        fprintf(fp4,"%lld",nobj);

        fprintf(fp4,"\n");
        fprintf(fp4,"%lld",numer_of_lines);
        fprintf(fp4,"\n");

        // fclose(fp4);
    for (lld rep = 0; rep < nobj; rep++)
    {
        // cin>>posx[rep]>>posy[rep]>>mass[rep];

        fscanf(fp, "%lf %lf %lf", &posx[rep],&posy[rep],&mass[rep]);
        velex[rep] = veley[rep] = 0;
    }

    if(DEBUG_MODE>=2)
    {
        for (lld rep = 0; rep < nobj; rep++)
        {
            // cin>>posx[rep]>>posy[rep]>>mass[rep];
    
            fprintf(fp3, "%lf %lf %lf\n", posx[rep],posy[rep],mass[rep]);
            // velex[rep] = veley[rep] = 0;
        }

    }
    


    float et;
    cudaEvent_t start, stop;
  
    clock_t startc, end;
    startc = clock();



    // fclose(fp);
    array_siez = nobj;
     lld array_bytes = array_siez*sizeof(double);
     lld threadsize = threads_in_block*sizeof(double);
    hoforecearrayx = (double*)malloc(threadsize);
    hoforecearrayy = (double*)malloc(threadsize);
    // hoforecearrayx = new double[threads_in_block];
    // hoforecearrayy = new double[threads_in_block];



    cudaMalloc((void **)&dposx,array_bytes );
    cudaMalloc((void **)&dposy,array_bytes );
    cudaMalloc((void **)&dmass,array_bytes );
    cudaMalloc((void **)&dvelex,array_bytes );
    cudaMalloc((void **)&dveley,array_bytes );
    cudaMalloc((void **)&forecearrayx,threadsize );
    cudaMalloc((void **)&forecearrayy,threadsize );










    for (lld timestep = 0; timestep*deltat < ToalTimeStep ; timestep++)
        {
            // fprintf(fp2,"%lf",timestep*deltat);
            // fprintf(fp2,"\n");
       
            for(lld i1 = 0 ;i1<nobj; i1++)
            {
                //each body send data to kernel 
                //a parallel compute vel

                //cpoy data to gpu

                cudaEventCreate(&start);
                cudaEventCreate(&stop);
                cudaEventRecord(start,0);
                cudaEventRecord(stop,0);
            

               

            cudaMemcpy(dposx,posx,array_bytes,cudaMemcpyHostToDevice);
            cudaMemcpy(dposy,posy,array_bytes,cudaMemcpyHostToDevice);
            cudaMemcpy(dmass,mass,array_bytes,cudaMemcpyHostToDevice);
            cudaMemcpy(dvelex,velex,array_bytes,cudaMemcpyHostToDevice);
            cudaMemcpy(dveley,veley,array_bytes,cudaMemcpyHostToDevice);
            // cudaMemcpy(dveley,veley,array_bytes,cudaMemcpyHostToDevice);
            // cudaMemcpy(dveley,veley,array_bytes,cudaMemcpyHostToDevice);




        kernelcomputenewvel<<<1,threads_in_block>>>(dmass,dvelex,dveley,dposx,dposy,forecearrayx,forecearrayy,i1,threads_in_block,nobj);
// __global__ void kernelcomputenewvel(double idx ,double masslist[], double velocitylistx[],double velocitylisty[], double postionx[],double postiony[],lld threadcunt,double Forcex1[],double Forcey2[] )
            

            // cudaMemcpy(hoforecearrayx,forecearrayx,array_bytes,cudaMemcpyDeviceToHost);
            // cudaMemcpy(hoforecearrayy,forecearrayy,array_bytes,cudaMemcpyDeviceToHost);


            cudaMemcpy(hoforecearrayx,forecearrayx,threadsize,cudaMemcpyDeviceToHost);
            cudaMemcpy(hoforecearrayy,forecearrayy,threadsize,cudaMemcpyDeviceToHost);

            cudaEventSynchronize(stop);
        cudaEventElapsedTime(&et, start, stop);
        cudaEventDestroy(start);
        cudaEventDestroy(stop);


            if(DEBUG_MODE>=1)
            {

                printf("\nGPU Time to generate kernel %lld:  %f  \n", timestep*nobj+i1,et);
            }
            average_kernel_time+=et; 

            Forcex1 = 0;Forcey2 = 0;
            for(int ij1=0;ij1<threads_in_block;ij1++)
            {
                Forcex1 += hoforecearrayx[ij1];
                Forcey2 += hoforecearrayy[ij1] ;
            }
           

            //once done update with new vector

            m1 = mass[i1];

            // accelerationx = Gconstant*Forcex1/m1;
            accelerationx = Gconstant*Forcex1;

            tempvx1 = velex[i1] + accelerationx*deltat;
            posx[i1] += (velex[i1]*deltat+(0.5*accelerationx*deltat*deltat) );
            velex[i1] = tempvx1;

    //moement in y direction 
            
            // accerlerationy = Gconstant*Forcey2/m1;
            accerlerationy = Gconstant*Forcey2;

            tempvy2 = veley[i1] + accerlerationy*deltat;
            posy[i1] += (veley[i1]*deltat+(0.5*accerlerationy*deltat*deltat) );
            veley[i1] = tempvy2;

             if(DEBUG_MODE>=2)
                {
                  cout<<" acceleration at each step of objects"<<accelerationx<<" "<< accerlerationy << endl;
                }




            }

            
            //write postion x ,position y


            for (lld xdata = 0; xdata < nobj; xdata++)
            {
    
                // printf("%lf,",posx[xdata]);
                // fprintf(fp2,"%lf ",posx[xdata]);
                fprintf(fp4,"%lf ",posx[xdata]);
    
            }
            fprintf(fp4,"\n");



            for (lld ydata = 0; ydata < nobj; ydata++)
            {
                // printf("%lf,",posx[ydata]);
                
                fprintf(fp4,"%lf ",posy[ydata]);
    
            }
            fprintf(fp4,"\n");




        }



        
    end = clock();

    printf("\naverage kernel time :  %f  \n", float((deltat*average_kernel_time) /(ToalTimeStep*nobj) ) );
    printf("\nTotal time:  %f  \n", float(end-startc) );











        cudaFree(dposx );
        cudaFree(dposy );
        cudaFree(dmass );
        cudaFree(dvelex );
        cudaFree(dveley );
        cudaFree(forecearrayx );
        cudaFree(forecearrayy );
    

        // cudaFree(d_in1);
        // cudaFree(d_in2);
        // cudaFree(d_out);
        free(hoforecearrayx);
        free(hoforecearrayy);
        // hoforecearrayx = new double[threads_in_block];
    // hoforecearrayy = new double[threads_in_block];

    fclose(fp);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);

        return 0;



}






