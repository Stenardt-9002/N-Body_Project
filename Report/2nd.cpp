// g++ -fopenmp 2nd.cpp -o run1 -lm
// OMP_NUM_THREADS=4 ./run1 

#include<omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include<iostream>
typedef long long int lld;
#define ToalTimeStep 10000
// #define ToalTimeStep 20

using namespace std;


// #define Gconstant 6.6674 * pow(10,-11)
double  Gconstant = 6.6674 * pow(10,-11);



int main(int argc, char const *argv[])
{

    lld n;
    int zemp;
    // cin>>n;
    // double *posx,*posy,*mass; //in kg
    // double *velex,*veley; //in m/s

    double posx[120000];
    double posy[120000];
    double mass[120000];
    double velex[120000];
    double veley[120000];


    double m1,newposxi1,newposyi1,m2,r1vector,r2vector,rdiffector;
    double Forcex1,Forcey2;
    double accelerationx,accerlerationy;
    double tempvx1,tempvy2;
    double deltat = 5;
    double taoltime1 = 0;
    // posx = new double[n];
    // posy = new double[n];
    // mass = new double[n];
    // velex = new double[n];
    // veley = new double[n];

    FILE *fp1,*fp2,*fp3,*fp4;

    fp1 = fopen("inputfile.txt","r");  
    fp4 = fopen("justcoordinates.txt","w");
    fp2 = fopen("outputtxt.txt","w");
    fp3 = fopen("outputtonlytimedone.txt","w");


    fscanf(fp1,"%lld",&n);
    n = 4000 ;//only 2000 objects taken 
    cout<<n;
    lld numoflines = ToalTimeStep/deltat;
    fprintf(fp4,"%lld",n);

    fprintf(fp4,"\n");
    fprintf(fp4,"%lld",numoflines);
    fprintf(fp4,"\n");



    for (lld rep = 0; rep < n; rep++)
        {
            // cin>>posx[rep]>>posy[rep]>>mass[rep];

            fscanf(fp1, "%lf %lf %lf", &posx[rep],&posy[rep],&mass[rep]);
            velex[rep] = veley[rep] = 0;
        }
    // for (int rep = 0; rep < n; rep++)
    //     {
    //         // cout<<posx[rep]<<" "<<posy[rep]<<" "<<mass[rep]<<endl;;
    //             fprintf(fp2,"%lf %lf ",posx[rep],posy[rep]);
    //             fprintf(fp2,"\n");
    //     }
        






     for (lld timestep = 0; timestep*deltat < ToalTimeStep ; timestep++)
    {
        fprintf(fp2,"%lf",timestep*deltat);
        fprintf(fp2,"\n");
        lld j1;

        double tbegin = omp_get_wtime();

        // #pragma omp for 
        #pragma omp parallel for schedule(static) private(j1,r1vector,r2vector,m1,m2,Forcex1,Forcey2,rdiffector,accelerationx,accerlerationy,velex,veley)
        for (lld i1 = 0; i1 < n; i1++)
        {
            //each body 
            // current status then update each body wrt static position of 

            m1 = mass[i1];
            newposxi1 = posx[i1];newposyi1 = posy[i1];
            Forcex1 = 0;
            Forcey2 = 0;
            //consider all bodies except i1 th to effect n(i1)
            for (j1 = 0; j1 < n ; j1++)
            {
                if (j1!=i1)
                {
                    m2 = mass[j1];
                    r1vector = posx[j1]-posx[i1] ; //wrt to current body
                    r2vector = posy[j1]-posy[i1] ; //wrt to current body
                    rdiffector  = r1vector*r1vector+r2vector*r2vector;
                    rdiffector = pow(rdiffector,1.5);

                    Forcex1 = Forcex1+ (((Gconstant*m1*m2)/rdiffector)*r1vector) ;
                    Forcey2 = Forcey2+ (((Gconstant*m1*m2)/rdiffector)*r2vector) ;



                }
                
            }
            // time difference delatat delta 
            
            accelerationx = Forcex1/m1;
            tempvx1 = velex[i1] + accelerationx*deltat;
            posx[i1] += (velex[i1]*deltat+(0.5*accelerationx*deltat*deltat) );
            velex[i1] = tempvx1;

    //moement in y direction 
            
            accerlerationy = Forcey2/m1;
            tempvy2 = veley[i1] + accerlerationy*deltat;
            posy[i1] += (veley[i1]*deltat+(0.5*accerlerationy*deltat*deltat) );
            veley[i1] = tempvy2;


        
        
            
        
        }

        // //debug mode 
        // cin>>zemp;
        // cout<<"FORCE1"<<Forcex1;
        // cout<<"FORCE2"<<Forcey2;
        // cin>>zemp;

          double wtime = omp_get_wtime() - tbegin;
          fprintf( fp3,"%lf", wtime );
          fprintf(fp3,"\n");
          taoltime1 +=wtime;





        //print data at each time step
        for (lld xdata = 0; xdata < n; xdata++)
        {

            // printf("%lf,",posx[xdata]);
            fprintf(fp2,"%lf ",posx[xdata]);
            fprintf(fp4,"%lf ",posx[xdata]);

        }
        fprintf(fp2,"\n");
        fprintf(fp4,"\n");


        for (lld ydata = 0; ydata < n; ydata++)
        {
            // printf("%lf,",posx[ydata]);
            
            fprintf(fp2,"%lf ",posy[ydata]);
            fprintf(fp4,"%lf ",posy[ydata]);

        }

        fprintf(fp2,"\n");
        fprintf(fp4,"\n");


        fprintf(fp2,"\n");




    }




















    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);



        cout<<"\nTime "<<taoltime1<<endl;




    // delete []posx;
    // delete []posy;
    // delete []mass;
    // delete []velex;
    // delete []veley;


    return 0;
}








