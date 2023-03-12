// g++ simulation.cpp -lGL -lGLU -lglut -lm -o pat1
// g++ 2second.cpp -lGL -lGLU -lglut -lm -o pat1

// ./pat1

#define GL_GLEXT_PROTOTYPES
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
const float DEG2RAD = 3.14159/180;

#include<iostream>
#include <vector>
using namespace std;

typedef long long int lld;

lld n;
lld numeroflines ;
vector<vector<double>>arrx1;
vector<vector<double>>arry2;
void drawsrray()
{


  // glClearColor (1.0, 1.0, 1.0, 1);
  glClearColor (0, 0, 0, 1);


  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glOrtho(-1000.0, 1000.0, -1000.0, 1000.0, -1000.0, 1000.0);
  glColor3f(1,1,1);
  for (lld varoflines = 0; varoflines < numeroflines; varoflines++)
  {
    glBegin(GL_POINTS);

    for (lld gotem = 0; gotem < n; gotem++)
    {
        glVertex2f(arrx1[varoflines][gotem],arry2[varoflines][gotem]);
    }

    for (lld curfew = 0; curfew < 10000000; curfew++)
    // for (lld curfew = 0; curfew < 10000; curfew++)
    {
      curfew++;curfew--;
    }
    glEnd();

    glFlush();
    // glClear();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


  }

  cout<<"\nEND REACHED";
}




int main(int argc, char **argv)
{

    // float a,b;
    // cin>>a>>b;
    // ELLIPS_A = a ; ELLIPS_B = b;
    FILE *fp4;
    fp4 = fopen("justcoordinates.txt","r");

    fscanf(fp4,"%lld",&n);
    fscanf(fp4,"%lld",&numeroflines);

    // cin>>n;
    // double arrayx1[120000];
    // double arrayy1[120000];
    vector<vector<double>>arx1(numeroflines);
    vector<vector<double>>ary2(numeroflines);
    lld var = 0;
    double temvar;
    while (var<numeroflines)
    {
      vector<double>tempx1(n) ;
      vector<double>tempy2(n) ;
      for (lld iw = 0; iw < n; iw++)
      {
        fscanf(fp4,"%lf",&temvar);
        tempx1[iw] = temvar;
      }
      for (lld iw = 0; iw < n; iw++)
      {
        fscanf(fp4,"%lf",&temvar);
        tempy2[iw] = temvar;
      }

      arx1[var] = tempx1;
      ary2[var] = tempy2;
      var++;
    }

    // for (long long int i1 = 0; i1 < n; i1++)
    // {
    //   cin>>arrayx1[i1]
    // }
    cout<<n<<numeroflines;

    fclose(fp4);
    arrx1 = arx1;
    arry2 = ary2;


    // for (lld gotem = 0; gotem < n; gotem++)
    // {
    //     cout<<arrx1[1][gotem]<<" "<<arry2[1][gotem]<<endl;
    // }


    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(1000,1000);

    glutInitWindowPosition(500,500);
    glutCreateWindow("Legacy OpenGL - Simulation ");

    // drawellipse_util(a,b);
    glutDisplayFunc(drawsrray);
    glutMainLoop();


    return 0;
}
