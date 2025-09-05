// PROGRAM POISSON 2D
//
// This program calculate the electric potential in a square zone, 0<=x,y<=L,
// with boundary conditions of Dirichlet type: V(0,y)=V1, V(L,y)=V2, V(x,0)=V3, V(x,L)=V4;
//   

# include <iostream>
# include <math.h>
# include <stdio.h>                      
# include <stdlib.h>

using namespace std;

int main()
{
 int npx=100,npy=npx;                 // Number of points along x and y axis, respectively.
 int it,itmax=10000;                  // Index of iteration and limit of iterations.
 double L=0.2,Ln;                     // Dimension of the square (m).
 double h,om=1.98;
 double V1=0.,V2=0.,V3=-100., V4=100.; 
 double Vmax=V4;
 double x[npx],y[npy];
 double v[npx][npy],vold[npx][npy];  // Two-dimensional Array for electric potential.
 double max,er_abs,tol=1e-6; 
 
 FILE *out1,*out2,*out3,*out4;

 out1=fopen("Potential_xy.txt","w");
 out2=fopen("Potential_x.txt","w");
 out3=fopen("Potential_y.txt","w");
 out4=fopen("Tol_vs_It.txt","w");


// Normalize Variables 

   V1=V1/Vmax;
   V2=V2/Vmax;
   V3=V3/Vmax;
   V4=V4/Vmax;

   Ln=L/L;

// To define grid points.

   h=Ln/(npx-1);          // paso espacial

   for(int i=0;i<npx;i++) 
   {
       x[i]=h*i;
   } 

   for(int j=0;j<npy;j++)
   {
       y[j]=h*j;
   }


// To initialize potential's matrix elements 

   for(int i=0;i<npx;i++)
   {
    for(int j=0;j<npy;j++)
    {
       v[i][j]=0.0;
    }
   }
 

//Discretization of boundary conditions

   for(int j=1;j<=npy-1;j++)
    {
      v[0][j]=V1;
      v[npx-1][j]=V2;
    }

   for(int i=0;i<npx;i++)
    {
      v[i][0]=V3;
      v[i][npy-1]=V4;
    }
 
// Start iterative process

   it=0;
 
   while(it<=itmax)
    {
      max=0.0;            // Clear maximun absolute error of last iteration.
      it=it+1;

            
      for(int i=0;i<npx;i++)
        {
          for(int j=0;j<npy;j++)
            {
              vold[i][j]=v[i][j];
            }
         }

 // Estimate the electric potential on points of the mesh.
 
      for(int i=1;i<npx-1;i++) 
        {
          for(int j=1;j<npy-1;j++)
           {
                     v[i][j]=0.25*(v[i+1][j]+v[i-1][j]+v[i][j+1]+v[i][j-1]);// (Gauss-Seidel Method)
                
                 //     v[i][j]= vold[i][j]+(om/4.)*(vold[i+1][j]+v[i-1][j]+ v[i][j-1]+vold[i][j+1]-4.*vold[i][j]); // (SOR)
           }
        }

 // Estimate maximun absolute error in the approximation.

      for(int i=1;i<npx-1;i++)
        {
          for(int j=1;j<npy-1;j++)
           {
             er_abs=fabs(v[i][j]-vold[i][j]);
            // er_abs=pow(er_abs,2);        // Warning
            // er_abs=sqrt(er_abs);
                        
             if(er_abs> max)
              {
                 max=er_abs;
              }
           }
         }

                 cout<<"it=";
                 cout<< it;
                 cout<<"  ";
                 cout<<"error=";
                 cout<< max;
                 cout<<"\n";

         fprintf(out4,"%i   %f\n",it,max);

            


      if(max<=tol) // To check the convergence to save results
 //       if (it==1000)
        {
          // Print two-dimensional array of potential
             
            for(int i=0;i<npx;i++)
              {
               for(int j=0;j<npy;j++)
                {
                    fprintf(out1,"%e",v[i][j]*Vmax); // Save results with the format specified.
                    fprintf(out1,"   ");

                }
                    fprintf(out1,"\n");
              }
                    fclose(out1);

           // Print table: V(x,L/2) vs x
           
             for(int i=0;i<npx;i++)
               {
                    fprintf(out2,"%10.3f",x[i]*L);
                    fprintf(out2,"   ");
                    fprintf(out2,"%e",v[i][npy/2]*Vmax);
                    fprintf(out2,"\n");
               }
                    fclose(out2);

            // Print table: V(L/2,y) vs y
           
             for(int j=0;j<npy;j++)
               {
                    fprintf(out3,"%10.3f",y[j]*L);
                    fprintf(out3,"   ");
                    fprintf(out3,"%e",v[npx/2][j]*Vmax);
                    fprintf(out3,"\n");
               }
                    fclose(out3);
                    

    	  exit(1);
         }
    }
     fclose(out4);

    if(it>=itmax)
      {
        cout<< " Maximun of iterations has been exceded...!";
        cout<< "\n";
      }

  return 0;
}

