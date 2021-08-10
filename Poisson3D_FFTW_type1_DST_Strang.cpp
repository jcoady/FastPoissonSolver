#include <stdio.h>
#include <math.h>
#include <cmath>
#include <fftw3.h>
#include <iostream>
#include <vector>

using namespace std;
int main() {

    int N1=8;
    int N2=16;
	int N3=24;

    double pi = 3.141592653589793;
    double L1  = 1.0;
    double dx = L1/(double)(N1+1);//+ instead of -1
    double L2= 2.0;
    double dy=L2/(double)(N2+1);
    double L3= 3.0;
    double dz=L3/(double)(N3+1);

    //double invL1s=1.0/(L1*L1);
    //double invL2s=1.0/(L2*L2);
    //double invL3s=1.0/(L3*L3);

    std::vector<double> in1(N1*N2*N3,0.0);
    std::vector<double> in2(N1*N2*N3,0.0);
    std::vector<double> out1(N1*N2*N3,0.0);
    std::vector<double> out2(N1*N2*N3,0.0);
    std::vector<double> X(N1,0.0);
    std::vector<double> Y(N2,0.0);
    std::vector<double> Z(N3,0.0);


    fftw_plan p, q;
    int i,j,k;
    p = fftw_plan_r2r_3d(N1,N2,N3, in1.data(), out1.data(), FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
    q = fftw_plan_r2r_3d(N1,N2,N3, in2.data(), out2.data(), FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);

    int l=0;

    for(i = 0;i <N1;i++){
        for(j = 0;j<N2;j++){
			for(k = 0;k<N3;k++){
				X[i] =dx+(double)i*dx ;      
				Y[j] =dy+ (double)j*dy ;
				Z[k] =dz+ (double)k*dz ;
				//in1[l]= sin(pi*X[i])*sin(pi*Y[j]) ; // row major ordering
				in1[l]=2*Z[k]*(L3-Z[k])*Y[j]*(L2-Y[j])+2*Z[k]*(L3-Z[k])*X[i]*(L1-X[i])+2*Y[j]*(L2-Y[j])*X[i]*(L1-X[i]);
				l=l+1;
			}
        }
    }

    fftw_execute(p);

    l=-1;
    for ( i = 0; i < N1; i++){   // f = g / ( kx² + ky² + kz**2 )  
        for( j = 0; j < N2; j++){
			for( k = 0; k < N3; k++){

				l=l+1;
				double fact=0;

				//fact=(2-2*cos((i+1)*pi/(N1+1)))*invL1s;

				//fact+= (2-2*cos((j+1)*pi/(N2+1)))*invL2s;

				//fact+= (2-2*cos((k+1)*pi/(N3+1)))*invL3s;

				fact=(2-2*cos((i+1)*pi/(N1+1)))/(dx*dx);

				fact+= (2-2*cos((j+1)*pi/(N2+1)))/(dy*dy);

				fact+= (2-2*cos((k+1)*pi/(N3+1)))/(dz*dz);

				in2[l] = out1[l]/fact;
			}
        }
    }

    fftw_execute(q);
    l=-1;
    double erl1 = 0.;
    for ( i = 0; i < N1; i++) {
        for( j = 0; j < N2; j++){
			for( k = 0; k < N3; k++){
				l=l+1;
				X[i] =dx+(double)i*dx ;      
				Y[j] =dy+ (double)j*dy ;
				Z[k] =dz+ (double)k*dz ;
				//double res=0.5/pi/pi*in1[l];
				double res=X[i]*(L1-X[i])*Y[j]*(L2-Y[j])*Z[k]*(L3-Z[k]);
				erl1 +=pow(fabs(res-  0.125*out2[l]/((double)(N1+1))/((double)(N2+1))/((double)(N3+1))),2); 
				printf("%3d %10.7g %10.7g\n", l, res, 0.125*out2[l]/((double)(N1+1))/((double)(N2+1))/((double)(N3+1)));
				//erl1 +=pow(fabs(res-  3*0.125*out2[l]/((double)(N1+1))/((double)(N2+1))/((double)(N3+1))/((double)(N1+1))/((double)(N2+1))/((double)(N3+1))),2); 
				//printf("%3d %10.5g %10.5g\n", l, res,  3*0.125*out2[l]/((double)(N1+1))/((double)(N2+1))/((double)(N3+1))/((double)(N1+1))/((double)(N2+1))/((double)(N3+1)));
			}
        }
    }
    erl1=erl1/((double)N1*N2*N3);
    cout<<"error=" <<sqrt(erl1) <<endl ;  
    fftw_destroy_plan(p); fftw_destroy_plan(q); fftw_cleanup();

    return 0;
}