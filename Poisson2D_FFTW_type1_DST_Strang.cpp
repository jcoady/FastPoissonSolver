#include <stdio.h>
#include <math.h>
#include <cmath>
#include <fftw3.h>
#include <iostream>
#include <vector>

using namespace std;
int main() {

    int N1=50;
    int N2=75;

    double pi = 3.141592653589793;
    double L1  = 2.0;
    double dx = L1/(double)(N1+1);//+ instead of -1
    double L2= 3.0;
    double dy=L2/(double)(N2+1);
	
    std::vector<double> in1(N1*N2,0.0);
    std::vector<double> in2(N1*N2,0.0);
    std::vector<double> out1(N1*N2,0.0);
    std::vector<double> out2(N1*N2,0.0);
    std::vector<double> X(N1,0.0);
    std::vector<double> Y(N2,0.0);


    fftw_plan p, q;
    int i,j;
    p = fftw_plan_r2r_2d(N1,N2, in1.data(), out1.data(), FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);
    q = fftw_plan_r2r_2d(N1,N2, in2.data(), out2.data(), FFTW_RODFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);

    int l=0;

	for(i = 0;i <N1;i++){
		for(j = 0;j<N2;j++){
           X[i] =dx+(double)i*dx ;      
            Y[j] =dy+ (double)j*dy ;
            in1[l]=2*Y[j]*(L2-Y[j])+2*X[i]*(L1-X[i]);
            l=l+1;
        }
    }

    fftw_execute(p);

    l=-1;
	for ( i = 0; i < N1; i++){   // f = g / ( kx² + ky² )  
		for( j = 0; j < N2; j++){

            l=l+1;
            double fact=0;

            fact=(2-2*cos((i+1)*pi/(N1+1)))/(dx*dx);

            fact+= (2-2*cos((j+1)*pi/(N2+1)))/(dy*dy);

            in2[l] = out1[l]/fact;

        }
    }

    fftw_execute(q);
    l=-1;
    double erl1 = 0.;
	for ( i = 0; i < N1; i++) {
		for( j = 0; j < N2; j++){
            l=l+1;
            X[i] =dx+(double)i*dx ;      
            Y[j] =dy+ (double)j*dy ;
            double res=X[i]*(L1-X[i])*Y[j]*(L2-Y[j]);
            erl1 +=pow(fabs(res-  0.25*out2[l]/((double)(N1+1))/((double)(N2+1))),2); 
            printf("%3d %10.15g %10.15g\n", l, res,  0.25*out2[l]/((double)(N1+1))/((double)(N2+1)));

        }
    }
    erl1=erl1/((double)N1*N2);
    cout<<"error=" <<sqrt(erl1) <<endl ;  
    fftw_destroy_plan(p); fftw_destroy_plan(q); fftw_cleanup();

    return 0;
}