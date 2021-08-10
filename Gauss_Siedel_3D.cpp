#include <stdio.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <vector>
#include "Field.h"

using Field = Field_<double>;

using namespace std;

int main() {

    int N1=10;
    int N2=18;
	int N3=26;

	Field phi(N1,N2,N3);			//potential
	
	unsigned max_solver_it = 5000;	//maximum number of solver iterations
	double tolerance = 1e-5;		//solver tolerance
	
    double pi = 3.141592653589793;
    double L1  = 1.0;
    double dx = L1/(double)(N1-1);//+ instead of -1
    double L2= 2.0;
    double dy=L2/(double)(N2-1);
    double L3= 3.0;
    double dz=L3/(double)(N3-1);

	//precompute 1/(dx^2)
    //double3 dh = world.getDh();
    double idx2 = 1.0/(dx*dx);
    double idy2 = 1.0/(dy*dy);
    double idz2 = 1.0/(dz*dz);

    double L2N=0;			//norm
    bool converged= false;
	
    std::vector<double> b(N1*N2*N3,0.0);
    std::vector<double> X(N1,0.0);
    std::vector<double> Y(N2,0.0);
    std::vector<double> Z(N3,0.0);

	int u=0;
	int l=0;
	for(int i = 0;i <N1;i++){
        for(int j = 0;j<N2;j++){
			for(int k = 0;k<N3;k++){
				//X[i] =dx+(double)i*dx ;      
				//Y[j] =dy+ (double)j*dy ;
				//Z[k] =dz+ (double)k*dz ;
				X[i] =(double)i*dx ;      
				Y[j] =(double)j*dy ;
				Z[k] =(double)k*dz ;
				u = i*N2*N3 + j*N3 + k;
				b[u]=2*Z[k]*(L3-Z[k])*Y[j]*(L2-Y[j])+2*Z[k]*(L3-Z[k])*X[i]*(L1-X[i])+2*Y[j]*(L2-Y[j])*X[i]*(L1-X[i]);
				if (u != l) cout << "ERROR u != l " << endl;
				l=l+1;
			}
        }
    }

    /*solve potential*/
    for (unsigned it=0;it<max_solver_it;it++)
    {
		 for (int i=1;i<N1-1;i++)
            for (int j=1;j<N2-1;j++)
                for (int k=1;k<N3-1;k++)
                {
					u = i*N2*N3 + j*N3 + k;
					//standard internal open node
					double phi_new = (b[u] +
									idx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
									idy2*(phi[i][j-1][k]+phi[i][j+1][k]) +
									idz2*(phi[i][j][k-1]+phi[i][j][k+1]))/(2*idx2+2*idy2+2*idz2);

					/*SOR*/
					phi[i][j][k] = phi[i][j][k] + 1.4*(phi_new-phi[i][j][k]);
				}

		 /*check for convergence*/
		 if (it%25==0)
		 {
			double sum = 0;
			for (int i=1;i<N1-1;i++)
				for (int j=1;j<N2-1;j++)
					for (int k=1;k<N3-1;k++)
					{
						u = i*N2*N3 + j*N3 + k;
						double R = -phi[i][j][k]*(2*idx2+2*idy2+2*idz2) +
									b[u] +
									idx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
									idy2*(phi[i][j-1][k]+phi[i][j+1][k]) +
									idz2*(phi[i][j][k-1]+phi[i][j][k+1]);

						sum += R*R;
					}

			L2N = sqrt(sum/(N1*N2*N3));
			if (L2N<tolerance) {converged=true;break;}
		}
    }

    if (!converged) cerr<<"GS failed to converge, L2="<<L2<<endl;
	
	
    l=-1;
    double erl1 = 0.;
    for (int i = 1; i < N1-1; i++) {
        for(int j = 1; j < N2-1; j++){
			for(int k = 1; k < N3-1; k++){
				l=l+1;
				//X[i] =dx+(double)i*dx ;      
				//Y[j] =dy+ (double)j*dy ;
				//Z[k] =dz+ (double)k*dz ;
				X[i] =(double)i*dx ;      
				Y[j] =(double)j*dy ;
				Z[k] =(double)k*dz ;
				//double res=0.5/pi/pi*in1[l];
				double res=X[i]*(L1-X[i])*Y[j]*(L2-Y[j])*Z[k]*(L3-Z[k]);
				erl1 +=pow(fabs(res -  phi[i][j][k]),2); 
				printf("%3d %10.7g %10.7g\n", l, res, phi[i][j][k]);
				//erl1 +=pow(fabs(res-  3*0.125*out2[l]/((double)(N1+1))/((double)(N2+1))/((double)(N3+1))/((double)(N1+1))/((double)(N2+1))/((double)(N3+1))),2); 
				//printf("%3d %10.5g %10.5g\n", l, res,  3*0.125*out2[l]/((double)(N1+1))/((double)(N2+1))/((double)(N3+1))/((double)(N1+1))/((double)(N2+1))/((double)(N3+1)));
			}
        }
    }
    erl1=erl1/((double)N1*N2*N3);
    cout<<"error=" <<sqrt(erl1) <<endl ;  

    return converged;


}