/* FDTD Polar Grid
   Solving Schrodinger Equation for a particle in 
   2D circular container with V = 0.
   
   dpsi/dt = (1/2) laplacian psi
   Boundary condition: psi(r=R, theta) = 0
   Periodic condition: psi(r, theta + 2pi) = psi(r, theta) 
   I Wayan Sudiarta

   Updated: 
   15/07/2021 - write simple code
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NR 50  // Number of grids for var r
#define NA 50  // var theta = alpha  
#define NT 100000

int main()
{
	int i, j, n;
	double psi[NR+1][NA+1];
	double psinew[NR+1][NA+1];
	double r[NR+1], r2[NR+1];
	double R;
	double pi = 2.0*asin(1.0);
	double psic, psicnew, temp;
	double dr, da, dt;
	double dtdr2, dtdr, dtda2;
    int k1, k2, k3, k4;	
	
	FILE *outfile;
	
	// indices for diff at center
	k1 = 0;
	k2 = (int)NA/4;
	k3 = 2*k2;
	k4 = 3*k2;
	
	R  = 1.0;
	dr = R/NR;
	dt = 0.01*dr*dr;

	da = 2.0*pi/NA;
	dtdr2 = 0.5*dt/(dr*dr);
	dtdr  = 0.5*dt/(2*dr);
	dtda2 = 0.5*dt/(da*da);
	
	// fill array r[i] and r2[i]
    for(i = 0; i <= NR; i++){
		temp  = i*dr;
        r[i]  = temp;
		r2[i] = temp*temp;
	}
	
    // STEP 1: Initialize //
    // initialize array to zero
	for(i = 0; i <= NR; i++){
		for(j = 0; j <= NA; j++){
            psi[i][j] = 0.0;
            psinew[i][j] = 0.0;			
	    }
	}
    // initial psi with random values
	srand(12345);
	for(i = 1; i < NR; i++){
		for(j = 0; j < NA; j++){
            psi[i][j] = rand()/(double)RAND_MAX;
	    }
	}
    // at the center, i = 0
    psic = rand()/RAND_MAX;    
	for(j = 0; j < NA; j++){
        psi[0][j] = psic;
	}

	outfile = fopen("initial.dat", "w");
	for(i = 0; i <= NR; i++){
		for(j = 0; j <= NA; j++){
			fprintf(outfile, "%lf  ", psi[i][j]);
		}
		fprintf(outfile, "\n");
	}


    // iteration
    for(n = 1; n <= NT; n++){
		
		//STEP 2: Update psi //
		for(i = 1; i < NR; i++){
			for(j = 0; j < NA; j++){
				psinew[i][j] = psi[i][j];
				psinew[i][j] += dtdr2*(psi[i-1][j] - 2*psi[i][j] + psi[i+1][j]);
				psinew[i][j] += dtdr*(psi[i+1][j] - psi[i-1][j])/r[i];
				psinew[i][j] += dtda2*(psi[i][j-1] - 2*psi[i][j] + psi[i][j+1])/r2[i];
			}
		}
		// index j = 0
		for(i = 1; i < NR; i++){
			psinew[i][0] = psi[i][0];
			psinew[i][0] += dtdr2*(psi[i-1][0] - 2*psi[i][0] + psi[i+1][0]);
			psinew[i][0] += dtdr*(psi[i+1][0] - psi[i-1][0])/r[i];
			psinew[i][0] += dtda2*(psi[i][NA-1] - 2*psi[i][0] + psi[i][1])/r2[i];
		}
		// at center i = 0
		psicnew = psic;
        psicnew	+= dtdr2*(psi[1][k1] - 2*psic + psi[1][k3]);	
        psicnew	+= dtdr2*(psi[1][k2] - 2*psic + psi[1][k4]);	
		for(j = 0; j <= NA; j++){
			psinew[0][j] = psicnew;
		}
		// periodic boundary psi[i][0] = psi[i][NR]
		for(i = 0; i < NR; i++){
			psinew[i][NA] = psinew[i][0];
		}

		// STEP 3: save for next iteration
		for(i = 0; i < NR; i++){
			for(j = 0; j <= NA; j++){
				psi[i][j] = psinew[i][j];
			}
		}
		// center
		psic = psicnew;
		
		// STEP 4: Compute Energy //
		
	}
	
	// STEP 5: Save Results
	outfile = fopen("results.dat", "w");
	for(i = 0; i <= NR; i++){
		for(j = 0; j <= NA; j++){
			fprintf(outfile, "%le ", psi[i][j]);
		}
		fprintf(outfile, "\n");
	}
	
    return 0;	
}