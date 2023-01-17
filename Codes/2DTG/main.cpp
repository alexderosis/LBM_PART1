#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
using namespace std;
///----------------------------------------------------------------------------------
/// LB variables
const bool plot_vtk = true;
const int np = 9;
const vector<int> cx = {0, 1, 0, -1,  0, 1, -1, -1,  1},
									cy = {0, 0, 1,  0, -1, 1,  1, -1, -1};
const vector<double> wf = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
const int nx = 32, ny = nx;
const double cs2 = 1./3., cs4 = 1./9., rho0 = 1., v0 = 0.01, Reynolds = 1000., kappa = 2.*M_PI/((double)nx), ni = v0*ny/Reynolds, T_ref = 1./(2.*kappa*kappa*ni), tau = ni*3.+0.5, omega = 1./tau, omega1 = 1.-omega;
const int nsteps = (int)(1*T_ref)+1, n_out = nsteps-1;
vector<double> f1(nx*ny*np, 0.), f2(nx*ny*np, 0.), u(nx*ny, 0.), v(nx*ny, 0.), rho(nx*ny, 0.), temp_pop(np, 0.);
double U, V, R, ftemp, kinetic_energy;
int newx, newy, id, idn;
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
void write_fluid_vtk(int time)
{
	/// Create filename
	stringstream output_filename;
	output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());

	/// Write VTK header
	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << nx << " " << ny << " 1" << "\n";
	output_file << "X_COORDINATES " << nx << " float\n";
	for(int i = 0; i < nx; ++i)
		output_file << i << " ";
	output_file << "\n";
	output_file << "Y_COORDINATES " << ny  << " float\n";
	for(int j = 0; j < ny ; ++j)
		output_file << j  << " ";
	output_file << "\n";
	output_file << "Z_COORDINATES " << 1 << " float\n";
	output_file << 0 << "\n";
	output_file << "POINT_DATA " << (nx) * (ny) << "\n";

	/// Write density difference
	output_file << "SCALARS density float 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
    {
      id = X*ny+Y;
      output_file << rho[id]-rho0<< "\n";
    }

	/// Write velocity
	output_file << "VECTORS velocity_vector float\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
    {
      id = X*ny+Y;
			output_file << u[id] << " " << v[id] << " 0\n";
    }

	/// Close file
	output_file.close();
}
///----------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
      id = x*ny+y;
			R = rho[id] = rho0/cs2*(1.-v0*v0/4./cs2*(cos(2.*kappa*x)+cos(2.*kappa*y)));
			U = u[id] = -v0*cos(kappa*x)*sin(kappa*y);
      V = v[id] = v0*sin(kappa*x)*cos(kappa*y);
			for(int k=0; k<np;k++)
        f1[id*np+k] = wf[k]*R*(1.+3.*(U*cx[k]+V*cy[k])+4.5*pow(U*cx[k]+V*cy[k],2)-1.5*(U*U+V*V));
		}
}
///----------------------------------------------------------------------------------------------------------------------------------
void algoLB()
{
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
      id = x*ny+y;
			U = V = R = 0.;
			for(int k=0; k<np; k++)
			{
				ftemp = temp_pop[k] = f1[id*np+k];
				R += ftemp;
				U += ftemp*cx[k];
				V += ftemp*cy[k];
			}
			rho[id] = R;
			U /= R;
			V /= R;
      u[id] = U;
			v[id] = V;
			for(int k=0; k<np; k++)
			{
				f1[id*np+k] = omega1*f1[id*np+k] + omega*wf[k]*R*(1.+3.*(U*cx[k]+V*cy[k])+4.5*pow(U*cx[k]+V*cy[k],2)-1.5*(U*U+V*V));
				newx = x+cx[k];
				newy = y+cy[k];
				if(x==0 || x==nx-1)
					newx = (newx+nx)%nx;
				if(y==0 || y==ny-1)
					newy = (newy+ny)%ny;
        idn = newx*ny+newy;
				f2[idn*np+k] = f1[id*np+k];
			}
		}
}
///----------------------------------------------------------------------------------------------------------------------------------
void compute_error()
{
	int time = nsteps-1;
	double u_an, v_an, vel_an, vel, num = 0., den = 0.;
  for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
    {
    	id = x*ny+y;
    	u_an = -v0*cos(kappa*x)*sin(kappa*y)*exp(-2.*kappa*kappa*time*ni);
    	v_an = v0*sin(kappa*x)*cos(kappa*y)*exp(-2.*kappa*kappa*time*ni);
    	vel_an = sqrt(u_an*u_an+v_an*v_an);
    	vel = sqrt(u[id]*u[id]+v[id]*v[id]);
    	num += pow(vel_an-vel,2);
    	den += pow(vel,2);
    }
   double error = sqrt(num/den);
   printf("Error = %e\n", error);

}
///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	//FILE *data = fopen("data.txt","wt");
	//system("mkdir vtk_fluid");
	initial_state();
	clock_t c_start = std::clock();
	for(int i=0; i<nsteps; i++)
  {
		algoLB();
    f1 = f2;
		if(plot_vtk==true && i%n_out==0)
			write_fluid_vtk(i);
		//printf("t %lf of %lf\n", (double)i/T_ref, (double)nsteps/T_ref);
    //if(check_mach==1)
      //goto labelA;
  }
  compute_error();
  //labelA:
  //clock_t c_end = std::clock();
  //double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
  //cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  //fclose(data);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
