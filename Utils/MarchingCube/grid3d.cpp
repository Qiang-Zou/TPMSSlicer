/*
	File: grid3d.cpp
	Programmed by: David Vega: dvega@uc.edu.ve
	August 2019
	August 2020
*/

#ifdef _MSC_VER
#pragma warning( disable : 4244 )
#endif

#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cstring>

#include "MC33.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

#ifndef GRD_orthogonal
void setIdentMat3x3d(double (*A)[3])
{
	for (double *d = A[0] + 8; --d != A[0];)
		d[0] = 0.0;
	for (int i = 0; i != 3; ++i)
		A[i][i] = 1.0;
}
#endif

//******************************************************************
const unsigned int* grid3d::get_N() {return N;}
const float* grid3d::get_L() {return L;}
const float* grid3d::get_Ang() {return Ang;}
const double* grid3d::get_r0() {return r0;}
const double* grid3d::get_d() {return d;}
GRD_data_type grid3d::get_grid_value(unsigned int i, unsigned int j, unsigned int k)
{
	if (F && i <= N[0] && j <= N[1] && k <= N[2])
		return F[k][j][i];
	return 0;
}
#ifndef GRD_orthogonal
const double (*grid3d::get__A())[3] {return _A;}
const double (*grid3d::get_A_())[3] {return A_;}
int grid3d::isnotorthogonal() {return nonortho;}
#endif
void grid3d::set_ratio_aspect(double rx, double ry, double rz)
{
	if (F)
	{
		d[0] = rx; d[1] = ry; d[2] = rz;
		for (int i = 0; i != 3; ++i)
			L[i] = N[i]*d[i];
	}
}

void grid3d::set_r0(double x, double y, double z)
{
	r0[0] = x; r0[1] = y; r0[2] = z;
}

void grid3d::set_Ang(float angle_bc, float angle_ca, float angle_ab)
{
	Ang[0] = angle_bc; Ang[1] = angle_ca; Ang[2] = angle_ab;
#ifndef GRD_orthogonal
	if (Ang[0] != 90 || Ang[1] != 90 || Ang[2] != 90)
	{
		nonortho = 1;
		double ca = cos(Ang[0]*(M_PI/180.0));
		double cb = cos(Ang[1]*(M_PI/180.0));
		double aux1 = Ang[2]*(M_PI/180.0);
		double sg = sin(aux1);
		double cg = cos(aux1);
		aux1 = ca - cb*cg;
		double aux2 = sqrt(sg*sg + 2*ca*cb*cg - ca*ca - cb*cb);
		_A[0][0] = A_[0][0] = 1.0;
		_A[0][1] = cg;
		_A[0][2] = cb;
		_A[1][1] = sg;
		A_[1][1] = cb = 1.0/sg;
		A_[0][1] = -cg*cb;
		_A[1][2] = aux1*cb;
		_A[2][2] = aux2*cb;
		aux2 = 1.0/aux2;
		A_[0][2] = (cg*aux1 - ca*sg*sg)*cb*aux2;
		A_[1][2] = -aux1*cb*aux2;
		A_[2][2] = sg*aux2;
		_A[1][0] = _A[2][0] = _A[2][1] = 0.0;
		A_[1][0] = A_[2][0] = A_[2][1] = 0.0;
	}
	else
	{
		nonortho = 0;
		setIdentMat3x3d(_A);
		setIdentMat3x3d(A_);
	}
#endif
}

void grid3d::set_grid_value(unsigned int i, unsigned int j, unsigned int k, GRD_data_type value)
{
	if (internal_data && F && i <= N[0] && j <= N[1] && k <= N[2])
		F[k][j][i] = value;
}
//******************************************************************
void grid3d::free_F()
{
	if (F)
	{
		for (unsigned int k = 0; k <= N[2]; ++k)
		{
			if (F[k] && internal_data)
			{
				for (unsigned int j = 0; j <= N[1]; ++j)
					free(F[k][j]);
			}
			else
			{
				k = N[2];
				break;
			}
			free(F[k]);
		}
		free(F);
		F = 0;
	}
}

grid3d::grid3d() : F(0) {}

grid3d::~grid3d()
{
	free_F();
}

int grid3d::alloc_F()
{
	unsigned int j, k;
	F = reinterpret_cast<GRD_data_type***>(malloc((N[2] + 1)*sizeof(GRD_data_type**)));
	if (!F)
		return -1;
	for (k = 0; k <= N[2]; ++k)
	{
		F[k] = reinterpret_cast<GRD_data_type**>(malloc((N[1] + 1)*sizeof(GRD_data_type*)));
		if (!F[k])
			return -1;
		for (j = 0; j <= N[1]; ++j)
		{
			F[k][j] = reinterpret_cast<GRD_data_type*>(malloc((N[0] + 1)*sizeof(GRD_data_type)));
			if (!F[k][j])
			{
				while (j)
					free(F[k][--j]);
				free(F[k]);
				F[k] = 0;
				return -1;
			}
		}
	}
	internal_data = 1;
	return 0;
}

int grid3d::set_grid_dimensions(unsigned int Nx, unsigned int Ny, unsigned int Nz)
{
	free_F();
	N[0] = Nx - 1;
	N[1] = Ny - 1;
	N[2] = Nz - 1;
	if (alloc_F())
	{
		free_F();
		return -1;
	}
	for (int i = 0; i != 3; ++i)
	{
		L[i] = N[i];
		d[i] = 1.0;
		r0[i] = 0.0;
#ifndef GRD_orthogonal
		Ang[i] = 90.0f;
#endif
	}
#ifndef GRD_orthogonal
	nonortho = 0;
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return 0;
}

int grid3d::set_data_pointer(unsigned int Nx, unsigned int Ny, unsigned int Nz, GRD_data_type* data)
{
	free_F();
	F = reinterpret_cast<GRD_data_type***>(malloc(Nz*sizeof(GRD_data_type**)));
	if (!F)
		return -1;
	N[0] = Nx - 1;
	N[1] = Ny - 1;
	N[2] = Nz - 1;
	for (unsigned int k = 0; k < Nz; ++k)
	{
		F[k] = reinterpret_cast<GRD_data_type**>(malloc(Ny*sizeof(GRD_data_type*)));
		if (!F[k])
		{
			while (k)
				free(F[--k]);
			free(F);
			F = 0;
			return -1;
		}
		for (unsigned int j = 0; j < Ny; ++j)
			F[k][j] = data + j*Nx;
		data += Ny*Nx;
	}
	for (int i = 0; i != 3; ++i)
	{
		L[i] = N[i];
		d[i] = 1.0;
		r0[i] = 0.0;
#ifndef GRD_orthogonal
		Ang[i] = 90.0f;
#endif
	}
#ifndef GRD_orthogonal
	nonortho = 0;
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	internal_data = 0;
	return 0;
}

//******************************************************************
/*
read_grd reads a filename file (the file must be a output *.grd file from the
DMol program), it returns a pointer to struct _GRD that contains all the grid
data.
*/
int grid3d::read_grd(const char *filename)
{
	char cs[32];
	unsigned int  i, j, k;
	int xi[3], order;

	ifstream in(filename);
	if (!in)
		return -1;
	free_F();
#ifndef GRD_orthogonal
	nonortho = 0;
#endif
	periodic = 0;
	in.getline(title, 159);
	in.ignore(60, '\n');
	in >> L[0] >> L[1] >> L[2] >> Ang[0] >> Ang[1] >> Ang[2] >> N[0] >> N[1] >> N[2];
	in >> order >> xi[0] >> cs >> xi[1] >> cs >> xi[2];
	in.ignore(20, '\n');
	if (N[0] < 2 || N[1] < 2 || N[2] < 2) return -1;
	if (order != 1 && order != 3) return -1;
	for (i = 0; i != 3; ++i)
	{
		d[i] = L[i]/N[i];
		r0[i] = xi[i]*d[i];
	}
	periodic = (xi[0] == 0)|((xi[1] == 0)<<1)|((xi[2] == 0)<<2);

#ifndef GRD_orthogonal
	if (Ang[0] != 90 || Ang[1] != 90 || Ang[2] != 90)
	{
		nonortho = 1;
		double ca = cos(Ang[0]*(M_PI/180.0));
		double cb = cos(Ang[1]*(M_PI/180.0));
		double aux1 = Ang[2]*(M_PI/180.0);
		double sg = sin(aux1);
		double cg = cos(aux1);
		aux1 = ca - cb*cg;
		double aux2 = sqrt(sg*sg + 2*ca*cb*cg - ca*ca - cb*cb);
		_A[0][0] = A_[0][0] = 1.0;
		_A[0][1] = cg;
		_A[0][2] = cb;
		_A[1][1] = sg;
		A_[1][1] = cb = 1.0/sg;
		A_[0][1] = -cg*cb;
		_A[1][2] = aux1*cb;
		_A[2][2] = aux2*cb;
		aux2 = 1.0/aux2;
		A_[0][2] = (cg*aux1 - ca*sg*sg)*cb*aux2;
		A_[1][2] = -aux1*cb*aux2;
		A_[2][2] = sg*aux2;
		_A[1][0] = _A[2][0] = _A[2][1] = 0.0;
		A_[1][0] = A_[2][0] = A_[2][1] = 0.0;
	}
	else
	{
		nonortho = 0;
		setIdentMat3x3d(_A);
		setIdentMat3x3d(A_);
	}
#endif
	if (alloc_F())
	{
		free_F();
		return -2;
	}
	for (k = 0; k <= N[2]; ++k)
		if (order == 1)
		{
			for (j = 0; j <= N[1]; ++j)
				for (i = 0; i <= N[0]; ++i)
				{
					in.getline(cs, 31);
					F[k][j][i] = stof(cs);
				}
		}
		else
		{
			for (i = 0; i <= N[0]; ++i)
				for (j = 0; j <= N[1]; ++j)
				{
					in.getline(cs, 31);
					F[k][j][i] = stof(cs);
				}
		}
	return (in.good()? 0: -4);
}

/*
internal binary format
*/
int grid3d::read_grd_binary(const char* filename)
{
	unsigned int i, j, k;
	ifstream in(filename, ios::binary);
	if (!in)
		return -1;
	free_F();
	//char identifier[5] = "_GRD";
#define grb_identifier 0x4452475f
	in.read(reinterpret_cast<char*>(&i),sizeof(int));
//	if (i != *(int*)identifier)
	if (i != grb_identifier)
		return -1;

	in.read(reinterpret_cast<char*>(&i),sizeof(int));
	if (i > 159)
		return -1;
	in.read(title,i*sizeof(char));
//	title[i] = '\0';
	in.read(reinterpret_cast<char*>(N),sizeof N);
	in.read(reinterpret_cast<char*>(L),sizeof L);
	in.read(reinterpret_cast<char*>(r0),sizeof r0);
	in.read(reinterpret_cast<char*>(d),sizeof d);
#ifndef GRD_orthogonal
	in.read(reinterpret_cast<char*>(&nonortho),sizeof(int));
	if (nonortho)
	{
		in.read(reinterpret_cast<char*>(Ang),3*sizeof(float));
		in.read(reinterpret_cast<char*>(_A),sizeof _A);
		in.read(reinterpret_cast<char*>(A_),sizeof A_);
		mult_Abf = _multA_bf;
	}
	else
	{
		setIdentMat3x3d(_A);
		setIdentMat3x3d(A_);
	}
#endif
	//periodic = 0;
	if (r0[0] == 0 && r0[1] == 0 && r0[2] == 0)
		periodic = 1;
	if (alloc_F())
	{
		free_F();
		return -2;
	}
	for (k = 0; k <= N[2]; ++k)
		for (j = 0; j <= N[1]; ++j)
			in.read(reinterpret_cast<char*>(F[k][j]),(N[0] + 1)*sizeof(GRD_data_type));
	return (in.good()? 0: -4);
}
//******************************************************************

/*
Reads a set of files that contain a slab of res*res scan data points, the data
points are read as unsigned short int (if order is different from 0, the bytes
of the unsigned short are exchanged). The filename must end with a number, and
the fuction read all files with end number greater or equal to filename.
(Some datasets: http://www.graphics.stanford.edu/data/voldata/voldata.html)
*/
int grid3d::read_scanfiles(const char *filename, unsigned int res, int order)
{
	string nm(filename);
	ifstream in;
	unsigned int i, j, l;
	int m, k = -1;
	unsigned short int n;
	free_F();
#ifndef GRD_orthogonal
	nonortho = 0;
#endif
	periodic = 0;
	r0[0] = r0[1] = r0[2] = 0.0;
	for (i = 0; i != 2; ++i)
	{
		d[i] = 1.0;
		N[i] = res - 1;
		L[i] = float(res - 1);
	}
	d[2] = 1.0;
	l = int(nm.size() - 1);
	while (nm[l] >= '0' && nm[l] <= '9')
		l--;
	m = stoi(nm.substr(++l));
	while (1)
	{
		nm.replace(l,6,to_string(m++));
		in.open(nm, ios::binary|ifstream::in);
		if (!in)
		{
			m = 0;
			break;
		}
		if (!((++k)&63))
		{
			void *pt = realloc(F,(k + 64)*sizeof(void*));
			if (!pt)
			{
				m = -2;
				break;
			}
			F = reinterpret_cast<GRD_data_type***>(pt);
		}
		F[k] = reinterpret_cast<GRD_data_type**>(malloc(res*sizeof(void*)));
		if (!F[k])
		{
			m = -2;
			break;
		}
		for (j = 0; j != res; ++j)
			F[k][j] = reinterpret_cast<GRD_data_type*>(malloc(res*sizeof(GRD_data_type)));
		if (!F[k][N[1]])
		{
			m = -2;
			break;
		}

		if (order)
			for (j = 0; j != res; ++j)
				for (i = 0; i != res; ++i)
				{
					in.read(reinterpret_cast<char*>(&n),sizeof(short int));
					F[k][j][i] = static_cast<unsigned short int>((n>>8)|(n<<8));
				}
		else
			for (j = 0; j != res; ++j)
				for (i = 0; i != res; ++i)
				{
					in.read(reinterpret_cast<char*>(&n),sizeof(short int));
					F[k][j][i] = n;
				}
		if (in.fail())
		{
			m = -4;
			break;
		}
		in.close();
	}
	N[2] = k;
	L[2] = float(k);
	internal_data = 1;
	if (m == -2)
		free_F();
	else
	{
		j = k>>1;
		for (i = 0; i != j; ++i)
			swap(F[i], F[k - i]);
	}
#ifndef GRD_orthogonal
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return m;
}

//******************************************************************
/*
read_raw_file reads a file that contains integer (8, 16 or 32 bits) or float (or double)
data points. byte is the number of bytes of the integer (1, 2 or 4), or of the float (4
or 8). If the data is big endian, byte must be negative (only for integers). The vector
n[3] contains the number of points in each dimension. The size of file must be
abs(byte)*n[0]*n[1]*n[2]. The function returns zero when succeeds.
*/
int grid3d::read_raw_file(const char *filename, unsigned int *n, int byte, int isfloat)
{
	unsigned int i, j, k;
	ifstream in(filename, ios::binary);
	if (!in)
		return -1;
	unsigned int ui = 0;
	if (isfloat)
	{
		if (byte != sizeof(float) && byte != sizeof(double))
			return -1;
	}
	else if (abs(byte) > 4 || abs(byte) == 3 || !byte)
		return -1;
	if (byte == -1) byte = 1;
	free_F();
#ifndef GRD_orthogonal
	nonortho = 0;
#endif
	periodic = 0;
	r0[0] = r0[1] = r0[2] = 0.0;
	for (int h = 0; h != 3; ++h)
	{
		d[h] = 1.0;
		N[h] = n[h] - 1;
		L[h] = float(N[h]);
	}
	if (alloc_F())
	{
		free_F();
		return -2;
	}
#if defined(integer_GRD)
	if (!isfloat && size_type_GRD == byte)
#else
	if (isfloat && size_type_GRD == byte)
#endif
	{
		byte *= n[0];
		for (k = 0; k != n[2]; ++k)
			for (j = 0; j != n[1]; ++j)
				in.read(reinterpret_cast<char*>(F[k][j]),byte);
	}
	else if (isfloat)
	{
		if (byte == 8)
		{
#if defined(integer_GRD) || size_type_GRD == 4
			double df;
			for (k = 0; k != n[2]; ++k)
				for (j = 0; j != n[1]; ++j)
					for (i = 0; i != n[0]; ++i)
					{
						in.read(reinterpret_cast<char*>(&df),byte);
						F[k][j][i] = GRD_data_type(df);
					}
#endif
		}
		else
		{
#if defined(integer_GRD) || size_type_GRD == 8
			float f;
			for (k = 0; k != n[2]; ++k)
				for (j = 0; j != n[1]; ++j)
					for (i = 0; i != n[0]; ++i)
					{
						in.read(reinterpret_cast<char*>(&f),byte);
						F[k][j][i] = GRD_data_type(f);
					}
#endif
		}
	}
	else if (byte < 0)
	{
		for (k = 0; k != n[2]; ++k)
			for (j = 0; j != n[1]; ++j)
				for (i = 0; i != n[0]; ++i)
				{
					in.read(reinterpret_cast<char*>(&ui), -byte);
					if (byte == -2)
						F[k][j][i] = GRD_data_type((ui>>8)|(ui<<8));
					else if (byte == -4)
						F[k][j][i] = GRD_data_type((ui>>24)|((ui>>8)&0x00f0)|((ui<<8)&0x0f00)|(ui<<24));
					else
						F[k][j][i] = GRD_data_type((ui>>16)|(ui&0x00f0)|(ui<<16));
				}
	}
	else
	{
		for (k = 0; k != n[2]; ++k)
			for (j = 0; j != n[1]; ++j)
				for (i = 0; i != n[0]; ++i)
				{
					in.read(reinterpret_cast<char*>(&ui),byte);
					F[k][j][i] = GRD_data_type(ui);
				}
	}
#ifndef GRD_orthogonal
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return (in.good()? 0: -4);
}

//******************************************************************
/*
Reads a dat file:
http://www.cg.tuwien.ac.at/research/vis/datasets/
*/
int grid3d::read_dat_file(const char *filename)
{
	unsigned short int n, nx, ny, nz;
	ifstream in(filename, ios::binary);
	if (!in)
		return -1;
	free_F();
#ifndef GRD_orthogonal
	nonortho = 0;
#endif
	periodic = 0;
	r0[0] = r0[1] = r0[2] = 0.0;
	for (int h = 0; h != 3; ++h)
		d[h] = 1.0;
	in.read(reinterpret_cast<char*>(&nx),sizeof(short int));
	in.read(reinterpret_cast<char*>(&ny),sizeof(short int));
	in.read(reinterpret_cast<char*>(&nz),sizeof(short int));
	N[0] = nx - 1;
	N[1] = ny - 1;
	N[2] = nz - 1;
	for (int h = 0; h != 3; ++h)
		L[h] = float(N[h]);
	if (alloc_F())
	{
		free_F();
		return -2;
	}
	while (nz--)
		for (unsigned int j = 0; j < ny; ++j)
			for (unsigned int i = 0; i < nx; ++i)
			{
				in.read(reinterpret_cast<char*>(&n),sizeof(short int));
				F[nz][j][i] = n;
			}
#ifndef GRD_orthogonal
	setIdentMat3x3d(_A);
	setIdentMat3x3d(A_);
#endif
	return (in.good()? 0: -4);
}

//******************************************************************
/*
Save all point values of the grid in binary format
*/
int grid3d::save_raw_file(const char *filename)
{
	if (!F)
		return -1;
	ofstream out(filename, ios::binary);
	if (!out)
		return -1;
	for (unsigned int k = 0; k <= N[2]; ++k)
		for (unsigned int j = 0; j <= N[1]; ++j)
			out.write(reinterpret_cast<char*>(F[k][j]),(N[0] + 1)*sizeof(GRD_data_type));
	return (out.good()? 0: -4);
}

