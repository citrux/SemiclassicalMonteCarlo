#define _USE_MATH_DEFINES
#include <math.h>
/*
	*
	* ���������, �������� ���������� ��������
	*
*/
struct __attribute__((target(mic:0))) Point
{
	double x;
	double y;
};
/*
	*
	* ���������, �������� ��������� ������
	*
*/
struct __attribute__((target(mic:0))) Params
{
	double beta;    // ��������� ������� ����������� ������ � ������ ��������
	double wla_max; // ��������� ��������-��������� �������������� (������������ ������)
	double wlo_max; // ��������� ��������-��������� �������������� (���������� ������)

	double Exc;     // ���������� ������������ �������������� ���� ����� ��� Ox
	double Ex;      // ��������� ���������� ������������ �������������� ���� ����� ��� ��
	double Eyc;     // ���������� ������������ �������������� ���� ����� ��� O�
	double Ey;      // ��������� ���������� ������������ �������������� ���� ����� ��� �y
	double H;       // ������������� ���������� ����, ����������������� � ����������� ������� 
	double vx0;     // ����������� � ��������� ��� ���������� �������� (1)
	double vy0;     // ����������� � ��������� ��� ���������� �������� (1)
	double wx;      // ������� �������������� ���� ������� �����, �������������� ����� ��� �� 
	double wy;		// ������� �������������� ���� ������� �����, �������������� ����� ��� ��
	double wy1;
	double wy2;
	double phi;		// ������������� ���� �����, �������������� ����� ��� ��, �� ��������� �
					// ������, �������������� ����� ��� ��

	double Ey1;
	double Ey2;

	int Nx;         // ���������� �����, �� ������� ����������� ������ ���� ��������� 
	int Ny;	
	double px_max;  // ������������ �������� ��������� �������������
	double py_max;
	double Ax;      // ���������� ������ ���������������, ���������������
	double Ay;		// ������ ���� ���������
	double Bx;
	double By;
	double Dx;
	double Dy;

	double Newton_abs_error; // ���������� ������ ���������� ����� ��������� ��� ������������� ������ �������
	int Newton_n_points;	 // ���������� ���������, �� ������� ������� �������, �� ������� ���� 
							 // ������� ��������� ������� �������

	int Simson_n;    // ���������� ���������, �� ������� ����������� ������� [0..2*Pi]
					 // ��� �������������� ������� ��������

	double dt;       // ��� �� ������� ��� ������������� (����� ��� ������ �����-�����)
	double all_time; // ����� ����� �������������

	int num_threads_openmp; // ���������� ����, �� ������� ������������������ ������ ����������� ���������
							// � ������ �������� �������� ��������� ���� 
	int n_part; // ���������� ������ n_part

	double T;		//�����������
	double Anorm;
	
	double max_prob;
};

/*
	*
	* ���������, �������� ��������� ���������� �������
	*
*/
struct Probability
{
	int flag_count_or_load; // ���� 1, ��������� �������� ������������ �� ������, ����� - ���������
	char* filename_opt;
	char* filename_ac;
	char* filename_res;
};

/*
	*
	* ���������, ������������, ��������� ��� ��������� ����������� ���������
	*
*/
struct Graphic
{
	double low; // ����������� �������� ����������, �� ������� ������ ������
	double high; // ������������ �������� ����������, �� ������� ������ ������
	double step; // ��� ��������� ����������, � ����������� �� ������� ������ ������
	int num_var; // ����� ����������, � ����������� �� ������� ������ ������
};



void set_init_params(Params *params);


/*
	*
	* ��������� ��� ��������������� ������� (� ���������� �����������)
	*
*/
double __attribute__((target(mic:0))) energy(double px, double py, Params* params);
/*
	*
	* ��������� ��� ��������������� ������� (� �������� �����������)
	*
*/
double __attribute__((target(mic:0))) energy_psi(double p, double psi, Params* params);
/*
	*
	* ����������� ������� �� ������ �������� (� �������� �����������)
	*
*/
double __attribute__((target(mic:0))) d_energy_psi(double p, double psi, Params* params);

/*
	*
	* ���������� ��������
	*
*/
double __attribute__((target(mic:0))) vx_fun(double px, double py, Params* params);
double __attribute__((target(mic:0))) vy_fun(double px, double py, Params* params);

/*
	*
	* ������ ����� ��������� ��������
	*
*/
double __attribute__((target(mic:0))) right_x(double px, double py, double t, Params* params);
double __attribute__((target(mic:0))) right_y(double px, double py, double t, Params* params);

/*
	*
	* ������� ������ ���� ���������
	*
*/
double __attribute__((target(mic:0))) pmax(double psi, Params* params);

/*
	*
	* �������, ���������� ������������ � ������ ���� ���������
	*
*/
Point __attribute__((target(mic:0))) ToFirstBand(Point p0, Params* params);

/*
	*
	* ���������� ����� � ��������� (k,m)
	*
*/
Point __attribute__((target(mic:0))) point_k_m(int k, int m, Params* params);

/*
	*
	* ������ ��������� ����� � ���������� ������������
	*
*/
void __attribute__((target(mic:0))) points_mas(double* px_mas, double* py_mas, Params* params);