#include "material_specific.h"
#include <omp.h>
/*
	*
	* ��������� ��������, ������ �� ������� ���� (���������� leftBound � rightBound - ��������
	* ������� �� �������� ���������������� ����������, ���� ��� �������� ����� ������ �����, �� �� ���������� ���� ������)
	*
*/
bool __attribute__((target(mic:0))) sign2(double leftBound, double rightBound);


/* 
	*
	* ������� ��� ���������� p1 �� ��������� energy(p1,psi) = energy_value �� ���������� (p1_low, p1_high),
	* �� ������� �������������� ���� ������
	* ������������ ��� ����������� ������ ������ ���� ��������� (��� ���������� ����������������� ������������)
	*
*/
double __attribute__((target(mic:0))) Newton_psi_energy(double p1_low, double p1_high, double psi, double energy_value, Params* params);

/*
	*
	* �������, ��������� ������ ������ ��������� energy(p1, psi)=energy_value �� ���������� (0, p1max).
	* ������������ ��� ���������� ����������������� ������������ 
	*
*/
double __attribute__((target(mic:0))) apply_Newton_psi_energy(double psi, bool &flag, double energy_value, double pmax, Params* params);


/*  
	* ���������� ����������� ��������� ��������� � ��������� p ��� ������
	* ���������� ������������ �� ������ ���������� ����������� ������������ ������ � ���������,
	* ����������� �� ���� ������ 
	* p - �������, ��� �������� ���� �����������. ������ ���� �������� � ������ ���� ���������.
	* px_max, py_max - ������� ������ ���� ��������� (�������������� ������������� ����)
	* Nx, Ny - ���������� ����� �� ��������������� ���� � �����
	* angle - ���� �������� ������� ��������� (����������� ������ ������� �������) ������������ 
	*         ��� ������������
	* points - ������ �������� �������� - ����� �����
	* Wer - ������ �������� �����������, ����������� �� ����� ��������� �����
*/
double __attribute__((target(mic:0))) getWer(double px, double py, double* px_mas, double* py_mas, double *Wer, Params* params);

/*
	*
	* ���������, �������� �� �������� ������� ������� ��������
	*
*/	
double __attribute__((target(mic:0))) simpson_function(double psi, double px, double py, Params* params);

/*
	*
	* ���������� ������ ��������
	*
*/
double  __attribute__((target(mic:0))) simpson(double px, double py, Params* params);

/*
	*
	* ���������� ���������� ����������� ���������
	*
*/
void __attribute__((target(mic:0))) full_probability_psi(double* px_mas, double *py_mas, double *res_mas, Params* params);


double __attribute__((target(mic:0))) distrib_function(double p, double psi, Params* params);


double A_norm(Params* params);


double ditrib_function_int_p(double phi, Params* params);

