#include "material_specific.h"
/* 
	*
	* Выражение для энергетического спектра (в декартовых координатах)
	*
*/
double  energy(double px, double py, Params* params)
{
	return (1-cos(px)*cos(py));
}
/* 
	*
	* Выражение для энергетического спектра (в полярных координатах)
	*
*/
double  energy_psi(double p, double psi, Params* params)
{
	return (1-cos(p*cos(psi))*cos(p*sin(psi)));
}

/*
	*
	* Производная энергии по модулю импульса (в полярных координатах)
	*
*/
double  d_energy_psi(double p, double psi, Params* params)
{
	return (sin(p*cos(psi))*cos(psi)*cos(p*sin(psi))+cos(p*cos(psi))*sin(p*sin(psi))*sin(psi));
}

/*
	*
	* Компоненты скорости
	*
*/
double  vx_fun(double px, double py, Params* params)
{
	return (*params).vx0*sin(px)*cos(py);
}

double  vy_fun(double px, double py, Params* params)
{
	return (*params).vy0*cos(px)*sin(py);
}

/*
	*
	* Правые части уравнений движения
	*
*/

double  right_x(double px, double py, double t, Params* params)
{
	//return ((*params).Exc + (*params).Ex*cos(t*(*params).wx)+(*params).H*vy_fun(px,py,params));
	return (*params).Exc;
}

double  right_y(double px, double py, double t, Params* params)
{
	//return ((*params).Eyc + (*params).Ey*cos(t*(*params).wy+(*params).phi) - (*params).H*vx_fun(px,py,params));
	return (*params).Ey1*cos(t*(*params).wy1)+(*params).Ey2*cos(t*(*params).wy2+(*params).phi);
}

/*
	*
	* Границы первой зоны Бриллюэна
	*
*/
double  pmax(double psi, Params* params)
{
	double res;
	if((psi >= 0) && (psi < M_PI_2))
		res = M_PI/(1 + tan(psi))/cos(psi);
	if(((psi >= M_PI_2) && (psi < M_PI)))
		res = - M_PI/(1 - tan(psi))/cos(psi);
	if(((psi >= M_PI) && (psi < M_PI*3/2)))
		res = - M_PI/(1 + tan(psi))/cos(psi);
	if((psi >= 3*M_PI/2) && (psi < 2*M_PI))
		res = M_PI/(1 - tan(psi))/cos(psi);
	return res;
}

/*
	*
	* Функция, приводящая квазиимпульс к первой зоне Бриллюэна
	*
*/
Point  ToFirstBand(Point p0, Params* params)
{
	double px = p0.x;
	double py = p0.y;
	double p = sqrt(px*px + py*py);
	double phi = atan2(py,px);
	double OK = M_PI/sqrt(2.0); // орт e1 (модуль этого орта)
	double OM = M_PI/sqrt(2.0); // орт e2 (модуль этого орта)
	double OK1 = p*cos(phi - M_PI/4); // проецируем исходный радиус-вектор на орт е1 
	double OM1 = p*sin(phi - M_PI/4); // проецируем исходный радиус-вектор на орт е2
	double pxi = floor(OK1/OK); // целая часть от деления проекции исходного радиус-вектора
								// на орт е1 на модуль этого вектора
	double pyi = floor(OM1/OM); // целая часть от деления проекции исходного радиус-вектора
								// на орт е2 на модуль этого вектора
	double pxf = OK1 - pxi*OK;  // дробная часть от деления проекции исходного радиус-вектора
								// на орт е1 на модуль этого вектора
	double pyf = OM1 - pyi*OM;  // дробная часть от деления проекции исходного радиус-вектора
								// на орт е1 на модуль этого вектора

	double pxN, pyN;
	if(fabs(pxi-floor(pxi/2)*2)<0.000001) // если pxi - четная
		pxN = pxf;
	else
		pxN = -(OK - pxf);

	if(fabs(pyi-floor(pyi/2)*2)<0.000001) // если pyi - четная
		pyN = pyf;
	else
		pyN = -(OM - pyf);
	double pxN1 = pxN/sqrt(2.0) - pyN/sqrt(2.0);
	double pyN1 = pxN/sqrt(2.0) + pyN/sqrt(2.0);
	Point p_res;
	p_res.x = pxN1;
	p_res.y = pyN1;
  	return p_res;
}

/*
	*
	* Возвращает точку с индексами (k,m)
	*
*/
Point  point_k_m(int k, int m, Params* params)
{
	double Ax = (*params).Ax; double Ay = (*params).Ay;
	double Bx = (*params).Bx; double By = (*params).By;
	double Dx = (*params).Dx; double Dy = (*params).Dy;

	double Kx = Ax + (Dx - Ax)*k/(*params).Nx;
	double Ky = Ay + (Dy - Ay)*k/(*params).Nx;
	double Mx = Ax + (Bx - Ax)*m/(*params).Ny;
	double My = Ay + (By - Ay)*m/(*params).Ny;
	double Ox = ((By - Ay)/(Bx - Ax)*Kx - (Dy - Ay)/(Dx - Ax)*Mx + My - Ky)/((By - Ay)/(Bx - Ax) - (Dy - Ay)/(Dx - Ax));
	double Oy = (By - Ay)/(Bx - Ax)*(Ox - Kx) + Ky;
	Point res;
	res.x = Ox;
	res.y = Oy;
	return res;
}

/*
	*
	* Массив координат точек в импульсном пространстве
	*
*/
void  points_mas(double* px_mas, double* py_mas, Params* params)
{
	
	int Nx = (*params).Nx;
	int Ny = (*params).Ny;
	Point temp_p;
	for(int i=0; i<Nx+1; i++)
		for(int j=0; j<Ny+1; j++)
		{
			temp_p = point_k_m(i, j, params);
			px_mas[j+i*(Ny+1)] = temp_p.x;
			py_mas[j+i*(Ny+1)] = temp_p.y;
		};
}

void set_init_params(Params *params)
{
	(*params).beta = 0;   
	(*params).wla_max = 0; 
	(*params).wlo_max = 0;
	(*params).Exc = 0; 
	(*params).Ex = 0;     
	(*params).Eyc = 0;     
	(*params).Ey  = 0; 
	(*params).Ey1 = 0;
	(*params).Ey2 = 0;
	
	(*params).H = 0;      
	(*params).vx0 = 0;     
	(*params).vy0 = 0;    
	(*params).wx = 0;  
    (*params).wy = 0;
	(*params).wy1 = 1.0;
	(*params).wy2 = 2.0;
	
	(*params).phi = 0;	
	(*params).Nx = 0;
	(*params).Ny = 0;
	(*params).px_max = 0;
	(*params).py_max = 0;
	(*params).Ax = 0;  
	(*params).Ay = 0;	
	(*params).Bx = 0;
	(*params).By = 0;
	(*params).Dx = 0;
	(*params).Dy = 0;
	(*params).Newton_abs_error = 0;
	(*params).Newton_n_points = 0;
	(*params).Simson_n = 0;
	(*params).dt = 0.0001;
	(*params).all_time = 0;
	(*params).num_threads_openmp = 0;
	(*params).n_part = 0;
	(*params).Anorm = 0;
	(*params).max_prob = 1;
}
