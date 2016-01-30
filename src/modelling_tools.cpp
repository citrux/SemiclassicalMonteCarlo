#include "modelling_tools.h"

double __attribute__((target(mic:0))) random_uniform(unsigned int& x1, unsigned int& y1, unsigned int& z1,  unsigned int& w1)
{
	unsigned int t1 = (x1 ^ (x1 << 11));
	x1 = y1;
	y1 = z1;
	z1 = w1;
	w1 = (w1 ^ (w1 >> 19)) ^ (t1 ^ (t1 >> 8));
	return ((double)w1) / ((double)UINT_MAX);
}

/*void RNGtest(double * mas, unsigned int * seed, int idx)
{
	unsigned int x1, y1, z1, w1; x1 = seed[idx]; y1 = 362436069; z1 = 521288629; w1 = 88675123;
	mas[idx] = random_uniform(x1, y1, z1, w1);
}*/

void __attribute__((target(mic:0))) runge(double &px, double &py, double t, Params* params)
{
	double px_1 = px;
	double py_1 = py;
	double kx1, kx2, kx3, kx4,
		  ky1, ky2, ky3, ky4;
	double dt = (*params).dt;

	kx1 = right_x(px_1, py_1, t, params);									ky1 = right_y(px_1, py_1, t, params);									
	kx2 = right_x(px_1 + dt/2 * kx1, py_1 + dt/2 * ky1, t + dt/2, params);	ky2 = right_y(px_1 + dt/2 * kx1, py_1 + dt/2 * ky1, t + dt/2, params);	
	kx3 = right_x(px_1 + dt/2 * kx2, py_1 + dt/2 * ky2, t + dt/2, params);	ky3 = right_y(px_1 + dt/2 * kx2, py_1 + dt/2 * ky2, t + dt/2, params);
	kx4 = right_x(px_1 + dt * kx3, py_1 + dt * ky3, t + dt, params);			ky4 = right_y(px_1 + dt * kx3, py_1 + dt * ky3, t + dt, params);

	px_1 += 1.0 / 6.0 * dt * (kx1 + 2 * kx2 + 2 * kx3 + kx4);
	py_1 += 1.0 / 6.0 * dt * (ky1 + 2 * ky2 + 2 * ky3 + ky4);
	px = px_1;
	py = py_1;
}

Point __attribute__((target(mic:0))) init_dist(unsigned int& x1, unsigned int& y1, unsigned int& z1, unsigned int& w1, Params* params)
{
	double psi;
	double p;
	double z;
	double f;
	Point point;
	do  
	{
		psi= 2*M_PI*random_uniform(x1, y1, z1, w1);
		p = pmax(psi, params)*random_uniform(x1, y1, z1, w1);
		z = random_uniform(x1, y1, z1, w1);
	    f = 1/((*params).Anorm)*distrib_function(p,psi,params);
	} while(z >= f);
	point.x = p*cos(psi);
	point.y = p*sin(psi);
	return point;
}

double Mean(double * arr, int count)
{
	double result = 0.0;

	for(int i = 0; i < count; i++)
	{
		result += arr[i];
	}

	return result / count;
}

double Mean(unsigned int * arr, int count)
{
	int result = 0;

	for(int i = 0; i < count; i++)
	{
		result += arr[i];
	}

	return ((double)result) / ((double)count);
}

double Std(double * arr, int count)
{
	double sum = 0.0;
	double mean = Mean(arr, count);

	for(int i = 0; i < count; i++)
	{
		sum += pow((arr[i] - mean), 2);
	}

	return sqrt(sum / (count - 1));
}

void __attribute__((target(mic:0))) jobKernel(
			   double * dev_average_value_x,
			   double * dev_average_value_y,
			   double * dev_average_time_array,
			   unsigned int  * dev_nAc,
			   unsigned int * dev_nOpt,

			   Params* params, double beta,
			   unsigned int rand_init_value, int idx,
			   
			   double *px_mas, double *py_mas,
			   double* res_ac, double* res_opt,
			   double* px_log, double* py_log,
			   int num_logs
			   )
{
	double* px_log_local = new double[num_logs];
	double* py_log_local = new double[num_logs];
	
	unsigned int x1, y1, z1, w1; x1 = rand_init_value; y1 = 362436069; z1 = 521288629; w1 = 88675123;// RNGInit
	
	Point p;
	bool flag;
	int iCount;
	double energy_value;
	double p_max;

	double vx = 0.0;
	double vy = 0.0;
	double value_x = 0.0;
	double value_y = 0.0;

	Point point = init_dist(x1, y1, z1, w1, params);
	double px = point.x;
	double py = point.y;
	double p1 = sqrt(px*px+py*py);
	double psi = atan2(py,px);
	
	double t = 0.0;
	double wsum = 0.0;
	double r1 = 0.0;
	double r2 = 0.0;
	double wla = 0.0;
	double wlo = 0.0;
	double r = - log(random_uniform(x1, y1, z1, w1));

	int n0 = 0;
	int nAc = 0;
	int nOpt = 0;

	Params params_opt, params_ac;
	params_opt = *params;
	params_opt.beta = beta;
	params_ac = *params;
	params_ac.beta = 0.0;

	int t_num = 0;

	while (t < (*params).all_time)
	{
		vx = vx_fun(px, py, params); // ��������� �������� ��������� ��������
		vy = vy_fun(px, py, params);
		

		value_x += vx * (*params).dt;
		value_y += vy * (*params).dt;

		runge(px, py, t, params); // ������ ��������� ��������

		/* �������� ������� � ���� */
		p.x = px;
		p.y = py;
		p = ToFirstBand(p, params);
		px = p.x;
		py = p.y;
        px_log_local[t_num] = px;
		py_log_local[t_num] = py;
		t += (*params).dt;
        
		// ��������� ����������� �������� � ���������� ��������� �� ������������ � ���������� �������
		wla = (*params).wla_max*getWer(p.x, p.y, px_mas, py_mas, res_ac, &params_ac);
		wlo = (*params).wlo_max*getWer(p.x, p.y, px_mas, py_mas, res_opt, &params_opt);
		
		r = -log(random_uniform(x1, y1, z1, w1));
		wsum += (wlo + wla) * (*params).dt;

		if(wsum > r)
		{
			n0++; // ���������� ������� ������ ����� ���������
			wsum = 0.0;
			r = random_uniform(x1, y1, z1, w1);
			if(wla > r * (wlo + wla))
			{
				nAc++; // ���������� ������� ��������� �� ������������ �������
				energy_value = energy_psi(sqrt(px*px+py*py),atan2(py,px),&params_ac);
				flag = false;
				iCount = 0;
				while(!flag)
				{
					psi = 2 * M_PI * random_uniform(x1, y1, z1, w1);   // ��������� ������� ����������� ���� �������������
					p_max = pmax(psi, &params_ac); // ������������ �������� ������ ������������� 
															  // � ����������� ���� psi 				
					p1 = apply_Newton_psi_energy(psi, flag, energy_value, p_max, &params_ac);
					// ���� p1 ����������, �� flag = true, � �� ��������� ��������� ���� ���������, ������� ������� �� �����
					
					// ���� �� 15 ������� �� ����� �������, ������� �� �����
					if(!flag)
						iCount++;
					if(iCount > 15 && !flag)
						flag = true;
				}
				px = p1*cos(psi); // ��������� ���������� ��������
				py = p1*sin(psi); // ��������� ���������� ��������
			}
			else
			{
				if((wlo>0.0001)&& (energy_psi(sqrt(px*px+py*py),atan2(py,px),&params_opt)>=params_opt.beta))
				{
					nOpt++; // ���������� ������� ��������� �� ���������� �������
					energy_value = energy_psi(sqrt(px*px+py*py),atan2(py,px),&params_opt)-params_opt.beta;
					flag = false;
					iCount = 0;

					while(!flag)
					{
						psi = 2 * M_PI * random_uniform(x1, y1, z1, w1);   // ��������� ������� ����������� ���� �������������
						p_max = pmax(psi, &params_opt); // ������������ �������� ������ ������������� 
																  // � ����������� ���� psi 				
						p1 = apply_Newton_psi_energy(psi, flag, energy_value, p_max, &params_opt);
						// ���� p1 ����������, �� flag = true, � �� ��������� ��������� ���� ���������, ������� ������� �� �����
						
						// ���� �� 15 ������� �� ����� �������, ������� �� �����
						if(!flag)
							iCount++;
						if(iCount > 15 && !flag)
							flag = true;
					}
					px = p1*cos(psi); // ��������� ���������� ��������
					py = p1*sin(psi); // ��������� ���������� ��������
				}				
			}
		}
		t_num++;
	}
	dev_average_value_x[idx] = value_x / t;
	dev_average_value_y[idx] = value_y / t;
	dev_average_time_array[idx] = t / n0;
	dev_nAc[idx] = nAc;
	dev_nOpt[idx] = nOpt;
	for(int i=0; i<num_logs; i++)
	{
		px_log[i+idx*num_logs] = px_log_local[i];
		py_log[i+idx*num_logs] = py_log_local[i];
	};
	delete [] px_log_local;
	delete [] py_log_local;
}

Result_one_point one_graphic_point(Params* params, double beta, double* px_mas, double* py_mas, double* WerOpt, double* WerAc, double var_value, char* filename_base)
{
	double* values_x = new double [(*params).n_part]; // ��������� ����������� ���� ����� �� (�� ���������� �� ��������)
	double* values_y = new double [(*params).n_part]; // ��������� ����������� ���� ����� �� (�� ���������� �� ��������)
	double*	av_time = new double [(*params).n_part];  // ������� ����� ���������� (�� ���������� �� ��������)
	unsigned int* nOpt = new unsigned int[(*params).n_part]; // ���������� ��������� �� ���������� �������
	unsigned int* nAc = new unsigned int[(*params).n_part];  // ���������� ��������� �� ������������ �������
		
	//�������-����
	int num_logs = (int)((*params).all_time / ((*params).dt))+1;
	double* px_log = new double[num_logs*((*params).n_part)];
	double* py_log = new double[num_logs*((*params).n_part)];
	
	// �������������� ��������� ��������� �����
	srand(time(NULL));
	unsigned int* seed = new unsigned int[(*params).n_part];
	for(int j = 0; j < (*params).n_part; j++)
	{
		seed[j] = ((unsigned int)rand()) % 100000000 + 100000000;
	};
	
	// ������ �������� ������������� �� ������������ Xeon Phi
    #pragma offload target(mic:0) inout(values_x,values_y,av_time,nAc,nOpt:length((*params).n_part)) in(seed:length((*params).n_part)) in(px_mas, py_mas, WerAc, WerOpt:length(((*params).Nx + 1) * ((*params).Ny + 1)))  in(beta) in(params) in(num_logs) inout(px_log, py_log: length(((*params).n_part) * num_logs))
	{
		omp_set_num_threads((*params).num_threads_openmp);
		#pragma omp parallel for
		for(int j = 0; j < (*params).n_part; j++)
		{
			jobKernel(  values_x,
						values_y,
						av_time,
						nAc,
						nOpt,
						params,
						beta,
						seed[j],
						j,
						px_mas,  py_mas,
						WerAc,	 WerOpt,
						px_log, py_log, 
						num_logs
					);
		};
	};
	
	// ���������� �������� ��������� ������������� � ������ ������ ������� ��� ������ �������
	cout << "Writing logs for " << var_value << " .. ";
	write_logs(filename_base, var_value, px_log, py_log, num_logs, *params);
	cout << "done " << endl;
	
	// �������� ��������� ������� ��� ����� ����� ������� � ���������
	Result_one_point result;
	result.result_value_mas_x = Mean(values_x, (*params).n_part);
	result.result_value_mas_y = Mean(values_y, (*params).n_part);
	result.std_values_mas_x = Std(values_x, (*params).n_part);
	result.std_values_mas_y = Std(values_y, (*params).n_part);
	result.result_av_time = Mean(av_time, (*params).n_part);
	result.result_nOpt = Mean(nOpt,(*params).n_part);
	result.result_nAc = Mean(nAc, (*params).n_part);
	
	// ������� ������� - ����
	delete [] px_log;
	delete [] py_log;
	
	// ������� �������, �������� �������� ��������, �������� ������� ���������� � �.�. �� ���������� �� ��������
	delete [] values_x;
	delete [] values_y;
	delete [] av_time;
	delete [] nOpt;
	delete [] nAc;

	return result;
}