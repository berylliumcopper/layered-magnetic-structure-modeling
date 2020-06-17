#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <random>
#include <fstream>
#include <string>

#include <thread>
#include <mutex>

using namespace std;

constexpr double pi{ 3.14159265358979323846264338327950288 };

constexpr double J{ 5.10 };
constexpr double J_intra{ 53.0 };
constexpr double K{ 1.58 };
//constexpr double T{ 8.617333262145e-5 * 22.0 / 4.607 / 5.7883818012e-5 / 1.000 };

template <int p>
class Status
{
public:
	Status()
	{
		array<double, p>theta0;
		array<double, p>phi0;
		for (int i = 0; i < p; i++)
		{
			theta0[i] = 0;
			phi0[i] = 0;
		}
		theta = theta0;
		phi = phi0;
	}
	Status(array<double, p>& theta0, array<double, p>& phi0)
	{
		theta = theta0;
		phi = phi0;
		for (int i = 0; i < p; i++)
		{
			phi[i] = phi0[i] - round(phi0[i] / (2 * pi)) * 2 * pi;
		}
	}

	double read_theta(int i) const
	{
		return theta.at(i);
	}

	double read_phi(int i) const
	{
		return phi.at(i);
	}

	double magnetization();
	double energy(double field);
	double energy_MF(double field, array<double, p>& const mag_x_MF, array<double, p>& const mag_y_MF, array<double, p>& const mag_z_MF, double rate);
	double old_energy_MF(double field, array<double, p>& const mag_x_MF, array<double, p>& const mag_y_MF, array<double, p>& const mag_z_MF, double rate);
	Status<p> shift(array<double, p>& x, array<double, p>& y);

private:
	array<double, p> theta;
	array<double, p> phi;
};

template<int p>
double Status<p>::energy(double field)
{
	double res = 0;
	for (int i = 0; i < p; i++)
	{
		res -= field * cos(theta.at(i));
		res -= 0.5 * K * (cos(theta.at(i))) * (cos(theta.at(i)));
	}
	for (int i = 1; i < p; i++)
	{
		res += 0.5 * J * (cos(theta.at(i - 1)) * cos(theta.at(i)) + sin(theta.at(i - 1)) * sin(theta.at(i)) * cos(phi.at(i) - phi.at(i - 1)));
	}
	return res;
}

template<int p>
double Status<p>::energy_MF(double field, array<double, p>& const mag_x_MF, array<double, p>& const mag_y_MF, array<double, p>& const mag_z_MF, double rate)
//mag_x, mag_z do not contain M_s
{
	double res = 0;
	for (int i = 0; i < p; i++)
	{
		res -= field * cos(theta.at(i));
		res -= 0.5 * K * (cos(theta.at(i))) * (cos(theta.at(i)));
		res -= 0.5 * J_intra * (mag_x_MF.at(i) * sin(theta.at(i)) * cos(phi.at(i)) + mag_y_MF.at(i) * sin(theta.at(i)) * sin(phi.at(i)) + mag_z_MF.at(i) * cos(theta.at(i))) * rate;
		//res += 0.5 * J_intra * (mag_x_MF.at(i) * sin(theta.at(i)) * cos(phi.at(i)) + mag_y_MF.at(i) * sin(theta.at(i)) * sin(phi.at(i)) + mag_z_MF.at(i) * cos(theta.at(i))) * (mag_x_MF.at(i) * sin(theta.at(i)) * cos(phi.at(i)) + mag_y_MF.at(i) * sin(theta.at(i)) * sin(phi.at(i)) + mag_z_MF.at(i) * cos(theta.at(i))) * rate;//
	}
	for (int i = 1; i < p; i++)
	{
		res += 0.5 * J * (cos(theta.at(i - 1)) * mag_z_MF.at(i) + sin(theta.at(i - 1)) * mag_x_MF.at(i) * cos(phi.at(i - 1)) + sin(theta.at(i - 1)) * mag_y_MF.at(i) * sin(phi.at(i - 1))) * rate;
		res += 0.5 * J * (cos(theta.at(i)) * mag_z_MF.at(i - 1) + sin(theta.at(i)) * mag_x_MF.at(i - 1) * cos(phi.at(i)) + sin(theta.at(i)) * mag_y_MF.at(i - 1) * sin(phi.at(i))) * rate;
		//res -= 0.5 * J * (mag_z_MF.at(i - 1) * mag_z_MF.at(i) + mag_x_MF.at(i - 1) * mag_x_MF.at(i) + mag_y_MF.at(i - 1) * mag_y_MF.at(i));
		//res += 0.5 * J * (cos(theta.at(i - 1)) * cos(theta.at(i)) + sin(theta.at(i - 1)) * sin(theta.at(i)) * cos(phi.at(i) - phi.at(i - 1)));
	}
	return res;
}

template<int p>
double Status<p>::old_energy_MF(double field, array<double, p>& const mag_x_MF, array<double, p>& const mag_y_MF, array<double, p>& const mag_z_MF, double rate)
//mag_x, mag_z do not contain M_s
{
	double res = 0;
	for (int i = 0; i < p; i++)
	{
		res -= field * cos(theta.at(i));
		res -= 0.5 * K * (cos(theta.at(i))) * (cos(theta.at(i)));
		res -= 0.5 * J_intra * (mag_x_MF.at(i) * sin(theta.at(i)) * cos(phi.at(i)) + mag_y_MF.at(i) * sin(theta.at(i)) * sin(phi.at(i)) + mag_z_MF.at(i) * cos(theta.at(i))) * rate;
	}
	for (int i = 1; i < p; i++)
	{
		//res += 0.5 * J * (cos(theta.at(i - 1)) * mag_z_MF.at(i) + sin(theta.at(i - 1)) * mag_x_MF.at(i) * cos(phi.at(i - 1)) + sin(theta.at(i - 1)) * mag_y_MF.at(i) * sin(phi.at(i - 1))) * rate;
		//res += 0.5 * J * (cos(theta.at(i)) * mag_z_MF.at(i - 1) + sin(theta.at(i)) * mag_x_MF.at(i - 1) * cos(phi.at(i)) + sin(theta.at(i)) * mag_y_MF.at(i - 1) * sin(phi.at(i))) * rate;
		//res -= 0.5 * J * (mag_z_MF.at(i - 1) * mag_z_MF.at(i) + mag_x_MF.at(i - 1) * mag_x_MF.at(i) + mag_y_MF.at(i - 1) * mag_y_MF.at(i));
		res += 0.5 * J * (cos(theta.at(i - 1)) * cos(theta.at(i)) + sin(theta.at(i - 1)) * sin(theta.at(i)) * cos(phi.at(i) - phi.at(i - 1)));
	}
	return res;
}

template<int p>
double Status<p>::magnetization()
{
	double res = 0;
	for (int i = 0; i < p; i++)
		res += cos(theta.at(i));
	return res;
}

template<int p>
Status<p> Status<p>::shift(array<double, p>& x, array<double, p>& y)
{
	array<double, p> n_theta;
	array<double, p> n_phi;
	for (int i = 0; i < p; i++)
	{
		n_theta.at(i) = acos(-sin(x.at(i)) * cos(y.at(i)) * sin(theta.at(i)) + cos(x.at(i)) * cos(theta.at(i)));
		n_phi.at(i) = atan2(sin(y.at(i)) * sin(x.at(i)), cos(theta.at(i)) * cos(y.at(i)) * sin(x.at(i)) + sin(theta.at(i)) * cos(x.at(i))) + phi.at(i);
	}
	return { n_theta,n_phi };
}

template<int p>
ostream& operator<<(ostream& os, const Status<p>& arg)
{
	for (int i = 0; i < p; i++)
	{
		os << arg.read_theta(i) << "\t" << arg.read_phi(i) << "\t";
	}
	os << "\n";
	return os;
}

template<int p>
class Result_M_C
{
public:
	array<double, p> mag_x;
	array<double, p> mag_y;
	array<double, p> mag_z;
	double ave_energy;
	Result_M_C()
	{
		array<double, p> mag_x0;
		array<double, p> mag_y0;
		array<double, p> mag_z0;
		for (int i = 0; i < p; i++)
		{
			mag_x0.at(i) = 0;
			mag_y0.at(i) = 0;
			mag_z0.at(i) = 0;
		}
		mag_x = mag_x0;
		mag_y = mag_x0;
		mag_z = mag_x0;
		ave_energy = 0.0;
	}
	Result_M_C(array<double, p>& const mag_x0, array<double, p>& const mag_y0, array<double, p>& const mag_z0, double ave_energy0)
	{
		mag_x = mag_x0;
		mag_y = mag_y0;
		mag_z = mag_z0;
		ave_energy = ave_energy0;
	}
};

template<int p>
class M_C
{
public:
	M_C(int num0, int heat0, double T0, double delta0, double field0, Status<p> temp0, Result_M_C<p> m_MF0, ofstream& os, double rate0)
	{
		if (num0 <= 0 || delta0 <= 0 || heat0 <= 0 || rate0 > 1.00000001 || rate0 < 0)
			throw runtime_error("M_C: invalid input");
		num = num0;
		heat = heat0;
		T = T0;
		delta = delta0;
		field = field0;
		buff = temp0;
		m_MF = m_MF0;
		rate = rate0;
		/*os << setprecision(15) << p << "\t" << num << "\t" << heat << "\t" << delta << "\t" << field << "\t" << rate << "\n" << buff << "---\n";
		for (int i = 0; i < p; i++)
		{
			os << m_MF.mag_x.at(i) << "\t" << m_MF.mag_y.at(i) << "\t" << m_MF.mag_z.at(i) << "\n";
		}
		os << "---\n";*/
		os.write((char*) & (*this), sizeof(*this));
	}

	Result_M_C<p> run(ofstream& os, std::mt19937_64& gen);
	Result_M_C<p> old_run(ofstream& os, std::mt19937_64& gen);

private:
	int num;//number of points
	int heat;
	double T;
	double delta;//the maxinum length of the status change
	double field;//magnetic field
	Result_M_C<p> m_MF;
	Status<p> curr;
	Status<p> buff;
	double rate;
};

template<int p>
Result_M_C<p> M_C<p>::run(ofstream& os, std::mt19937_64& gen)
{
	/*std::random_device rd;
	std::mt19937_64 gen(rd());*/
	std::uniform_real_distribution<> dis_x(0.0, delta);
	std::uniform_real_distribution<> dis_y(0.0, 2.0 * pi);
	std::uniform_real_distribution<> dis_a(0.0, 1.0);
	array<double, p>sum_x, sum_y, sum_z;
	for (int i = 0; i < p; i++)
	{
		sum_x.at(i) = 0.0;
		sum_y.at(i) = 0.0;
		sum_z.at(i) = 0.0;
	}
	double sum_energy = 0.0;
	double buff_energy = buff.energy_MF(field, m_MF.mag_x, m_MF.mag_y, m_MF.mag_z, rate);
	for (int i = 0; i < heat + num; i++)
	{
		//step 1: generate the random vectors
		array<double, p>dx;
		array<double, p>dy;
		for (int j = 0; j < p; j++)
		{
			dx.at(j) = dis_x(gen);
			dy.at(p - 1 - j) = dis_y(gen);
		}
		//step 2: calculate the energy of the new status
		curr = buff.shift(dx, dy);
		double curr_energy = curr.energy_MF(field, m_MF.mag_x, m_MF.mag_y, m_MF.mag_z, rate);
		//step 3: judge the transition
		if (curr_energy <= buff_energy)//probability of transition P=1
		{
			buff = curr;
			buff_energy = curr_energy;
		}
		else
		{
			double prob = exp((buff_energy - curr_energy) / T);
			double judge = dis_a(gen);
			if (judge < prob)
			{
				buff = curr;
				buff_energy = curr_energy;
			}
		}
		//step 4: calculate the magnetization
		if (i >= heat)
		{
			for (int j = 0; j < p; j++)
			{
				sum_x.at(j) += sin(buff.read_theta(j)) * cos(buff.read_phi(j));
				sum_y.at(j) += sin(buff.read_theta(j)) * sin(buff.read_phi(j));
				sum_z.at(j) += cos(buff.read_theta(j));
			}
		}
		sum_energy += buff_energy;
		//step 5: write the original data
		//os << buff;
	}
	for (int i = 0; i < p; i++)
	{
		sum_x.at(i) = sum_x.at(i) / num;//
		sum_y.at(i) = sum_y.at(i) / num;//
		sum_z.at(i) = sum_z.at(i) / num;//
	}
	sum_energy = sum_energy / num;
	for (int i = 0; i < p; i++)
	{
		sum_energy += 0.25 * J_intra * (sum_x.at(i) * sum_x.at(i) + sum_y.at(i) * sum_y.at(i) + sum_z.at(i) * sum_z.at(i));
	}
	for (int i = 1; i < p; i++)
	{
		sum_energy -= 0.5 * J * (sum_x.at(i - 1) * sum_x.at(i) + sum_y.at(i - 1) * sum_y.at(i) + sum_z.at(i - 1) * sum_z.at(i));
	}
	return { sum_x,sum_y,sum_z,sum_energy };
}

template<int p>
Result_M_C<p> M_C<p>::old_run(ofstream& os, std::mt19937_64& gen)
{
	/*std::random_device rd;
	std::mt19937_64 gen(rd());*/
	std::uniform_real_distribution<> dis_x(0.0, delta);
	std::uniform_real_distribution<> dis_y(0.0, 2.0 * pi);
	std::uniform_real_distribution<> dis_a(0.0, 1.0);
	array<double, p>sum_x, sum_y, sum_z;
	for (int i = 0; i < p; i++)
	{
		sum_x.at(i) = 0.0;
		sum_y.at(i) = 0.0;
		sum_z.at(i) = 0.0;
	}
	double buff_energy = buff.old_energy_MF(field, m_MF.mag_x, m_MF.mag_y, m_MF.mag_z, rate);
	for (int i = 0; i < heat + num; i++)
	{
		//step 1: generate the random vectors
		array<double, p>dx;
		array<double, p>dy;
		for (int j = 0; j < p; j++)
		{
			dx.at(j) = dis_x(gen);
			dy.at(p - 1 - j) = dis_y(gen);
		}
		//step 2: calculate the energy of the new status
		curr = buff.shift(dx, dy);
		double curr_energy = curr.old_energy_MF(field, m_MF.mag_x, m_MF.mag_y, m_MF.mag_z, rate);
		//step 3: judge the transition
		if (curr_energy <= buff_energy)//probability of transition P=1
		{
			buff = curr;
			buff_energy = curr_energy;
		}
		else
		{
			double prob = exp((buff_energy - curr_energy) / T);
			double judge = dis_a(gen);
			if (judge < prob)
			{
				buff = curr;
				buff_energy = curr_energy;
			}
		}
		//step 4: calculate the magnetization
		if (i >= heat)
		{
			for (int j = 0; j < p; j++)
			{
				sum_x.at(j) += sin(buff.read_theta(j)) * cos(buff.read_phi(j));
				sum_y.at(j) += sin(buff.read_theta(j)) * sin(buff.read_phi(j));
				sum_z.at(j) += cos(buff.read_theta(j));
			}
		}
		//step 5: write the original data
		//os << buff;
	}
	for (int i = 0; i < p; i++)
	{
		sum_x.at(i) = sum_x.at(i) / num;
		sum_y.at(i) = sum_y.at(i) / num;
		sum_z.at(i) = sum_z.at(i) / num;
	}
	return { sum_x,sum_y,sum_z,0.0 };
}

constexpr int N{ 2 };

int counte{ 0 };

void sys_run(const char* ofile, const char* ofile_M_H, double temperature0)
{
	std::random_device rd;
	std::mt19937_64 gen(rd());
	double temperature = temperature0 * 8.617333262145e-5 / 4.607 / 5.7883818012e-5 / 1.000;
	array<double, N>mx = { 0.0, 0.0 };
	array<double, N>my = { 0.0, 0.0 };
	array<double, N>mz = { 0.0, 0.0 };
	Result_M_C<N> mm{ mx,my,mz,0.0 };
	array<double, N>theta_s{ 0.0, 0.0 };
	array<double, N>phi_s{ 0.0, 0.0 };
	Status<N> ss{ theta_s,phi_s };
	ofstream os;
	ofstream os_M_H;
	//const char* ofile{ "E:\\PKU\\YuGroup\\MBT even-odd effect\\T\\test\\test31.dat" };
	os.open(ofile, ios::out | ios::binary);
	os_M_H.open(ofile_M_H, ios::out);
	if (!os)
	{
		throw runtime_error("file wrong");
	}
	if(!os_M_H)
	{
		os.close();
		throw runtime_error("file wrong");
	}
	os_M_H << setprecision(15);
	for (int mag_i = 239; mag_i >= 0; mag_i--)
	{//
		double sm = 0;
		double se = 0;
		for (int iter = 0; iter < 64 * 16; iter++)
		{
			double rate = 1.0;
			//if (iter < 20)
				//rate = (iter) * 0.05;
			M_C<N> au{ 4096*2,2048,temperature,0.8,0.05 * mag_i,ss,mm,os,rate };
			if (iter < 64)
				mm = au.old_run(os,gen);
			else
				mm = au.run(os,gen);
			for (int i = 0; i < N; i++)
			{
				//cout << mm.mag_x.at(i) << "\t" << mm.mag_y.at(i) << "\t" << mm.mag_z.at(i) << "\n";
			}
			//cout << "\n";
			for (int i = 0; i < N; i++)
			{
				theta_s.at(i) = acos(mm.mag_z.at(i) / sqrt(mm.mag_x.at(i) * mm.mag_x.at(i) + mm.mag_y.at(i) * mm.mag_y.at(i) + mm.mag_z.at(i) * mm.mag_z.at(i)));
				phi_s.at(i) = atan2(mm.mag_y.at(i), mm.mag_z.at(i));
			}
		}
		bool judge = false;
		array<double, N>s_each;
		do {
			sm = 0.0;
			se = 0.0;
			for (int i = 0; i < N; i++)
				s_each.at(i) = 0.0;
			judge = false;
			for (int iter = 0; iter < 32; iter++)
			{
				double rate = 1.0;
				//if (iter < 20)
					//rate = (iter) * 0.05;
				M_C<N> au{ 8192 * 16,2048,temperature,0.5,0.05 * mag_i,ss,mm,os,rate };
				mm = au.run(os,gen);
				for (int i = 0; i < N; i++)
				{
					//cout << mm.mag_x.at(i) << "\t" << mm.mag_y.at(i) << "\t" << mm.mag_z.at(i) << "\n";
				}
				//cout << "\n";
				for (int i = 0; i < N; i++)
				{
					theta_s.at(i) = acos(mm.mag_z.at(i) / sqrt(mm.mag_x.at(i) * mm.mag_x.at(i) + mm.mag_y.at(i) * mm.mag_y.at(i) + mm.mag_z.at(i) * mm.mag_z.at(i)));
					phi_s.at(i) = atan2(mm.mag_y.at(i), mm.mag_z.at(i));
				}
				if (iter >= 16)
				{
					double gr = 0.0;
					for (int i = 0; i < N; i++) 
					{
						gr += mm.mag_z.at(i);
						s_each.at(i) += mm.mag_z.at(i) / 16.0;
					}
					sm += gr / 16.0;
					if (iter == 31 && abs(gr - sm) >= 0.01 * N)
					{
						cout << "warning!\n";
						judge = true;
					}
					se += mm.ave_energy / 16.0;
				}
			}
		} while (judge);
		/*double sm = 0;
		for (int i = 0; i < N; i++)
			sm += mm.mag_z.at(i);*/
		os_M_H << 0.05 * mag_i << "\t" << sm << "\t" << se;
		for (int i = 0; i < N; i++)
			os_M_H << "\t" << s_each.at(i);
		os_M_H << "\n";
		cout << temperature0 << "\t" << 0.05 * mag_i << "\n";
	}//
	for (int mag_i = 0; mag_i < 240; mag_i++)
	{//
		double sm = 0;
		double se = 0;
		for (int iter = 0; iter < 64 * 16; iter++)
		{
			double rate = 1.0;
			//if (iter < 20)
				//rate = (iter) * 0.05;
			M_C<N> au{ 4096*2,2048,temperature,0.8,0.05 * mag_i,ss,mm,os,rate };
			if (iter < 64)
				mm = au.old_run(os,gen);
			else
				mm = au.run(os,gen);
			for (int i = 0; i < N; i++)
			{
				//cout << mm.mag_x.at(i) << "\t" << mm.mag_y.at(i) << "\t" << mm.mag_z.at(i) << "\n";
			}
			//cout << "\n";
			for (int i = 0; i < N; i++)
			{
				theta_s.at(i) = acos(mm.mag_z.at(i) / sqrt(mm.mag_x.at(i) * mm.mag_x.at(i) + mm.mag_y.at(i) * mm.mag_y.at(i) + mm.mag_z.at(i) * mm.mag_z.at(i)));
				phi_s.at(i) = atan2(mm.mag_y.at(i), mm.mag_z.at(i));
			}
		}
		bool judge = false;
		array<double, N>s_each;
		do {
			sm = 0.0;
			se = 0.0;
			for (int i = 0; i < N; i++)
				s_each.at(i) = 0.0;
			judge = false;
			for (int iter = 0; iter < 32; iter++)
			{
				double rate = 1.0;
				//if (iter < 20)
					//rate = (iter) * 0.05;
				M_C<N> au{ 8192 * 16,2048,temperature,0.5,0.05 * mag_i,ss,mm,os,rate };
				mm = au.run(os,gen);
				for (int i = 0; i < N; i++)
				{
					//cout << mm.mag_x.at(i) << "\t" << mm.mag_y.at(i) << "\t" << mm.mag_z.at(i) << "\n";
				}
				//cout << "\n";
				for (int i = 0; i < N; i++)
				{
					theta_s.at(i) = acos(mm.mag_z.at(i) / sqrt(mm.mag_x.at(i) * mm.mag_x.at(i) + mm.mag_y.at(i) * mm.mag_y.at(i) + mm.mag_z.at(i) * mm.mag_z.at(i)));
					phi_s.at(i) = atan2(mm.mag_y.at(i), mm.mag_z.at(i));
				}
				if (iter >= 16)
				{
					double gr = 0.0;
					for (int i = 0; i < N; i++)
					{
						gr += mm.mag_z.at(i);
						s_each.at(i) += mm.mag_z.at(i) / 16.0;
					}
					sm += gr / 16.0;
					if (iter == 31 && abs(gr - sm) >= 0.01 * N)
					{
						cout << "warning!\n";
						judge = true;
					}
					se += mm.ave_energy / 16.0;
				}
			}
		} while (judge);
		/*double sm = 0;
		for (int i = 0; i < N; i++)
			sm += mm.mag_z.at(i);*/
		os_M_H << 0.05 * mag_i << "\t" << sm << "\t" << se;
		for (int i = 0; i < N; i++)
			os_M_H << "\t" << s_each.at(i);
		os_M_H << "\n";
		cout << temperature0 << "\t" << 0.05 * mag_i << "\n";
	}//
	os.close();
	os_M_H.close();
	cout << temperature0 << "\tfinished." << counte << "\n";
	counte++;
}

int main()
{
	try 
	{
		string str;
		int x1;
		cout << "Enter a path:\n";
		cin >> str;
		cout << "Temperature 0:\n";
		cin >> x1;
		string s1 = str + to_string(N) + "-" + to_string(x1) + "_0.dat";
		string s1_M_H = str + to_string(N) + "-" + to_string(x1) + "_0-MH.txt";
		string s2 = str + to_string(N) + "-" + to_string(x1) + "_5.dat";
		string s2_M_H = str + to_string(N) + "-" + to_string(x1) + "_5-MH.txt";
		string s3 = str + to_string(N) + "-" + to_string(x1 + 1) + "_0.dat";
		string s3_M_H = str + to_string(N) + "-" + to_string(x1 + 1) + "_0-MH.txt";
		string s4 = str + to_string(N) + "-" + to_string(x1 + 1) + "_5.dat";
		string s4_M_H = str + to_string(N) + "-" + to_string(x1 + 1) + "_5-MH.txt";
		thread t1{ sys_run, s1.c_str(), s1_M_H.c_str(), x1 * 1.0 };
		thread t2{ sys_run, s2.c_str(), s2_M_H.c_str(), x1 * 1.0 + 0.5 };
		thread t3{ sys_run, s3.c_str(), s3_M_H.c_str(), x1 * 1.0 + 1.0 };
		thread t4{ sys_run, s4.c_str(), s4_M_H.c_str(), x1 * 1.0 + 1.5 };
		t1.detach();
		t2.detach();
		t3.detach();
		t4.detach();
		system("pause");
	}
	/*catch (runtime_error& e)
	{
		cerr << e.what();
	}*/
	catch (...)
	{
		cerr << "unknown error";
	}
	return 0;
}
