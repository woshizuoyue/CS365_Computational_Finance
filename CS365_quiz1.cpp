// CS365_quiz1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

class Bond
{
public: 

	double Face;
	double issue; 
	double maturity;
	int cpnFreq, numCpnPeriods;

	vector<double> cpnAmt;
	vector<double> cpnDate;

public:

	Bond(double F, double issue_date, int num_periods, int freq, const vector<double> &c)
	{
		if (F < 0)
		{
			cout << "the Face must be positive";
			exit(1);
		}

		Face = F;

		if (freq < 1)
		{
			cout << "the Frequency at least one per year";

			exit(1);
		}

		cpnFreq = freq;

		if (num_periods < 1)
		{
			cout << "at least one coupon period";
			exit(1);
		}

		numCpnPeriods = num_periods;

		issue = issue_date;

		maturity = issue + (double)numCpnPeriods / cpnFreq;

		cpnAmt.resize(numCpnPeriods);

		cpnDate.resize(numCpnPeriods);
		
		//question!!!

		double temp = issue + 1.0 /cpnFreq;

		cpnDate[0] = temp;

		for (int i = 1; i < numCpnPeriods; i++)
		{
			cpnDate[i] = cpnDate[i - 1] + 1 / (double)cpnFreq;
		}

		setCoupons(c);
	}
	
	~Bond() {}

	void setFlatCoupons(double c)
	{
		if (c < 0.0)c = 0.0;
		fill(cpnAmt.begin(), cpnAmt.end(), c);
	}

	void setCoupons(const vector<double> &c)
	{
		// for loop to set each coupon amount with coupon date;
		
		for (int i = 0; i < cpnAmt.size(); i++)
		{
			if (i >= c.size())cpnAmt[i] = c.back();

			else {

				if (c[i] >= 0)
				{
					cpnAmt[i] = c[i];
				}
				else
				{
					cpnAmt[i] = 0;
				}
			}
		}
			
	}

	double FairValue(double t0, double y)const
	{
		double B = 0;
		double dummy1 = 0;
		double dummy2 = 0;
		FV_duration(t0, y, B, dummy1, dummy2);

		return B;
	}

	int FV_duration(double t0, double y, double &B, double &Mac_dur, double &mod_dur)const
	{
		B = 0;
		Mac_dur = 0;
		mod_dur = 0;

		y = y * 0.01;

		const double tol = 1.0e-6; 

		if (t0 < issue || t0 >= maturity)return 1; // fail;

		for (int i = 0; i < numCpnPeriods; i++)
		{
			double time = cpnDate[i] - t0;

			if (time >= tol)
			{
				if (i == numCpnPeriods - 1)
				{
					B = B + (Face + cpnAmt[i] / cpnFreq) / pow((1 + y / cpnFreq), time*cpnFreq);

					Mac_dur = Mac_dur + time * (Face + cpnAmt[i] / cpnFreq) / pow((1 + y / cpnFreq), time*cpnFreq);

				}

				else {

					B = B + (cpnAmt[i] / cpnFreq) / pow((1 + y / cpnFreq), cpnFreq*time);
					Mac_dur = Mac_dur + time * (cpnAmt[i] / cpnFreq) / pow((1 + y / cpnFreq), time*cpnFreq);
				}

			}
		}

		Mac_dur  = Mac_dur / B;
		mod_dur = Mac_dur / (1 + y / cpnFreq);

		return 0;

	}
	
};

int yield(double &y, int &num_iter, const Bond &bond, double B_target, double t0, double tol = 1.0e-4, int max_iter = 100)
{
	y = 0;
	num_iter = 0;

	if (B_target <= 0.0 || t0 < bond.issue || t0 >= bond.maturity) return 1;

	double y_low = 0.0;

	double B_y_low = bond.FairValue(t0, y_low);
	
	double diff_B_y_low = B_y_low - B_target;

	
	if (abs(diff_B_y_low) <= tol)
	{
		y = y_low;

		return 0;
	}
	

	double y_high = 100.00;

	double B_y_high = bond.FairValue(t0, y_high);

	double differ_B_y_high = B_y_high - B_target;
	
	if (abs(differ_B_y_high) <= tol)
	{
		y = y_high;

		return 0;
	}
	
	if (diff_B_y_low*differ_B_y_high > 0)
	{
		y = 0;

		return 1;
	}
	
	// main bisection iteration loop

	for (num_iter = 1; num_iter < max_iter; ++num_iter)
	{
		
		y = (y_low + y_high) / 2.0;

		double diff_B = bond.FairValue(t0, y) - B_target;

		if (abs(diff_B) <= tol)return 0;

		if (diff_B*diff_B_y_low > 0) y_low = y;
		else
			y_high = y;

		if (abs(y_high - y_low) <= tol)
		{
			return 0;
		}

	}

	y = 0;

	return 1;
}


int main()
{

	double F = 100.0;
	int f = 2;
	vector<double> c(1, 4.0);
	int num_cpn1 = 4;
	int num_cpn2 = 8;
	double issue1 = 0.0;
	double issue2 = -0.56;
	double t0 = 0.0;
	double y = 5.0;
	double FV = 0.0;
	double Mac_dur = 0.0;
	double mod_dur = 0.0;

	Bond bond1(F, issue1, num_cpn1, f, c);
	Bond bond2(F, issue2, num_cpn2, f, c);

	bond1.FV_duration(t0, y, FV, Mac_dur, mod_dur);
	cout << "#1a: " << FV << "   " << Mac_dur << "   " << mod_dur << endl;

	bond2.FV_duration(t0, y, FV, Mac_dur, mod_dur);
	cout << "#2a: " << FV << "   " << Mac_dur << "   " << mod_dur << endl;

	vector<double> cpn(100, 4.0);
	bond1.setCoupons(cpn);
	bond2.setCoupons(cpn);

	bond1.FV_duration(t0, y, FV, Mac_dur, mod_dur);
	cout << "#1b: " << FV << "   " << Mac_dur << "   " << mod_dur << endl;

	bond2.FV_duration(t0, y, FV, Mac_dur, mod_dur);
	cout << "#2b: " << FV << "   " << Mac_dur << "   " << mod_dur << endl;


	return 0;
}


