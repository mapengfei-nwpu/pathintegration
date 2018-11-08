#include"structs.h"
#include<cmath>
double RHS(double M10, double M01, double M11, double M20, double M02, double TM, int MARK)
{
	//double ALPHA = 1.0, BETA = 0.3, YETA = 0.1, FF = 0.2, OMIGA = 1.2, SEGMA2 = 0.01;
	double GAMA = 0.2, BETA = 0.0, EPSON = 0.3, OMIGA1 = 1.0, OMIGA = 1.6, SIGMA1 = 1.0, SIGMA2 = 0.4;
	double ret;
	if (MARK == 1)
		ret = M01;
	if (MARK == 2)
		ret = -(GAMA*BETA / OMIGA1 / OMIGA1)*(3.0*M01*M02 - 2.0*M01*M01*M01) - OMIGA1*OMIGA1*EPSON*(3.0*M10*M20 - 2.0*M10*M10*M10) - GAMA*M01 - OMIGA1*OMIGA1*M10 + OMIGA1*SIGMA1*sin(OMIGA*TM);
	//ret = -YETA*M01 - ALPHA*M10 - 3.0*BETA*M10*M20 + 2.0*BETA*M10*M10*M10 + FF*cos(OMIGA*TM);

	if (MARK == 3)
		ret = M02 + (GAMA*BETA / OMIGA1 / OMIGA1)*(2.0*M01*M01*M01*M10 - 3.0*M02*M11) + OMIGA1*OMIGA1*EPSON*(2.0*M10*M10*M10*M10 - 3.0*M20*M20) - GAMA*M11 - OMIGA1*OMIGA1*M20 + OMIGA1*SIGMA1*M10*sin(OMIGA*TM);
	//ret = M02 - YETA*M11 - ALPHA*M20 - 3.0*BETA*M20*M20 + 2.0*BETA*M10*M10*M10*M10 + M10*FF*cos(OMIGA*TM);
	if (MARK == 4)
		ret = 2.0*M11;
	if (MARK == 5)
		ret = (GAMA*BETA / OMIGA1 / OMIGA1)*(4.0*M01*M01*M01*M01 - 6.0*M02*M02) + OMIGA1*OMIGA1*EPSON*(4.0*M10*M10*M10*M01 - 6.0*M20*M11) - 2.0*GAMA*M02 - 2.0*OMIGA1*OMIGA1*M11 + 2.0*OMIGA1*SIGMA1*M01*sin(OMIGA*TM) + OMIGA1*OMIGA1*SIGMA2*SIGMA2;
	//ret = -2.0*YETA*M02 - 2.0*ALPHA*M11 - 6.0*BETA*M20*M11 + 4.0*BETA*M10*M10*M10*M01 + SEGMA2 + 2.0*M01*FF*cos(OMIGA*TM);
	return ret;
}
void oderk1(Moment *moments, int n, double t, double dt) {
	for (int i = 0; i < n; i++) {
		double m10 = moments[i]._1;
		double m01 = moments[i]._2;
		double m11 = moments[i]._3;
		double m20 = moments[i]._4;
		double m02 = moments[i]._5;
		moments[i]._1 = m10 + dt*RHS(m10, m01, m11, m20, m02, t, 1);
		moments[i]._2 = m01 + dt*RHS(m10, m01, m11, m20, m02, t, 2);
		moments[i]._3 = m11 + dt*RHS(m10, m01, m11, m20, m02, t, 3);
		moments[i]._4 = m20 + dt*RHS(m10, m01, m11, m20, m02, t, 4);
		moments[i]._5 = m02 + dt*RHS(m10, m01, m11, m20, m02, t, 5);
	}
}
void oderk4(Moment *moments, int n, double t, double dt) {

	double M10_1, M01_1, M11_1, M20_1, M02_1;
	double M10_2, M01_2, M11_2, M20_2, M02_2;
	double M10_3, M01_3, M11_3, M20_3, M02_3;
	for (int i = 0; i < n; i++) {
		double M10 = moments[i]._1;
		double M01 = moments[i]._2;
		double M11 = moments[i]._3;
		double M20 = moments[i]._4;
		double M02 = moments[i]._5;
		double DT0 = dt / 100;
		double TM0 = t - dt;
		while (abs(TM0 - t) > 1.0e-3) {
			//STEP 1
			M10_1 = M10 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 1);
			M01_1 = M01 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 2);
			M11_1 = M11 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 3);
			M20_1 = M20 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 4);
			M02_1 = M02 + DT0*RHS(M10, M01, M11, M20, M02, TM0, 5);
			//STEP 2
			M10_2 = 0.75*M10 + 0.25*M10_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 1);
			M01_2 = 0.75*M01 + 0.25*M01_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 2);
			M11_2 = 0.75*M11 + 0.25*M11_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 3);
			M20_2 = 0.75*M20 + 0.25*M20_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 4);
			M02_2 = 0.75*M02 + 0.25*M02_1 + 0.25*DT0*RHS(M10_1, M01_1, M11_1, M20_1, M02_1, TM0, 5);
			//STEP 3
			M10_3 = M10 / 3.0 + 2.0*M10_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 1) / 3.0;
			M01_3 = M01 / 3.0 + 2.0*M01_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 2) / 3.0;
			M11_3 = M11 / 3.0 + 2.0*M11_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 3) / 3.0;
			M20_3 = M20 / 3.0 + 2.0*M20_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 4) / 3.0;
			M02_3 = M02 / 3.0 + 2.0*M02_2 / 3.0 + 2.0*DT0*RHS(M10_2, M01_2, M11_2, M20_2, M02_2, TM0, 5) / 3.0;
			//UPDATE
			M10 = M10_3;
			M01 = M01_3;
			M11 = M11_3;
			M20 = M20_3;
			M02 = M02_3;
			//LOOP CONDITIONS
			TM0 = TM0 + DT0;
			//cout << i << " ";
		}
		moments[i]._1 = M10;
		moments[i]._2 = M01;
		moments[i]._3 = M11;
		moments[i]._4 = M20;
		moments[i]._5 = M02;
	}
}