#include "Matrix.h"

void main()
{
	setlocale(LC_CTYPE, "Russian");
	//Matrix<float> A;
	Matrix<double> A;
	A.MatInput();
	A.ChangeToTight();
	A.GilbertMat();
	A.ChangeToProf();
	cout << "Gauss" << endl;
	A.Gauss();
	cout << "LLt" << endl;
	A.LLt();
	A.calc_y();
	A.calc_x();
}
