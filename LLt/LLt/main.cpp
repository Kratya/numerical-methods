#include "Matrix.h"

void main()
{
	setlocale(LC_CTYPE, "Russian");
	Matrix<float> A;
	//Matrix<double> A;
	A.MatInput();
	A.ChangeToTight();
	A.GilbertMat();
	A.Gauss();
	A.LLt();
	A.calc_y();
	A.calc_x();
}