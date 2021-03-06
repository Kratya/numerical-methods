#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

template <class T>
class Matrix
{
private:
	int N, N_al;
	vector<T> di;
	vector<T> al;
	vector<int> ia;
	vector<T> b;
	//vector<T> x;
	//vector<T> y;
	vector<T> f;
	vector<vector<T>> tightMat;
public:
	void MatInput();
	void calc_x();
	void calc_y();
	void LLt();
	void ChangeToTight();
	void GilbertVec();
	void GilbertMat();
	void getCol();
	void ChangeToProf();
	void Gauss();
};

template <class T>
void Matrix<T>::MatInput()
{
	ifstream fsize, fdi, fal, fia, fb;
	ofstream fout;
	T tal, tdi, tb;
	int tia;
	cout << "Open files..." << endl;
	fsize.open("fsize.txt");
	fdi.open("fdi.txt");
	fal.open("fal.txt");
	fia.open("fia.txt");
	fb.open("fb.txt");
	fout.open("out.txt");
	fout.precision(5);
	cout << "Done." << endl;

	fsize >> N;
	cout << "N =  " << N << endl;
	ia.resize(N + 1);
	b.resize(N);
	di.resize(N);
	int j = 0;

	j = 0;
	while (fia >> tia)
	{
		ia[j] = tia;
		j++;
	}
	
	j = 0;
	while (fdi >> tdi)
	{
		di[j] = tdi;
		j++;
	}
	
	j = 0;
	while (fb >> tb)
	{
		b[j] = tb;
		j++;
	}

	al.resize(ia[N]);
	j = 0;
	while (fal >> tal)
	{
		al[j] = tal;
		j++;
	}
	
	//x.resize(N);
	//y.resize(N);
	f.resize(N);

	fout << "The size of matrix is: " << N << "x" << N << endl;
	fout << "Elements of DI: [";
	for (int i = 0; i < di.size(); i++)
		fout << di[i] << " ";
	fout << "]" << endl;
	fout << "Elements of AL: [";
	for (int i = 0; i < al.size(); i++)
		fout << al[i] << " ";
	fout << "]" << endl;
	fout << "Elements of IA: [";
	for (int i = 0; i < ia.size(); i++)
		fout << ia[i] << " ";
	fout << "]" << endl;
	fout << "Elements of B: [";
	for (int i = 0; i < b.size(); i++)
		fout << b[i] << " ";
	fout << "]" << endl;

	fsize.close();
	fdi.close();
	fal.close();
	fia.close();
	fb.close();
	fout.close();
}

template <class T>
void Matrix<T>::LLt()
{
	T sumal = 0;
	di[0] = sqrt(di[0]);
	for (int i = 1; i < N; i++)
	{
		int i_first = ia[i];
		int i_last = ia[i + 1];
		int nCol = i - (i_last - i_first);
		T sumdi = 0;
		//int nColbuf = nCol;
		for (int k = i_first; k < i_last; k++, nCol++)
		{
			int j_first = ia[nCol];
			int j_last = ia[nCol + 1];
			int nRow = nCol - (j_last - j_first);
			int ti = i_first;
			int tj = j_first;
			/*
			int ki = i_first;
			int kj = j_first;
			sumal = 0;

			if (nRow <= nCol) //<=
				ki += nRow - nCol;
			else
				kj += nCol - nRow;

			*/
			int ki = k - i_first;
			int kj = j_last - j_first;
			sumal = 0;

			if (kj > ki)
				tj -= ki - kj;
			else
			{
				nRow += ki - kj;
				ti += ki - kj;
			}

			for (; ti < k; ti++, tj++)
			{
				sumal += al[ti] * al[tj];
				//cout << al[ki] << " * " << al[kj] << endl;
			}
			al[k] = (al[k] - sumal) / di[nCol];
			sumdi += al[k] * al[k];
		}
		di[i] = sqrt(di[i] - sumdi);
	}
}

template <class T>
void Matrix<T>::calc_y()
{
	ofstream fout;
	//fout.precision(15);
	fout.open("out.txt", ios::app);
	fout << fixed;
	fout << setprecision(15);
	fout << "Calc Y by LLT method:" << endl;
	for (int i = 0; i < N; i++)
	{
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);
		T sum = 0;

		for (int k = i0; k < i1; k++, j++)
			sum += al[k] * b[j];

		b[i] = (b[i] - sum) / di[i];
	}
	cout << "y:" << endl;
	fout << "y = [";

	for (int i = 0; i < N; i++)
	{
		cout << b[i] << endl;
		fout << b[i] << " ";
	}

	fout << "]" << endl;
	fout.close();
	cout << endl;
}

template <class T>
void Matrix<T>::calc_x()
{
	ofstream fout;
	//fout.precision(15);
	fout.open("out.txt", ios::app);
	fout << fixed;
	fout << setprecision(15);
	fout << "Calc X by LLT method:" << endl;
	for (int i = N - 1; i >= 0; i--)
	{
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);
		T xi = b[i] / di[i];
		for (int k = i0; k < i1; k++, j++)
			b[j] -= al[k] * xi;
		b[i] = xi;
	}
	cout << "Writing x..." << endl;
	fout << "x = [";
	for (int i = 0; i < N; i++)
	{
		fout << b[i] << " ";
		cout << b[i] << endl;
	}
	fout << "]" << endl;
	fout.close();
	cout << "Done." << endl;
}

template <class T>
void Matrix<T>::ChangeToTight()
{
	ofstream fout;
	//fout.precision(15);
	fout.open("out.txt", ios::app);
	fout << fixed;
	fout << setprecision(15);
	fout << "Tight matrix:" << endl;
	tightMat.resize(N, vector<T>(N));
	for (auto& i : tightMat)
		fill(i.begin(), i.end(), 0);

	for (int i = 0; i < N; i++)
	{
		int j = i - (ia[i + 1] - ia[i]);
		for (int k = ia[i]; k < ia[i + 1]; k++, j++)
			tightMat[i][j] = tightMat[j][i] = al[k];
		tightMat[i][i] = di[i];
	}
	cout << "Tight matrix:" << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			fout << tightMat[i][j] << " ";
			cout << tightMat[i][j] << " ";
		}
		fout << endl;
		cout << endl;
	}
	fout << endl;
	cout << endl;
	fout.close();
}

template <class T>
void Matrix<T>::GilbertVec()
{
	for (int i = 0; i < N; i++)
	{
		double sum = 0;
		for (int k = 0; k < N; k++)
			sum += tightMat[i][k] * (k + 1);
		b[i] = sum;
	}
}

template <class T>
void Matrix<T>::GilbertMat() // ���������� ������� ���������
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			tightMat[i][j] = (double)1.0 / (i + j + 1);
	}
	GilbertVec();
}

template <class T>
void Matrix<T>::getCol()
{
	int col = 0;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < i; j++)
		{
			if (tightMat[i][j] != 0 || tightMat[j][i] != 0)
				col++;
		}
	N_al = col;
}

template <class T>
void Matrix<T>::ChangeToProf() // ������� ������� � ���������� ������
{
	getCol();
	al.resize(N_al);
	int s = 0;
	for (int i = 0; i < N; i++)
	{
		di[i] = tightMat[i][i];
		ia[i] = s;
		for (int j = 0; j < i; j++)
		{
			if (tightMat[i][j])
			{
				al[s] = tightMat[i][j];
				s++;
			}
		}
	}
	ia[N] = s;

	cout << "al: ";
	for (int i = 0; i < N_al; i++)
		cout << al[i] << " ";
	cout << endl;

	cout << "di: ";
	for (int i = 0; i < N; i++)
		cout << di[i] << " ";
	cout << endl;

	cout << "ia: ";
	for (int i = 0; i < N + 1; i++)
		cout << ia[i] << " ";
	cout << endl;
}

template <class T>
void Matrix<T>::Gauss() // ����� ���c�� �������
{
	ofstream fout;
	//fout.precision(15);
	fout << fixed;
	fout << setprecision(15);
	fout.open("out.txt", ios::app);
	fout << "Calc by Gauss method:" << endl;
	for (int i = 0; i < N; i++)
		f[i] = b[i];

	for (int i = 1; i < N; i++)
		for (int j = i; j < N; j++)
		{
			T m = tightMat[j][i - 1] / tightMat[i - 1][i - 1];
			for (int k = 0; k < N; k++)
				tightMat[j][k] = tightMat[j][k] - m * tightMat[i - 1][k];
			f[j] = f[j] - m * f[i - 1];
		}

	for (int k = N - 1; k >= 0; k--)
	{
		T buf = 0;
		for (int j = k + 1; j < N; j++)
			buf += tightMat[k][j] * f[j];

		f[k] = f[k] - buf;
		f[k] = f[k] / tightMat[k][k];

	}

	cout << endl;
	fout << "x = [";
	for (int i = 0; i < N; i++)
	{
		fout << f[i] << " ";
		cout << f[i] << endl;
	}
	fout << "]" << endl;
	fout.close();
}