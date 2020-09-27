#pragma once
#include <iostream>
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
	vector<T> x;
	vector<T> y;
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
	T tal, tdi, tb;
	int tia;
	cout << "Open files..." << endl;
	fsize.open("fsize.txt");
	fdi.open("fdi.txt");
	fal.open("fal.txt");
	fia.open("fia.txt");
	fb.open("fb.txt");
	cout << "Done." << endl;

	fsize >> N;
	cout << N << endl;
	/*
	for (fal >> tal; !fal.eof(); fal >> tal)
		al.push_back(tal);
	cout << "al=" << al.size() << endl;
	for (fia >> tia; !fia.eof(); fia >> tia)
		ia.push_back(tia);
	cout << "ia=" << ia.size() << endl;
	for (fdi >> tdi; !fdi.eof(); fdi >> tdi)
		di.push_back(tdi);
	cout << "di=" << di.size() << endl;
	for (fb >> tb; !fb.eof(); fb >> tb)
		b.push_back(tb);
	cout << "b=" << b.size() << endl;
	*/
	for (int i = 0; i < 3; i++)
	{
		fal >> tal;
		al.push_back(tal);
		cout << "al=" << tal << endl;
	}

	for (int i = 0; i < 4; i++)
	{
		fia >> tia;
		ia.push_back(tia);
		cout << "ia=" << tia << endl;
	}

	for (int i = 0; i < 3; i++)
	{
		fdi >> tdi;
		di.push_back(tdi);
		cout << "di=" << tdi << endl;
	}

	for (int i = 0; i < 3; i++)
	{
		fb >> tb;
		b.push_back(tb);
		cout << "b=" << tb << endl;
	}

	x.resize(N);
	y.resize(N);
	f.resize(N);

	fsize.close();
	fdi.close();
	fal.close();
	fia.close();
	fb.close();
}

template <class T>
void Matrix<T>::LLt()
{
	T sumal = 0;
	di[0] = sqrt(di[0]);
	for (int i = 1; i < N; i++)
	{
		cout << "------------------" << endl;
		cout << "i = " << i << endl;
		int i_first = ia[i];
		int i_last = ia[i + 1];
		int nCol = i - (i_last - i_first);
		T sumdi = 0;
		for (int k = i_first; k < i_last; k++, nCol++)
		{
			int j_first = ia[nCol];
			int j_last = ia[nCol + 1];
			int nRow = nCol - (j_last - j_first);
			int ki = i_first;
			int kj = j_first;
			sumal = 0;

			if (nRow <= nCol)
				ki += nRow - nCol;
			else
				kj += nCol - nRow;

			for (; kj < j_last; ki++, kj++)
			{
				sumal += al[ki] * al[kj];
				cout << al[ki] << " * " << al[kj] << endl;
			}
			cout << "sumal = " << sumal << endl;
			cout << "a = " << al[k] << endl;
			al[k] = (al[k] - sumal) / di[nCol];
			sumdi += al[k] * al[k];
		}
		di[i] = sqrt(di[i] - sumdi);
	}
}

template <class T>
void Matrix<T>::calc_y()
{
	for (int i = 0; i < N; i++)
	{
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);
		T sum = 0;

		for (int k = i0; k < i1; k++, j++)
			sum += al[k] * y[j];

		y[i] = (b[i] - sum) / di[i];
	}
	cout << "y:" << endl;

	for (int i = 0; i < N; i++)
		cout << y[i] << endl;

	cout << endl;
}

template <class T>
void Matrix<T>::calc_x()
{
	for (int i = N - 1; i >= 0; i--)
	{
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);
		T xi = y[i] / di[i];
		for (int k = i0; k < i1; k++, j++)
			y[j] -= al[k] * xi;
		x[i] = xi;
	}
	fstream out_x;
	out_x.open("output_x.txt", 'w');
	cout << "Writing x..." << endl;
	for (int i = 0; i < N; i++)
	{
		out_x << x[i] << endl;
		cout << x[i] << endl;
	}
	cout << "Done." << endl;
}

template <class T>
void Matrix<T>::ChangeToTight()
{
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
			cout << tightMat[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}

template <class T>
void Matrix<T>::GilbertVec()
{
	for (int i = 0; i < N; i++)
	{
		T sum = 0;
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
			tightMat[i][j] = (T)1.0 / (i + j + 1);
	}
	GilbertVec();
	//output
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

	for (int i = 0; i < N; i++)
		cout << f[i] << endl;
}