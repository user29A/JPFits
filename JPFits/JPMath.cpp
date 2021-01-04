/*Copyright 2017 Joseph Edwin Postma

Contact email: joepostma@live.ca

This file is part of JPMath.

JPMath is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

JPMath is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with JPMath. If not, see http://www.gnu.org/licenses/. */


#include "stdafx.h"
#include "JPFITS.h"

JPFITS::JPMath::PointD::PointD(double x, double y, double value)
{
	POINTX = x;
	POINTY = y;
	POINTVAL = value;
}

inline double JPFITS::JPMath::PointD::DistanceTo(JPMath::PointD^ other_point)
{
	double dx = POINTX - other_point->X, dy = POINTY - other_point->Y;
	return Math::Sqrt(dx * dx + dy * dy);
}

inline double JPFITS::JPMath::PointD::ISLEFT(JPMath::PointD^ P0, JPMath::PointD^ P1, JPMath::PointD^ P2)
{
	return ((P1->X - P0->X) * (P2->Y - P0->Y) - (P2->X - P0->X) * (P1->Y - P0->Y));
}

bool JPFITS::JPMath::PointD::PointInPoly(JPMath::PointD^ P, array<JPMath::PointD^>^ V, int n)
{
	int    wn = 0;    // the  winding number counter

   // loop through all edges of the polygon
	for (int i = 0; i < n; i++) {   // edge from V[i] to  V[i+1]
		if (V[i]->Y <= P->Y) {          // start y <= P.y
			if (V[i + 1]->Y > P->Y)      // an upward crossing
				if (ISLEFT(V[i], V[i + 1], P) > 0)  // P left of  edge
					++wn;            // have  a valid up intersect
		}
		else {                        // start y > P.y (no test needed)
			if (V[i + 1]->Y <= P->Y)     // a downward crossing
				if (ISLEFT(V[i], V[i + 1], P) < 0)  // P right of  edge
					--wn;            // have  a valid down intersect
		}
	}
	
	if (wn == 0)
		return false;
	return true;
}

void JPFITS::JPMath::PointD::PolygonInteriorPointsRegion(array<bool, 2>^ region, array<JPMath::PointD^>^ polygon, int Xmin, int Ymin, int Xmax, int Ymax)
{
	JPMath::PointD^ P;
	#pragma omp parallel for private(P)
	for (int x = Xmin; x <= Xmax; x++)
		for (int y = Ymin; y <= Ymax; y++)
		{
			P = gcnew JPMath::PointD(double(x), double(y), 0);

			int wn = 0;
			for (int i = 0; i < polygon->Length - 1; i++)
				if (polygon[i]->Y <= P->Y)
				{
					if (polygon[i + 1]->Y > P->Y)
						if (ISLEFT(polygon[i], polygon[i + 1], P) > 0)
							wn++;
				}
				else
					if (polygon[i + 1]->Y <= P->Y)
						if (ISLEFT(polygon[i], polygon[i + 1], P) < 0)
							wn--;

			if (wn != 0)
				region[x, y] = true;
		}
}

JPFITS::JPMath::Triangle::Triangle(JPMath::PointD^ point0, JPMath::PointD^ point1, JPMath::PointD^ point2)
{
	POINTS = gcnew array<JPMath::PointD^>(3) { point0, point1, point2 };
	SORTTRIANGLE();
	MAKEVERTEXANGLES();
	MAKEFIELDVECTORS();
	TOTALVERTEXPOINTVALUESUM = POINTS[0]->Value + POINTS[1]->Value + POINTS[2]->Value;
}

JPFITS::JPMath::Triangle::Triangle(array<JPMath::PointD^>^ points)
{
	POINTS = gcnew array<JPMath::PointD^>(3) { points[0], points[1], points[2] };
	SORTTRIANGLE();
	MAKEVERTEXANGLES();
	MAKEFIELDVECTORS();
	TOTALVERTEXPOINTVALUESUM = POINTS[0]->Value + POINTS[1]->Value + POINTS[2]->Value;
}

inline void JPFITS::JPMath::Triangle::SORTTRIANGLE()
{
	array<JPMath::PointD^>^ points = gcnew array<JPMath::PointD^>(3) { POINTS[0], POINTS[1], POINTS[2] };
	double D01 = points[0]->DistanceTo(points[1]);
	double D02 = points[0]->DistanceTo(points[2]);
	double D12 = points[1]->DistanceTo(points[2]);
	SIDELENGTHS = gcnew array<double>(3) { D01, D02, D12 };
	array<int>^ distseq = gcnew array<int>(3) { 0, 1, 2 };
	Array::Sort(SIDELENGTHS, distseq);

	int common0;

	switch (distseq[0])
	{
	case 0:
		if (distseq[1] == 1)
			common0 = 0;
		else
			common0 = 1;
		break;
	case 1:
		if (distseq[1] == 0)
			common0 = 0;
		else
			common0 = 2;
		break;
	case 2:
		if (distseq[1] == 1)
			common0 = 2;
		else
			common0 = 1;
		break;
	}

	int common1, common2;

	switch (common0)
	{
	case 0:
		if (distseq[0] == 0)
		{
			common1 = 1;
			common2 = 2;
		}
		else
		{
			common1 = 2;
			common2 = 1;
		}
		break;
	case 1:
		if (distseq[0] == 0)
		{
			common1 = 0;
			common2 = 2;
		}
		else
		{
			common1 = 2;
			common2 = 0;
		}
		break;
	case 2:
		if (distseq[0] == 1)
		{
			common1 = 0;
			common2 = 1;
		}
		else
		{
			common1 = 1;
			common2 = 0;
		}
		break;
	}

	POINTS = gcnew array<JPMath::PointD^>(3) { points[common0], points[common1], points[common2] };
}

inline void JPFITS::JPMath::Triangle::MAKEVERTEXANGLES()
{
	VERTEXANGLES = gcnew array<double>(3);
	double D01 = SIDELENGTHS[0];// POINTS[0]->DistanceTo(POINTS[1]);
	double D02 = SIDELENGTHS[1];// POINTS[0]->DistanceTo(POINTS[2]);
	double D12 = SIDELENGTHS[2];// POINTS[1]->DistanceTo(POINTS[2]);
	//c2 = a2 + b2 − 2ab cos(C)
	VERTEXANGLES[0] = Math::Acos(-(D12 * D12 - D01 * D01 - D02 * D02) / (2 * D01 * D02));
	VERTEXANGLES[1] = Math::Acos(-(D02 * D02 - D01 * D01 - D12 * D12) / (2 * D01 * D12));
	VERTEXANGLES[2] = Math::Acos(-(D01 * D01 - D02 * D02 - D12 * D12) / (2 * D02 * D12));
}

inline void JPFITS::JPMath::Triangle::MAKEFIELDVECTORS()
{
	FIELDVECTOR = gcnew JPMath::PointD(POINTS[2]->X - POINTS[1]->X, POINTS[2]->Y - POINTS[1]->Y, SIDELENGTHS[2]);
	FIELDVECTORRADANGLE = Math::Atan2(FIELDVECTOR->Y, FIELDVECTOR->X);
}

unsigned __int64 JPFITS::JPMath::Factorial(unsigned __int64 N)
{
	if (N == 0)
		return 1;
	return N * Factorial(N - 1);
}

double JPFITS::JPMath::Area_Triangle(array<double>^ x, array<double>^ y)
{
	double area = 0.5;
	area *= y[0] * (x[1] - x[2]) + y[1] * (x[2] - x[0]) + y[2] * (x[0] - x[1]);
	return Math::Abs(area);
}

void JPFITS::JPMath::RADecSexToDegree(String^ ra, String^ dec, String^ separator, double& RA, double& DEC)
{
	if (separator == nullptr || separator == "")//determine the separator
		for (int i = 0; i < ra->Length; i++)
		{
			if (IsNumeric(ra[i].ToString()))
				continue;
			separator = ra[i].ToString();
			break;
		}

	double h, m, s, d, am, as;
	int lastind = 0;
	int nextind = ra->IndexOf(separator, lastind);
	h = Convert::ToDouble(ra->Substring(lastind, nextind - lastind));
	lastind = nextind + 1;
	nextind = ra->IndexOf(separator, lastind);
	m = Convert::ToDouble(ra->Substring(lastind, nextind - lastind));
	lastind = nextind + 1;
	s = Convert::ToDouble(ra->Substring(lastind));

	lastind = 0;
	nextind = dec->IndexOf(separator, lastind);
	d = Convert::ToDouble(dec->Substring(lastind, nextind - lastind));
	lastind = nextind + 1;
	nextind = dec->IndexOf(separator, lastind);
	am = Convert::ToDouble(dec->Substring(lastind, nextind - lastind));
	lastind = nextind + 1;
	as = Convert::ToDouble(dec->Substring(lastind));

	RA = h / 24 * 360 + m / 60 / 24 * 360 + s / 3600 / 24 * 360;
	DEC = Math::Abs(d) + am / 60 + as / 3600;
	if (d < 0)
		DEC *= -1;
}

void JPFITS::JPMath::RADecDegreeToSex(double RA_deg, double DEC_deg, String^ separator, String^ &ra_sex, String^ &dec_sex)
{
	double h = Math::Floor(RA_deg / 360 * 24);
	double m = Math::Floor((RA_deg / 360 * 24 - h) * 60);
	double s = Math::Round((RA_deg / 360 * 24 - h - m / 60) * 3600, 4);

	ra_sex = h.ToString() + separator + m.ToString() + separator + s.ToString();

	String^ sign = "+";
	if (DEC_deg < 0)
	{
		sign = "-";
		DEC_deg = Math::Abs(DEC_deg);
	}

	double d = Math::Floor(DEC_deg);
	double am = Math::Floor((DEC_deg - d) * 60);
	double as = Math::Round((DEC_deg - d - am / 60) * 3600, 4);

	dec_sex = sign + d.ToString() + separator + am.ToString() + separator + as.ToString();
}

array<double, 2>^ JPFITS::JPMath::Pad(array<double, 2>^ data, array<int>^ padding, bool do_parallel)
{
	int width = data->GetLength(0), height = data->GetLength(1);
	array<double, 2>^ result = gcnew array<double, 2>( width + padding[0] + padding[1], height + padding[2] + padding[3]);

	#pragma omp parallel for if (do_parallel)
	for (int x = 0; x < width; x++)
		for (int y = 0; y < height; y++)
			result[x + padding[0], y + padding[2]] = data[x, y];

	return result;
}

array<double, 2>^ JPFITS::JPMath::Crop(array<double, 2>^ data, array<int>^ cropping, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(cropping[1] - cropping[0] + 1, cropping[3] - cropping[2] + 1);
	int width = result->GetLength(0), height = result->GetLength(1);

	#pragma omp parallel for if (do_parallel)
	for (int x = 0; x < width; x++)
		for (int y = 0; y < height; y++)
			result[x, y] = data[x + cropping[0] - 1, y + cropping[2] - 1];

	return result;
}

array<double>^ JPFITS::JPMath::Bin(array<double>^ data, int Nx)
{
	if (Nx > data->GetLength(0))
		Nx = data->GetLength(0);
	int Lx = (int)System::Math::Floor((double)data->GetLength(0) / (double)Nx);

	array<double>^ result = gcnew array<double>(Lx);

	for (int i = 0; i < Lx; i++)
	{
		double s = 0;
		for (int k = i * Nx; k < i*Nx + Nx; k++)
			s = s + data[k];

		result[i] = s;
	}

	return result;
}

array<double, 2>^ JPFITS::JPMath::Bin(cli::array<double, 2> ^data, int Nx, int Ny, bool do_parallel)//remainders are dropped
{
	if (Nx > data->GetLength(0))
		Nx = data->GetLength(0);
	if (Ny > data->GetLength(1))
		Ny = data->GetLength(1);

	int Lx = (int)System::Math::Floor((double)data->GetLength(0) / (double)Nx);
	int Ly = (int)System::Math::Floor((double)data->GetLength(1) / (double)Ny);

	array<double, 2>^ result = gcnew array<double, 2>(Lx, Ly);
	double s = 0;

	#pragma omp parallel for if (do_parallel) private(s)
	for (int i = 0; i < Lx; i++)
		for (int j = 0; j < Ly; j++)
		{
			s = 0;
			for (int k = i * Nx; k < i*Nx + Nx; k++)
				for (int l = j * Ny; l < j*Ny + Ny; l++)
					s = s + data[k, l];

			result[i, j] = s;
		}

	return result;
}

array<unsigned int, 2>^ JPFITS::JPMath::Bin(array<unsigned int, 2> ^data, int Nx, int Ny, bool do_parallel)//remainders are dropped
{
	if (Nx > data->GetLength(0))
		Nx = data->GetLength(0);
	if (Ny > data->GetLength(1))
		Ny = data->GetLength(1);

	int Lx = (int)System::Math::Floor((double)data->GetLength(0) / (double)Nx);
	int Ly = (int)System::Math::Floor((double)data->GetLength(1) / (double)Ny);
	array<unsigned int, 2>^ result = gcnew array<unsigned int, 2>(Lx, Ly);
	unsigned int s = 0;

	#pragma omp parallel for if (do_parallel) private(s)
	for (int i = 0; i < Lx; i++)
		for (int j = 0; j < Ly; j++)
		{
			s = 0;
			for (int k = i * Nx; k < i*Nx + Nx; k++)
				for (int l = j * Ny; l < j*Ny + Ny; l++)
					s = s + data[k, l];

			result[i, j] = s;
		}

	return result;
}

array<double>^ JPFITS::JPMath::Min(array<double, 2> ^data, int dim, array<int>^ &indices, bool do_parallel)
{
	array<double>^ result = gcnew array<double>(data->GetLength(dim));
	//throw gcnew Exception("Dimension not x (0) or y (1).");

	double min = System::Double::MaxValue;

	if (dim == 0)//collapses array horizontally, i.e., makes a 'vertical' vector
	{
		#pragma omp parallel for if (do_parallel) private(min)
		for (int j = 0; j < data->GetLength(1); j++)
		{
			min = System::Double::MaxValue;
			for (int i = 0; i < data->GetLength(0); i++)
				if (data[i, j] < min)
				{
					min = data[i, j];
					result[j] = min;
					indices[j] = i;
				}
		}
	}

	if (dim == 1)//collpases array vertically...i.e., makes a 'horizontal' vector
	{
		#pragma omp parallel for if (do_parallel) private(min)
		for (int i = 0; i < data->GetLength(0); i++)
		{
			min = System::Double::MaxValue;
			for (int j = 0; j < data->GetLength(1); j++)
				if (data[i, j] < min)
				{
					min = data[i, j];
					result[i] = min;
					indices[i] = j;
				}
		}
	}

	return result;
}

double JPFITS::JPMath::Min(array<double, 2> ^data, int& x, int& y, bool do_parallel)
{
	double min = System::Double::MaxValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			if (data[i, j] < min)
			{
				#pragma omp critical
				{
					if (data[i, j] < min)
					{
						min = data[i, j];
						x = i;
						y = j;
					}
				}
			}

	return min;
}

double JPFITS::JPMath::Min(array<double, 2> ^data, bool do_parallel)
{
	double min = System::Double::MaxValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			if (data[i, j] < min)
			{
				#pragma omp critical
				{
					if (data[i, j] < min)
						min = data[i, j];
				}
			}

	return min;
}

double JPFITS::JPMath::Min(array<double> ^data, int& index, bool do_parallel)
{
	double min = System::Double::MaxValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->Length; i++)
		if (data[i] < min)
		{
			#pragma omp critical
			{
				if (data[i] < min)
				{
					min = data[i];
					index = i;
				}
			}
		}

	return min;
}

double JPFITS::JPMath::Min(array<double> ^data, bool do_parallel)
{
	double min = System::Double::MaxValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->Length; i++)
		if (data[i] < min)
		{
			#pragma omp critical
			{
				if (data[i] < min)
					min = data[i];
			}
		}

	return min;
}

void JPFITS::JPMath::MinSetIndexes(array<double>^ data, array<int>^ indexes, bool do_parallel)
{
	int datalength = data->Length;
	int n = indexes->Length;
	array<double>^ mins = gcnew array<double>(n);
	
	array<int>^ dataindexes = gcnew array<int>(datalength);
	for (int i = 0; i < datalength; i++)
		dataindexes[i] = i;

	int k = 0, l = 0;
	while (k < n)
	{
		if (!Double::IsNaN(data[l]))
		{
			mins[k] = data[l];
			indexes[k] = l;
			k++;
		}
		l++;
	}

	double tempmin;
	int tempind;

	#pragma omp parallel for if (do_parallel) private(tempmin, tempind)
	for (int i = l; i < datalength; i++)
		if (!Double::IsNaN(data[i]))
			for (int j = 0; j < n; j++)
				if (data[i] <= mins[j])
				{
					#pragma omp critical
					{
						if (data[i] <= mins[j])
						{
							tempmin = mins[j];
							tempind = indexes[j];
							mins[j] = data[i];
							indexes[j] = dataindexes[i];
							data[i] = tempmin;
							dataindexes[i] = tempind;
						}
					}
				}
}

void JPFITS::JPMath::MinMax(array<double, 2>^ data, double &min, double &max, bool do_parallel)
{
	min = System::Double::MaxValue;
	max = System::Double::MinValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
		{
			if (data[i, j] < min)
			{
				#pragma omp critical
				{
					if (data[i, j] < min)
						min = data[i, j];
				}
			}

			if (data[i, j] > max)
			{
				#pragma omp critical
				{
					if (data[i, j] > max)
						max = data[i, j];
				}
			}
		}
}

void JPFITS::JPMath::MinMax(array<double>^ data, double &min, double &max, bool do_parallel)
{
	min = System::Double::MaxValue;
	max = System::Double::MinValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
	{
		if (data[i] < min)
		{
			#pragma omp critical
			{
				if (data[i] < min)
					min = data[i];
			}
		}

		if (data[i] > max)
		{
			#pragma omp critical
			{
				if (data[i] > max)
					max = data[i];
			}
		}
	}
}

array<double>^ JPFITS::JPMath::Max(array<double, 2>^ data, int dim, array<int>^ &indices, bool do_parallel)
{
	array<double>^ result = gcnew array<double>(data->GetLength(dim));
	//throw gcnew Exception("Dimension not x (0) or y (1).");

	double max = System::Double::MinValue;

	if (dim == 0)//collapses array horizontally, i.e., makes a 'vertical' vector
	{
		#pragma omp parallel for if (do_parallel) private(max)
		for (int j = 0; j < data->GetLength(1); j++)
		{
			max = System::Double::MinValue;
			for (int i = 0; i < data->GetLength(0); i++)
				if (data[i, j] > max)
				{
					max = data[i, j];
					result[j] = max;
					indices[j] = i;
				}
		}
	}

	if (dim == 1)//collpases array vertically...i.e., makes a 'horizontal' vector
	{
		#pragma omp parallel for if (do_parallel) private(max)
		for (int i = 0; i < data->GetLength(0); i++)
		{
			max = System::Double::MinValue;
			for (int j = 0; j < data->GetLength(1); j++)
				if (data[i, j] > max)
				{
					max = data[i, j];
					result[i] = max;
					indices[i] = j;
				}
		}
	}

	return result;
}

double JPFITS::JPMath::Max(array<double, 2>^ data, int& x, int& y, bool do_parallel)
{
	double max = System::Double::MinValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			if (data[i, j] > max)
			{
				#pragma omp critical
				{
					if (data[i, j] > max)
					{
						max = data[i, j];
						x = i;
						y = j;
					}
				}
			}

	return max;
}

double JPFITS::JPMath::Max(array<double, 2>^ data, bool do_parallel)
{
	double max = System::Double::MinValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			if (data[i, j] > max)
			{
				#pragma omp critical
				{
					if (data[i, j] > max)
						max = data[i, j];
				}
			}

	return max;
}

double JPFITS::JPMath::Max(array<double> ^data, int& index, bool do_parallel)
{
	double max = System::Double::MinValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->Length; i++)
		if (data[i] > max)
		{
			#pragma omp critical
			{
				if (data[i] > max)
				{
					max = data[i];
					index = i;
				}
			}
		}

	return max;
}

double JPFITS::JPMath::Max(array<double> ^data, bool do_parallel)
{
	double max = System::Double::MinValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->Length; i++)
		if (data[i] > max)
		{
			#pragma omp critical
			{
				if (data[i] > max)
					max = data[i];
			}
		}

	return max;
}

double JPFITS::JPMath::Max(cli::array<double> ^data, int startIndex, int endIndex, int &maxIndex, bool do_parallel)
{
	if (startIndex < 0)
		startIndex = 0;
	if (endIndex > data->Length)
		endIndex = data->Length - 1;

	double max = System::Double::MinValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = startIndex; i < endIndex; i++)
		if (data[i] > max)
		{
			#pragma omp critical
			{
				if (data[i] > max)
				{
					max = data[i];
					maxIndex = i;
				}
			}
		}

	return max;
}

unsigned int JPFITS::JPMath::Max(array<unsigned int, 2>^ data, int& x, int& y, bool do_parallel)
{
	unsigned int max = System::UInt32::MinValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			if (data[i, j] > max)
			{
				#pragma omp critical
				{
					if (data[i, j] > max)
					{
						max = data[i, j];
						x = i;
						y = j;
					}
				}
			}

	return max;
}

double JPFITS::JPMath::Min(array<double> ^data, int startIndex, int endIndex, int &minIndex, bool do_parallel)
{
	if (startIndex < 0)
		startIndex = 0;
	if (endIndex > data->Length)
		endIndex = data->Length - 1;

	double min = System::Double::MaxValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = startIndex; i < endIndex; i++)
		if (data[i] < min)
		{
			#pragma omp critical
			{
				if (data[i] < min)
				{
					min = data[i];
					minIndex = i;
				}
			}
		}

	return min;
}

int JPFITS::JPMath::Max(array<int> ^data, int &index, bool do_parallel)
{
	int max = System::Int32::MinValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->Length; i++)
		if (data[i] > max)
		{
			#pragma omp critical
			{
				if (data[i] > max)
				{
					max = data[i];
					index = i;
				}
			}
		}

	return max;
}

int JPFITS::JPMath::Max(array<int> ^data, bool do_parallel)
{
	int max = System::Int32::MinValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->Length; i++)
		if (data[i] > max)
		{
			#pragma omp critical
			{
				if (data[i] > max)
					max = data[i];
			}
		}

	return max;
}

array<double, 1>^ JPFITS::JPMath::Stdv(array<double, 2>^ data, int dim, bool do_parallel)
{
	int d = 1;
	if (dim == 1)
		d = 0;

	array<double>^ mean = JPMath::Mean(data, dim, true);

	array<double, 1>^ result = gcnew array<double, 1>(data->GetLength(d));
	double dimOppL = (double)data->GetLength(d);
	double std = 0;

	if (dim == 0)//collapses array horizontally, i.e., makes a 'vertical' vector
	{
		#pragma omp parallel for if (do_parallel) private(std)
		for (int j = 0; j < data->GetLength(1); j++)
		{
			std = 0;
			for (int i = 0; i < data->GetLength(0); i++)
				std += (data[i, j] - mean[j]) * (data[i, j] - mean[j]);

			result[j] = System::Math::Sqrt(std / (data->GetLength(0) - 1.0));
		}
	}

	if (dim == 1)//collpases array vertically, i.e., makes a 'horizontal' vector
	{
		#pragma omp parallel for if (do_parallel) private(std)
		for (int i = 0; i < data->GetLength(0); i++)
		{
			std = 0;
			for (int j = 0; j < data->GetLength(1); j++)
				std += (data[i, j] - mean[i]) * (data[i, j] - mean[i]);

			result[i] = System::Math::Sqrt(std / (data->GetLength(1) - 1.0));
		}
	}

	return result;
}

double JPFITS::JPMath::Stdv(cli::array<double, 2> ^data, bool do_parallel)
{
	double mean = 0;
	double std = 0;
	double N = double(data->Length);
	#pragma omp parallel for if (do_parallel) reduction(+:mean)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			mean += data[i, j];

	mean = mean / N;

	#pragma omp parallel for if (do_parallel) reduction(+:std)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			std += (data[i, j] - mean)*(data[i, j] - mean);

	std = System::Math::Sqrt(std / (N - 1.0));

	return std;
}

double JPFITS::JPMath::Stdv(cli::array<double, 2> ^data, double known_mean, bool do_parallel)
{
	double std = 0;
	double N = double(data->Length);

	#pragma omp parallel for if (do_parallel) reduction(+:std)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			std += (data[i, j] - known_mean)*(data[i, j] - known_mean);

	std = System::Math::Sqrt(std / (N - 1.0));

	return std;
}


double JPFITS::JPMath::Stdv(cli::array<double, 1> ^data, bool do_parallel)
{
	double mean = 0;
	double std = 0;
	double N = double(data->Length);
	#pragma omp parallel for if (do_parallel) reduction(+:mean)
	for (int i = 0; i < (int)N; i++)
		mean += data[i];

	mean = mean / N;

	#pragma omp parallel for if (do_parallel) reduction(+:std)
	for (int i = 0; i < (int)N; i++)
		std += (data[i] - mean)*(data[i] - mean);

	std = System::Math::Sqrt(std / (N - 1.0));

	return std;
}

double JPFITS::JPMath::Stdv(cli::array<double, 1> ^data, double known_mean, bool do_parallel)
{
	double std = 0;
	double N = double(data->Length);

	#pragma omp parallel for if (do_parallel) reduction(+:std)
	for (int i = 0; i < (int)N; i++)
		std += (data[i] - known_mean)*(data[i] - known_mean);

	std = System::Math::Sqrt(std / (N - 1.0));

	return std;
}

array<int, 2>^ JPFITS::JPMath::Find(array<double, 2>^ data, double val, String^ style, bool do_parallel)
{
	ArrayList^ ptslist = gcnew ArrayList();

	int method = 0;
	if (style->Equals("<"))
		method = 0;
	else if (style->Equals("<="))
		method = 1;
	else if (style->Equals("==") || style->Equals("="))
		method = 2;
	else if (style->Equals(">="))
		method = 3;
	else if (style->Equals(">"))
		method = 4;
	else if (style->Equals("!="))
		method = 5;
	else
		throw gcnew Exception("Error:  Search style '" + style + "' not meaningful.");

	if (Double::IsNaN(val))
		method = 6;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
		{
			switch (method)
			{
				case 0://<
				{
					if (data[i, j] < val)
					{
						#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}
					}
					break;
				}
				case 1://<=
				{
					if (data[i, j] <= val)
					{
						#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}
					}
					break;
				}
				case 2://== || =
				{
					if (data[i, j] == val)
					{
						#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}
					}
					break;
				}
				case 3://>=
				{
					if (data[i, j] >= val)
					{
						#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}
					}
					break;
				}
				case 4://>
				{
					if (data[i, j] > val)
					{
						#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}
					}
					break;
				}
				case 5://!=
				{
					if (data[i, j] != val)
					{
						#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}
					}
					break;
				}
				case 6://= NaN
				{
					if (Double::IsNaN(data[i, j]))
					{
						#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}
					}
					break;
				}
			}
		}

	int N = int(double(ptslist->Count) / 2);
	array<int, 2>^ result = gcnew array<int, 2>(N, 2);

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < N; i++)
	{
		result[i, 0] = System::Convert::ToInt32(ptslist[i * 2]);//x
		result[i, 1] = System::Convert::ToInt32(ptslist[i * 2 + 1]);//y
	}

	return result;
}

array<int, 1>^ JPFITS::JPMath::Find(array<double, 1>^ data, double val, String^ style)//style can be <, <=, ==, >=, >
{
	ArrayList^ ptslist = gcnew ArrayList();

	int method = 0;
	if (style->Equals("<"))
		method = 0;
	else if (style->Equals("<="))
		method = 1;
	else if (style->Equals("==") || style->Equals("="))
		method = 2;
	else if (style->Equals(">="))
		method = 3;
	else if (style->Equals(">"))
		method = 4;
	else if (style->Equals("!="))
		method = 5;
	else
		throw gcnew Exception("Error:  Search style '" + style + "' not meaningful.");

	for (int i = 0; i < data->Length; i++)
	{
		switch (method)
		{
		case 0://<
		{
			if (data[i] < val)
				ptslist->Add(i);
			break;
		}
		case 1://<=
		{
			if (data[i] <= val)
				ptslist->Add(i);
			break;
		}
		case 2://==
		{
			if (data[i] == val)
				ptslist->Add(i);
			break;
		}
		case 3://>=
		{
			if (data[i] >= val)
				ptslist->Add(i);
			break;
		}
		case 4://>
		{
			if (data[i] > val)
				ptslist->Add(i);
			break;
		}
		case 5://!=
		{
			if (data[i] != val)
				ptslist->Add(i);
			break;
		}
		}
	}

	int N = ptslist->Count;
	array<int, 1>^ result = gcnew array<int, 1>(N);

	for (int i = 0; i < N; i++)
		result[i] = System::Convert::ToInt32(ptslist[i]);

	return result;
}

array<int, 1>^ JPFITS::JPMath::Find(array<double, 1>^ data, double val, String^ style, int startindex)//style can be <, <=, ==, >=, >, !=
{
	if (startindex < 0)
		startindex = 0;

	ArrayList^ ptslist = gcnew ArrayList();

	int method = 0;
	if (style->Equals("<"))
		method = 0;
	else if (style->Equals("<="))
		method = 1;
	else if (style->Equals("==") || style->Equals("="))
		method = 2;
	else if (style->Equals(">="))
		method = 3;
	else if (style->Equals(">"))
		method = 4;
	else if (style->Equals("!="))
		method = 5;
	else
		throw gcnew Exception("Error:  Search style '" + style + "' not meaningful.");

	for (int i = startindex; i < data->Length; i++)
	{
		switch (method)
		{
		case 0://<
		{
			if (data[i] < val)
				ptslist->Add(i);
			break;
		}
		case 1://<=
		{
			if (data[i] <= val)
				ptslist->Add(i);
			break;
		}
		case 2://==
		{
			if (data[i] == val)
				ptslist->Add(i);
			break;
		}
		case 3://>=
		{
			if (data[i] >= val)
				ptslist->Add(i);
			break;
		}
		case 4://>
		{
			if (data[i] > val)
				ptslist->Add(i);
			break;
		}
		case 5://!=
		{
			if (data[i] != val)
				ptslist->Add(i);
			break;
		}
		}
	}

	int N = ptslist->Count;
	array<int, 1>^ result = gcnew array<int, 1>(N);

	for (int i = 0; i < N; i++)
		result[i] = System::Convert::ToInt32(ptslist[i]);

	return result;
}

array<int>^ JPFITS::JPMath::Find(cli::array<double, 1> ^data, double val, System::String ^style, int startindex, int endindex)
{
	if (startindex < 0)
		startindex = 0;
	if (endindex > data->Length)
		endindex = data->Length - 1;

	ArrayList^ ptslist = gcnew ArrayList();

	int method = 0;
	if (style->Equals("<"))
		method = 0;
	else if (style->Equals("<="))
		method = 1;
	else if (style->Equals("==") || style->Equals("="))
		method = 2;
	else if (style->Equals(">="))
		method = 3;
	else if (style->Equals(">"))
		method = 4;
	else if (style->Equals("!="))
		method = 5;
	else
		throw gcnew Exception("Error:  Search style '" + style + "' not meaningful.");

	for (int i = startindex; i <= endindex; i++)
	{
		switch (method)
		{
		case 0://<
		{
			if (data[i] < val)
				ptslist->Add(i);
			break;
		}
		case 1://<=
		{
			if (data[i] <= val)
				ptslist->Add(i);
			break;
		}
		case 2://==
		{
			if (data[i] == val)
				ptslist->Add(i);
			break;
		}
		case 3://>=
		{
			if (data[i] >= val)
				ptslist->Add(i);
			break;
		}
		case 4://>
		{
			if (data[i] > val)
				ptslist->Add(i);
			break;
		}
		case 5://!=
		{
			if (data[i] != val)
				ptslist->Add(i);
			break;
		}
		}
	}

	int N = ptslist->Count;
	array<int, 1>^ result = gcnew array<int, 1>(N);

	for (int i = 0; i < N; i++)
		result[i] = System::Convert::ToInt32(ptslist[i]);

	return result;
}

int JPFITS::JPMath::Find(array<double, 1>^ data, double val, String^ style, bool return_first_true_last_false)
{
	int method = 0;
	if (style == "<")
		method = 0;
	else if (style == "<=")
		method = 1;
	else if (style == "==" || style == "=")
		method = 2;
	else if (style == ">=")
		method = 3;
	else if (style == ">")
		method = 4;
	else if (style == "!=")
		method = 5;
	else
		throw gcnew Exception("Error:  Search style '" + style + "' not meaningful.");

	int index = -1;

	if (return_first_true_last_false)
		for (int i = 0; i < data->Length; i++)
		{
			switch (method)
			{
			case 0://<
			{
				if (data[i] < val)
					index = i;
				break;
			}
			case 1://<=
			{
				if (data[i] <= val)
					index = i;
				break;
			}
			case 2://==
			{
				if (data[i] == val)
					index = i;
				break;
			}
			case 3://>=
			{
				if (data[i] >= val)
					index = i;
				break;
			}
			case 4://>
			{
				if (data[i] > val)
					index = i;
				break;
			}
			case 5://!=
			{
				if (data[i] != val)
					index = i;
				break;
			}
			}
			if (index != -1)
				break;
		}
	else
		for (int i = data->Length - 1; i >= 0; i--)
		{
			switch (method)
			{
			case 0://<
			{
				if (data[i] < val)
					index = i;
				break;
			}
			case 1://<=
			{
				if (data[i] <= val)
					index = i;
				break;
			}
			case 2://==
			{
				if (data[i] == val)
					index = i;
				break;
			}
			case 3://>=
			{
				if (data[i] >= val)
					index = i;
				break;
			}
			case 4://>
			{
				if (data[i] > val)
					index = i;
				break;
			}
			case 5://!=
			{
				if (data[i] != val)
					index = i;
				break;
			}
			}
			if (index != -1)
				break;
		}

	return index;
}

array<double, 2>^ JPFITS::JPMath::Replace(array<double, 2>^ data, array<int, 2>^ coords, double val, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	Array::Copy(data, result, data->GetLength(0)*data->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < coords->GetLength(0); i++)
		result[coords[i, 0], coords[i, 1]] = val;

	return result;
}

array<double>^ JPFITS::JPMath::Replace(array<double, 1>^ data, array<int, 1>^ coords, double val)
{
	array<double>^ result = data;
	for (int i = 0; i < coords->Length; i++)
		result[coords[i]] = val;

	return result;
}

double JPFITS::JPMath::Sum(array<double, 1> ^data, bool do_parallel)
{
	double S = 0;
	#pragma omp parallel for if (do_parallel) reduction(+:S)
	for (int i = 0; i < data->Length; i++)
		S += data[i];

	return S;
}

int JPFITS::JPMath::Sum(array<int, 1> ^data, bool do_parallel)
{
	int S = 0;
	#pragma omp parallel for if (do_parallel) reduction(+:S)
	for (int i = 0; i < data->Length; i++)
		S += data[i];

	return S;
}

double JPFITS::JPMath::Sum(array<double, 2>^ data, bool do_parallel)
{
	double S = 0;
	#pragma omp parallel for if (do_parallel) reduction(+:S)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			S += data[i, j];

	return S;
}

int JPFITS::JPMath::Sum(array<int, 2>^ data, bool do_parallel)
{
	int S = 0;
	#pragma omp parallel for if (do_parallel) reduction(+:S)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			S += data[i, j];

	return S;
}

array<double, 1>^ JPFITS::JPMath::Sum(array<double, 2>^ data, int dim, bool do_parallel)//sum allong a dimension
{
	int d = 1;
	if (dim == 1)
		d = 0;

	array<double, 1>^ result = gcnew array<double, 1>(data->GetLength(d));
	double S;

	if (dim == 0)//collapses array horizontally, i.e., makes a 'vertical' vector
	{
		#pragma omp parallel for if(do_parallel) private(S)
		for (int j = 0; j < data->GetLength(1); j++)
		{
			S = 0;
			for (int i = 0; i < data->GetLength(0); i++)
				S += data[i, j];

			result[j] = S;
		}
	}

	if (dim == 1)//collpases array vertically...i.e., makes a 'horizontal' vector
	{
		#pragma omp parallel for if(do_parallel) private(S)
		for (int i = 0; i < data->GetLength(0); i++)
		{
			S = 0;
			for (int j = 0; j < data->GetLength(1); j++)
				S += data[i, j];

			result[i] = S;
		}
	}

	return result;
}

array<int, 1>^ JPFITS::JPMath::Sum(array<int, 2>^ data, int dim, bool do_parallel)//sum allong a dimension
{
	int d = 1;
	if (dim == 1)
		d = 0;

	array<int, 1>^ result = gcnew array<int, 1>(data->GetLength(d));
	int S;

	if (dim == 0)//collapses array horizontally, i.e., makes a 'vertical' vector
	{
		#pragma omp parallel for if(do_parallel) private(S)
		for (int j = 0; j < data->GetLength(1); j++)
		{
			S = 0;
			for (int i = 0; i < data->GetLength(0); i++)
				S += data[i, j];

			result[j] = S;
		}
	}

	if (dim == 1)//collpases array vertically...i.e., makes a 'horizontal' vector
	{
		#pragma omp parallel for if(do_parallel) private(S)
		for (int i = 0; i < data->GetLength(0); i++)
		{
			S = 0;
			for (int j = 0; j < data->GetLength(1); j++)
				S += data[i, j];

			result[i] = S;
		}
	}

	return result;
}

double JPFITS::JPMath::Mean(array<double, 1> ^data, bool do_parallel)
{
	double M = 0;
	#pragma omp parallel for if (do_parallel) reduction(+:M)
	for (int i = 0; i < data->Length; i++)
		M += data[i];

	M = M / (double)data->Length;
	return M;
}

double JPFITS::JPMath::Mean(array<double, 2>^ data, bool do_parallel)
{
	double M = 0;
	#pragma omp parallel for if (do_parallel) reduction(+:M)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			M += data[i, j];

	M = M / (double)data->Length;
	return M;
}


array<double, 1>^ JPFITS::JPMath::Mean(array<double, 2>^ data, int dim, bool do_parallel)
{
	int d = 1;
	if (dim == 1)
		d = 0;

	array<double, 1>^ result = gcnew array<double, 1>(data->GetLength(d));
	double dimOppL = (double)data->GetLength(d);
	double M;

	if (dim == 0)//collapses array horizontally, i.e., makes a 'vertical' vector
	{
		#pragma omp parallel for if (do_parallel) private(M)
		for (int j = 0; j < data->GetLength(1); j++)
		{
			M = 0;
			for (int i = 0; i < data->GetLength(0); i++)
				M += data[i, j];

			result[j] = M / dimOppL;
		}
	}

	if (dim == 1)//collpases array vertically, i.e., makes a 'horizontal' vector
	{
		#pragma omp parallel for if (do_parallel) private(M)
		for (int i = 0; i < data->GetLength(0); i++)
		{
			M = 0;
			for (int j = 0; j < data->GetLength(1); j++)
				M += data[i, j];

			result[i] = M / dimOppL;
		}
	}

	return result;
}

double JPFITS::JPMath::Mean_RobustClipped(array<double>^ data, double sigma)
{
	array<double>^ clipper = data;
	double m = Mean(clipper, true);
	double s = Stdv(clipper, true);

	array<int>^ pts = Find(Abs(VectorSubScalar(clipper, m, true), true), sigma*s, ">");

	while (pts->Length > 0)
	{
		clipper = Replace(clipper, pts, Median(clipper));
		m = Mean(clipper, true);
		s = Stdv(clipper, true);
		pts = Find(Abs(VectorSubScalar(clipper, m, true), true), sigma*s, ">");
	}
	return m;
}

double JPFITS::JPMath::Mean_RobustClipped(array<double, 2>^ data, double sigma)
{
	array<double, 2>^ clipper = data;
	double m = Mean(clipper, true);
	double s = Stdv(clipper, true);

	array<int, 2>^ pts = Find(Abs(MatrixSubScalar(clipper, m, true), true), sigma*s, ">", true);

	while (pts->Length > 0)
	{
		clipper = Replace(clipper, pts, Median(clipper), true);
		m = Mean(clipper, true);
		s = Stdv(clipper, true);
		pts = Find(Abs(MatrixSubScalar(clipper, m, true), true), sigma*s, ">", true);
	}
	return m;
}

inline double JPFITS::JPMath::Median(array<double>^ data)
{
	int len = data->Length;
	array<double>^ lindata = gcnew array<double>(len);

	::Array::Copy(data, lindata, len);

	pin_ptr<double> p = &lindata[0];

	return MedianSTD(p, len);
}

inline double JPFITS::JPMath::Median(array<double, 2>^ data)
{
	int lenx = data->GetLength(0);
	int leny = data->GetLength(1);
	array<double, 2>^ lindata = gcnew array<double, 2>(lenx, leny);

	::Array::Copy(data, lindata, lenx*leny);

	pin_ptr<double> p = &lindata[0, 0];

	return MedianSTD(p, lenx*leny);
}

inline double JPFITS::JPMath::MedianSTD(double *np, int len)
{
	int n = (len - 1) / 2;

	std::nth_element(np, np + n, np + len);

	double m = np[n];

	if (len & 1)
		return m;

	std::nth_element(np, np + ++n, np + len);

	return (m + np[n]) / 2;
}

array<double, 1>^ JPFITS::JPMath::CosineBell(int length)
{
	array<double>^ result = gcnew array<double>(length);
	double L = (double)length;

	for (int i = 0; i < int(L); i++)
		result[i] = -(1 + Math::Cos((double(i) / (L - 1)) * 2 * Math::PI))*.5 + 1;

	return result;
}

array<double, 1>^ JPFITS::JPMath::Hanning(array<double, 1>^ data)
{
	array<double, 1>^ result = gcnew array<double, 1>(data->Length);
	double L = (double)data->Length;
	double bell;
	for (int i = 0; i < int(L); i++)
	{
		bell = -(1 + Math::Cos((double(i) / (L - 1)) * 2 * Math::PI))*.5 + 1;
		result[i] = data[i] * bell;
	}

	return result;
}

array<double, 2>^ JPFITS::JPMath::Hanning(array<double, 2>^ data, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	double HL = (double)data->GetLength(0);
	double VL = (double)data->GetLength(1);
	double bell;

	if (HL > 1 && VL > 1)
	{
		#pragma omp parallel for if (do_parallel) private(bell)
		for (int x = 0; x < data->GetLength(0); x++)
			for (int y = 0; y < data->GetLength(1); y++)
			{
				bell = -(1 + Math::Cos((double(x) / (HL - 1)) * 2 * Math::PI))*.5 + 1;
				bell *= (-(1 + Math::Cos((double(y) / (VL - 1)) * 2 * Math::PI))*.5 + 1);
				result[x, y] = data[x, y] * bell;
			}
	}

	if (HL == 1)
	{
		#pragma omp parallel for if (do_parallel) private(bell)
		for (int y = 0; y < data->GetLength(1); y++)
		{
			bell = (-(1 + Math::Cos((double(y) / (VL - 1)) * 2 * Math::PI))*.5 + 1);
			result[0, y] = data[0, y] * bell;
		}
	}

	if (VL == 1)
	{
		#pragma omp parallel for if (do_parallel) private(bell)
		for (int x = 0; x < data->GetLength(0); x++)
		{
			bell = -(1 + Math::Cos((double(x) / (HL - 1)) * 2 * Math::PI))*.5 + 1;
			result[x, 0] = data[x, 0] * bell;
		}
	}

	return result;
}

array<double>^ JPFITS::JPMath::XCorr(array<double>^ reference, array<double>^ relative, array<int>^ lags, bool do_parallel)
{
	int L = reference->Length;
	array<double>^ result = gcnew array<double>(2 * L - 1);

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < 2 * L - 1; i++)
	{
		if (i < L)
			for (int j = 0; j < i + 1; j++)
				result[i] += reference[L - 1 - i + j] * relative[j];
		else
			for (int j = 0; j < 2 * L - i - 1; j++)
				result[i] += reference[j] * relative[i - L + j + 1];

		lags[i] = -L + 1 + i;
	}

	return result;
}

void JPFITS::JPMath::XCorrImageLagShifts(array<double, 2>^ reference, array<double, 2>^ COMPARISON, bool autoDeBias_refX, bool autoDeBias_refY, bool autoDeBias_COMX, bool autoDeBias_COMY, bool autoHanning_ref, bool autoHanning_COM, double& xshift, double& yshift, bool do_parallel)
{
	array<double, 2>^ refref = gcnew array<double, 2>(reference->GetLength(0), reference->GetLength(1));
	array<double, 2>^ COMCOM = gcnew array<double, 2>(COMPARISON->GetLength(0), COMPARISON->GetLength(1));
	Array::Copy(reference, refref, reference->LongLength);
	Array::Copy(COMPARISON, COMCOM, COMCOM->LongLength);

	if (autoDeBias_refX)
		refref = DeGradient(refref, 0, do_parallel);
	if (autoDeBias_refY)
		refref = DeGradient(refref, 1, do_parallel);
	if (autoDeBias_COMX)
		COMCOM = DeGradient(COMCOM, 0, do_parallel);
	if (autoDeBias_COMY)
		COMCOM = DeGradient(COMCOM, 1, do_parallel);
	if (autoHanning_ref)
		refref = JPFITS::JPMath::Hanning(refref, do_parallel);
	if (autoHanning_COM)
		COMCOM = JPFITS::JPMath::Hanning(COMCOM, do_parallel);

	array<double>^ Href = JPFITS::JPMath::Sum(refref, 1, do_parallel);
	array<double>^ Vref = JPFITS::JPMath::Sum(refref, 0, do_parallel);
	array<double>^ HCOM = JPFITS::JPMath::Sum(COMCOM, 1, do_parallel);
	array<double>^ VCOM = JPFITS::JPMath::Sum(COMCOM, 0, do_parallel);
	double meanref = JPFITS::JPMath::Mean(Href, do_parallel);
	double meanCOM = JPFITS::JPMath::Mean(HCOM, do_parallel);
	Href = JPFITS::JPMath::VectorSubScalar(Href, meanref, do_parallel);
	Vref = JPFITS::JPMath::VectorSubScalar(Vref, meanref, do_parallel);
	HCOM = JPFITS::JPMath::VectorSubScalar(HCOM, meanCOM, do_parallel);
	VCOM = JPFITS::JPMath::VectorSubScalar(VCOM, meanCOM, do_parallel);

	array<int>^ Hcorr_Lag = gcnew array<int>(2 * HCOM->Length - 1);
	array<int>^ Vcorr_Lag = gcnew array<int>(2 * VCOM->Length - 1);
	array<double>^ Hcorr_Amp = JPFITS::JPMath::XCorr(Href, HCOM, Hcorr_Lag, do_parallel);
	array<double>^ Vcorr_Amp = JPFITS::JPMath::XCorr(Vref, VCOM, Vcorr_Lag, do_parallel);

	int maxHAmpindex;
	double maxH = JPFITS::JPMath::Max(Hcorr_Amp, maxHAmpindex, do_parallel);
	int xmax = Hcorr_Lag[maxHAmpindex];
	int maxVAmpindex;
	double maxV = JPFITS::JPMath::Max(Vcorr_Amp, maxVAmpindex, do_parallel);
	int ymax = Vcorr_Lag[maxVAmpindex];

	array<double>^ Hx = gcnew array<double>(3) { double(xmax - 1), double(xmax), double(xmax + 1) };
	array<double>^ Hy = gcnew array<double>(3) { Hcorr_Amp[maxHAmpindex - 1], Hcorr_Amp[maxHAmpindex], Hcorr_Amp[maxHAmpindex + 1] };
	array<double>^ Vx = gcnew array<double>(3) { double(ymax - 1), double(ymax), double(ymax + 1) };
	array<double>^ Vy = gcnew array<double>(3) { Vcorr_Amp[maxVAmpindex - 1], Vcorr_Amp[maxVAmpindex], Vcorr_Amp[maxVAmpindex + 1] };

	//array<double>^ shift = gcnew array<double>(2);//return value, x,y
	xshift = JPFITS::JPMath::QuadFit3PtsCenterPos(Hx, Hy);
	yshift = JPFITS::JPMath::QuadFit3PtsCenterPos(Vx, Vy);

	//return shift;
}

void JPFITS::JPMath::XCorrImageLagShifts(array<double>^ referenceX, array<double>^ referenceY, array<double, 2>^ COMPARISON, bool autoDeBias_COMX, bool autoDeBias_COMY, bool autoHanning_COM, double& xshift, double& yshift, bool do_parallel)
{
	array<double, 2>^ COMCOM = gcnew array<double, 2>(COMPARISON->GetLength(0), COMPARISON->GetLength(1));
	Array::Copy(COMPARISON, COMCOM, COMCOM->LongLength);

	if (autoDeBias_COMX)
		COMCOM = DeGradient(COMCOM, 0, do_parallel);
	if (autoDeBias_COMY)
		COMCOM = DeGradient(COMCOM, 1, do_parallel);
	if (autoHanning_COM)
		COMCOM = JPFITS::JPMath::Hanning(COMCOM, do_parallel);

	array<double>^ HCOM = JPFITS::JPMath::Sum(COMCOM, 1, do_parallel);
	array<double>^ VCOM = JPFITS::JPMath::Sum(COMCOM, 0, do_parallel);
	double meanCOM = JPFITS::JPMath::Mean(HCOM, do_parallel);
	HCOM = JPFITS::JPMath::VectorSubScalar(HCOM, meanCOM, do_parallel);
	VCOM = JPFITS::JPMath::VectorSubScalar(VCOM, meanCOM, do_parallel);

	array<int>^ Hcorr_Lag = gcnew array<int>(2 * HCOM->Length - 1);
	array<double>^ Hcorr_Amp = JPFITS::JPMath::XCorr(referenceX, HCOM, Hcorr_Lag, do_parallel);
	int maxHAmpindex = 0;
	double maxH = JPFITS::JPMath::Max(Hcorr_Amp, maxHAmpindex, do_parallel);
	int xmax = Hcorr_Lag[maxHAmpindex];

	array<int>^ Vcorr_Lag = gcnew array<int>(2 * VCOM->Length - 1);
	array<double>^ Vcorr_Amp = JPFITS::JPMath::XCorr(referenceY, VCOM, Vcorr_Lag, do_parallel);
	int maxVAmpindex = 0;
	double maxV = JPFITS::JPMath::Max(Vcorr_Amp, maxVAmpindex, do_parallel);
	int ymax = Vcorr_Lag[maxVAmpindex];

	array<double>^ Hx = gcnew array<double>(3) { double(xmax - 1), double(xmax), double(xmax + 1) };
	array<double>^ Hy = gcnew array<double>(3) { Hcorr_Amp[maxHAmpindex - 1], Hcorr_Amp[maxHAmpindex], Hcorr_Amp[maxHAmpindex + 1] };
	array<double>^ Vx = gcnew array<double>(3) { double(ymax - 1), double(ymax), double(ymax + 1) };
	array<double>^ Vy = gcnew array<double>(3) { Vcorr_Amp[maxVAmpindex - 1], Vcorr_Amp[maxVAmpindex], Vcorr_Amp[maxVAmpindex + 1] };

	//array<double>^ shift = gcnew array<double>(2);//return value, x,y
	xshift = JPFITS::JPMath::QuadFit3PtsCenterPos(Hx, Hy);
	yshift = JPFITS::JPMath::QuadFit3PtsCenterPos(Vx, Vy);

	//return shift;
}

double JPFITS::JPMath::QuadFit3PtsCenterPos(array<double, 1> ^x, array<double, 1> ^y)
{
	double x1 = x[0];
	double x2 = x[1];
	double x3 = x[2];
	double y1 = y[0];
	double y2 = y[1];
	double y3 = y[2];

	double det = x1 * x1*(x2 - x3) - x2 * x2*(x1 - x3) + x3 * x3*(x1 - x2);
	double a1 = (y1*(x2 - x3) - y2 * (x1 - x3) + y3 * (x1 - x2)) / det;
	double b1 = (x1*x1*(y2 - y3) - x2 * x2*(y1 - y3) + x3 * x3*(y1 - y2)) / det;
	double c1 = (x1*x1*(x2*y3 - x3 * y2) - x2 * x2*(x1*y3 - x3 * y1) + x3 * x3*(x1*y2 - x2 * y1)) / det;

	double X0 = -0.5*(b1 / a1);

	return X0;
}

array<double>^ JPFITS::JPMath::QuadFit3PtsParams(array<double, 1> ^x, array<double, 1> ^y)
{
	array<double>^ result = gcnew array<double>(3);
	double x1 = x[0];
	double x2 = x[1];
	double x3 = x[2];
	double y1 = y[0];
	double y2 = y[1];
	double y3 = y[2];

	double det = x1 * x1*(x2 - x3) - x2 * x2*(x1 - x3) + x3 * x3*(x1 - x2);
	result[0] = (y1*(x2 - x3) - y2 * (x1 - x3) + y3 * (x1 - x2)) / det;
	result[1] = (x1*x1*(y2 - y3) - x2 * x2*(y1 - y3) + x3 * x3*(y1 - y2)) / det;
	result[2] = (x1*x1*(x2*y3 - x3 * y2) - x2 * x2*(x1*y3 - x3 * y1) + x3 * x3*(x1*y2 - x2 * y1)) / det;

	return result;
}

double JPFITS::JPMath::VectorDotProdVector(array<double, 1> ^vector1, array<double, 1> ^vector2, bool do_parallel)
{
	double result = 0;
	#pragma omp parallel for if (do_parallel) reduction(+:result)
	for (int i = 0; i < vector1->Length; i++)
		result += vector1[i] * vector2[i];

	return result;
}

array<double, 1>^ JPFITS::JPMath::VectorDivScalar(array<double, 1> ^vector, double scalar, bool do_parallel)
{
	array<double, 1>^ result = gcnew array<double, 1>(vector->Length);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < vector->Length; i++)
		result[i] = vector[i] / scalar;

	return result;
}

array<double, 1>^ JPFITS::JPMath::VectorDivVector(array<double, 1> ^vector1, array<double, 1> ^vector2, bool do_parallel)
{
	array<double, 1>^ result = gcnew array<double, 1>(vector1->Length);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < vector1->Length; i++)
		result[i] = vector1[i] / vector2[i];

	return result;
}

array<double, 2>^ JPFITS::JPMath::MatrixDivScalar(array<double, 2> ^matrix, double scalar, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(matrix->GetLength(0), matrix->GetLength(1));
	scalar = 1 / scalar;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < matrix->GetLength(0); i++)
		for (int j = 0; j < matrix->GetLength(1); j++)
			result[i, j] = matrix[i, j] * scalar;

	return result;
}

array<double, 2>^ JPFITS::JPMath::MatrixDivMatrix(array<double, 2> ^matrix1, array<double, 2> ^matrix2, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(matrix1->GetLength(0), matrix1->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < matrix1->GetLength(0); i++)
		for (int j = 0; j < matrix1->GetLength(1); j++)
			result[i, j] = matrix1[i, j] / matrix2[i, j];

	return result;
}

array<double, 1>^ JPFITS::JPMath::VectorMultScalar(array<double, 1> ^vector, double scalar, bool do_parallel)
{
	array<double, 1>^ result = gcnew array<double, 1>(vector->Length);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < vector->Length; i++)
		result[i] = vector[i] * scalar;

	return result;
}

array<double, 1>^ JPFITS::JPMath::VectorMultVector(array<double, 1> ^vector1, array<double, 1> ^vector2, bool do_parallel)
{
	array<double, 1>^ result = gcnew array<double, 1>(vector1->Length);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < vector1->Length; i++)
		result[i] = vector1[i] * vector2[i];

	return result;
}

array<double, 2>^ JPFITS::JPMath::MatrixMultScalar(array<double, 2> ^matrix, double scalar, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(matrix->GetLength(0), matrix->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < matrix->GetLength(0); i++)
		for (int j = 0; j < matrix->GetLength(1); j++)
			result[i, j] = matrix[i, j] * scalar;

	return result;
}

array<double, 2>^ JPFITS::JPMath::MatrixMultMatrix(array<double, 2> ^matrix1, array<double, 2> ^matrix2, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(matrix1->GetLength(0), matrix1->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < matrix1->GetLength(0); i++)
		for (int j = 0; j < matrix1->GetLength(1); j++)
			result[i, j] = matrix1[i, j] * matrix2[i, j];

	return result;
}

array<double, 1>^ JPFITS::JPMath::VectorAddScalar(array<double, 1> ^vector, double scalar, bool do_parallel)
{
	array<double, 1>^ result = gcnew array<double, 1>(vector->Length);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < vector->Length; i++)
		result[i] = vector[i] + scalar;

	return result;
}

array<double, 1>^ JPFITS::JPMath::VectorAddVector(array<double, 1> ^vector1, array<double, 1> ^vector2, bool do_parallel)
{
	array<double, 1>^ result = gcnew array<double, 1>(vector1->Length);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < vector1->Length; i++)
		result[i] = vector1[i] + vector2[i];
	return result;
}

array<double, 2>^ JPFITS::JPMath::MatrixAddScalar(array<double, 2> ^matrix, double scalar, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(matrix->GetLength(0), matrix->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < matrix->GetLength(0); i++)
		for (int j = 0; j < matrix->GetLength(1); j++)
			result[i, j] = matrix[i, j] + scalar;

	return result;
}

array<double, 2>^ JPFITS::JPMath::MatrixAddMatrix(array<double, 2> ^matrix1, array<double, 2> ^matrix2, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(matrix1->GetLength(0), matrix1->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < matrix1->GetLength(0); i++)
		for (int j = 0; j < matrix1->GetLength(1); j++)
			result[i, j] = matrix1[i, j] + matrix2[i, j];

	return result;
}

array<double, 1>^ JPFITS::JPMath::VectorSubScalar(array<double, 1> ^vector, double scalar, bool do_parallel)
{
	array<double, 1>^ result = gcnew array<double, 1>(vector->Length);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < vector->Length; i++)
		result[i] = vector[i] - scalar;

	return result;
}

array<double, 1>^ JPFITS::JPMath::VectorSubVector(array<double, 1> ^vector1, array<double, 1> ^vector2, bool do_parallel)
{
	array<double, 1>^ result = gcnew array<double, 1>(vector1->Length);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < vector1->Length; i++)
		result[i] = vector1[i] - vector2[i];

	return result;
}

array<double, 2>^ JPFITS::JPMath::MatrixSubScalar(array<double, 2> ^matrix, double scalar, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(matrix->GetLength(0), matrix->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < matrix->GetLength(0); i++)
		for (int j = 0; j < matrix->GetLength(1); j++)
			result[i, j] = matrix[i, j] - scalar;

	return result;
}

array<double, 2>^ JPFITS::JPMath::MatrixSubMatrix(array<double, 2> ^matrix1, array<double, 2> ^matrix2, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(matrix1->GetLength(0), matrix1->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < matrix1->GetLength(0); i++)
		for (int j = 0; j < matrix1->GetLength(1); j++)
			result[i, j] = matrix1[i, j] - matrix2[i, j];

	return result;
}

array<double, 1>^ JPFITS::JPMath::Abs(cli::array<double, 1> ^data, bool do_parallel)
{
	array<double, 1>^ result = gcnew array<double, 1>(data->Length);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->Length; i++)
		result[i] = System::Math::Abs(data[i]);

	return result;
}

array<double, 2>^ JPFITS::JPMath::Abs(cli::array<double, 2> ^data, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			result[i, j] = System::Math::Abs(data[i, j]);

	return result;
}

array<double, 2>^ JPFITS::JPMath::Round(cli::array<double, 2> ^data, int digits, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			result[i, j] = System::Math::Round(data[i, j], digits);

	return result;
}

array<double, 2>^ JPFITS::JPMath::Floor(cli::array<double, 2> ^data, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			result[i, j] = System::Math::Floor(data[i, j]);

	return result;
}

array<double, 2>^ JPFITS::JPMath::Floor(cli::array<double, 2> ^data, double clip_floor, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
		{
			if (data[i, j] < clip_floor)
				result[i, j] = clip_floor;
			else
				result[i, j] = data[i, j];
		}

	return result;
}

array<double, 2>^ JPFITS::JPMath::Ceil(cli::array<double, 2> ^data, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			result[i, j] = System::Math::Ceiling(data[i, j]);

	return result;
}

array<double, 2>^ JPFITS::JPMath::Power(cli::array<double, 2> ^data, double exponent, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			result[i, j] = System::Math::Pow(data[i, j], exponent);

	return result;
}

array<double, 2>^ JPFITS::JPMath::Sqrt(cli::array<double, 2> ^data, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			result[i, j] = System::Math::Sqrt(data[i, j]);

	return result;
}

array<double, 2>^ JPFITS::JPMath::Log(array<double, 2>^ data, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
		{
			if (data[i, j] <= 0)
				result[i, j] = 0;
			else
				result[i, j] = System::Math::Log10(data[i, j]);
		}

	return result;
}

array<double, 2>^ JPFITS::JPMath::Log(array<double, 2>^ data, double base, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
		{
			if (data[i, j] <= 0)
				result[i, j] = 0;
			else
				result[i, j] = System::Math::Log(data[i, j], base);
		}

	return result;
}

array<double, 2>^ JPFITS::JPMath::Ln(array<double, 2>^ data, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
		{
			if (data[i, j] <= 0)
				result[i, j] = 0;
			else
				result[i, j] = System::Math::Log(data[i, j]);
		}

	return result;
}

array<double, 2>^ JPFITS::JPMath::Exp(array<double, 2>^ data, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			result[i, j] = System::Math::Exp(data[i, j]);

	return result;
}

array<double, 2>^ JPFITS::JPMath::Exp(array<double, 2>^ data, double base, bool do_parallel)
{
	array<double, 2>^ result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			result[i, j] = System::Math::Pow(base, data[i, j]);

	return result;
}

bool JPFITS::JPMath::IsEven(int x)
{
	if (System::Math::IEEERemainder(double(x), 2) == 0)
		return true;
	else
		return false;
}

inline bool JPFITS::JPMath::IsNumeric(String^ x)
{
	try
	{
		Convert::ToDouble(x);
		return true;
	}
	catch (...)
	{
		return false;
	}
}

inline bool JPFITS::JPMath::IsInteger(double x)
{
	if (Math::Truncate(x) == x)
		return true;
	else
		return false;
}

/*array<double,2>^ JPFITS::JPMath::MatrixMedianFilter(array<double,2>^ data, int filter_width)
{
	//width must always be an odd valued number, for median filtering geometry to make sense
	if (IsEven(filter_width) || filter_width < 3)
		throw gcnew Exception("Width '" + filter_width + "' must be odd-valued and greater than 3.");

	int HW = (filter_width -1)/2;
	array<double,2>^ result = gcnew array<double,2>(data->GetLength(0),data->GetLength(1));
	Array::Copy(data, result, data->LongLength);

	#pragma omp parallel for
	for (int j = HW; j < data->GetLength(1)-HW; j++)
		for (int i = HW; i < data->GetLength(0)-HW; i++)
		{
			//make a subarray around this i,j point, compute the median, replace this i,j value with median
			array<double, 2>^ subim = gcnew array<double, 2>(filter_width, filter_width);
			int iii = -1;
			for (int ii = i-HW; ii <= i+HW; ii++)
			{
				iii++;
				int jjj = -1;
				for (int jj = j-HW; jj <= j+HW; jj++)
				{
					jjj++;
					subim[iii,jjj] = data[ii,jj];
				}
			}
			//done creating subarray around i,j point; now get median and assign it to i,j point
			result[i,j] = Median(subim);
		}

	return result;
}*/

array<double, 2>^ JPFITS::JPMath::MedianFilter(array<double, 2>^ data, int kernelHalfWidth, bool do_parallel)
{
	int szx = data->GetLength(0);
	int szy = data->GetLength(1);
	array<double, 2>^ result = gcnew array<double, 2>(szx, szy);
	szx -= kernelHalfWidth;
	szy -= kernelHalfWidth;
	int Nkernpix = (kernelHalfWidth * 2 + 1) * (kernelHalfWidth * 2 + 1);

	#pragma omp parallel for if (do_parallel)
	for (int x = kernelHalfWidth; x < szx; x++)
	{
		array<double>^ kernel = gcnew array<double>(Nkernpix);
		int kxmax = x + kernelHalfWidth;
		for (int y = kernelHalfWidth; y < szy; y++)
		{
			int i = 0;
			int kx = x - kernelHalfWidth;
			int kymax = y + kernelHalfWidth;
			for (kx; kx <= kxmax; kx++)
			{
				int ky = y - kernelHalfWidth;
				for (ky; ky <= kymax; ky++)
				{
					kernel[i] = data[kx, ky];
					i++;
				}
			}
			result[x, y] = Median(kernel);
		}
	}

	return result;
}

array<double, 2>^ JPFITS::JPMath::MatrixConvolveMatrix(array<double, 2>^ primary, array<double, 2>^ kernel, bool do_parallel)
{
	int KFWX = kernel->GetLength(0);
	int KFWY = kernel->GetLength(1);
	int KHWX = (KFWX - 1) / 2;
	int KHWY = (KFWY - 1) / 2;
	int PFWX = primary->GetLength(0);
	int PFWY = primary->GetLength(1);
	int PMX = PFWX - KHWX, PMY = PFWY - KHWY;

	array<double, 2>^ result = gcnew array<double, 2>(PFWX, PFWY);

	#pragma omp parallel for if (do_parallel)
	for (int x = KHWX; x < PMX; x++)
		for (int y = KHWY; y < PMY; y++)
		{
			int ky = y - KHWY;
			for (int xx = 0; xx < KFWX; xx++)
			{
				int kx = x - KHWX + xx;
				for (int yy = 0; yy < KFWY; yy++)
					result[x, y] += kernel[xx, yy] * primary[kx, ky + yy];
			}
		}

	return result;
}

array<double, 2>^ JPFITS::JPMath::Gaussian(double Amplitude, double FWHM, int HalfWidth, bool do_parallel)
{
	double sig = FWHM / (2 * System::Math::Sqrt(2 * System::Math::Log(2)));
	int width = HalfWidth * 2 + 1;
	array<double, 2>^ result = gcnew array<double, 2>(width, width);
	double twosigsq = 2 * sig * sig;

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < width; i++)
	{
		double imid = double(i) - HalfWidth;
		for (int j = 0; j < width; j++)
		{
			double jmid = double(j) - HalfWidth;
			result[i, j] = Amplitude * System::Math::Exp(-(imid * imid + jmid * jmid) / twosigsq);
		}
	}

	return result;
}

array<double, 2>^ JPFITS::JPMath::Sersic(double Effective_Radius_Re_PIXELS, double Effective_Radius_Io, double SersicFactor_n, int CalculationRadius_NRe)
{
	int width = int(CalculationRadius_NRe * Effective_Radius_Re_PIXELS) * 2 + 1;
	array<double, 2>^ result = gcnew array<double, 2>(width, width);
	int hw = (width - 1) / 2;

	/*#pragma omp parallel for
	for (int x = -hw; x <= hw; x++)
		for (int y = -hw; y <= hw; y++)
		{
			double r_sq = x * x + y * y;
			result[x, y] = 0;
		}*/

	return result;
}

array<double>^ JPFITS::JPMath::Smooth(array<double> ^data, int kernelsize, bool do_parallel)
{
	int rem = 0;
	Math::DivRem(kernelsize, 2, rem);
	if (rem == 0 || kernelsize == 1)
	return data;

	array<double>^ result = gcnew array<double>(data->Length);
	int kernelHW = (kernelsize - 1) / 2;

	double dkernelsize = (double)kernelsize;
	double val = 0;
	int jstart = 0;
	int jend = 0;
	int kern;
	int L = data->Length;

	#pragma omp parallel for if (do_parallel) private(val, kern)
	for (int i = 0; i < L; i++)
	{
		val = 0;

		kern = i - kernelHW;
		if (kern < 0)
		{
			kern = kern + kernelHW;
			for (int j = i - kern; j <= i + kern; j++)
				val += data[j];
			result[i] = val / (2 * kern + 1);
			continue;
		}

		kern = i + kernelHW;
		if (kern > L - 1)
		{
			kern = -(kern - kernelHW - L + 1);
			for (int j = i - kern; j <= i + kern; j++)
				val += data[j];
			result[i] = val / (2 * kern + 1);
			continue;
		}

		for (int j = i - kernelHW; j <= i + kernelHW; j++)
			val += data[j];
		result[i] = val / dkernelsize;
	}

	return result;
}

array<double, 2>^ JPFITS::JPMath::ShiftArrayInt(array<double, 2> ^data, int xshift, int yshift, bool do_parallel)
{
	int width = data->GetLength(0);
	int height = data->GetLength(1);
	array<double, 2>^ arr = gcnew array<double, 2>(width, height);

	if (xshift == 0 && yshift == 0)
		return data;

	int istart = 0;
	int iend = width;
	//if shift is +ve, then it is shifting right, so need blank columns at left,
	if (xshift > 0)
	{
		istart = xshift;
		iend = width;
	}
	//if shift is -ve, then it is shifting left, so need blank columns at right,
	if (xshift < 0)
	{
		istart = 0;
		iend = width + xshift;
	}

	int jstart = 0;
	int jend = height;
	//if shift is +ve, then it is shifting down, so need blank columns at top,
	if (yshift > 0)
	{
		jstart = yshift;
		jend = height;
	}
	//if shift is -ve, then it is shifting up, so need blank columns at bottom,
	if (yshift < 0)
	{
		jstart = 0;
		jend = height + yshift;
	}

	#pragma omp parallel for if (do_parallel)
	for (int i = istart; i < iend; i++)
		for (int j = jstart; j < jend; j++)
			arr[i, j] = data[i - xshift, j - yshift];

	return arr;
}

inline double JPFITS::JPMath::InterpolateBiLinear(array<double, 2>^ data, int width, int height, double x, double y)
{
	int xoldfloor = (int)Math::Floor(x);
	int yoldfloor = (int)Math::Floor(y);
	if (xoldfloor < 0 || yoldfloor < 0)
		return 0;

	int xoldciel = (int)Math::Ceiling(x);
	int yoldciel = (int)Math::Ceiling(y);
	if (xoldciel >= width || yoldciel >= height)
		return 0;

	double deltaX = x - (double)xoldfloor;
	double onemdeltaX = 1 - deltaX;
	double deltaY = y - (double)yoldfloor;

	//linearly interpolate horizontally between top neighbours
	double pixTop = onemdeltaX * data[xoldfloor, yoldfloor] + deltaX * data[xoldciel, yoldfloor];

	// linearly interpolate horizontally between bottom neighbours
	double pixBottom = onemdeltaX * data[xoldfloor, yoldciel] + deltaX * data[xoldciel, yoldciel];
	
	// linearly interpolate vertically between top and bottom interpolated results
	return (1 - deltaY) * pixTop + deltaY * pixBottom;
}

inline double JPFITS::JPMath::InterpolateLanczos(array<double, 2>^ data, int width, int height, double x, double y, int n)
{
	double xoldfloor = Math::Floor(x);
	if (xoldfloor <= n || xoldfloor >= width - n)
		return 0;
	double xoldfloorminisx = xoldfloor - x;

	double yoldfloor = Math::Floor(y);
	if (yoldfloor <= n || yoldfloor >= height - n)
		return 0;
	double yoldfloorminisy = yoldfloor - y;

	double w = 0, val = 0, Lx = 0, Ly = 0;
	for (int i = -n + 1; i <= n; i++)
	{
		Lx = Lanczos(double(i) + xoldfloorminisx, n);
		for (int j = -n + 1; j <= n; j++)
		{
			Ly = Lx * Lanczos(double(j) + yoldfloorminisy, n);
			w += Ly;
			val += data[(int)xoldfloor + i, (int)yoldfloor + j] * Ly;
		}
	}
	return val / w;
}

array<double, 2>^ JPFITS::JPMath::RotateShiftArray(array<double, 2>^ data, double radians, double x_center, double y_center, String^ style, int xshift, int yshift, bool do_parallel)
{
	int width = data->GetLength(0);
	int height = data->GetLength(1);
	double xmid = double(width) / 2;
	double ymid = double(height) / 2;
	if (x_center != Double::MaxValue)
		xmid = x_center;
	if (y_center != Double::MaxValue)
		ymid = y_center;
	array<double, 2>^ arr = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1));

	if (style->ToLower() == "nearest")
	{
		int xold = 0, yold = 0;
		#pragma omp parallel for if (do_parallel) private(xold, yold)
		for (int x = 0; x < width; x++)
			for (int y = 0; y < height; y++)
			{
				xold = (int)Math::Round(((double)x - xmid - (double)xshift) * Math::Cos(radians) - ((double)y - ymid - (double)yshift) * Math::Sin(radians) + xmid);
				yold = (int)Math::Round(((double)x - xmid - (double)xshift) * Math::Sin(radians) + ((double)y - ymid - (double)yshift) * Math::Cos(radians) + ymid);
				if (xold >= 0 && xold < width && yold >= 0 && yold < height)
					arr[x, y] = data[xold, yold];
			}
	}

	if (style->ToLower() == "bilinear")
	{
		double xold, yold;
		#pragma omp parallel for if (do_parallel) private(xold, yold)
		for (int x = 0; x < width; x++)
			for (int y = 0; y < height; y++)
			{
				xold = ((double)x - xmid - (double)xshift) * Math::Cos(radians) - ((double)y - ymid - (double)yshift) * Math::Sin(radians) + xmid;
				yold = ((double)x - xmid - (double)xshift) * Math::Sin(radians) + ((double)y - ymid - (double)yshift) * Math::Cos(radians) + ymid;
				arr[x, y] = InterpolateBiLinear(data, width, height, xold, yold);
			}
	}

	if (style->Contains("lanc"))
	{
		String^ sn = style->Substring(style->Length - 1);
		if (!JPMath::IsNumeric(sn))
			throw gcnew Exception("Lanczos order " + sn + " indeterminable. ");
		int n = Convert::ToInt32(sn);
		if (n < 3 || n > 5)
			throw gcnew Exception("Lanczos order " + n + " not computable. Allowed n = 3, 4, 5. ");

		double xold, yold;
		#pragma omp parallel for if (do_parallel) private(xold, yold)
		for (int x = n; x < width - n; x++)
			for (int y = n; y < height - n; y++)
			{
				xold = ((double)x - xmid - (double)xshift) * Math::Cos(radians) - ((double)y - ymid - (double)yshift) * Math::Sin(radians) + xmid;
				yold = ((double)x - xmid - (double)xshift) * Math::Sin(radians) + ((double)y - ymid - (double)yshift) * Math::Cos(radians) + ymid;
				arr[x, y] = InterpolateLanczos(data, width, height, xold, yold, n);
			}
	}

	return arr;
}

inline double JPFITS::JPMath::Lanczos(double x, int n)
{
	if (Math::Abs(x) > n)
		return 0;
	else
		return sinc_n(x) * sinc_n(x / double(n));
}

//normalized sinc function
inline double JPFITS::JPMath::sinc_n(double x)
{
	if (x == 0)
		return 1;
	else
		return Math::Sin(Math::PI * x) / Math::PI / x;
}

//un-normalized sinc function
inline double JPFITS::JPMath::sinc_un(double x) 
{
	if (x == 0)
		return 1;
	else
		return Math::Sin(x) / x;
}

array<double, 2>^ JPFITS::JPMath::Histogram_IntegerStep(array<double> ^values, int step)
{
	Array::Sort(values);
	double min = Math::Floor(values[0]);
	double max = Math::Ceiling(values[values->Length - 1]);

	int NDivs = int(Math::Floor((max - min + step) / step));
	array<double>^ posts = gcnew array<double>(NDivs + 1);
	int PL = posts->Length;
	for (int i = 0; i < PL; i++)
		posts[i] = min + double(step*i);

	array<double>^ histogram = gcnew array<double>(NDivs);

	int VL = values->Length;
	int ind = 0;
	for (int j = 0; j < VL; j++)
	{
	recheck:;
		if (values[j] >= posts[ind] && values[j] < posts[ind + 1])
			histogram[ind]++;
		else
		{
			ind++;
			goto recheck;
		}
	}

	array<double, 2>^ result = gcnew array<double, 2>(NDivs, 2);
	for (int i = 0; i < NDivs; i++)
	{
		result[i, 0] = posts[i];
		result[i, 1] = histogram[i];
	}

	return result;
}

array<double, 2>^ JPFITS::JPMath::Histogram_IntegerDivisions(array<double> ^values, int NDivs)
{
	Array::Sort(values);
	array<double>^ posts = gcnew array<double>(NDivs + 1);//histogram subsection bounds
	double min = Math::Floor(values[0]);
	double max = Math::Ceiling(values[values->Length - 1]);
	double step = (max - min) / (double)NDivs;

	for (int i = 0; i <= NDivs; i++)
		posts[i] = min + step * double(i);

	array<double>^ histogram = gcnew array<double>(NDivs);

	int ind = 0;
	int skipped = 0;
	for (int j = 0; j < values->Length; j++)
	{
	recheck:;
		if (values[j] < posts[0])
		{
			skipped++;
			continue;
		}
		if (values[j] >= posts[NDivs])
		{
			skipped++;
			continue;
		}
		if (values[j] >= posts[ind] && values[j] < posts[ind + 1])
		{
			histogram[ind]++;
		}
		else
		{
			ind++;
			goto recheck;
		}
	}

	array<double, 2>^ result = gcnew array<double, 2>(NDivs, 2);
	for (int i = 0; i < NDivs; i++)
	{
		result[i, 0] = posts[i];
		result[i, 1] = histogram[i];
	}

	return result;
}

array<double, 2>^ JPFITS::JPMath::DeGradient(cli::array<double, 2> ^data, int dim, bool do_parallel)
{
	int width = data->GetLength(0);
	int height = data->GetLength(1);
	array<double, 2>^ result = gcnew array<double, 2>(width, height);
	double med;

	if (dim == 0)
	{
		#pragma omp parallel for if (do_parallel) private(med)
		for (int ii = 0; ii < width; ii++)
		{
			array<double>^ column = gcnew array<double>(height);
			for (int jj = 0; jj < height; jj++)
			{
				column[jj] = data[ii, jj];
				result[ii, jj] = data[ii, jj];
			}

			med = Median(column);
			for (int kk = 0; kk < height; kk++)
				result[ii, kk] -= med;
		}
	}

	if (dim == 1)
	{
		#pragma omp parallel for if (do_parallel) private(med)
		for (int ii = 0; ii < height; ii++)
		{
			array<double>^ row = gcnew array<double>(width);
			for (int jj = 0; jj < width; jj++)
			{
				row[jj] = data[jj, ii];
				result[jj, ii] = data[jj, ii];
			}

			med = Median(row);
			for (int jj = 0; jj < width; jj++)
				result[jj, ii] -= med;
		}
	}

	return result;
}

double JPFITS::JPMath::aTanAbsoluteAngle(double x, double y)
{
	/*if (x >= 0)
		if (y >= 0)
			return Math::Atan(Math::Abs(y) / Math::Abs(x));
		else
			return -Math::Atan(Math::Abs(y) / Math::Abs(x));

	if (x < 0)
		if(y >= 0)
			return Math::PI - Math::Atan(Math::Abs(y) / Math::Abs(x));
		else
			return -(Math::PI - Math::Atan(Math::Abs(y) / Math::Abs(x)));

	return Double::NaN;*/

	return Math::Atan2(y, x);
}

array<double>^ JPFITS::JPMath::Interpolate1d(array<double>^ xdata, array<double>^ ydata, array<double>^ xinterp, String^ style, bool do_parallel)
{
	array<double>^ result = gcnew array<double>(xinterp->Length);

	if (style != "linear" && style != "cubic" && style != "mono" && style != "akima" && style != "catmullrom")
	{
		throw gcnew Exception("Interpolation style '" + style + "' not recognized.");
		return nullptr;
	}

	alglib::spline1dinterpolant^ sp = gcnew alglib::spline1dinterpolant();

	if (style == "linear")
		alglib::spline1dbuildlinear(xdata, ydata, sp);

	if (style == "cubic")
		alglib::spline1dbuildcubic(xdata, ydata, sp);

	if (style == "mono")
		alglib::spline1dbuildmonotone(xdata, ydata, sp);

	if (style == "akima")
		alglib::spline1dbuildakima(xdata, ydata, sp);

	if (style == "catmullrom")
		alglib::spline1dbuildcatmullrom(xdata, ydata, sp);

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < xinterp->Length; i++)
		result[i] = alglib::spline1dcalc(sp, xinterp[i]);

	return result;
}

array<double, 2>^ JPFITS::JPMath::Interpolate2d(array<double>^ xdata, array<double>^ ydata, array<double, 2>^ surfdata, int xinterpdelta_inv, int yinterpdelta_inv, array<double>^ xinterp, array<double>^ yinterp, bool do_parallel)
{
	int xw = surfdata->GetLength(0);
	int yh = surfdata->GetLength(1);
	array<double, 2>^ result = gcnew array<double, 2>(xw*xinterpdelta_inv, yh*yinterpdelta_inv);

	array<double>^ surfdata_vec = gcnew array<double>(surfdata->Length);
	for (int j = 0; j < yh; j++)
		for (int i = 0; i < xw; i++)
			surfdata_vec[j*xw + i] = surfdata[i, j];

	if (xdata == nullptr)
	{
		xdata = gcnew array<double>(xw);
		for (int i = 0; i < xw; i++)
			xdata[i] = double(i);
	}
	if (ydata == nullptr)
	{
		ydata = gcnew array<double>(yh);
		for (int i = 0; i < yh; i++)
			ydata[i] = double(i);
	}

	alglib::spline2dinterpolant^ s;
	alglib::spline2dbuildbicubicv(xdata, xw, ydata, yh, surfdata_vec, 1, s);

	if (xinterp == nullptr)
		xinterp = gcnew array<double>(xw*xinterpdelta_inv);
	for (int i = 0; i < xw*xinterpdelta_inv; i++)
		xinterp[i] = xdata[0] + double(i) / double(xinterpdelta_inv);
	if (yinterp == nullptr)
		yinterp = gcnew array<double>(yh*yinterpdelta_inv);
	for (int j = 0; j < yh*yinterpdelta_inv; j++)
		yinterp[j] = ydata[0] + double(j) / double(yinterpdelta_inv);

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < xw*xinterpdelta_inv; i++)
		for (int j = 0; j < yh*yinterpdelta_inv; j++)
			result[i, j] = alglib::spline2dcalc(s, xinterp[i], yinterp[j]);

	return result;
}

void JPFITS::JPMath::Gaussian1d(array<double>^ xdata, array<double>^ &G, array<double>^ p)
{
	int xw = G->Length;
	double xhw = double(xw - 1) / 2.0;

	if (xdata == nullptr)
	{
		xdata = gcnew array<double>(xw);
		for (int i = 0; i < xw; i++)
			xdata[i] = double(i) - xhw;
	}

	for (int i = 0; i < xw; i++)
		G[i] = p[0] * Math::Exp(-((xdata[i] - p[1])*(xdata[i] - p[1])) / (2 * p[2] * p[2])) + p[3];
}

void JPFITS::JPMath::Fit_Gaussian1d(array<double>^ xdata, array<double>^ Gdata, array<double>^ &p, array<double>^ p_lbnd, array<double>^ p_ubnd, array<double>^ p_err, array<double>^ fit_residuals)
{
	int xw = Gdata->Length;
	double xhw = double(xw - 1) / 2.0;

	if (xdata == nullptr)
	{
		xdata = gcnew array<double>(xw);
		for (int i = 0; i < xw; i++)
			xdata[i] = double(i) - xhw;
	}

	array<double, 2>^ x = gcnew array<double, 2>(xw, 1);
	for (int i = 0; i < xw; i++)
		x[i, 0] = xdata[i];

	alglib::ndimensional_pfunc^ pf = gcnew alglib::ndimensional_pfunc(alglib_Gauss_1d);
	alglib::ndimensional_pgrad ^ pg = gcnew alglib::ndimensional_pgrad(alglib_Gauss_1d_grad);
	alglib::ndimensional_rep^ rep;
	Object^ obj;
	alglib::lsfitstate^ state;
	alglib::lsfitreport^ report;
	double diffstep = 0.0001;
	double epsx = 0.000001;
	int maxits = 0;
	int info;
	array<double>^ scale;

	double min = JPFITS::JPMath::Min(Gdata, false);
	double max = JPFITS::JPMath::Max(Gdata, false);
	double amp = max - min;
	if (amp == 0)
		amp = 1;
	double x0 = xdata[(int)xhw];
	if (x0 == 0)
		x0 = 1;
	double bias = min;
	if (bias == 0)
		bias = 1;

	scale = gcnew array<double>(4) { amp, x0, 2, bias };

	if (p[0] == 0 && p[1] == 0 && p[2] == 0 && p[3] == 0)
	{
		p[0] = amp;
		p[1] = x0;
		p[2] = 2;
		p[3] = bias;
	}

	if (p_lbnd == nullptr)
		p_lbnd = gcnew array<double>(4) { 0, xdata[0], (double)xw / 50, 0 };
	if (p_ubnd == nullptr)
		p_ubnd = gcnew array<double>(4) { 2 * amp, xdata[xw - 1], (double)xw, max };

	alglib::lsfitcreatefg(x, Gdata, p, false, state);
	alglib::lsfitsetcond(state, epsx, maxits);
	alglib::lsfitsetscale(state, scale);
	alglib::lsfitsetbc(state, p_lbnd, p_ubnd);
	alglib::lsfitfit(state, pf, pg, rep, obj);
	alglib::lsfitresults(state, info, p, report);

	if (p_err != nullptr)
		for (int i = 0; i < p_err->Length; i++)
			p_err[i] = report->errpar[i];

	if (fit_residuals != nullptr)
	{
		double val;
		array<double>^ X = gcnew array<double>(1);
		for (int i = 0; i < xw; i++)
		{
			X[0] = xdata[i];
			alglib_Gauss_1d(p, X, val, nullptr);
			fit_residuals[i] = Gdata[i] - val;
		}
	}
}

void JPFITS::JPMath::alglib_Gauss_1d(array<double>^ p, array<double>^ x, double %val, Object^ obj)
{
	val = p[0] * Math::Exp(-((x[0] - p[1])*(x[0] - p[1])) / (2 * p[2] * p[2])) + p[3];
}

void JPFITS::JPMath::alglib_Gauss_1d_grad(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj)
{
	val = p[0] * Math::Exp(-((x[0] - p[1])*(x[0] - p[1])) / (2 * p[2] * p[2])) + p[3];
	grad[0] = Math::Exp(-((x[0] - p[1])*(x[0] - p[1])) / (2 * p[2] * p[2]));
	grad[1] = -(p[0] * Math::Exp(-(p[1] - x[0])*(p[1] - x[0]) / (2 * p[2] * p[2]))*(2 * p[1] - 2 * x[0])) / (2 * p[2] * p[2]);
	grad[2] = (p[0] * Math::Exp(-(p[1] - x[0])*(p[1] - x[0]) / (2 * p[2] * p[2]))*(p[1] - x[0])*(p[1] - x[0])) / (p[2] * p[2] * p[2]);
	grad[3] = 1;
}

void JPFITS::JPMath::Gaussian2d(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ &G, array<double>^ p, bool do_parallel)
{
	int xw = G->GetLength(0);
	int yh = G->GetLength(1);

	if (p->Length == 5)
	{
		#pragma omp parallel for if (do_parallel)
		for (int x = 0; x < xw; x++)
			for (int y = 0; y < yh; y++)
				G[x, y] = p[0] * Math::Exp(-(((double)xdata[x] - p[1])*((double)xdata[x] - p[1]) + ((double)ydata[y] - p[2])*((double)ydata[y] - p[2])) / (2 * p[3] * p[3])) + p[4];
	}

	if (p->Length == 7)
	{
		#pragma omp parallel for if (do_parallel)
		for (int x = 0; x < xw; x++)
			for (int y = 0; y < yh; y++)
				G[x, y] = p[0] * Math::Exp(-Math::Pow(((double)xdata[x] - p[1])*Math::Cos(p[3]) + ((double)ydata[y] - p[2])*Math::Sin(p[3]), 2) / (2 * p[4] * p[4]) - Math::Pow(-((double)xdata[x] - p[1])*Math::Sin(p[3]) + ((double)ydata[y] - p[2])*Math::Cos(p[3]), 2) / (2 * p[5] * p[5])) + p[6];
	}
}

void JPFITS::JPMath::Fit_Gaussian2d(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ Gdata, array<double>^ &p, array<double>^ p_LB, array<double>^ p_UB, array<double>^ &p_err, array<double, 2>^ &fit_residuals)
{
	int N = Gdata->Length;
	array<double, 2>^ x = gcnew array<double, 2>(N, 2);
	array<double>^ y = gcnew array<double>(N);
	int xw = xdata->Length;
	int yh = ydata->Length;

	int i = 0;
	for (int xaxis = 0; xaxis < xw; xaxis++)
		for (int yaxis = 0; yaxis < yh; yaxis++)
		{
			x[i, 0] = (double)xdata[xaxis];
			x[i, 1] = (double)ydata[yaxis];
			y[i] = Gdata[xaxis, yaxis];
			i++;
		}

	alglib::ndimensional_pfunc^ pf = gcnew alglib::ndimensional_pfunc(alglib_Gauss_2d);
	alglib::ndimensional_pgrad ^ pg = gcnew alglib::ndimensional_pgrad(alglib_Gauss_2d_grad);
	alglib::ndimensional_rep^ rep;
	Object^ obj;
	alglib::lsfitstate^ state;
	alglib::lsfitreport^ report;
	double epsx = 0.00001;
	int maxits = 5000;//automatic stopping conditions w epsx & maxits = 0
	int info;
	array<double>^ scale;

	double min = JPFITS::JPMath::Min(Gdata, false);
	double max = JPFITS::JPMath::Max(Gdata, false);
	double amp = max - min;
	if (amp == 0)
		amp = 1;
	double x0 = (double)xdata[xw / 2];
	if (x0 == 0)
		x0 = 1;
	double y0 = (double)ydata[yh / 2];
	if (y0 == 0)
		y0 = 1;
	double bias = min;
	if (bias == 0)
		bias = 1;

	if (p->Length == 5)
	{
		scale = gcnew array<double>(5) { amp, x0, y0, 2, bias };
		if (p_LB == nullptr || p_LB->Length == 0)
			p_LB = gcnew array<double>(5) { 0, (double)xdata[0], (double)ydata[0], 0.01, min - amp };
		if (p_UB == nullptr || p_UB->Length == 0)
			p_UB = gcnew array<double>(5) { 2 * amp, (double)xdata[xw - 1], (double)ydata[yh - 1], (double)xw * 5, max };
		if (p[0] == 0 && p[1] == 0 && p[2] == 0 && p[3] == 0 && p[4] == 0)
		{
			p[0] = amp;
			p[1] = x0;
			p[2] = y0;
			p[3] = 2;
			p[4] = bias;
		}
	}

	if (p->Length == 7)
	{
		scale = gcnew array<double>(7) { amp, x0, y0, 1, 2, 2, bias };
		if (p_LB == nullptr || p_LB->Length == 0)
			p_LB = gcnew array<double>(7) { 0, (double)xdata[0], (double)ydata[0], -Math::PI, 0.01, 0.01, min - amp };
		if (p_UB == nullptr || p_UB->Length == 0)
			p_UB = gcnew array<double>(7) { 2 * amp, (double)xdata[xw - 1], (double)ydata[yh - 1], Math::PI, (double)xw * 5, (double)yh * 5, max };
		if (p[0] == 0 && p[1] == 0 && p[2] == 0 && p[3] == 0 && p[4] == 0 && p[5] == 0 && p[6] == 0)
		{
			p[0] = amp;
			p[1] = x0;
			p[2] = y0;
			p[3] = 0;
			p[4] = 2;
			p[5] = 2;
			p[6] = bias;
		}
	}

	alglib::lsfitcreatefg(x, y, p, false, state);
	alglib::lsfitsetcond(state, epsx, maxits);
	alglib::lsfitsetscale(state, scale);
	alglib::lsfitsetbc(state, p_LB, p_UB);
	alglib::lsfitfit(state, pf, pg, rep, obj);
	alglib::lsfitresults(state, info, p, report);

	if (p_err->Length != 0)
		for (int i = 0; i < p_err->Length; i++)
			p_err[i] = report->errpar[i];

	if (fit_residuals->Length != 0)
	{
		double val;
		array<double>^ X = gcnew array<double>(2);
		for (int i = 0; i < xw; i++)
			for (int j = 0; j < xw; j++)
			{
				X[0] = (double)xdata[i];
				X[1] = (double)ydata[j];
				alglib_Gauss_2d(p, X, val, nullptr);
				fit_residuals[i, j] = Gdata[i, j] - val;
			}
	}
}

void JPFITS::JPMath::alglib_Gauss_2d(array<double>^ p, array<double>^ x, double %val, Object^ obj)
{
	if (p->Length == 5)
		val = p[0] * Math::Exp(-((x[0] - p[1])*(x[0] - p[1]) + (x[1] - p[2])*(x[1] - p[2])) / (2 * p[3] * p[3])) + p[4];
	if (p->Length == 7)
		val = p[0] * Math::Exp(-Math::Pow((x[0] - p[1])*Math::Cos(p[3]) + (x[1] - p[2])*Math::Sin(p[3]), 2) / (2 * p[4] * p[4]) - Math::Pow(-(x[0] - p[1])*Math::Sin(p[3]) + (x[1] - p[2])*Math::Cos(p[3]), 2) / (2 * p[5] * p[5])) + p[6];
}

void JPFITS::JPMath::alglib_Gauss_2d_grad(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj)
{
	if (p->Length == 5)
	{
		val = p[0] * Math::Exp(-((x[0] - p[1])*(x[0] - p[1]) + (x[1] - p[2])*(x[1] - p[2])) / (2 * p[3] * p[3])) + p[4];
		grad[0] = Math::Exp(-((p[1] - x[0])*(p[1] - x[0]) + (p[2] - x[1])*(p[2] - x[1])) / (2 * p[3] * p[3]));
		grad[1] = -(p[0] * Math::Exp(-((p[1] - x[0])*(p[1] - x[0]) + (p[2] - x[1])*(p[2] - x[1])) / (2 * p[3] * p[3]))*(2 * p[1] - 2 * x[0])) / (2 * p[3] * p[3]);
		grad[2] = -(p[0] * Math::Exp(-((p[1] - x[0])*(p[1] - x[0]) + (p[2] - x[1])*(p[2] - x[1])) / (2 * p[3] * p[3]))*(2 * p[2] - 2 * x[1])) / (2 * p[3] * p[3]);
		grad[3] = (p[0] * Math::Exp(-((p[1] - x[0])*(p[1] - x[0]) + (p[2] - x[1])*(p[2] - x[1])) / (2 * p[3] * p[3]))*((p[1] - x[0])*(p[1] - x[0]) + (p[2] - x[1])*(p[2] - x[1]))) / (p[3] * p[3] * p[3]);
		grad[4] = 1;
	}
	if (p->Length == 7)
	{
		val = p[0] * Math::Exp(-Math::Pow((x[0] - p[1])*Math::Cos(p[3]) + (x[1] - p[2])*Math::Sin(p[3]), 2) / (2 * p[4] * p[4]) - Math::Pow(-(x[0] - p[1])*Math::Sin(p[3]) + (x[1] - p[2])*Math::Cos(p[3]), 2) / (2 * p[5] * p[5])) + p[6];
		grad[0] = Math::Exp(-Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (2 * p[4] * p[4]) - Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (2 * p[5] * p[5]));
		grad[1] = -p[0] * Math::Exp(-Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (2 * p[4] * p[4]) - Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (2 * p[5] * p[5]))*((Math::Cos(p[3])*(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]))) / (p[4] * p[4]) - (Math::Sin(p[3])*(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]))) / (p[5] * p[5]));
		grad[2] = -p[0] * Math::Exp(-Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (2 * p[4] * p[4]) - Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (2 * p[5] * p[5]))*((Math::Cos(p[3])*(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]))) / (p[5] * p[5]) + (Math::Sin(p[3])*(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]))) / (p[4] * p[4]));
		grad[3] = -p[0] * Math::Exp(-Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (2 * p[4] * p[4]) - Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (2 * p[5] * p[5]))*(((Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]))*(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]))) / (p[4] * p[4]) - ((Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]))*(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]))) / (p[5] * p[5]));
		grad[4] = (p[0] * Math::Exp(-Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (2 * p[4] * p[4]) - Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (2 * p[5] * p[5]))*Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2)) / (p[4] * p[4] * p[4]);
		grad[5] = (p[0] * Math::Exp(-Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (2 * p[4] * p[4]) - Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (2 * p[5] * p[5]))*Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2)) / (p[5] * p[5] * p[5]);
		grad[6] = 1;
	}
}

void JPFITS::JPMath::Fit_Gaussian2d_Compound(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ Gdata, array<double, 2>^ &p, array<double, 2>^ p_LB, array<double, 2>^ p_UB, array<double, 2>^ &p_err, array<double, 2>^ &fit_residuals)
{
	int N = Gdata->Length;
	array<double, 2>^ x = gcnew array<double, 2>(N, 2);
	array<double>^ y = gcnew array<double>(N);
	int xw = xdata->Length;
	int yh = ydata->Length;
	int func = p->GetLength(0);
	int count = p->GetLength(1);

	int i = 0;
	for (int xaxis = 0; xaxis < xw; xaxis++)
		for (int yaxis = 0; yaxis < yh; yaxis++)
		{
			x[i, 0] = (double)xdata[xaxis];
			x[i, 1] = (double)ydata[yaxis];
			y[i] = Gdata[xaxis, yaxis];
			i++;
		}

	alglib::ndimensional_pfunc^ pf = gcnew alglib::ndimensional_pfunc(alglib_Gauss_2d_compound);
	alglib::ndimensional_pgrad ^ pg = gcnew alglib::ndimensional_pgrad(alglib_Gauss_2d_compound_grad);
	alglib::ndimensional_rep^ rep;
	array<int>^ fcobj = gcnew array<int>(2) { func, count };
	alglib::lsfitstate^ state;
	alglib::lsfitreport^ report;
	double epsx = 0;// 0.00001;
	int maxits = 0;// 5000;//automatic stopping conditions w epsx & maxits = 0
	int info;
	array<double>^ P = gcnew array<double>(p->Length - count + 1);//remove repeated background, but keep one
	array<double>^ PLB = gcnew array<double>(p->Length - count + 1);//remove repeated background, but keep one
	array<double>^ PUB = gcnew array<double>(p->Length - count + 1);//remove repeated background, but keep one
	array<double>^ scale = gcnew array<double>(p->Length - count + 1);//remove repeated background, but keep one

	int c = 0;
	for (int i = 0; i < count; i++)
		for (int j = 0; j < func - 1; j++)
		{
			P[c] = p[j, i];
			PLB[c] = p_LB[j, i];
			PUB[c] = p_UB[j, i];
			scale[c] = P[c];
			if (scale[c] == 0)
				scale[c] = 1;
			c++;
		}
	P[P->Length - 1] = p[func - 1, count - 1];//background
	PLB[P->Length - 1] = p_LB[func - 1, count - 1];//background
	PUB[P->Length - 1] = p_UB[func - 1, count - 1];//background
	scale[P->Length - 1] = P[P->Length - 1];//background
	if (scale[P->Length - 1] == 0)
		scale[P->Length - 1] = 1;

	alglib::lsfitcreatefg(x, y, P, false, state);
	alglib::lsfitsetcond(state, epsx, maxits);
	alglib::lsfitsetscale(state, scale);
	alglib::lsfitsetbc(state, PLB, PUB);
	alglib::lsfitfit(state, pf, pg, rep, (Object^)fcobj);
	alglib::lsfitresults(state, info, P, report);

	for (int i = 0; i < count; i++)
		for (int j = 0; j < func - 1; j++)
			p[j, i] = P[j + i * (func - 1)];
	for (int i = 0; i < count; i++)
		p[func - 1, i] = P[P->Length - 1];//background

	if (p_err->Length != 0)
	{
		for (int i = 0; i < count; i++)
			for (int j = 0; j < func - 1; j++)
				p_err[j, i] = report->errpar[j + i * (func - 1)];
		for (int i = 0; i < count; i++)
			p_err[func - 1, i] = report->errpar[P->Length - 1];//background
	}

	if (fit_residuals->Length != 0)
	{
		double val;
		array<double>^ X = gcnew array<double>(2);
		for (int i = 0; i < xw; i++)
			for (int j = 0; j < yh; j++)
			{
				X[0] = (double)xdata[i];
				X[1] = (double)ydata[j];
				alglib_Gauss_2d_compound(P, X, val, fcobj);
				fit_residuals[i, j] = Gdata[i, j] - val;
			}
	}
}

void JPFITS::JPMath::alglib_Gauss_2d_compound(array<double>^ p, array<double>^ x, double %val, Object^ obj)
{
	int func = ((array<int>^)obj)[0];
	int count = ((array<int>^)obj)[1];

	val = p[p->Length - 1];//background added in once

	if (func == 5)
		for (int i = 0; i < count; i++)
		{
			int index = i * (func - 1);
			val += p[0 + index] * Math::Exp(-((x[0] - p[1 + index])*(x[0] - p[1 + index]) + (x[1] - p[2 + index])*(x[1] - p[2 + index])) / (2 * p[3 + index] * p[3 + index]));
		}
	else
		for (int i = 0; i < count; i++)
		{
			int index = i * (func - 1);
			val += p[0 + index] * Math::Exp(-Math::Pow((x[0] - p[1 + index])*Math::Cos(p[3 + index]) + (x[1] - p[2 + index])*Math::Sin(p[3 + index]), 2) / (2 * p[4 + index] * p[4 + index]) - Math::Pow(-(x[0] - p[1 + index])*Math::Sin(p[3 + index]) + (x[1] - p[2 + index])*Math::Cos(p[3 + index]), 2) / (2 * p[5 + index] * p[5 + index]));
		}
}

void JPFITS::JPMath::alglib_Gauss_2d_compound_grad(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj)
{
	int func = ((array<int>^)obj)[0];
	int count = ((array<int>^)obj)[1];

	val = p[p->Length - 1];//background added in once
	grad[p->Length - 1] = 1;//background added in once

	if (func == 5)
		for (int i = 0; i < count; i++)
		{
			int index = i * (func - 1);
			val += p[0 + index] * Math::Exp(-((x[0] - p[1 + index])*(x[0] - p[1 + index]) + (x[1] - p[2 + index])*(x[1] - p[2 + index])) / (2 * p[3 + index] * p[3 + index]));
			grad[0 + index] = Math::Exp(-((p[1 + index] - x[0])*(p[1 + index] - x[0]) + (p[2 + index] - x[1])*(p[2 + index] - x[1])) / (2 * p[3 + index] * p[3 + index]));
			grad[1 + index] = -(p[0 + index] * Math::Exp(-((p[1 + index] - x[0])*(p[1 + index] - x[0]) + (p[2 + index] - x[1])*(p[2 + index] - x[1])) / (2 * p[3 + index] * p[3 + index]))*(2 * p[1 + index] - 2 * x[0])) / (2 * p[3 + index] * p[3 + index]);
			grad[2 + index] = -(p[0 + index] * Math::Exp(-((p[1 + index] - x[0])*(p[1 + index] - x[0]) + (p[2 + index] - x[1])*(p[2 + index] - x[1])) / (2 * p[3 + index] * p[3 + index]))*(2 * p[2 + index] - 2 * x[1])) / (2 * p[3 + index] * p[3 + index]);
			grad[3 + index] = (p[0 + index] * Math::Exp(-((p[1 + index] - x[0])*(p[1 + index] - x[0]) + (p[2 + index] - x[1])*(p[2 + index] - x[1])) / (2 * p[3 + index] * p[3 + index]))*((p[1 + index] - x[0])*(p[1 + index] - x[0]) + (p[2 + index] - x[1])*(p[2 + index] - x[1]))) / (p[3 + index] * p[3 + index] * p[3 + index]);
		}
	else
		for (int i = 0; i < count; i++)
		{
			int index = i * (func - 1);
			val += p[0 + index] * Math::Exp(-Math::Pow((x[0] - p[1 + index])*Math::Cos(p[3 + index]) + (x[1] - p[2 + index])*Math::Sin(p[3 + index]), 2) / (2 * p[4 + index] * p[4 + index]) - Math::Pow(-(x[0] - p[1 + index])*Math::Sin(p[3 + index]) + (x[1] - p[2 + index])*Math::Cos(p[3 + index]), 2) / (2 * p[5 + index] * p[5 + index]));
			grad[0 + index] = Math::Exp(-Math::Pow(Math::Cos(p[3 + index])*(p[1 + index] - x[0]) + Math::Sin(p[3 + index])*(p[2 + index] - x[1]), 2) / (2 * p[4 + index] * p[4 + index]) - Math::Pow(Math::Cos(p[3 + index])*(p[2 + index] - x[1]) - Math::Sin(p[3 + index])*(p[1 + index] - x[0]), 2) / (2 * p[5 + index] * p[5 + index]));
			grad[1 + index] = -p[0 + index] * Math::Exp(-Math::Pow(Math::Cos(p[3 + index])*(p[1 + index] - x[0]) + Math::Sin(p[3 + index])*(p[2 + index] - x[1]), 2) / (2 * p[4 + index] * p[4 + index]) - Math::Pow(Math::Cos(p[3 + index])*(p[2 + index] - x[1]) - Math::Sin(p[3 + index])*(p[1 + index] - x[0]), 2) / (2 * p[5 + index] * p[5 + index]))*((Math::Cos(p[3 + index])*(Math::Cos(p[3 + index])*(p[1 + index] - x[0]) + Math::Sin(p[3 + index])*(p[2 + index] - x[1]))) / (p[4 + index] * p[4 + index]) - (Math::Sin(p[3 + index])*(Math::Cos(p[3 + index])*(p[2 + index] - x[1]) - Math::Sin(p[3 + index])*(p[1 + index] - x[0]))) / (p[5 + index] * p[5 + index]));
			grad[2 + index] = -p[0 + index] * Math::Exp(-Math::Pow(Math::Cos(p[3 + index])*(p[1 + index] - x[0]) + Math::Sin(p[3 + index])*(p[2 + index] - x[1]), 2) / (2 * p[4 + index] * p[4 + index]) - Math::Pow(Math::Cos(p[3 + index])*(p[2 + index] - x[1]) - Math::Sin(p[3 + index])*(p[1 + index] - x[0]), 2) / (2 * p[5 + index] * p[5 + index]))*((Math::Cos(p[3 + index])*(Math::Cos(p[3 + index])*(p[2 + index] - x[1]) - Math::Sin(p[3 + index])*(p[1 + index] - x[0]))) / (p[5 + index] * p[5 + index]) + (Math::Sin(p[3 + index])*(Math::Cos(p[3 + index])*(p[1 + index] - x[0]) + Math::Sin(p[3 + index])*(p[2 + index] - x[1]))) / (p[4 + index] * p[4 + index]));
			grad[3 + index] = -p[0 + index] * Math::Exp(-Math::Pow(Math::Cos(p[3 + index])*(p[1 + index] - x[0]) + Math::Sin(p[3 + index])*(p[2 + index] - x[1]), 2) / (2 * p[4 + index] * p[4 + index]) - Math::Pow(Math::Cos(p[3 + index])*(p[2 + index] - x[1]) - Math::Sin(p[3 + index])*(p[1 + index] - x[0]), 2) / (2 * p[5 + index] * p[5 + index]))*(((Math::Cos(p[3 + index])*(p[1 + index] - x[0]) + Math::Sin(p[3 + index])*(p[2 + index] - x[1]))*(Math::Cos(p[3 + index])*(p[2 + index] - x[1]) - Math::Sin(p[3 + index])*(p[1 + index] - x[0]))) / (p[4 + index] * p[4 + index]) - ((Math::Cos(p[3 + index])*(p[1 + index] - x[0]) + Math::Sin(p[3 + index])*(p[2 + index] - x[1]))*(Math::Cos(p[3 + index])*(p[2 + index] - x[1]) - Math::Sin(p[3 + index])*(p[1 + index] - x[0]))) / (p[5 + index] * p[5 + index]));
			grad[4 + index] = (p[0 + index] * Math::Exp(-Math::Pow(Math::Cos(p[3 + index])*(p[1 + index] - x[0]) + Math::Sin(p[3 + index])*(p[2 + index] - x[1]), 2) / (2 * p[4 + index] * p[4 + index]) - Math::Pow(Math::Cos(p[3 + index])*(p[2 + index] - x[1]) - Math::Sin(p[3 + index])*(p[1 + index] - x[0]), 2) / (2 * p[5 + index] * p[5 + index]))*Math::Pow(Math::Cos(p[3 + index])*(p[1 + index] - x[0]) + Math::Sin(p[3 + index])*(p[2 + index] - x[1]), 2)) / (p[4 + index] * p[4 + index] * p[4 + index]);
			grad[5 + index] = (p[0 + index] * Math::Exp(-Math::Pow(Math::Cos(p[3 + index])*(p[1 + index] - x[0]) + Math::Sin(p[3 + index])*(p[2 + index] - x[1]), 2) / (2 * p[4 + index] * p[4 + index]) - Math::Pow(Math::Cos(p[3 + index])*(p[2 + index] - x[1]) - Math::Sin(p[3 + index])*(p[1 + index] - x[0]), 2) / (2 * p[5 + index] * p[5 + index]))*Math::Pow(Math::Cos(p[3 + index])*(p[2 + index] - x[1]) - Math::Sin(p[3 + index])*(p[1 + index] - x[0]), 2)) / (p[5 + index] * p[5 + index] * p[5 + index]);
		}
}

array<double>^ JPFITS::JPMath::COG(array<double, 2>^ ROI, int N_last_fit_pts, array<double>^& N_points_COG, double& background_signal_per_pix, double& source_signal)
{
	if (ROI->GetLength(0) != ROI->GetLength(1))
	{
		throw gcnew Exception("Error: ROI array must be square.");
		return nullptr;
	}
	if (ROI->GetLength(0) < 5)
	{
		throw gcnew Exception("Error: Region of interest SubWindow must be at least 5x5 pixels...");
		return nullptr;
	}
	if (JPMath::IsEven(ROI->GetLength(0)) || JPMath::IsEven(ROI->GetLength(1)))
	{
		throw gcnew Exception("Error: ROI array not odd-size.");
		return nullptr;
	}

	int HalfWidth = (ROI->GetLength(0) - 1) / 2;
	N_points_COG = gcnew array<double>(HalfWidth + 1);
	array<double>^ cog = gcnew array<double>(HalfWidth + 1);

	//#pragma omp parallel for - cant do this because of pass by ref to N_points_COG[i] = diskmask; in loop
	for (int i = 0; i <= HalfWidth; i++)
	{
		double data = 0;
		double diskmask = 0;
		for (int ii = -i; ii <= i; ii++)
		{
			for (int jj = -i; jj <= i; jj++)
			{
				if (ii * ii + jj * jj > i * i)
					continue;

				data += ROI[HalfWidth + jj, HalfWidth + ii];
				diskmask++;
			}
		}
		cog[i] = data;
		N_points_COG[i] = diskmask;
	}

	if (N_last_fit_pts > cog->Length)
		N_last_fit_pts = cog->Length;

	array<double>^ xdata = gcnew array<double>(N_last_fit_pts);
	array<double>^ ydata = gcnew array<double>(N_last_fit_pts);
	for (int i = cog->Length - N_last_fit_pts; i < cog->Length; i++)
	{
		xdata[i - (cog->Length - N_last_fit_pts)] = N_points_COG[i];
		ydata[i - (cog->Length - N_last_fit_pts)] = cog[i];
	}
	array<double>^ param = gcnew array<double>(2);
	JPMath::Fit_Poly1d(xdata, ydata, 1, true, param);

	source_signal = param[0];
	background_signal_per_pix = param[1];

	return cog;
}

void JPFITS::JPMath::Radial_Profile_Normalized(array<double, 2>^ Mdata, array<int>^ XSUBRANGE, array<int>^ YSUBRANGE, double axisscale, array<double>^& radial_x, array<double>^& radial_y)
{
	int x, y;
	double center_val = JPMath::Max(Mdata, x, y, false);
	if (x != (Mdata->GetLength(0) - 1) / 2 || y != (Mdata->GetLength(1) - 1) / 2)
	{
		throw gcnew Exception("Error: Mdata array either not odd-size or maximum of array not at the center of the array. Cannot solve radial fit");
		return;
	}
	if (Mdata->GetLength(0) != Mdata->GetLength(1))
	{
		throw gcnew Exception("Error: Mdata array must be square.");
		return;
	}
	if (Mdata->GetLength(0) < 5)
	{
		throw gcnew Exception("Error: Region of interest SubWindow must be at least 5x5 pixels...");
		return;
	}

	array<double, 2>^ SUBIMAGE_radplot = Mdata;

	int SUBIMAGE_HWX = (Mdata->GetLength(0) - 1) / 2;
	int SUBIMAGE_HWY = (Mdata->GetLength(1) - 1) / 2;
	double X0 = (double)XSUBRANGE[SUBIMAGE_HWX];
	double Y0 = (double)YSUBRANGE[SUBIMAGE_HWY];
	int interp_delta = 1;

	array<double>^ XSUBRANGE_radplot;
	array<double>^ YSUBRANGE_radplot;

	if (SUBIMAGE_HWX > 15)
	{
		XSUBRANGE_radplot = gcnew array<double>(XSUBRANGE->Length);
		YSUBRANGE_radplot = gcnew array<double>(YSUBRANGE->Length);
		for (int i = 0; i < XSUBRANGE->Length; i++)
			XSUBRANGE_radplot[i] = double(XSUBRANGE[i]);
		for (int i = 0; i < YSUBRANGE->Length; i++)
			YSUBRANGE_radplot[i] = double(YSUBRANGE[i]);
	}
	else
	{
		interp_delta = 5;
		XSUBRANGE_radplot = gcnew array<double>(XSUBRANGE->Length * interp_delta);
		YSUBRANGE_radplot = gcnew array<double>(YSUBRANGE->Length * interp_delta);
		array<double>^ xx = gcnew array<double>(XSUBRANGE->Length);
		array<double>^ yy = gcnew array<double>(YSUBRANGE->Length);
		for (int i = 0; i < XSUBRANGE->Length; i++)
			xx[i] = double(XSUBRANGE[i]);
		for (int i = 0; i < YSUBRANGE->Length; i++)
			yy[i] = double(YSUBRANGE[i]);

		SUBIMAGE_radplot = JPMath::Interpolate2d(xx, yy, SUBIMAGE_radplot, interp_delta, interp_delta, XSUBRANGE_radplot, YSUBRANGE_radplot, true);
	}

	ArrayList^ distances_LIST = gcnew ArrayList();
	ArrayList^ values_LIST = gcnew ArrayList();
	double SUBIMAGE_HWX_sq = double(SUBIMAGE_HWX * SUBIMAGE_HWX);

	for (int x = 0; x < XSUBRANGE_radplot->Length; x++)
	{
		double dx_sq = (XSUBRANGE_radplot[x] - X0);
		dx_sq *= dx_sq;
		for (int y = 0; y < YSUBRANGE_radplot->Length; y++)
		{
			double dy = YSUBRANGE_radplot[y] - Y0;
			double d_sq = dx_sq + dy * dy;
			if (d_sq > SUBIMAGE_HWX_sq)
				continue;

			distances_LIST->Add(d_sq);
			values_LIST->Add(SUBIMAGE_radplot[x, y]);
		}
	}

	array<double>^ distances_sq = gcnew array<double>(distances_LIST->Count);
	array<double>^ values = gcnew array<double>(distances_LIST->Count);
	for (int q = 0; q < distances_sq->Length; q++)
	{
		distances_sq[q] = double(distances_LIST[q]);
		values[q] = (double)values_LIST[q];
	}
	Array::Sort(distances_sq, values);
	values = JPMath::VectorDivScalar(values, center_val, false);//normalize to max count for radial profile plot

	ArrayList^ r_binnedlist = gcnew ArrayList();
	ArrayList^ v_binnedlist = gcnew ArrayList();
	double d0 = distances_sq[0];
	for (int i = 0; i < distances_sq->Length; i++)
	{
		int dcounter = 0;
		double val = 0;
		while ((i + dcounter) < distances_sq->Length && d0 == distances_sq[i + dcounter])
		{
			val += values[i + dcounter];
			dcounter++;
		}
		r_binnedlist->Add(d0);
		v_binnedlist->Add(val / (double(dcounter)));
		if ((i + dcounter) < distances_sq->Length)
			d0 = distances_sq[i + dcounter];
		i += dcounter - 1;
	}

	array<double>^ r_binned = gcnew array<double>(r_binnedlist->Count);
	array<double>^ v_binned = gcnew array<double>(r_binnedlist->Count);
	for (int q = 0; q < r_binned->Length; q++)
	{
		r_binned[q] = Math::Sqrt((double)r_binnedlist[q]);
		v_binned[q] = (double)v_binnedlist[q];
	}
	if (axisscale != 0 && axisscale != 1 && axisscale > 0)
		for (int q = 0; q < r_binned->Length; q++)
			r_binned[q] *= axisscale;

	radial_x = r_binned;
	radial_y = v_binned;
}

void JPFITS::JPMath::Radial_Profile_Normalized_Fit_Moffat(array<double>^ radial_x, array<double>^ radial_y, array<double>^& p, double& FWHM, array<double>^ &interp_radial_x, array<double>^ & interp_radial_y)
{
	array<double>^ PFit = gcnew array<double>(5) { 1, 0, 1, 1, 0 };
	array<double>^ PFitL = gcnew array<double>(5) { 1, 0, radial_x[radial_x->Length - 1] / 50, radial_x[radial_x->Length - 1] / 50, 0 };
	array<double>^ PFitU = gcnew array<double>(5) { 1, 0, radial_x[radial_x->Length - 1], radial_x[radial_x->Length - 1], 0 };
	array<double>^ resids;
	JPMath::Fit_Moffat1d(radial_x, radial_y, PFit, PFitL, PFitU, nullptr, resids);
	FWHM = 2 * PFit[2] * Math::Sqrt(Math::Pow(2, 1 / PFit[3]) - 1);
	p = PFit;

	interp_radial_x = gcnew array<double>(radial_x->Length * 10);
	interp_radial_y = gcnew array<double>(radial_x->Length * 10);
	double step = (radial_x[radial_x->Length - 1] - radial_x[0]) / interp_radial_x->Length;
	for (int i = 0; i < interp_radial_x->Length; i++)
		interp_radial_x[i] = radial_x[0] + double(i) * step;
	JPMath::Moffat1d(interp_radial_x, interp_radial_y, PFit);
}

void JPFITS::JPMath::Moffat1d(array<double>^ xdata, array<double>^ &M, array<double>^ p)
{
	int xw = M->Length;
	double xhw = double(xw - 1) / 2.0;

	if (xdata == nullptr)
	{
		xdata = gcnew array<double>(xw);
		for (int i = 0; i < xdata->GetLength(0); i++)
			xdata[i] = double(i) - xhw;
	}

	for (int i = 0; i < xw; i++)
		M[i] = p[0] * Math::Pow(1.0 + ((xdata[i] - p[1])*(xdata[i] - p[1])) / p[2] / p[2], -p[3]) + p[4];
}

void JPFITS::JPMath::Fit_Moffat1d(array<double>^ xdata, array<double>^ Mdata, array<double>^ &p, array<double>^ p_lbnd, array<double>^ p_ubnd, array<double>^ p_err, array<double>^ fit_residuals)
{
	int xw = Mdata->Length;
	double xhw = double(xw - 1) / 2.0;

	if (xdata == nullptr)
	{
		xdata = gcnew array<double>(xw);
		for (int i = 0; i < xdata->GetLength(0); i++)
			xdata[i] = double(i) - xhw;
	}

	array<double, 2>^ x = gcnew array<double, 2>(xw, 1);
	for (int i = 0; i < xw; i++)
		x[i, 0] = xdata[i];

	alglib::ndimensional_pfunc^ pf = gcnew alglib::ndimensional_pfunc(alglib_Moffat_1d);
	alglib::ndimensional_pgrad ^ pg = gcnew alglib::ndimensional_pgrad(alglib_Moffat_1d_grad);
	alglib::ndimensional_rep^ rep;
	Object^ obj;
	alglib::lsfitstate^ state;
	alglib::lsfitreport^ report;
	double diffstep = 0.0001;
	double epsx = 0.000001;
	int maxits = 0;
	int info;
	array<double>^ scale;

	double min = JPFITS::JPMath::Min(Mdata, false);
	double max = JPFITS::JPMath::Max(Mdata, false);
	double amp = max - min;
	if (amp == 0)
		amp = 1;
	double x0 = xdata[(int)xhw];
	if (x0 == 0)
		x0 = 1;
	double bias = min;
	if (bias == 0)
		bias = 1;

	scale = gcnew array<double>(5) { amp, x0, 2, 2, bias };

	if (p[0] == 0 && p[1] == 0 && p[2] == 0 && p[3] == 0 && p[4] == 0)
	{
		p[0] = amp;
		p[1] = x0;
		p[2] = 2;
		p[3] = 2;
		p[4] = bias;
	}

	if (p_lbnd == nullptr)
		p_lbnd = gcnew array<double>(5) { 0, xdata[0], (double)xw / 500, (double)xw / 500, 0 };
	if (p_ubnd == nullptr)
		p_ubnd = gcnew array<double>(5) { 2 * amp, xdata[xw - 1], (double)xw * 5, (double)xw * 5, max };

	alglib::lsfitcreatefg(x, Mdata, p, false, state);
	alglib::lsfitsetcond(state, epsx, maxits);
	alglib::lsfitsetscale(state, scale);
	alglib::lsfitsetbc(state, p_lbnd, p_ubnd);
	alglib::lsfitfit(state, pf, pg, rep, obj);
	alglib::lsfitresults(state, info, p, report);

	if (p_err != nullptr)
		for (int i = 0; i < p_err->Length; i++)
			p_err[i] = report->errpar[i];

	if (fit_residuals != nullptr)
	{
		double val;
		array<double>^ X = gcnew array<double>(1);
		for (int i = 0; i < xw; i++)
		{
			X[0] = xdata[i];
			alglib_Moffat_1d(p, X, val, nullptr);
			fit_residuals[i] = Mdata[i] - val;
		}
	}
}

void JPFITS::JPMath::alglib_Moffat_1d(array<double>^ p, array<double>^ x, double %val, Object^ obj)
{
	val = p[0] * Math::Pow(1.0 + ((x[0] - p[1])*(x[0] - p[1])) / p[2] / p[2], -p[3]) + p[4];
}

void JPFITS::JPMath::alglib_Moffat_1d_grad(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj)
{
	val = p[0] * Math::Pow(1.0 + ((x[0] - p[1])*(x[0] - p[1])) / p[2] / p[2], -p[3]) + p[4];
	grad[0] = Math::Pow((p[1] - x[0])*(p[1] - x[0]) / (p[2] * p[2]) + 1, -p[3]);
	grad[1] = -(p[0] * p[3] * (2 * p[1] - 2 * x[0])) / (p[2] * p[2] * Math::Pow((p[1] - x[0])*(p[1] - x[0]) / (p[2] * p[2]) + 1, (p[3] + 1)));
	grad[2] = (2 * p[0] * p[3] * (p[1] - x[0])*(p[1] - x[0])) / (p[2] * p[2] * p[2] * Math::Pow((p[1] - x[0])*(p[1] - x[0]) / (p[2] * p[2]) + 1, (p[3] + 1)));
	grad[3] = -(p[0] * Math::Log((p[1] - x[0])*(p[1] - x[0]) / (p[2] * p[2]) + 1)) / Math::Pow((p[1] - x[0])*(p[1] - x[0]) / (p[2] * p[2]) + 1, p[3]);
	grad[4] = 1;
}

void JPFITS::JPMath::Moffat2d(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ &M, array<double>^ p, bool do_parallel)
{
	int xw = M->GetLength(0);
	int yh = M->GetLength(1);

	if (p->Length == 6)
	{
		#pragma omp parallel for if (do_parallel)
		for (int i = 0; i < xw; i++)
			for (int j = 0; j < yh; j++)
				M[i, j] = p[0] * Math::Pow(1.0 + (((double)xdata[i] - p[1])*((double)xdata[i] - p[1]) + ((double)ydata[j] - p[2])*((double)ydata[j] - p[2])) / p[3] / p[3], -p[4]) + p[5];
	}

	if (p->Length == 8)
	{
		#pragma omp parallel for if (do_parallel)
		for (int i = 0; i < xw; i++)
			for (int j = 0; j < yh; j++)
				M[i, j] = p[0] * Math::Pow(1.0 + (Math::Pow(((double)xdata[i] - p[1])*Math::Cos(p[3]) + ((double)ydata[j] - p[2])*Math::Sin(p[3]), 2)) / (p[4] * p[4]) + (Math::Pow(-((double)xdata[i] - p[1])*Math::Sin(p[3]) + ((double)ydata[j] - p[2])*Math::Cos(p[3]), 2)) / (p[5] * p[5]), -p[6]) + p[7];
	}
}

void JPFITS::JPMath::Fit_Moffat2d(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ Mdata, array<double>^ &p, array<double>^ p_LB, array<double>^ p_UB, array<double>^ &p_err, array<double, 2>^ &fit_residuals)
{
	int N = Mdata->Length;
	array<double, 2>^ x = gcnew array<double, 2>(N, 2);
	array<double>^ y = gcnew array<double>(N);
	int xw = Mdata->GetLength(0);
	int yh = Mdata->GetLength(1);

	int i = 0;
	for (int xaxis = 0; xaxis < xw; xaxis++)
		for (int yaxis = 0; yaxis < yh; yaxis++)
		{
			x[i, 0] = (double)xdata[xaxis];
			x[i, 1] = (double)ydata[yaxis];
			y[i] = Mdata[xaxis, yaxis];
			i++;
		}

	alglib::ndimensional_pfunc^ pf = gcnew alglib::ndimensional_pfunc(alglib_Moffat_2d);
	alglib::ndimensional_pgrad ^ pg = gcnew alglib::ndimensional_pgrad(alglib_Moffat_2d_grad);
	alglib::ndimensional_rep^ rep;
	Object^ obj;
	alglib::lsfitstate^ state;
	alglib::lsfitreport^ report;
	double epsx = 0.00001;
	int maxits = 5000;//automatic stopping conditions w epsx & maxits = 0
	int info;
	array<double>^ scale;

	double min = JPFITS::JPMath::Min(Mdata, false);
	double max = JPFITS::JPMath::Max(Mdata, false);
	double amp = max - min;
	if (amp == 0)
		amp = 1;
	double x0 = (double)xdata[xw / 2];
	if (x0 == 0)
		x0 = 1;
	double y0 = (double)ydata[yh / 2];
	if (y0 == 0)
		y0 = 1;
	double bias = min;
	if (bias == 0)
		bias = 1;

	if (p->Length == 6)
	{
		scale = gcnew array<double>(6) { amp, x0, y0, 2, 3, bias };
		if (p_LB == nullptr || p_LB->Length == 0)
			p_LB = gcnew array<double>(6) { 0, (double)xdata[0], (double)ydata[0], (double)xw / 500, 1.01, min - amp };
		if (p_UB == nullptr || p_UB->Length == 0)
			p_UB = gcnew array<double>(6) { 2 * amp, (double)xdata[xw - 1], (double)ydata[yh - 1], (double)xw * 5, (double)xw * 5, max };
		if (p[0] == 0 && p[1] == 0 && p[2] == 0 && p[3] == 0 && p[4] == 0 && p[5] == 0)
		{
			p[0] = amp;
			p[1] = x0;
			p[2] = y0;
			p[3] = 2;
			p[4] = 4;
			p[5] = bias;
		}
	}
	if (p->Length == 8)
	{
		scale = gcnew array<double>(8) { amp, x0, y0, 1, 2, 2, 3, bias };
		if (p_LB == nullptr || p_LB->Length == 0)
			p_LB = gcnew array<double>(8) { 0, (double)xdata[0], (double)ydata[0], -Math::PI, (double)xw / 500, (double)yh / 500, 1.01, min - amp };
		if (p_UB == nullptr || p_UB->Length == 0)
			p_UB = gcnew array<double>(8) { 2 * amp, (double)xdata[xw - 1], (double)ydata[yh - 1], Math::PI, (double)xw * 5, (double)xw * 5, (double)xw * 5, max };
		if (p[0] == 0 && p[1] == 0 && p[2] == 0 && p[3] == 0 && p[4] == 0 && p[5] == 0 && p[6] == 0 && p[7] == 0)
		{
			p[0] = amp;
			p[1] = x0;
			p[2] = y0;
			p[3] = 0;
			p[4] = 2;
			p[5] = 2;
			p[6] = 4;
			p[7] = bias;
		}
	}

	alglib::lsfitcreatefg(x, y, p, false, state);
	alglib::lsfitsetcond(state, epsx, maxits);
	alglib::lsfitsetscale(state, scale);
	alglib::lsfitsetbc(state, p_LB, p_UB);
	alglib::lsfitfit(state, pf, pg, rep, obj);
	alglib::lsfitresults(state, info, p, report);

	if (p_err->Length != 0)
		for (int i = 0; i < p_err->Length; i++)
			p_err[i] = report->errpar[i];

	if (fit_residuals->Length != 0)
	{
		double val;
		array<double>^ X = gcnew array<double>(2);
		for (int i = 0; i < xw; i++)
			for (int j = 0; j < xw; j++)
			{
				X[0] = (double)xdata[i];
				X[1] = (double)ydata[j];
				alglib_Moffat_2d(p, X, val, nullptr);
				fit_residuals[i, j] = Mdata[i, j] - val;
			}
	}
}

void JPFITS::JPMath::alglib_Moffat_2d(array<double>^ p, array<double>^ x, double %val, Object^ obj)
{
	if (p->Length == 6)
		val = p[0] * Math::Pow(1.0 + ((x[0] - p[1])*(x[0] - p[1]) + (x[1] - p[2])*(x[1] - p[2])) / p[3] / p[3], -p[4]) + p[5];
	if (p->Length == 8)
		val = p[0] * Math::Pow(1.0 + (Math::Pow((x[0] - p[1])*Math::Cos(p[3]) + (x[1] - p[2])*Math::Sin(p[3]), 2)) / (p[4] * p[4]) + (Math::Pow(-(x[0] - p[1])*Math::Sin(p[3]) + (x[1] - p[2])*Math::Cos(p[3]), 2)) / (p[5] * p[5]), -p[6]) + p[7];
}

void JPFITS::JPMath::alglib_Moffat_2d_grad(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj)
{
	if (p->Length == 6)
	{
		val = p[0] * Math::Pow(1.0 + ((x[0] - p[1])*(x[0] - p[1]) + (x[1] - p[2])*(x[1] - p[2])) / p[3] / p[3], -p[4]) + p[5];
		grad[0] = Math::Pow(((p[1] - x[0])*(p[1] - x[0]) + (p[2] - x[1])*(p[2] - x[1])) / (p[3] * p[3]) + 1, -p[4]);
		grad[1] = -(p[0] * p[4] * (2 * p[1] - 2 * x[0])) / (p[3] * p[3] * Math::Pow(((p[1] - x[0])*(p[1] - x[0]) + (p[2] - x[1])*(p[2] - x[1])) / (p[3] * p[3]) + 1, (p[4] + 1)));
		grad[2] = -(p[0] * p[4] * (2 * p[2] - 2 * x[1])) / (p[3] * p[3] * Math::Pow(((p[1] - x[0])*(p[1] - x[0]) + (p[2] - x[1])*(p[2] - x[1])) / (p[3] * p[3]) + 1, (p[4] + 1)));
		grad[3] = (2 * p[0] * p[4] * ((p[1] - x[0])*(p[1] - x[0]) + (p[2] - x[1])*(p[2] - x[1]))) / (p[3] * p[3] * p[3] * Math::Pow(((p[1] - x[0])*(p[1] - x[0]) + (p[2] - x[1])*(p[2] - x[1])) / (p[3] * p[3]) + 1, (p[4] + 1)));
		grad[4] = -(p[0] * Math::Log(((p[1] - x[0])*(p[1] - x[0]) + (p[2] - x[1])*(p[2] - x[1])) / (p[3] * p[3]) + 1)) / Math::Pow(((p[1] - x[0])*(p[1] - x[0]) + (p[2] - x[1])*(p[2] - x[1])) / (p[3] * p[3]) + 1, p[4]);
		grad[5] = 1;
	}
	if (p->Length == 8)
	{
		val = p[0] * Math::Pow(1.0 + (Math::Pow((x[0] - p[1])*Math::Cos(p[3]) + (x[1] - p[2])*Math::Sin(p[3]), 2)) / (p[4] * p[4]) + (Math::Pow(-(x[0] - p[1])*Math::Sin(p[3]) + (x[1] - p[2])*Math::Cos(p[3]), 2)) / (p[5] * p[5]), -p[6]) + p[7];
		grad[0] = Math::Pow(Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (p[4] * p[4]) + Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (p[5] * p[5]) + 1, -p[6]);
		grad[1] = -(p[0] * p[6] * ((2 * Math::Cos(p[3])*(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]))) / (p[4] * p[4]) - (2 * Math::Sin(p[3])*(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]))) / (p[5] * p[5]))) / Math::Pow(Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (p[4] * p[4]) + Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (p[5] * p[5]) + 1, (p[6] + 1));
		grad[2] = -(p[0] * p[6] * ((2 * Math::Cos(p[3])*(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]))) / (p[5] * p[5]) + (2 * Math::Sin(p[3])*(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]))) / (p[4] * p[4]))) / Math::Pow(Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (p[4] * p[4]) + Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (p[5] * p[5]) + 1, (p[6] + 1));
		grad[3] = -(p[0] * p[6] * ((2 * (Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]))*(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]))) / (p[4] * p[4]) - (2 * (Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]))*(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]))) / (p[5] * p[5]))) / Math::Pow(Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (p[4] * p[4]) + Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (p[5] * p[5]) + 1, (p[6] + 1));
		grad[4] = (2 * p[0] * p[6] * Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2)) / (p[4] * p[4] * p[4] * Math::Pow(Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (p[4] * p[4]) + Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (p[5] * p[5]) + 1, (p[6] + 1)));
		grad[5] = (2 * p[0] * p[6] * Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2)) / (p[5] * p[5] * p[5] * Math::Pow(Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (p[4] * p[4]) + Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (p[5] * p[5]) + 1, (p[6] + 1)));
		grad[6] = -(p[0] * Math::Log(Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (p[4] * p[4]) + Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (p[5] * p[5]) + 1)) / Math::Pow(Math::Pow(Math::Cos(p[3])*(p[1] - x[0]) + Math::Sin(p[3])*(p[2] - x[1]), 2) / (p[4] * p[4]) + Math::Pow(Math::Cos(p[3])*(p[2] - x[1]) - Math::Sin(p[3])*(p[1] - x[0]), 2) / (p[5] * p[5]) + 1, p[6]);
		grad[7] = 1;
	}
}

void JPFITS::JPMath::Fit_Moffat2d_Compound(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ Mdata, array<double, 2>^ &p, array<double, 2>^ p_LB, array<double, 2>^ p_UB, array<double, 2>^ &p_err, array<double, 2>^ &fit_residuals)
{
	int N = Mdata->Length;
	array<double, 2>^ x = gcnew array<double, 2>(N, 2);
	array<double>^ y = gcnew array<double>(N);
	int xw = xdata->Length;
	int yh = ydata->Length;
	int func = p->GetLength(0);
	int count = p->GetLength(1);

	int i = 0;
	for (int xaxis = 0; xaxis < xw; xaxis++)
		for (int yaxis = 0; yaxis < yh; yaxis++)
		{
			x[i, 0] = (double)xdata[xaxis];
			x[i, 1] = (double)ydata[yaxis];
			y[i] = Mdata[xaxis, yaxis];
			i++;
		}

	alglib::ndimensional_pfunc^ pf = gcnew alglib::ndimensional_pfunc(alglib_Moffat_2d_compound);
	alglib::ndimensional_pgrad ^ pg = gcnew alglib::ndimensional_pgrad(alglib_Moffat_2d_compound_grad);
	alglib::ndimensional_rep^ rep;
	array<int>^ fcobj = gcnew array<int>(2) { func, count };
	alglib::lsfitstate^ state;
	alglib::lsfitreport^ report;
	double epsx = 0;// 0.00001;
	int maxits = 0;// 5000;//automatic stopping conditions w epsx & maxits = 0
	int info;
	array<double>^ P = gcnew array<double>(p->Length - count + 1);//remove repeated background, but keep one
	array<double>^ PLB = gcnew array<double>(p->Length - count + 1);//remove repeated background, but keep one
	array<double>^ PUB = gcnew array<double>(p->Length - count + 1);//remove repeated background, but keep one
	array<double>^ scale = gcnew array<double>(p->Length - count + 1);//remove repeated background, but keep one

	int c = 0;
	for (int i = 0; i < count; i++)
		for (int j = 0; j < func - 1; j++)
		{
			P[c] = p[j, i];
			PLB[c] = p_LB[j, i];
			PUB[c] = p_UB[j, i];
			scale[c] = P[c];
			if (scale[c] <= 0)
				scale[c] = 1;
			c++;
		}
	P[P->Length - 1] = p[func - 1, count - 1];//background
	PLB[P->Length - 1] = p_LB[func - 1, count - 1];//background
	PUB[P->Length - 1] = p_UB[func - 1, count - 1];//background
	scale[P->Length - 1] = P[P->Length - 1];//background
	if (scale[P->Length - 1] <= 0)
		scale[P->Length - 1] = 1;

	alglib::lsfitcreatefg(x, y, P, false, state);
	alglib::lsfitsetcond(state, epsx, maxits);
	alglib::lsfitsetscale(state, scale);
	alglib::lsfitsetbc(state, PLB, PUB);
	alglib::lsfitfit(state, pf, pg, rep, (Object^)fcobj);
	alglib::lsfitresults(state, info, P, report);

	for (int i = 0; i < count; i++)
		for (int j = 0; j < func - 1; j++)
			p[j, i] = P[j + i * (func - 1)];
	for (int i = 0; i < count; i++)
		p[func - 1, i] = P[P->Length - 1];//background

	if (p_err->Length != 0)
	{
		for (int i = 0; i < count; i++)
			for (int j = 0; j < func - 1; j++)
				p_err[j, i] = report->errpar[j + i * (func - 1)];
		for (int i = 0; i < count; i++)
			p_err[func - 1, i] = report->errpar[P->Length - 1];//background
	}

	if (fit_residuals->Length != 0)
	{
		double val;
		array<double>^ X = gcnew array<double>(2);
		for (int i = 0; i < xw; i++)
			for (int j = 0; j < yh; j++)
			{
				X[0] = (double)xdata[i];
				X[1] = (double)ydata[j];
				alglib_Moffat_2d_compound(P, X, val, fcobj);
				fit_residuals[i, j] = Mdata[i, j] - val;
			}
	}
}

void JPFITS::JPMath::alglib_Moffat_2d_compound(array<double>^ p, array<double>^ x, double %val, Object^ obj)
{
	int func = ((array<int>^)obj)[0];
	int count = ((array<int>^)obj)[1];

	val = p[p->Length - 1];//background added in once

	if (func == 6)
		for (int i = 0; i < count; i++)
		{
			int index = i * (func - 1);
			val += p[0 + index] * Math::Pow(1.0 + ((x[0] - p[1 + index])*(x[0] - p[1 + index]) + (x[1] - p[2 + index])*(x[1] - p[2 + index])) / p[3 + index] / p[3 + index], -p[4 + index]);
		}
	if (func == 8)
		for (int i = 0; i < count; i++)
		{
			int index = i * (func - 1);
			val += p[0 + index] * Math::Pow(1.0 + (Math::Pow((x[0] - p[1 + index])*Math::Cos(p[3 + index]) + (x[1] - p[2 + index])*Math::Sin(p[3 + index]), 2)) / (p[4 + index] * p[4 + index]) + (Math::Pow(-(x[0] - p[1 + index])*Math::Sin(p[3 + index]) + (x[1] - p[2 + index])*Math::Cos(p[3 + index]), 2)) / (p[5 + index] * p[5 + index]), -p[6 + index]);
		}
}

void JPFITS::JPMath::alglib_Moffat_2d_compound_grad(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj)
{
	int func = ((array<int>^)obj)[0];
	int count = ((array<int>^)obj)[1];

	val = p[p->Length - 1];//background added in once
	grad[p->Length - 1] = 1;//background added in once

	if (func == 6)
	{
		for (int i = 0; i < count; i++)
		{
			int index = i * (func - 1);
			val += p[index + 0] * Math::Pow(1.0 + ((x[0] - p[index + 1])*(x[0] - p[index + 1]) + (x[1] - p[index + 2])*(x[1] - p[index + 2])) / p[index + 3] / p[index + 3], -p[index + 4]);
			grad[0 + index] = Math::Pow(((p[index + 1] - x[0])*(p[index + 1] - x[0]) + (p[index + 2] - x[1])*(p[index + 2] - x[1])) / (p[index + 3] * p[index + 3]) + 1, -p[index + 4]);
			grad[1 + index] = -(p[index + 0] * p[index + 4] * (2 * p[index + 1] - 2 * x[0])) / (p[index + 3] * p[index + 3] * Math::Pow(((p[index + 1] - x[0])*(p[index + 1] - x[0]) + (p[index + 2] - x[1])*(p[index + 2] - x[1])) / (p[index + 3] * p[index + 3]) + 1, (p[index + 4] + 1)));
			grad[2 + index] = -(p[index + 0] * p[index + 4] * (2 * p[index + 2] - 2 * x[1])) / (p[index + 3] * p[index + 3] * Math::Pow(((p[index + 1] - x[0])*(p[index + 1] - x[0]) + (p[index + 2] - x[1])*(p[index + 2] - x[1])) / (p[index + 3] * p[index + 3]) + 1, (p[index + 4] + 1)));
			grad[3 + index] = (2 * p[index + 0] * p[index + 4] * ((p[index + 1] - x[0])*(p[index + 1] - x[0]) + (p[index + 2] - x[1])*(p[index + 2] - x[1]))) / (p[index + 3] * p[index + 3] * p[index + 3] * Math::Pow(((p[index + 1] - x[0])*(p[index + 1] - x[0]) + (p[index + 2] - x[1])*(p[index + 2] - x[1])) / (p[index + 3] * p[index + 3]) + 1, (p[index + 4] + 1)));
			grad[4 + index] = -(p[index + 0] * Math::Log(((p[index + 1] - x[0])*(p[index + 1] - x[0]) + (p[index + 2] - x[1])*(p[index + 2] - x[1])) / (p[index + 3] * p[index + 3]) + 1)) / Math::Pow(((p[index + 1] - x[0])*(p[index + 1] - x[0]) + (p[index + 2] - x[1])*(p[index + 2] - x[1])) / (p[index + 3] * p[index + 3]) + 1, p[index + 4]);
		}
	}
	if (func == 8)
	{
		for (int i = 0; i < count; i++)
		{
			int index = i * (func - 1);
			val += p[index + 0] * Math::Pow(1.0 + (Math::Pow((x[0] - p[index + 1])*Math::Cos(p[index + 3]) + (x[1] - p[index + 2])*Math::Sin(p[index + 3]), 2)) / (p[index + 4] * p[index + 4]) + (Math::Pow(-(x[0] - p[index + 1])*Math::Sin(p[index + 3]) + (x[1] - p[index + 2])*Math::Cos(p[index + 3]), 2)) / (p[index + 5] * p[index + 5]), -p[index + 6]);
			grad[0 + index] = Math::Pow(Math::Pow(Math::Cos(p[index + 3])*(p[index + 1] - x[0]) + Math::Sin(p[index + 3])*(p[index + 2] - x[1]), 2) / (p[index + 4] * p[index + 4]) + Math::Pow(Math::Cos(p[index + 3])*(p[index + 2] - x[1]) - Math::Sin(p[index + 3])*(p[index + 1] - x[0]), 2) / (p[index + 5] * p[index + 5]) + 1, -p[index + 6]);
			grad[1 + index] = -(p[index + 0] * p[index + 6] * ((2 * Math::Cos(p[index + 3])*(Math::Cos(p[index + 3])*(p[index + 1] - x[0]) + Math::Sin(p[index + 3])*(p[index + 2] - x[1]))) / (p[index + 4] * p[index + 4]) - (2 * Math::Sin(p[index + 3])*(Math::Cos(p[index + 3])*(p[index + 2] - x[1]) - Math::Sin(p[index + 3])*(p[index + 1] - x[0]))) / (p[index + 5] * p[index + 5]))) / Math::Pow(Math::Pow(Math::Cos(p[index + 3])*(p[index + 1] - x[0]) + Math::Sin(p[index + 3])*(p[index + 2] - x[1]), 2) / (p[index + 4] * p[index + 4]) + Math::Pow(Math::Cos(p[index + 3])*(p[index + 2] - x[1]) - Math::Sin(p[index + 3])*(p[index + 1] - x[0]), 2) / (p[index + 5] * p[index + 5]) + 1, (p[index + 6] + 1));
			grad[2 + index] = -(p[index + 0] * p[index + 6] * ((2 * Math::Cos(p[index + 3])*(Math::Cos(p[index + 3])*(p[index + 2] - x[1]) - Math::Sin(p[index + 3])*(p[index + 1] - x[0]))) / (p[index + 5] * p[index + 5]) + (2 * Math::Sin(p[index + 3])*(Math::Cos(p[index + 3])*(p[index + 1] - x[0]) + Math::Sin(p[index + 3])*(p[index + 2] - x[1]))) / (p[index + 4] * p[index + 4]))) / Math::Pow(Math::Pow(Math::Cos(p[index + 3])*(p[index + 1] - x[0]) + Math::Sin(p[index + 3])*(p[index + 2] - x[1]), 2) / (p[index + 4] * p[index + 4]) + Math::Pow(Math::Cos(p[index + 3])*(p[index + 2] - x[1]) - Math::Sin(p[index + 3])*(p[index + 1] - x[0]), 2) / (p[index + 5] * p[index + 5]) + 1, (p[index + 6] + 1));
			grad[3 + index] = -(p[index + 0] * p[index + 6] * ((2 * (Math::Cos(p[index + 3])*(p[index + 1] - x[0]) + Math::Sin(p[index + 3])*(p[index + 2] - x[1]))*(Math::Cos(p[index + 3])*(p[index + 2] - x[1]) - Math::Sin(p[index + 3])*(p[index + 1] - x[0]))) / (p[index + 4] * p[index + 4]) - (2 * (Math::Cos(p[index + 3])*(p[index + 1] - x[0]) + Math::Sin(p[index + 3])*(p[index + 2] - x[1]))*(Math::Cos(p[index + 3])*(p[index + 2] - x[1]) - Math::Sin(p[index + 3])*(p[index + 1] - x[0]))) / (p[index + 5] * p[index + 5]))) / Math::Pow(Math::Pow(Math::Cos(p[index + 3])*(p[index + 1] - x[0]) + Math::Sin(p[index + 3])*(p[index + 2] - x[1]), 2) / (p[index + 4] * p[index + 4]) + Math::Pow(Math::Cos(p[index + 3])*(p[index + 2] - x[1]) - Math::Sin(p[index + 3])*(p[index + 1] - x[0]), 2) / (p[index + 5] * p[index + 5]) + 1, (p[index + 6] + 1));
			grad[4 + index] = (2 * p[index + 0] * p[index + 6] * Math::Pow(Math::Cos(p[index + 3])*(p[index + 1] - x[0]) + Math::Sin(p[index + 3])*(p[index + 2] - x[1]), 2)) / (p[index + 4] * p[index + 4] * p[index + 4] * Math::Pow(Math::Pow(Math::Cos(p[index + 3])*(p[index + 1] - x[0]) + Math::Sin(p[index + 3])*(p[index + 2] - x[1]), 2) / (p[index + 4] * p[index + 4]) + Math::Pow(Math::Cos(p[index + 3])*(p[index + 2] - x[1]) - Math::Sin(p[index + 3])*(p[index + 1] - x[0]), 2) / (p[index + 5] * p[index + 5]) + 1, (p[index + 6] + 1)));
			grad[5 + index] = (2 * p[index + 0] * p[index + 6] * Math::Pow(Math::Cos(p[index + 3])*(p[index + 2] - x[1]) - Math::Sin(p[index + 3])*(p[index + 1] - x[0]), 2)) / (p[index + 5] * p[index + 5] * p[index + 5] * Math::Pow(Math::Pow(Math::Cos(p[index + 3])*(p[index + 1] - x[0]) + Math::Sin(p[index + 3])*(p[index + 2] - x[1]), 2) / (p[index + 4] * p[index + 4]) + Math::Pow(Math::Cos(p[index + 3])*(p[index + 2] - x[1]) - Math::Sin(p[index + 3])*(p[index + 1] - x[0]), 2) / (p[index + 5] * p[index + 5]) + 1, (p[index + 6] + 1)));
			grad[6 + index] = -(p[index + 0] * Math::Log(Math::Pow(Math::Cos(p[index + 3])*(p[index + 1] - x[0]) + Math::Sin(p[index + 3])*(p[index + 2] - x[1]), 2) / (p[index + 4] * p[index + 4]) + Math::Pow(Math::Cos(p[index + 3])*(p[index + 2] - x[1]) - Math::Sin(p[index + 3])*(p[index + 1] - x[0]), 2) / (p[index + 5] * p[index + 5]) + 1)) / Math::Pow(Math::Pow(Math::Cos(p[index + 3])*(p[index + 1] - x[0]) + Math::Sin(p[index + 3])*(p[index + 2] - x[1]), 2) / (p[index + 4] * p[index + 4]) + Math::Pow(Math::Cos(p[index + 3])*(p[index + 2] - x[1]) - Math::Sin(p[index + 3])*(p[index + 1] - x[0]), 2) / (p[index + 5] * p[index + 5]) + 1, p[index + 6]);
		}
	}
}

/*void JPFITS::JPMath::Transform2d(array<double>^ x0, array<double>^ y0, array<double>^ p, array<double>^ &xp, array<double>^ &yp)
{
	if (p->Length == 4)
		for (int i = 0; i < x0->Length; i++)
		{
			xp[i] = p[0] * (Math::Cos(p[1]) * (x0[i] - p[2]) - Math::Sin(p[1]) * (y0[i] - p[3]));
			yp[i] = p[0] * (Math::Sin(p[1]) * (x0[i] - p[2]) + Math::Cos(p[1]) * (y0[i] - p[3]));

		}

	if (p->Length == 6)
		for (int i = 0; i < x0->Length; i++)
		{
			xp[i] = p[0] * (x0[i] - p[4]) + p[1] * (y0[i] - p[5]);
			yp[i] = p[2] * (x0[i] - p[4]) + p[3] * (y0[i] - p[5]);
		}
}*/

void JPFITS::JPMath::Fit_WCSTransform2d(array<double>^ x_intrmdt, array<double>^ y_intrmdt, array<double>^ x_pix, array<double>^ y_pix, array<double>^ &p, array<double>^ p_lbnd, array<double>^ p_ubnd, array<double>^ p_scale)
{
	double epsx = 0.00000000001;
	int maxits = 0;
	alglib::minlmstate^ state;
	alglib::minlmreport^ report;
	alglib::ndimensional_fvec^ minfunc = gcnew alglib::ndimensional_fvec(alglib_WCSTransform2d);
	//alglib::ndimensional_jac^ jac = gcnew alglib::ndimensional_jac(alglib_Transform2d_jac);
	array<Object^>^ objj = gcnew array<Object^>(4);
	objj[0] = (Object^)x_pix;
	objj[1] = (Object^)y_pix;
	objj[2] = (Object^)x_intrmdt;
	objj[3] = (Object^)y_intrmdt;
	Object^ obj = (Object^)objj;

	//alglib::minlmcreatevj(1, p, state);
	alglib::minlmcreatev(1, p, 0.00001, state);
	alglib::minlmsetcond(state, epsx, maxits);
	alglib::minlmsetscale(state, p_scale);
	alglib::minlmsetbc(state, p_lbnd, p_ubnd);
	//alglib::minlmoptimize(state, minfunc, jac, nullptr, obj);
	alglib::minlmoptimize(state, minfunc, nullptr, obj);
	alglib::minlmresults(state, p, report);

	/*if (p_err != nullptr)
		for (int i = 0; i < p_err->Length; i++)
			p_err[i] = report->errerrpar[i];*/
}

void JPFITS::JPMath::alglib_WCSTransform2d(array<double>^ p, array<double>^ f, Object^ obj)
{
	array<Object^>^ objj = (array<Object^>^)obj;
	array<double>^ x_pix = (array<double>^)objj[0];
	array<double>^ y_pix = (array<double>^)objj[1];
	array<double>^ x_intrmdt = (array<double>^)objj[2];
	array<double>^ y_intrmdt = (array<double>^)objj[3];
	f[0] = 0;
	double xprime, yprime;

	if (p->Length == 4)
	{
		for (int i = 0; i < x_intrmdt->Length; i++)
		{
			xprime = p[0] * (Math::Cos(p[1]) * (x_pix[i] - p[2]) - Math::Sin(p[1]) * (y_pix[i] - p[3]));
			yprime = p[0] * (Math::Sin(p[1]) * (x_pix[i] - p[2]) + Math::Cos(p[1]) * (y_pix[i] - p[3]));
			xprime -= x_intrmdt[i];
			yprime -= y_intrmdt[i];
			f[0] += xprime * xprime + yprime * yprime;
			//f[0] += (xprime - x_intrmdt[i])*(xprime - x_intrmdt[i]) + (yprime - y_intrmdt[i])*(yprime - y_intrmdt[i]);
		}
	}

	if (p->Length == 6)
	{
		for (int i = 0; i < x_pix->Length; i++)
		{
			xprime = p[0] * (x_pix[i] - p[4]) + p[1] * (y_pix[i] - p[5]);
			yprime = p[2] * (x_pix[i] - p[4]) + p[3] * (y_pix[i] - p[5]);
			xprime -= x_intrmdt[i];
			yprime -= y_intrmdt[i];
			f[0] += xprime * xprime + yprime * yprime;
			//f[0] += (xprime - xref[i])*(xprime - xref[i]) + (yprime - yref[i])*(yprime - yref[i]);
		}
	}
}

void JPFITS::JPMath::Fit_GeneralTransform2d(array<double>^ x_ref, array<double>^ y_ref, array<double>^ x_tran, array<double>^ y_tran, array<double>^ &p, array<double>^ p_lbnd, array<double>^ p_ubnd, array<double>^ p_scale)
{
	double epsx = 0.00000000001;
	int maxits = 0;
	alglib::minlmstate^ state;
	alglib::minlmreport^ report;
	alglib::ndimensional_fvec^ minfunc = gcnew alglib::ndimensional_fvec(alglib_GeneralTransform2d);
	array<Object^>^ objj = gcnew array<Object^>(4);
	objj[0] = (Object^)x_tran;
	objj[1] = (Object^)y_tran;
	objj[2] = (Object^)x_ref;
	objj[3] = (Object^)y_ref;
	Object^ obj = (Object^)objj;

	alglib::minlmcreatev(1, p, 0.00001, state);
	alglib::minlmsetcond(state, epsx, maxits);
	alglib::minlmsetscale(state, p_scale);
	alglib::minlmsetbc(state, p_lbnd, p_ubnd);
	alglib::minlmoptimize(state, minfunc, nullptr, obj);
	alglib::minlmresults(state, p, report);
}

void JPFITS::JPMath::alglib_GeneralTransform2d(array<double>^ p, array<double>^ f, Object^ obj)
{
	array<Object^>^ objj = (array<Object^>^)obj;
	array<double>^ x_tran = (array<double>^)objj[0];
	array<double>^ y_tran = (array<double>^)objj[1];
	array<double>^ x_ref = (array<double>^)objj[2];
	array<double>^ y_ref = (array<double>^)objj[3];
	f[0] = 0;
	double xprime, yprime;

	if (p->Length == 6)
	{
		for (int i = 0; i < x_ref->Length; i++)
		{
			xprime = p[0] * (Math::Cos(p[1]) * (x_tran[i] - p[2]) - Math::Sin(p[1]) * (y_tran[i] - p[3])) + p[2] + p[4];
			yprime = p[0] * (Math::Sin(p[1]) * (x_tran[i] - p[2]) + Math::Cos(p[1]) * (y_tran[i] - p[3])) + p[3] + p[5];
			xprime -= x_ref[i];
			yprime -= y_ref[i];
			f[0] += xprime * xprime + yprime * yprime;
			//f[0] += (xprime - x_intrmdt[i])*(xprime - x_intrmdt[i]) + (yprime - y_intrmdt[i])*(yprime - y_intrmdt[i]);
		}
	}

	if (p->Length == 8)
	{
		for (int i = 0; i < x_tran->Length; i++)
		{
			xprime = p[0] * (x_tran[i] - p[4]) + p[1] * (y_tran[i] - p[5]) + p[4] + p[6];
			yprime = p[2] * (x_tran[i] - p[4]) + p[3] * (y_tran[i] - p[5]) + p[5] + p[7];
			xprime -= x_ref[i];
			yprime -= y_ref[i];
			f[0] += xprime * xprime + yprime * yprime;
			//f[0] += (xprime - xref[i])*(xprime - xref[i]) + (yprime - yref[i])*(yprime - yref[i]);
		}
	}
}

/*void JPFITS::JPMath::alglib_Transform2d_jac(array<double>^ p, array<double>^ f, array<double, 2>^ jac, Object^ obj)
{

}*/

void JPFITS::JPMath::Fit_Poly1d(array<double>^ xdata, array<double>^ ydata, int poly_degree, bool robust, array<double>^ &poly_coeffs)
{
	array<double>^ weights = gcnew array<double>(xdata->Length);
	for (int i = 0; i < xdata->Length; i++)
		weights[i] = 1.0;

	array<double>^ xc = gcnew array<double>(0);
	array<double>^ yc = gcnew array<double>(0);
	array<int>^ dc = gcnew array<int>(0);

	int m = poly_degree + 1;
	int info;
	alglib::barycentricinterpolant^ p;
	alglib::polynomialfitreport^ rep;

	alglib::polynomialfitwc(xdata, ydata, weights, xc, yc, dc, m, info, p, rep);
	alglib::polynomialbar2pow(p, poly_coeffs);

	if (!robust || xdata->Length <= 2)
		return;

	//if robust then determine some weights via the residuals and recalculate solution a few times until some convergence criteria is found
	int iteration_count = 0;
	double sigma = rep->rmserror, yfit, rmsrat = Double::MaxValue;

	while (Math::Abs(rmsrat - 1) > 0.0000001 && iteration_count < 50)
	{
		rmsrat = rep->rmserror;//get the previous rms

		for (int i = 0; i < xdata->Length; i++)
		{
			yfit = alglib::barycentriccalc(p, xdata[i]);
			weights[i] = Math::Exp(-(ydata[i] - yfit)*(ydata[i] - yfit) / (2 * sigma*sigma));
			weights[i] *= weights[i];
		}

		alglib::polynomialfitwc(xdata, ydata, weights, xc, yc, dc, m, info, p, rep);

		sigma = rep->rmserror;
		rmsrat /= rep->rmserror;
		iteration_count++;
	}
	alglib::polynomialbar2pow(p, poly_coeffs);
}

