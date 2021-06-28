/*Copyright 2021 Joseph Edwin Postma

Contact email: joepostma@live.ca

This file is part of JPFITS.

JPFITS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

JPFITS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

See http://www.gnu.org/licenses/. */


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

array<double, 2>^ JPFITS::JPMath::Excise(array<double, 2>^ data, bool column, int X0, int halfWidth, bool do_parallel)
{
	array<double, 2>^ result;
	array<int>^ xrange;
	array<int>^ yrange;
	if (column)
	{
		result = gcnew array<double, 2>(data->GetLength(0) - (halfWidth * 2 + 1), data->GetLength(1));
		yrange = gcnew array<int>(result->GetLength(1));
		for (int i = 0; i < yrange->Length; i++)
			yrange[i] = i;
		
		xrange = gcnew array<int>(result->GetLength(0));
		for (int i = 0; i < data->GetLength(0); i++)
			if (i < X0 - halfWidth)
				xrange[i] = i;
			else if (i > X0 + halfWidth)
				xrange[i - (halfWidth * 2 + 1)] = i;
	}
	else
	{
		result = gcnew array<double, 2>(data->GetLength(0), data->GetLength(1) - (halfWidth * 2 + 1));
		xrange = gcnew array<int>(result->GetLength(0));
		for (int i = 0; i < xrange->Length; i++)
			xrange[i] = i;

		yrange = gcnew array<int>(result->GetLength(1));
		for (int i = 0; i < data->GetLength(1); i++)
			if (i < X0 - halfWidth)
				yrange[i] = i;
			else if (i > X0 + halfWidth)
				yrange[i - (halfWidth * 2 + 1)] = i;
	}

	#pragma omp parallel for if (do_parallel)
	for (int x = 0; x < result->GetLength(0); x++)
		for (int y = 0; y < result->GetLength(1); y++)
			result[x, y] = data[xrange[x], yrange[y]];

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

double JPFITS::JPMath::Max(array<double> ^data, int startIndex, int endIndex, int &maxIndex, bool do_parallel)
{
	if (startIndex < 0)
		startIndex = 0;
	if (endIndex > data->Length)
		endIndex = data->Length - 1;

	double max = System::Double::MinValue;

	#pragma omp parallel for if (do_parallel)
	for (int i = startIndex; i <= endIndex; i++)
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
	for (int i = startIndex; i <= endIndex; i++)
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

	if (Double::IsPositiveInfinity(val))
		method = 7;
	
	if (Double::IsNegativeInfinity(val))
		method = 8;

	array<bool, 2>^ finds = gcnew array<bool, 2>(data->GetLength(0), data->GetLength(1));
	int Nfinds = 0;

	#pragma omp parallel for if (do_parallel) reduction (+: Nfinds)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
		{
			switch (method)
			{
				case 0://<
				{
					if (data[i, j] < val)
					{
						finds[i, j] = true;
						Nfinds++;
						/*#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}*/
					}
					break;
				}
				case 1://<=
				{
					if (data[i, j] <= val)
					{
						finds[i, j] = true;
						Nfinds++;
						/*#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}*/
					}
					break;
				}
				case 2://== || =
				{
					if (data[i, j] == val)
					{
						finds[i, j] = true;
						Nfinds++;
						/*#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}*/
					}
					break;
				}
				case 3://>=
				{
					if (data[i, j] >= val)
					{
						finds[i, j] = true;
						Nfinds++;
						/*#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}*/
					}
					break;
				}
				case 4://>
				{
					if (data[i, j] > val)
					{
						finds[i, j] = true;
						Nfinds++;
						/*#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}*/
					}
					break;
				}
				case 5://!=
				{
					if (data[i, j] != val)
					{
						finds[i, j] = true;
						Nfinds++;
						/*#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}*/
					}
					break;
				}
				case 6://= NaN
				{
					if (Double::IsNaN(data[i, j]))
					{
						finds[i, j] = true;
						Nfinds++;
						/*#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}*/
					}
					break;
				}
				case 7://= +ve inf
				{
					if (Double::IsPositiveInfinity(data[i, j]))
					{
						finds[i, j] = true;
						Nfinds++;
						/*#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}*/
					}
					break;
				}
				case 8://= -ve inf
				{
					if (Double::IsNegativeInfinity(data[i, j]))
					{
						finds[i, j] = true;
						Nfinds++;
						/*#pragma omp critical
						{
							ptslist->Add(i);
							ptslist->Add(j);
						}*/
					}
					break;
				}
			}
		}

	/*int N = int(double(ptslist->Count) / 2);
	array<int, 2>^ result = gcnew array<int, 2>(N, 2);	

	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < N; i++)
	{
		result[i, 0] = System::Convert::ToInt32(ptslist[i * 2]);//x
		result[i, 1] = System::Convert::ToInt32(ptslist[i * 2 + 1]);//y
	}*/

	array<int, 2>^ result = gcnew array<int, 2>(Nfinds, 2);
	Nfinds = 0;
	//#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < data->GetLength(0); i++)
		for (int j = 0; j < data->GetLength(1); j++)
			if (finds[i, j])
			{
				//#pragma omp critical
				{
					result[Nfinds, 0] = i;//x
					result[Nfinds, 1] = j;//y
					Nfinds++;
				}
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

/*double JPFITS::JPMath::Sum(array<double, 1> ^data, bool do_parallel)
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
}*/

double JPFITS::JPMath::Sum(Object^ vectorOrArray, bool do_parallel)
{
	TypeCode type = Type::GetTypeCode((((Array^)vectorOrArray)->GetType())->GetElementType());
	int rank = ((Array^)vectorOrArray)->Rank;
	double res = 0;
	int naxis0, naxis1;
	if (rank == 1)
		naxis0 = ((Array^)vectorOrArray)->Length;
	else
	{
		naxis0 = ((Array^)vectorOrArray)->GetLength(0);
		naxis1 = ((Array^)vectorOrArray)->GetLength(1);
	}

	switch (type)
	{
		case TypeCode::Double:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					res += ((array<double>^)vectorOrArray)[x];
			}
			else //rank == 2
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					for (int y = 0; y < naxis1; y++)
						res += ((array<double, 2>^)vectorOrArray)[x, y];
			}
			return res;
		}

		case TypeCode::Single:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					res += (double)((array<float>^)vectorOrArray)[x];
			}
			else //rank == 2
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					for (int y = 0; y < naxis1; y++)
						res += (double)((array<float, 2>^)vectorOrArray)[x, y];
			}
			return res;
		}

		case TypeCode::UInt64:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					res += (double)((array<unsigned __int64>^)vectorOrArray)[x];
			}
			else //rank == 2
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					for (int y = 0; y < naxis1; y++)
						res += (double)((array<unsigned __int64, 2>^)vectorOrArray)[x, y];
			}
			return res;
		}

		case TypeCode::Int64:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					res += (double)((array<__int64>^)vectorOrArray)[x];
			}
			else //rank == 2
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					for (int y = 0; y < naxis1; y++)
						res += (double)((array<__int64, 2>^)vectorOrArray)[x, y];
			}
			return res;
		}

		case TypeCode::UInt32:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					res += (double)((array<unsigned __int32>^)vectorOrArray)[x];
			}
			else //rank == 2
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					for (int y = 0; y < naxis1; y++)
						res += (double)((array<unsigned __int32, 2>^)vectorOrArray)[x, y];
			}
			return res;
		}

		case TypeCode::Int32:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					res += (double)((array<__int32>^)vectorOrArray)[x];
			}
			else //rank == 2
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					for (int y = 0; y < naxis1; y++)
						res += (double)((array<__int32, 2>^)vectorOrArray)[x, y];
			}
			return res;
		}

		case TypeCode::UInt16:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					res += (double)((array<unsigned __int16>^)vectorOrArray)[x];
			}
			else //rank == 2
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					for (int y = 0; y < naxis1; y++)
						res += (double)((array<unsigned __int16, 2>^)vectorOrArray)[x, y];
			}
			return res;
		}

		case TypeCode::Int16:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					res += (double)((array<__int16>^)vectorOrArray)[x];
			}
			else //rank == 2
			{
				#pragma omp parallel for if (do_parallel) reduction(+:res)
				for (int x = 0; x < naxis0; x++)
					for (int y = 0; y < naxis1; y++)
						res += (double)((array<__int16, 2>^)vectorOrArray)[x, y];
			}
			return res;
		}

		default:
		{
			throw gcnew Exception("Typecode '" + type.ToString() + "' not supported for Sum.");
			return 0;
		}
	}
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

/*double JPFITS::JPMath::Mean(array<double, 1> ^data, bool do_parallel)
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
}*/

double JPFITS::JPMath::Mean(Object^ vectorOrArray, bool do_parallel)
{
	return JPMath::Sum(vectorOrArray, do_parallel) / double(((Array^)vectorOrArray)->Length);
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
	return (System::Math::IEEERemainder(double(x), 2) == 0);

	/*if (System::Math::IEEERemainder(double(x), 2) == 0)
		return true;
	else
		return false;*/
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
		throw gcnew Exception("Interpolation style '" + style + "' not recognized.");

	void *x = 0;
	alglib::spline1dinterpolant^ sp = gcnew alglib::spline1dinterpolant(x);

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
	double epsx = 0.000000001;
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

void JPFITS::JPMath::Fit_PointSource_Compound(String^ model_name, String^ minimization_type, array<int>^ xdata, array<int>^ ydata, array<double, 2>^ source, array<double>^ xpositions, array<double>^ ypositions, double position_radius, array<double, 2>^ &params, array<double, 2>^ &p_err, array<double, 2>^ &fit_residuals, String^ &termination_msg)
{
	if (params->GetLength(1) != xpositions->Length)
	{
		throw gcnew Exception("Parameter array params not consistent with the number of source xpositions indicated.");
		return;
	}

	try
	{
		alglib::ndimensional_grad^ grad;

		array<double>^ XDATA = gcnew array<double>(xdata->Length);
		array<double>^ YDATA = gcnew array<double>(ydata->Length);
		int func = params->GetLength(0);
		int count = params->GetLength(1);
		array<int>^ fcobj = gcnew array<int>(2) { func, count };
		array<double>^ P0I = gcnew array<double>(params->Length - count + 1);//remove repeated background, but keep one
		array<double>^ PLB = gcnew array<double>(params->Length - count + 1);//remove repeated background, but keep one
		array<double>^ PUB = gcnew array<double>(params->Length - count + 1);//remove repeated background, but keep one
		array<double>^ scale = gcnew array<double>(params->Length - count + 1);//remove repeated background, but keep one

		for (int i = 0; i < XDATA->Length; i++)
			XDATA[i] = double(xdata[i]);
		for (int i = 0; i < YDATA->Length; i++)
			YDATA[i] = double(ydata[i]);

		double min, max;
		JPMath::MinMax(source, min, max, false);
		double bias = min;
		if (bias < 1)
			bias = 1;
		P0I[P0I->Length - 1] = bias;//background
		PLB[P0I->Length - 1] = 0;//background
		PUB[P0I->Length - 1] = max / 4;//background
		scale[P0I->Length - 1] = bias;//background

		int funcindex;
		if (model_name == "Gaussian")
		{
			if (func == 5)
			{
				for (int i = 0; i < count; i++)
				{
					funcindex = i * (func - 1);

					//amplitude
					P0I[funcindex] = source[(int)(xpositions[i] - XDATA[0] + 0.5), (int)(ypositions[i] - YDATA[0] + 0.5)] - bias;
					PLB[funcindex] = P0I[funcindex] / 3;
					PUB[funcindex] = 2 * P0I[funcindex];
					scale[funcindex] = P0I[funcindex];

					//x0
					PLB[funcindex + 1] = xpositions[i] - position_radius;
					P0I[funcindex + 1] = xpositions[i];
					PUB[funcindex + 1] = xpositions[i] + position_radius;
					scale[funcindex + 1] = P0I[funcindex + 1];

					//y0
					PLB[funcindex + 2] = ypositions[i] - position_radius;
					P0I[funcindex + 2] = ypositions[i];
					PUB[funcindex + 2] = ypositions[i] + position_radius;
					scale[funcindex + 2] = P0I[funcindex + 2];

					//sigma
					PLB[funcindex + 3] = 1e-1;
					P0I[funcindex + 3] = 2;
					PUB[funcindex + 3] = 10;
					scale[funcindex + 3] = P0I[funcindex + 3];
				}
			}
			else if (func == 7)
			{
				for (int i = 0; i < count; i++)
				{
					funcindex = i * (func - 1);

					//amplitude
					P0I[funcindex] = source[(int)(xpositions[i] - XDATA[0]), (int)(ypositions[i] - YDATA[0])] - bias;
					PLB[funcindex] = P0I[funcindex] / 3;
					PUB[funcindex] = 2 * P0I[funcindex];
					scale[funcindex] = P0I[funcindex];

					//x0
					PLB[funcindex + 1] = xpositions[i] - position_radius;
					P0I[funcindex + 1] = xpositions[i];
					PUB[funcindex + 1] = xpositions[i] + position_radius;
					scale[funcindex + 1] = P0I[funcindex + 1];

					//y0
					PLB[funcindex + 2] = ypositions[i] - position_radius;
					P0I[funcindex + 2] = ypositions[i];
					PUB[funcindex + 2] = ypositions[i] + position_radius;
					scale[funcindex + 2] = P0I[funcindex + 2];

					//phi
					PLB[funcindex + 3] = -Math::PI;
					P0I[funcindex + 3] = 0;
					PUB[funcindex + 3] = Math::PI;
					scale[funcindex + 3] = 1;

					//sigmaX
					PLB[funcindex + 4] = 1e-1;
					P0I[funcindex + 4] = 2;
					PUB[funcindex + 4] = 10;
					scale[funcindex + 4] = P0I[funcindex + 4];

					//sigmaY
					PLB[funcindex + 5] = 1e-1;
					P0I[funcindex + 5] = 2;
					PUB[funcindex + 5] = 10;
					scale[funcindex + 5] = P0I[funcindex + 5];
				}
			}
			else
			{
				throw gcnew Exception("Parameter length does not correspond to either circular (5 params) or elliptical (7 params) Gaussian; params length = " + func);
				return;
			}

			if (minimization_type->ToLower() == "ls")
				grad = gcnew alglib::ndimensional_grad(alglib_Gauss_2d_LM_LS_grad_compound);
			else if (minimization_type->ToLower() == "chisq")
				grad = gcnew alglib::ndimensional_grad(alglib_Gauss_2d_LM_LS_CHISQ_grad_compound);
			else if (minimization_type->ToLower() == "robust")
				grad = gcnew alglib::ndimensional_grad(alglib_Gauss_2d_LM_LS_ROBUST_grad_compound);
			else if (minimization_type->ToLower() == "cstat")
				grad = gcnew alglib::ndimensional_grad(alglib_Gauss_2d_LM_LS_CSTAT_grad_compound);
			else
			{
				throw gcnew Exception("Fit Type not recognized: '" + minimization_type);
				return;
			}
		}
		else if (model_name == "Moffat")
		{
			if (func == 6)
			{
				for (int i = 0; i < count; i++)
				{
					funcindex = i * (func - 1);

					//amplitude
					P0I[funcindex] = source[(int)(xpositions[i] - XDATA[0]), (int)(ypositions[i] - YDATA[0])] - bias;
					PLB[funcindex] = P0I[funcindex] / 3;
					PUB[funcindex] = 2 * P0I[funcindex];
					scale[funcindex] = P0I[funcindex];

					//x0
					PLB[funcindex + 1] = xpositions[i] - position_radius;
					P0I[funcindex + 1] = xpositions[i];
					PUB[funcindex + 1] = xpositions[i] + position_radius;
					scale[funcindex + 1] = P0I[funcindex + 1];

					//y0
					PLB[funcindex + 2] = ypositions[i] - position_radius;
					P0I[funcindex + 2] = ypositions[i];
					PUB[funcindex + 2] = ypositions[i] + position_radius;
					scale[funcindex + 2] = P0I[funcindex + 2];

					//alpha
					PLB[funcindex + 3] = 1e-1;
					P0I[funcindex + 3] = 2;
					PUB[funcindex + 3] = 10;
					scale[funcindex + 3] = P0I[funcindex + 3];

					//beta
					PLB[funcindex + 4] = 1 + 1e-6;
					P0I[funcindex + 4] = 2;
					PUB[funcindex + 4] = 10;
					scale[funcindex + 4] = P0I[funcindex + 4];
				}
			}
			else if (func == 8)
			{
				for (int i = 0; i < count; i++)
				{
					funcindex = i * (func - 1);

					//amplitude
					P0I[funcindex] = source[(int)(xpositions[i] - XDATA[0]), (int)(ypositions[i] - YDATA[0])] - bias;
					PLB[funcindex] = P0I[funcindex] / 3;					
					PUB[funcindex] = 2 * P0I[funcindex];
					scale[funcindex] = P0I[funcindex];

					//x0
					PLB[funcindex + 1] = xpositions[i] - position_radius;
					P0I[funcindex + 1] = xpositions[i];
					PUB[funcindex + 1] = xpositions[i] + position_radius;
					scale[funcindex + 1] = P0I[funcindex + 1];

					//y0
					PLB[funcindex + 2] = ypositions[i] - position_radius;
					P0I[funcindex + 2] = ypositions[i];
					PUB[funcindex + 2] = ypositions[i] + position_radius;
					scale[funcindex + 2] = P0I[funcindex + 2];

					//phi
					PLB[funcindex + 3] = -Math::PI;
					P0I[funcindex + 3] = 0;
					PUB[funcindex + 3] = Math::PI;
					scale[funcindex + 3] = 1;

					//alphaX
					PLB[funcindex + 4] = 1e-1;
					P0I[funcindex + 4] = 2;
					PUB[funcindex + 4] = 10;
					scale[funcindex + 4] = P0I[funcindex + 4];

					//alphaY
					PLB[funcindex + 5] = 1e-1;
					P0I[funcindex + 5] = 2;
					PUB[funcindex + 5] = 10;
					scale[funcindex + 5] = P0I[funcindex + 5];

					//beta
					PLB[funcindex + 6] = 1 + 1e-6;
					P0I[funcindex + 6] = 2;
					PUB[funcindex + 6] = 10;
					scale[funcindex + 6] = P0I[funcindex + 6];
				}
			}
			else
			{
				throw gcnew Exception("Parameter length does not correspond to either circular (6 params) or elliptical (8 params) Moffat; params length = " + func);
				return;
			}

			if (minimization_type->ToLower() == "ls")
				/*grad = gcnew alglib::ndimensional_grad(alglib_Moffat_2d_LM_LS_compound_grad)*/;
			else if (minimization_type->ToLower() == "chisq")
				/*grad = gcnew alglib::ndimensional_grad(alglib_Moffat_2d_LM_LS_compound_CHISQ_grad)*/;
			else if (minimization_type->ToLower() == "robust")
				/*grad = gcnew alglib::ndimensional_grad(alglib_Moffat_2d_LM_LS_compound_ROBUST_grad)*/;
			else if (minimization_type->ToLower() == "cstat")
				/*grad = gcnew alglib::ndimensional_grad(alglib_Moffat_2d_LM_LS_compound_CSTAT_grad)*/;
			else
			{
				throw gcnew Exception("Fit Type not recognized: '" + minimization_type);
				return;
			}
		}
		else
		{
			throw gcnew Exception("Fit Model not recognized: '" + model_name);
			return;
		}

		double epsx = 1e-6;
		double epsg = 0;
		double epsf = 0;
		int maxits = 2000;
		array<Object^>^ arrays = gcnew array<Object^>(4);
		arrays[0] = source;
		arrays[1] = XDATA;
		arrays[2] = YDATA;
		arrays[3] = fcobj;

		alglib::minbcstate^ bcstate;
		alglib::minbccreate(P0I, bcstate);
		alglib::minbcsetcond(bcstate, epsg, epsf, epsx, maxits);
		alglib::minbcsetbc(bcstate, PLB, PUB);
		alglib::minbcsetscale(bcstate, scale);
		alglib::minbcoptimize(bcstate, grad, nullptr, arrays, alglib::parallel);
		alglib::minbcreport^ report;
		alglib::minbcresults(bcstate, P0I, report);

		switch (report->terminationtype)
		{
			case -8:
			{
				termination_msg = "Internal integrity control detected infinite or NAN values in function or gradient. Abnormal termination signalled.";
				break;
			}
			case -3:
			{
				termination_msg = "Inconsistent constraints.";
				break;
			}
			case 1:
			{
				termination_msg = "Relative function improvement is no more than EpsF.";
				break;
			}
			case 2:
			{
				termination_msg = "Scaled step is no more than EpsX.";
				break;
			}
			case 4:
			{
				termination_msg = "Scaled gradient norm is no more than EpsG.";
				break;
			}
			case 5:
			{
				termination_msg = "MaxIts steps was taken.";
				break;
			}
			case 7:
			{
				termination_msg = "Stopping conditions are too stringent, further improvement is impossible, X contains best point found so far.";
				break;
			}
			case 8:
			{
				termination_msg = "Terminated by user.";
				break;
			}
		}

		for (int i = 0; i < count; i++)
			for (int j = 0; j < func - 1; j++)
				params[j, i] = P0I[j + i * (func - 1)];
		for (int i = 0; i < count; i++)
			params[func - 1, i] = P0I[P0I->Length - 1];//background

		/*if (p_err->Length != 0)
		{
			if (model_name == "Gaussian")
				p_err = Gauss_2D_param_err(params, XDATA, YDATA, source);
			else
				p_err = Moffat_2D_param_err(params, XDATA, YDATA, source);
		}*/

		if (fit_residuals->Length != 0)
		{
			double val;
			array<double>^ X = gcnew array<double>(2);
			for (int i = 0; i < XDATA->Length; i++)
				for (int j = 0; j < YDATA->Length; j++)
				{
					X[0] = XDATA[i];
					X[1] = YDATA[j];
					alglib_Gauss_2d_compound(P0I, X, val, fcobj);
					fit_residuals[i, j] = source[i, j] - val;
				}
		}
	}
	catch (Exception^ e)
	{
		MessageBox::Show(e->Data + "	" + e->InnerException + "	" + e->Message + "	" + e->Source + "	" + e->StackTrace + "	" + e->TargetSite);
	}
}

void JPFITS::JPMath::Fit_PointSource(String^ model_name, String^ minimization_type, array<int>^ xdata, array<int>^ ydata, array<double, 2>^ source, array<double>^ &params, array<double>^ params_LB, array<double>^ params_UB, array<double>^ &p_err, array<double, 2>^ &fit_residuals, String^ &termination_msg)
{
	alglib::ndimensional_grad^ grad;
	array<double>^ scale;

	array<double>^ XDATA = gcnew array<double>(xdata->Length);
	array<double>^ YDATA = gcnew array<double>(ydata->Length);

	for (int i = 0; i < XDATA->Length; i++)
		XDATA[i] = double(xdata[i]);
	for (int i = 0; i < YDATA->Length; i++)
		YDATA[i] = double(ydata[i]);
		
	int maxx, maxy;
	double max = JPMath::Max(source, maxx, maxy, false);
	double min = JPMath::Min(source, false);
	double amp = max - min;
	if (amp == 0)
		amp = 1;
	double x0 = XDATA[maxx];
	if (x0 == 0)
		x0 = 1;
	double y0 = YDATA[maxy];
	if (y0 == 0)
		y0 = 1;
	double bias = min;
	if (bias < 1)
		bias = 1;	

	if (model_name == "Gaussian")
	{
		if (params->Length == 5)
		{
			scale = gcnew array<double>(5) { amp, x0, y0, 2, bias };
			if (params_LB == nullptr || params_LB->Length == 0)
				params_LB = gcnew array<double>(5) { amp / 3, XDATA[0], YDATA[0], 0.1, 0 };
			if (params_UB == nullptr || params_UB->Length == 0)
				params_UB = gcnew array<double>(5) { 2 * amp, XDATA[XDATA->Length - 1], YDATA[YDATA->Length - 1], 10, max / 4 };
			if (params[0] == 0 && params[1] == 0 && params[2] == 0 && params[3] == 0 && params[4] == 0)
			{
				params[0] = amp;
				params[1] = x0;
				params[2] = y0;
				params[3] = 2;
				params[4] = bias;
			}
		}
		else if (params->Length == 7)
		{
			scale = gcnew array<double>(7) { amp, x0, y0, 1, 2, 2, bias };
			if (params_LB == nullptr || params_LB->Length == 0)
				params_LB = gcnew array<double>(7) { amp / 3, XDATA[0], YDATA[0], -Math::PI, 0.1, 0.1, 0 };
			if (params_UB == nullptr || params_UB->Length == 0)
				params_UB = gcnew array<double>(7) { 2 * amp, XDATA[XDATA->Length - 1], YDATA[YDATA->Length - 1], Math::PI, 10, 10, max / 4 };
			if (params[0] == 0 && params[1] == 0 && params[2] == 0 && params[3] == 0 && params[4] == 0 && params[5] == 0 && params[6] == 0)
			{
				params[0] = amp;
				params[1] = x0;
				params[2] = y0;
				params[3] = 0;
				params[4] = 2;
				params[5] = 2;
				params[6] = bias;
			}
		}
		else
		{
			throw gcnew Exception("Parameter length does not correspond to either circular (5 params) or elliptical (7 params) Gaussian; params length = " + params->Length);
			return;
		}

		if (minimization_type->ToLower() == "ls")
			grad = gcnew alglib::ndimensional_grad(alglib_Gauss_2d_LM_LS_grad);
		else if (minimization_type->ToLower() == "chisq")
			grad = gcnew alglib::ndimensional_grad(alglib_Gauss_2d_LM_LS_CHISQ_grad);
		else if (minimization_type->ToLower() == "robust")
			grad = gcnew alglib::ndimensional_grad(alglib_Gauss_2d_LM_LS_ROBUST_grad);
		else if (minimization_type->ToLower() == "cstat")
			grad = gcnew alglib::ndimensional_grad(alglib_Gauss_2d_LM_LS_CSTAT_grad);
		else
		{
			throw gcnew Exception("Fit Type not recognized: '" + minimization_type);
			return;
		}
	}
	else if (model_name == "Moffat")
	{
		if (params->Length == 6)
		{
			scale = gcnew array<double>(6) { amp, x0, y0, 2, 2, bias };
			if (params_LB == nullptr || params_LB->Length == 0)
				params_LB = gcnew array<double>(6) { amp / 3, XDATA[0], YDATA[0], 0.1, 1.000001, 0 };
			if (params_UB == nullptr || params_UB->Length == 0)
				params_UB = gcnew array<double>(6) { 2 * amp, XDATA[XDATA->Length - 1], YDATA[YDATA->Length - 1], 10, 10, max / 4 };
			if (params[0] == 0 && params[1] == 0 && params[2] == 0 && params[3] == 0 && params[4] == 0 && params[5] == 0)
			{
				params[0] = amp;
				params[1] = x0;
				params[2] = y0;
				params[3] = 2;
				params[4] = 2;
				params[5] = bias;
			}
		}
		else if (params->Length == 8)
		{
			scale = gcnew array<double>(8) { amp, x0, y0, 1, 2, 2, 2, bias };
			if (params_LB == nullptr || params_LB->Length == 0)
				params_LB = gcnew array<double>(8) { amp / 3, XDATA[0], YDATA[0], -Math::PI, 0.1, 0.1, 1.000001, 0 };
			if (params_UB == nullptr || params_UB->Length == 0)
				params_UB = gcnew array<double>(8) { 2 * amp, XDATA[XDATA->Length - 1], YDATA[YDATA->Length - 1], Math::PI, 10, 10, 10, max / 4 };
			if (params[0] == 0 && params[1] == 0 && params[2] == 0 && params[3] == 0 && params[4] == 0 && params[5] == 0 && params[6] == 0 && params[7] == 0)
			{
				params[0] = amp;
				params[1] = x0;
				params[2] = y0;
				params[3] = 0;
				params[4] = 2;
				params[5] = 2;
				params[6] = 2;
				params[7] = bias;
			}
		}
		else
		{
			throw gcnew Exception("Parameter length does not correspond to either circular (6 params) or elliptical (8 params) Moffat; params length = " + params->Length);
			return;
		}

		if (minimization_type->ToLower() == "ls")
			grad = gcnew alglib::ndimensional_grad(alglib_Moffat_2d_LM_LS_grad);
		else if (minimization_type->ToLower() == "chisq")
			grad = gcnew alglib::ndimensional_grad(alglib_Moffat_2d_LM_LS_CHISQ_grad);
		else if (minimization_type->ToLower() == "robust")
			grad = gcnew alglib::ndimensional_grad(alglib_Moffat_2d_LM_LS_ROBUST_grad);
		else if (minimization_type->ToLower() == "cstat")
			grad = gcnew alglib::ndimensional_grad(alglib_Moffat_2d_LM_LS_CSTAT_grad);
		else
		{
			throw gcnew Exception("Fit Type not recognized: '" + minimization_type);
			return;
		}
	}
	else
	{
		throw gcnew Exception("Fit Model not recognized: '" + model_name);
		return;
	}

	double epsx = 1e-6;
	double epsg = 0;
	double epsf = 0;
	int maxits = 1000;
	array<Object^>^ arrays = gcnew array<Object^>(3);
	arrays[0] = source;
	arrays[1] = XDATA;
	arrays[2] = YDATA;

	alglib::minbcstate^ bcstate;
	alglib::minbccreate(params, bcstate);	
	alglib::minbcsetcond(bcstate, epsg, epsf, epsx, maxits);
	alglib::minbcsetbc(bcstate, params_LB, params_UB);
	alglib::minbcsetscale(bcstate, scale);	
	alglib::minbcoptimize(bcstate, grad, nullptr, arrays, alglib::parallel);
	alglib::minbcreport^ report;
	alglib::minbcresults(bcstate, params, report);

	switch (report->terminationtype)
	{
		case -8:
		{
			termination_msg = "Internal integrity control detected infinite or NAN values in function or gradient. Abnormal termination signalled.";
			break;
		}
		case -3:
		{
			termination_msg = "Inconsistent constraints.";
			break;
		}
		case 1:
		{
			termination_msg = "Relative function improvement is no more than EpsF.";
			break;
		}
		case 2:
		{
			termination_msg = "Scaled step is no more than EpsX.";
			break;
		}
		case 4:
		{
			termination_msg = "Scaled gradient norm is no more than EpsG.";
			break;
		}
		case 5:
		{
			termination_msg = "MaxIts steps was taken.";
			break;
		}
		case 7:
		{
			termination_msg = "Stopping conditions are too stringent, further improvement is impossible, X contains best point found so far.";
			break;
		}
		case 8:
		{
			termination_msg = "Terminated by user.";
			break;
		}
	}

	if (p_err->Length != 0)
	{
		if (model_name == "Gaussian")
			p_err = Gauss_2D_param_err(params, XDATA, YDATA, source);
		else
			p_err = Moffat_2D_param_err(params, XDATA, YDATA, source);
	}

	if (fit_residuals->Length != 0)
	{
		array<double, 2>^ model = gcnew array<double, 2>(XDATA->Length, YDATA->Length);
		if (model_name->ToLower()->Contains("gaus"))
			Gaussian2d(xdata, ydata, model, params, false);
		else
			Moffat2d(xdata, ydata, model, params, false);

		for (int i = 0; i < XDATA->Length; i++)
			for (int j = 0; j < YDATA->Length; j++)
				fit_residuals[i, j] = source[i, j] - model[i, j];
	}
}

void JPFITS::JPMath::alglib_Gauss_2d_LM_LS_grad(array<double>^ p, double %f, array<double>^ grad, Object^ obj)
{
	f = 0;
	for (int i = 0; i < grad->Length; i++)
		grad[i] = 0;
	array<double, 2>^ Z = (array<double, 2>^)(((array<Object^>^)obj)[0]);
	array<double>^ x = (array<double>^)(((array<Object^>^)obj)[1]);
	array<double>^ y = (array<double>^)(((array<Object^>^)obj)[2]);

	if (p->Length == 5)
	{
		double p3sq = p[3] * p[3];
		double twop3sq = 2 * p3sq;
		double p3cu = p3sq * p[3];
		double twop1 = 2 * p[1];
		double twop2 = 2 * p[2];
		double dx, dxsq, dy, dxsqpdysq, expres, exparg, p0expres, a;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			dxsq = dx * dx;

			for (int j = 0; j < y->Length; j++)
			{
				dy = p[2] - y[j];
				dxsqpdysq = dxsq + dy * dy;
				exparg = -dxsqpdysq / twop3sq;
				expres = Math::Exp(exparg);
				p0expres = p[0] * expres;
				a = p[4] - Z[i, j] + p0expres;

				f += a * a;

				grad[0] += 2 * expres*a;
				grad[1] += -(p0expres*(twop1 - 2 * x[i])*a) / p3sq;
				grad[2] += -(p0expres*(twop2 - 2 * y[j])*a) / p3sq;
				grad[3] += (2 * p0expres*dxsqpdysq*a) / p3cu;
				grad[4] += 2 * a;
			}
		}
		return;
	}

	if (p->Length == 7)
	{
		double cosp3 = Math::Cos(p[3]);
		double sinp3 = Math::Sin(p[3]);
		double p4sq = p[4] * p[4];
		double twop4sq = 2 * p4sq;
		double p5sq = p[5] * p[5];
		double twop5sq = 2 * p5sq;
		double p4cu = p4sq * p[4];
		double p5cu = p5sq * p[5];
		double dx, cosp3dx, sinp3dx, dy, sinp3dy, a, cosp3dy, b, c, exparg, asq, bsq, expexparg, p0expexparg, twop0expexparg;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			cosp3dx = cosp3 * dx;
			sinp3dx = sinp3 * dx;

			for (int j = 0; j < y->Length; j++)
			{
				dy = p[2] - y[j];
				sinp3dy = sinp3 * dy;
				a = cosp3dx + sinp3dy;
				cosp3dy = cosp3 * dy;
				b = cosp3dy - sinp3dx;
				exparg = -a * a / twop4sq - b * b / twop5sq;
				asq = a * a;
				bsq = b * b;
				expexparg = Math::Exp(exparg);
				p0expexparg = p[0] * expexparg;
				c = p[6] - Z[i, j] + p0expexparg;
				twop0expexparg = 2 * p0expexparg;

				f += c * c;

				grad[0] += 2 * expexparg*c;
				grad[1] += -twop0expexparg * ((cosp3*a) / p4sq - (sinp3*b) / p5sq)*c;
				grad[2] += -twop0expexparg * ((cosp3*b) / p5sq + (sinp3*a) / p4sq)*c;
				grad[3] += -twop0expexparg * ((a*b) / p4sq - (a*b) / p5sq)*c;
				grad[4] += (twop0expexparg*asq * c) / p4cu;
				grad[5] += (twop0expexparg*bsq * c) / p5cu;
				grad[6] += 2 * p[6] - 2 * Z[i, j] + twop0expexparg;
			}
		}
	}
}

void JPFITS::JPMath::alglib_Gauss_2d_LM_LS_grad_compound(array<double>^ p, double %f, array<double>^ grad, Object^ obj)
{
	f = 0;
	for (int i = 0; i < grad->Length; i++)
		grad[i] = 0;
	array<double>^ g = gcnew array<double>(grad->Length);
	array<double, 2>^ Z = (array<double, 2>^)(((array<Object^>^)obj)[0]);
	array<double>^ x = (array<double>^)(((array<Object^>^)obj)[1]);
	array<double>^ y = (array<double>^)(((array<Object^>^)obj)[2]);
	int func = ((array<int>^)(((array<Object^>^)obj)[3]))[0];
	int count = ((array<int>^)(((array<Object^>^)obj)[3]))[1];
	int pindex;

	if (func == 5)
	{
		array<double>^ p3sq = gcnew array<double>(count);// = p[3] * p[3];
		array<double>^ twop3sq = gcnew array<double>(count);// = 2 * p3sq;
		array<double>^ p3cu = gcnew array<double>(count);// = p3sq * p[3];
		array<double>^ twop1 = gcnew array<double>(count); // = 2 * p[1];
		array<double>^ twop2 = gcnew array<double>(count);// = 2 * p[2];
		for (int k = 0; k < count; k++)
		{
			pindex = k * (func - 1);
			p3sq[k] = p[pindex + 3] * p[pindex + 3];
			twop3sq[k] = 2 * p3sq[k];
			p3cu[k] = p3sq[k] * p[pindex + 3];
			twop1[k] = 2 * p[pindex + 1];
			twop2[k] = 2 * p[pindex + 2];
		}

		double dy, dxsqpdysq, expres, exparg, p0expres, a;
		array<double>^ dxsq = gcnew array<double>(count);

		for (int i = 0; i < x->Length; i++)
		{
			for (int k = 0; k < count; k++)
			{
				pindex = k * (func - 1);

				dxsq[k] = p[pindex + 1] - x[i];
				dxsq[k] *= dxsq[k];
			}

			for (int j = 0; j < y->Length; j++)
			{
				p0expres = 0;
				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					dy = p[pindex + 2] - y[j];
					dxsqpdysq = dxsq[k] + dy * dy;
					exparg = -dxsqpdysq / twop3sq[k];
					expres = Math::Exp(exparg);

					p0expres += p[pindex + 0] * expres;
					g[pindex + 0] = expres;
					g[pindex + 1] = -((twop1[k] - 2 * x[i])) / p3sq[k] * p[pindex + 0] * expres;
					g[pindex + 2] = -((twop2[k] - 2 * y[j])) / p3sq[k] * p[pindex + 0] * expres;
					g[pindex + 3] = (dxsqpdysq) / p3cu[k] * p[pindex + 0] * expres;
				}

				a = p[p->Length - 1] - Z[i, j] + p0expres;

				f += a * a;
				grad[grad->Length - 1] += 2 * a;
				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					grad[pindex + 0] += 2 * g[pindex + 0] * a;
					grad[pindex + 1] += g[pindex + 1] * a;
					grad[pindex + 2] += g[pindex + 2] * a;
					grad[pindex + 3] += 2 * g[pindex + 3] * a;
				}
			}
		}
		return;
	}

	if (func == 7)
	{
		array<double>^ cosp3 = gcnew array<double>(count);// = Math::Cos(p[3]);
		array<double>^ sinp3 = gcnew array<double>(count);// = Math::Sin(p[3]);
		array<double>^ p4sq = gcnew array<double>(count);// = p[4] * p[4];
		array<double>^ twop4sq = gcnew array<double>(count);// = 2 * p4sq;
		array<double>^ p5sq = gcnew array<double>(count);// = p[5] * p[5];
		array<double>^ twop5sq = gcnew array<double>(count);// = 2 * p5sq;
		array<double>^ p4cu = gcnew array<double>(count);// = p4sq * p[4];
		array<double>^ p5cu = gcnew array<double>(count);// = p5sq * p[5];
		for (int k = 0; k < count; k++)
		{
			pindex = k * (func - 1);
			cosp3[k] = Math::Cos(p[pindex + 3]);
			sinp3[k] = Math::Sin(p[pindex + 3]);
			p4sq[k] = p[pindex + 4] * p[pindex + 4];
			twop4sq[k] = 2 * p4sq[k];
			p5sq[k] = p[pindex + 5] * p[pindex + 5];
			twop5sq[k] = 2 * p5sq[k];
			p4cu[k] = p4sq[k] * p[pindex + 4];
			p5cu[k] = p5sq[k] * p[pindex + 5];
		}

		double dx, cosp3dx, sinp3dx, dy, sinp3dy, a, cosp3dy, b, c, exparg, asq, bsq, expexparg, p0expexparg, twop0expexparg;

		for (int i = 0; i < x->Length; i++)
		{
			for (int j = 0; j < y->Length; j++)
			{
				p0expexparg = 0;
				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					dx = p[pindex + 1] - x[i];
					cosp3dx = cosp3[k] * dx;
					sinp3dx = sinp3[k] * dx;
					dy = p[pindex + 2] - y[j];
					sinp3dy = sinp3[k] * dy;
					a = cosp3dx + sinp3dy;
					cosp3dy = cosp3[k] * dy;
					b = cosp3dy - sinp3dx;
					exparg = -a * a / twop4sq[k] - b * b / twop5sq[k];
					asq = a * a;
					bsq = b * b;
					expexparg = Math::Exp(exparg);
					p0expexparg += p[pindex + 0] * expexparg;

					g[pindex + 0] = 2 * expexparg;
					g[pindex + 1] = ((cosp3[k] *a) / p4sq[k] - (sinp3[k] *b) / p5sq[k]) * 2 * p[pindex + 0] * expexparg;
					g[pindex + 2] = ((cosp3[k] *b) / p5sq[k] + (sinp3[k] *a) / p4sq[k]) * 2 * p[pindex + 0] * expexparg;
					g[pindex + 3] = ((a*b) / p4sq[k] - (a*b) / p5sq[k]) * 2 * p[pindex + 0] * expexparg;
					g[pindex + 4] = (asq) / p4cu[k] * 2 * p[pindex + 0] * expexparg;
					g[pindex + 5] = (bsq) / p5cu[k] * 2 * p[pindex + 0] * expexparg;
				}

				twop0expexparg = 2 * p0expexparg;
				c = p[p->Length - 1] - Z[i, j] + p0expexparg;

				f += c * c;
				grad[grad->Length - 1] += 2 * c;

				for (int k = 0; k < count; k++)				
				{
					pindex = k * (func - 1);

					grad[pindex + 0] += g[pindex + 0] * c;
					grad[pindex + 1] += - g[pindex + 1] * c;
					grad[pindex + 2] += - g[pindex + 2] * c;
					grad[pindex + 3] += - g[pindex + 3] * c;
					grad[pindex + 4] += (c) * g[pindex + 4];
					grad[pindex + 5] += (c) * g[pindex + 5];
				}				
			}
		}
	}
}

void JPFITS::JPMath::alglib_Gauss_2d_LM_LS_CHISQ_grad(array<double>^ p, double %f, array<double>^ grad, Object^ obj)
{
	f = 0;
	for (int i = 0; i < grad->Length; i++)
		grad[i] = 0;
	array<double, 2>^ Z = (array<double, 2>^)(((array<Object^>^)obj)[0]);
	array<double>^ x = (array<double>^)(((array<Object^>^)obj)[1]);
	array<double>^ y = (array<double>^)(((array<Object^>^)obj)[2]);
	double den;

	if (p->Length == 5)
	{
		double p3sq = p[3] * p[3];
		double p3cu = p[3] * p3sq;
		double twop3sq = 2 * p3sq;
		double twop1 = 2 * p[1];
		double twop2 = 2 * p[2];

		for (int i = 0; i < x->Length; i++)
		{
			double dx = p[1] - x[i];
			double dxsq = dx * dx;
			double twoxi = 2 * x[i];

			for (int j = 0; j < y->Length; j++)
			{
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Z[i, j];

				double dy = p[2] - y[j];
				double dysq = dy * dy;
				double exparg = (dxsq + dysq) / twop3sq;
				double expmexparg = Math::Exp(-exparg);
				double p0expmexparg = p[0] * expmexparg;
				double a = p[4] - Z[i, j] + p0expmexparg;

				f += a * a / den;

				grad[0] += (2 * expmexparg*a) / den;
				grad[1] += -(p0expmexparg*(twop1 - twoxi)*a) / (den*p3sq);
				grad[2] += -(p0expmexparg*(twop2 - 2 * y[j])*a) / (den*p3sq);
				grad[3] += (2 * p0expmexparg*(dxsq + dysq)*a) / (den*p3cu);
				grad[4] += (2 * p[4] - 2 * Z[i, j] + 2 * p0expmexparg) / den;
			}
		}
		return;
	}

	if (p->Length == 7)
	{
		double cosp3 = Math::Cos(p[3]);
		double sinp3 = Math::Sin(p[3]);
		double p4sq = p[4] * p[4];
		double p5sq = p[5] * p[5];
		double twop4sq = 2 * p4sq;
		double twop5sq = 2 * p5sq;
		double p4cu = p4sq * p[4];
		double p5cu = p5sq * p[5];
		double twop6 = 2 * p[6];
		double dx, cosp3dx, dy, sinp3dy, cosp3dxsinp3dy, cosp3dxsinp3dysq, cosp3dysinp3dx, cosp3dysinp3dxsq, exparg, expres, a, p0expres, twop0expres;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			cosp3dx = cosp3 * dx;

			for (int j = 0; j < y->Length; j++)
			{
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Z[i, j];

				dy = p[2] - y[j];
				sinp3dy = sinp3 * dy;
				cosp3dxsinp3dy = cosp3dx + sinp3dy;
				cosp3dxsinp3dysq = cosp3dxsinp3dy * cosp3dxsinp3dy;
				cosp3dysinp3dx = cosp3 * dy - sinp3 * dx;
				cosp3dysinp3dxsq = cosp3dysinp3dx * cosp3dysinp3dx;
				exparg = -cosp3dxsinp3dysq / twop4sq - cosp3dysinp3dxsq / twop5sq;
				expres = Math::Exp(exparg);
				a = p[6] - Z[i, j] + p[0] * expres;
				p0expres = p[0] * expres;
				twop0expres = 2 * p0expres;

				f += a * a / den;

				grad[0] += (2 * expres*a) / den;
				grad[1] += -(twop0expres*((cosp3*cosp3dxsinp3dy) / p4sq - (sinp3*cosp3dysinp3dx) / p5sq)*a) / den;
				grad[2] += -(twop0expres*((cosp3*cosp3dysinp3dx) / p5sq + (sinp3*cosp3dxsinp3dy) / p4sq)*a) / den;
				grad[3] += -(twop0expres*((cosp3dxsinp3dy*cosp3dysinp3dx) / p4sq - (cosp3dxsinp3dy*cosp3dysinp3dx) / p5sq)*a) / den;
				grad[4] += (twop0expres*cosp3dxsinp3dysq * a) / (den*p4cu);
				grad[5] += (twop0expres*cosp3dysinp3dxsq * a) / (den*p5cu);
				grad[6] += (twop6 - 2 * Z[i, j] + twop0expres) / den;
			}
		}
	}
}

void JPFITS::JPMath::alglib_Gauss_2d_LM_LS_CHISQ_grad_compound(array<double>^ p, double% f, array<double>^ grad, Object^ obj)
{
	f = 0;
	for (int i = 0; i < grad->Length; i++)
		grad[i] = 0;
	array<double>^ g = gcnew array<double>(grad->Length);
	array<double, 2>^ Z = (array<double, 2>^)(((array<Object^>^)obj)[0]);
	array<double>^ x = (array<double>^)(((array<Object^>^)obj)[1]);
	array<double>^ y = (array<double>^)(((array<Object^>^)obj)[2]);
	int func = ((array<int>^)(((array<Object^>^)obj)[3]))[0];
	int count = ((array<int>^)(((array<Object^>^)obj)[3]))[1];
	int pindex;

	if (func == 5)
	{
		array<double>^ p3sq = gcnew array<double>(count);// = p[3] * p[3];
		array<double>^ twop3sq = gcnew array<double>(count);// = 2 * p3sq;
		array<double>^ p3cu = gcnew array<double>(count);// = p3sq * p[3];
		array<double>^ twop1 = gcnew array<double>(count); // = 2 * p[1];
		array<double>^ twop2 = gcnew array<double>(count);// = 2 * p[2];
		for (int k = 0; k < count; k++)
		{
			pindex = k * (func - 1);
			p3sq[k] = p[pindex + 3] * p[pindex + 3];
			twop3sq[k] = 2 * p3sq[k];
			p3cu[k] = p3sq[k] * p[pindex + 3];
			twop1[k] = 2 * p[pindex + 1];
			twop2[k] = 2 * p[pindex + 2];
		}

		double dy, dxsqpdysq, expres, exparg, p0expres, a, den;
		array<double>^ dxsq = gcnew array<double>(count);

		for (int i = 0; i < x->Length; i++)
		{
			for (int k = 0; k < count; k++)
			{
				pindex = k * (func - 1);

				dxsq[k] = p[pindex + 1] - x[i];
				dxsq[k] *= dxsq[k];
			}

			for (int j = 0; j < y->Length; j++)
			{
				p0expres = 0;
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Z[i, j];

				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					dy = p[pindex + 2] - y[j];
					dxsqpdysq = dxsq[k] + dy * dy;
					exparg = -dxsqpdysq / twop3sq[k];
					expres = Math::Exp(exparg);

					p0expres += p[pindex + 0] * expres;
					g[pindex + 0] = expres;
					g[pindex + 1] = -((twop1[k] - 2 * x[i])) / p3sq[k] * p[pindex + 0] * expres;
					g[pindex + 2] = -((twop2[k] - 2 * y[j])) / p3sq[k] * p[pindex + 0] * expres;
					g[pindex + 3] = (dxsqpdysq) / p3cu[k] * p[pindex + 0] * expres;
				}

				a = p[p->Length - 1] - Z[i, j] + p0expres;

				f += a * a / den;
				grad[grad->Length - 1] += 2 * a / den;
				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					grad[pindex + 0] += 2 * g[pindex + 0] * a / den;
					grad[pindex + 1] += g[pindex + 1] * a / den;
					grad[pindex + 2] += g[pindex + 2] * a / den;
					grad[pindex + 3] += 2 * g[pindex + 3] * a / den;
				}
			}
		}
		return;
	}

	if (func == 7)
	{
		array<double>^ cosp3 = gcnew array<double>(count);// = Math::Cos(p[3]);
		array<double>^ sinp3 = gcnew array<double>(count);// = Math::Sin(p[3]);
		array<double>^ p4sq = gcnew array<double>(count);// = p[4] * p[4];
		array<double>^ twop4sq = gcnew array<double>(count);// = 2 * p4sq;
		array<double>^ p5sq = gcnew array<double>(count);// = p[5] * p[5];
		array<double>^ twop5sq = gcnew array<double>(count);// = 2 * p5sq;
		array<double>^ p4cu = gcnew array<double>(count);// = p4sq * p[4];
		array<double>^ p5cu = gcnew array<double>(count);// = p5sq * p[5];
		for (int k = 0; k < count; k++)
		{
			pindex = k * (func - 1);
			cosp3[k] = Math::Cos(p[pindex + 3]);
			sinp3[k] = Math::Sin(p[pindex + 3]);
			p4sq[k] = p[pindex + 4] * p[pindex + 4];
			twop4sq[k] = 2 * p4sq[k];
			p5sq[k] = p[pindex + 5] * p[pindex + 5];
			twop5sq[k] = 2 * p5sq[k];
			p4cu[k] = p4sq[k] * p[pindex + 4];
			p5cu[k] = p5sq[k] * p[pindex + 5];
		}

		double dx, cosp3dx, sinp3dx, dy, sinp3dy, a, cosp3dy, b, c, exparg, asq, bsq, expexparg, p0expexparg, twop0expexparg, den;

		for (int i = 0; i < x->Length; i++)
		{
			for (int j = 0; j < y->Length; j++)
			{
				p0expexparg = 0;
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Z[i, j];

				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					dx = p[pindex + 1] - x[i];
					cosp3dx = cosp3[k] * dx;
					sinp3dx = sinp3[k] * dx;
					dy = p[pindex + 2] - y[j];
					sinp3dy = sinp3[k] * dy;
					a = cosp3dx + sinp3dy;
					cosp3dy = cosp3[k] * dy;
					b = cosp3dy - sinp3dx;
					exparg = -a * a / twop4sq[k] - b * b / twop5sq[k];
					asq = a * a;
					bsq = b * b;
					expexparg = Math::Exp(exparg);
					p0expexparg += p[pindex + 0] * expexparg;

					g[pindex + 0] = 2 * expexparg;
					g[pindex + 1] = ((cosp3[k] * a) / p4sq[k] - (sinp3[k] * b) / p5sq[k]) * 2 * p[pindex + 0] * expexparg;
					g[pindex + 2] = ((cosp3[k] * b) / p5sq[k] + (sinp3[k] * a) / p4sq[k]) * 2 * p[pindex + 0] * expexparg;
					g[pindex + 3] = ((a * b) / p4sq[k] - (a * b) / p5sq[k]) * 2 * p[pindex + 0] * expexparg;
					g[pindex + 4] = (asq) / p4cu[k] * 2 * p[pindex + 0] * expexparg;
					g[pindex + 5] = (bsq) / p5cu[k] * 2 * p[pindex + 0] * expexparg;
				}

				twop0expexparg = 2 * p0expexparg;
				c = p[p->Length - 1] - Z[i, j] + p0expexparg;

				f += c * c / den;
				grad[grad->Length - 1] += 2 * c / den;

				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					grad[pindex + 0] += g[pindex + 0] * c / den;
					grad[pindex + 1] += -g[pindex + 1] * c / den;
					grad[pindex + 2] += -g[pindex + 2] * c / den;
					grad[pindex + 3] += -g[pindex + 3] * c / den;
					grad[pindex + 4] += (c)*g[pindex + 4] / den;
					grad[pindex + 5] += (c)*g[pindex + 5] / den;
				}
			}
		}
	}
}

void JPFITS::JPMath::alglib_Gauss_2d_LM_LS_ROBUST_grad(array<double>^ p, double %f, array<double>^ grad, Object^ obj)
{
	f = 0;
	for (int i = 0; i < grad->Length; i++)
		grad[i] = 0;
	array<double, 2>^ Z = (array<double, 2>^)(((array<Object^>^)obj)[0]);
	array<double>^ x = (array<double>^)(((array<Object^>^)obj)[1]);
	array<double>^ y = (array<double>^)(((array<Object^>^)obj)[2]);
	double den;

	if (p->Length == 5)
	{
		double p3sq = p[3] * p[3];
		double twop3sq = 2 * p3sq;
		double p3cu = p3sq * p[3];
		double twop1 = 2 * p[1];
		double twop2 = 2 * p[2];
		double dx, dxsq, dy, dysq, dxsqdysq, exparg, expres, p0expres, a;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			dxsq = dx * dx;

			for (int j = 0; j < y->Length; j++)
			{
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Math::Sqrt(Z[i, j]);

				dy = p[2] - y[j];
				dysq = dy * dy;
				dxsqdysq = dxsq + dysq;
				exparg = -dxsqdysq / twop3sq;
				expres = Math::Exp(exparg);
				p0expres = p[0] * expres;
				a = p[4] - Z[i, j] + p0expres;

				f += a * a / den;

				grad[0] += (2 * expres*a) / den;
				grad[1] += -(p0expres*(twop1 - 2 * x[i])*a) / (den*p3sq);
				grad[2] += -(p0expres*(twop2 - 2 * y[j])*a) / (den*p3sq);
				grad[3] += (2 * p0expres*dxsqdysq*a) / (den*p3cu);
				grad[4] += (2 * p[4] - 2 * Z[i, j] + 2 * p0expres) / den;
			}
		}
		return;
	}

	if (p->Length == 7)
	{
		double cosp3 = Math::Cos(p[3]);
		double sinp3 = Math::Sin(p[3]);
		double p4sq = p[4] * p[4];
		double twop4sq = 2 * p4sq;
		double p5sq = p[5] * p[5];
		double twop5sq = 2 * p5sq;
		double p4cu = p[4] * p4sq;
		double p5cu = p[5] * p5sq;
		double dx, cosp3dx, dy, sinp3dy, a, b, asq, bsq, exparg, expres, c, p0expres, twop0expres;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			cosp3dx = cosp3 * dx;

			for (int j = 0; j < y->Length; j++)
			{
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Math::Sqrt(Z[i, j]);

				dy = p[2] - y[j];
				sinp3dy = sinp3 * dy;
				a = cosp3dx + sinp3dy;
				b = cosp3 * dy - sinp3 * dx;
				asq = a * a;
				bsq = b * b;
				exparg = -asq / twop4sq - bsq / twop5sq;
				expres = Math::Exp(exparg);
				c = p[6] - Z[i, j] + p[0] * expres;
				p0expres = p[0] * expres;
				twop0expres = 2 * p0expres;

				f += c * c / den;

				grad[0] += (2 * expres*c) / den;
				grad[1] += -(twop0expres*((cosp3*a) / p4sq - (sinp3*b) / p5sq)*c) / den;
				grad[2] += -(twop0expres*((cosp3*b) / p5sq + (sinp3*a) / p4sq)*c) / den;
				grad[3] += -(twop0expres*((a*b) / p4sq - (a*b) / p5sq)*c) / den;
				grad[4] += (twop0expres*asq * c) / (den*p4cu);
				grad[5] += (twop0expres*bsq * c) / (den*p5cu);
				grad[6] += (2 * p[6] - 2 * Z[i, j] + twop0expres) / den;
			}
		}
	}
}

void JPFITS::JPMath::alglib_Gauss_2d_LM_LS_ROBUST_grad_compound(array<double>^ p, double% f, array<double>^ grad, Object^ obj)
{
	f = 0;
	for (int i = 0; i < grad->Length; i++)
		grad[i] = 0;
	array<double>^ g = gcnew array<double>(grad->Length);
	array<double, 2>^ Z = (array<double, 2>^)(((array<Object^>^)obj)[0]);
	array<double>^ x = (array<double>^)(((array<Object^>^)obj)[1]);
	array<double>^ y = (array<double>^)(((array<Object^>^)obj)[2]);
	int func = ((array<int>^)(((array<Object^>^)obj)[3]))[0];
	int count = ((array<int>^)(((array<Object^>^)obj)[3]))[1];
	int pindex;

	if (func == 5)
	{
		array<double>^ p3sq = gcnew array<double>(count);// = p[3] * p[3];
		array<double>^ twop3sq = gcnew array<double>(count);// = 2 * p3sq;
		array<double>^ p3cu = gcnew array<double>(count);// = p3sq * p[3];
		array<double>^ twop1 = gcnew array<double>(count); // = 2 * p[1];
		array<double>^ twop2 = gcnew array<double>(count);// = 2 * p[2];
		for (int k = 0; k < count; k++)
		{
			pindex = k * (func - 1);
			p3sq[k] = p[pindex + 3] * p[pindex + 3];
			twop3sq[k] = 2 * p3sq[k];
			p3cu[k] = p3sq[k] * p[pindex + 3];
			twop1[k] = 2 * p[pindex + 1];
			twop2[k] = 2 * p[pindex + 2];
		}

		double dy, dxsqpdysq, expres, exparg, p0expres, a, den;
		array<double>^ dxsq = gcnew array<double>(count);

		for (int i = 0; i < x->Length; i++)
		{
			for (int k = 0; k < count; k++)
			{
				pindex = k * (func - 1);

				dxsq[k] = p[pindex + 1] - x[i];
				dxsq[k] *= dxsq[k];
			}

			for (int j = 0; j < y->Length; j++)
			{				
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Math::Sqrt(Z[i, j]);

				p0expres = 0;
				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					dy = p[pindex + 2] - y[j];
					dxsqpdysq = dxsq[k] + dy * dy;
					exparg = -dxsqpdysq / twop3sq[k];
					expres = Math::Exp(exparg);

					p0expres += p[pindex + 0] * expres;
					g[pindex + 0] = expres;
					g[pindex + 1] = -((twop1[k] - 2 * x[i])) / p3sq[k] * p[pindex + 0] * expres;
					g[pindex + 2] = -((twop2[k] - 2 * y[j])) / p3sq[k] * p[pindex + 0] * expres;
					g[pindex + 3] = (dxsqpdysq) / p3cu[k] * p[pindex + 0] * expres;
				}

				a = p[p->Length - 1] - Z[i, j] + p0expres;

				f += a * a / den;
				grad[grad->Length - 1] += 2 * a / den;
				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					grad[pindex + 0] += 2 * g[pindex + 0] * a / den;
					grad[pindex + 1] += g[pindex + 1] * a / den;
					grad[pindex + 2] += g[pindex + 2] * a / den;
					grad[pindex + 3] += 2 * g[pindex + 3] * a / den;
				}
			}
		}
		return;
	}

	if (func == 7)
	{
		array<double>^ cosp3 = gcnew array<double>(count);// = Math::Cos(p[3]);
		array<double>^ sinp3 = gcnew array<double>(count);// = Math::Sin(p[3]);
		array<double>^ p4sq = gcnew array<double>(count);// = p[4] * p[4];
		array<double>^ twop4sq = gcnew array<double>(count);// = 2 * p4sq;
		array<double>^ p5sq = gcnew array<double>(count);// = p[5] * p[5];
		array<double>^ twop5sq = gcnew array<double>(count);// = 2 * p5sq;
		array<double>^ p4cu = gcnew array<double>(count);// = p4sq * p[4];
		array<double>^ p5cu = gcnew array<double>(count);// = p5sq * p[5];
		for (int k = 0; k < count; k++)
		{
			pindex = k * (func - 1);
			cosp3[k] = Math::Cos(p[pindex + 3]);
			sinp3[k] = Math::Sin(p[pindex + 3]);
			p4sq[k] = p[pindex + 4] * p[pindex + 4];
			twop4sq[k] = 2 * p4sq[k];
			p5sq[k] = p[pindex + 5] * p[pindex + 5];
			twop5sq[k] = 2 * p5sq[k];
			p4cu[k] = p4sq[k] * p[pindex + 4];
			p5cu[k] = p5sq[k] * p[pindex + 5];
		}

		double dx, cosp3dx, sinp3dx, dy, sinp3dy, a, cosp3dy, b, c, exparg, asq, bsq, expexparg, p0expexparg, twop0expexparg, den;

		for (int i = 0; i < x->Length; i++)
		{
			for (int j = 0; j < y->Length; j++)
			{
				p0expexparg = 0;
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Math::Sqrt(Z[i, j]);

				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					dx = p[pindex + 1] - x[i];
					cosp3dx = cosp3[k] * dx;
					sinp3dx = sinp3[k] * dx;
					dy = p[pindex + 2] - y[j];
					sinp3dy = sinp3[k] * dy;
					a = cosp3dx + sinp3dy;
					cosp3dy = cosp3[k] * dy;
					b = cosp3dy - sinp3dx;
					exparg = -a * a / twop4sq[k] - b * b / twop5sq[k];
					asq = a * a;
					bsq = b * b;
					expexparg = Math::Exp(exparg);
					p0expexparg += p[pindex + 0] * expexparg;

					g[pindex + 0] = 2 * expexparg;
					g[pindex + 1] = ((cosp3[k] * a) / p4sq[k] - (sinp3[k] * b) / p5sq[k]) * 2 * p[pindex + 0] * expexparg;
					g[pindex + 2] = ((cosp3[k] * b) / p5sq[k] + (sinp3[k] * a) / p4sq[k]) * 2 * p[pindex + 0] * expexparg;
					g[pindex + 3] = ((a * b) / p4sq[k] - (a * b) / p5sq[k]) * 2 * p[pindex + 0] * expexparg;
					g[pindex + 4] = (asq) / p4cu[k] * 2 * p[pindex + 0] * expexparg;
					g[pindex + 5] = (bsq) / p5cu[k] * 2 * p[pindex + 0] * expexparg;
				}

				twop0expexparg = 2 * p0expexparg;
				c = p[p->Length - 1] - Z[i, j] + p0expexparg;

				f += c * c / den;
				grad[grad->Length - 1] += 2 * c / den;

				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					grad[pindex + 0] += g[pindex + 0] * c / den;
					grad[pindex + 1] += -g[pindex + 1] * c / den;
					grad[pindex + 2] += -g[pindex + 2] * c / den;
					grad[pindex + 3] += -g[pindex + 3] * c / den;
					grad[pindex + 4] += (c)*g[pindex + 4] / den;
					grad[pindex + 5] += (c)*g[pindex + 5] / den;
				}
			}
		}
	}
}

void JPFITS::JPMath::alglib_Gauss_2d_LM_LS_CSTAT_grad(array<double>^ p, double %f, array<double>^ grad, Object^ obj)
{
	f = 0;
	for (int i = 0; i < grad->Length; i++)
		grad[i] = 0;
	array<double, 2>^ Z = (array<double, 2>^)(((array<Object^>^)obj)[0]);
	array<double>^ x = (array<double>^)(((array<Object^>^)obj)[1]);
	array<double>^ y = (array<double>^)(((array<Object^>^)obj)[2]);

	if (p->Length == 5)
	{
		double p3sq = p[3] * p[3];
		double twop3sq = 2 * p3sq;
		double p3cu = p3sq * p[3];
		double twop1 = 2 * p[1];
		double twop2 = 2 * p[2];
		double dx, dxsq, twoxi, dy, dysq, dxsqdysq, exparg, expres, p0expres, a;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			dxsq = dx * dx;
			twoxi = 2 * x[i];

			for (int j = 0; j < y->Length; j++)
			{
				dy = p[2] - y[j];
				dysq = dy * dy;
				dxsqdysq = dxsq + dysq;
				exparg = -dxsqdysq / twop3sq;
				expres = Math::Exp(exparg);
				p0expres = p[0] * expres;
				a = p[4] + p0expres;

				f += 2 * p[4] - 2 * Z[i, j] * Math::Log(a) + 2 * p0expres;

				grad[0] += 2 * expres - (2 * Z[i, j] * expres) / (a);
				grad[1] += (Z[i, j] * p0expres*(twop1 - twoxi)) / (p3sq * (a)) - (p0expres*(twop1 - twoxi)) / p3sq;
				grad[2] += (Z[i, j] * p0expres*(twop2 - 2 * y[j])) / (p3sq * (a)) - (p0expres*(twop2 - 2 * y[j])) / p3sq;
				grad[3] += (2 * p0expres*dxsqdysq) / p3cu - (2 * Z[i, j] * p0expres*dxsqdysq) / (p3cu * (a));
				grad[4] += 2 - (2 * Z[i, j]) / (a);
			}
		}
		return;
	}

	if (p->Length == 7)
	{
		double p4sq = p[4] * p[4];
		double p5sq = p[5] * p[5];
		double cosp3 = Math::Cos(p[3]);
		double sinp3 = Math::Sin(p[3]);
		double twop4sq = 2 * p4sq;
		double twop5sq = 2 * p5sq;
		double p4cu = p4sq * p[4];
		double p5cu = p5sq * p[5]; 
		double dy, dx, cosp3dx, sinp3dx, sinp3dy, a, b, c, d, g, ab, h, cosp3dy, asq, bsq, exparg, expres, p0expres, twop0expres, twoZijp0expres, twoZij, cosp3a, cosp3b, sinp3b, sinp3a;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			cosp3dx = cosp3 * dx;
			sinp3dx = sinp3 * dx;			

			for (int j = 0; j < y->Length; j++)
			{
				dy = p[2] - y[j];
				sinp3dy = sinp3 * dy;
				a = cosp3dx + sinp3dy;
				cosp3dy = cosp3 * dy;
				b = cosp3dy - sinp3dx;
				asq = a * a;
				bsq = b * b;
				exparg = -asq / twop4sq - bsq / twop5sq;
				expres = Math::Exp(exparg);
				c = p[6] + p[0] * expres;
				p0expres = p[0] * expres;
				twoZijp0expres = 2 * Z[i, j] * p0expres;
				twop0expres = 2 * p0expres;
				twoZij = 2 * Z[i, j];
				cosp3a = cosp3 * a;
				cosp3b = cosp3 * b;
				sinp3b = sinp3 * b;
				sinp3a = sinp3 * a;
				d = cosp3a / p4sq - sinp3b / p5sq;
				g = cosp3b / p5sq + sinp3a / p4sq;
				ab = a * b;
				h = ab / p4sq - ab / p5sq;

				f += 2 * p[6] - twoZij * Math::Log(c) + twop0expres;

				grad[0] += 2 * expres - (twoZij*expres) / (c);
				grad[1] += (twoZijp0expres*d) / (c)-twop0expres * d;
				grad[2] += (twoZijp0expres*g) / (c)-twop0expres * g;
				grad[3] += (twoZijp0expres*h) / (c)-twop0expres * h;
				grad[4] += (twop0expres*asq) / p4cu - (twoZijp0expres*asq) / (p4cu * (c));
				grad[5] += (twop0expres*bsq) / p5cu - (twoZijp0expres*bsq) / (p5cu * (c));
				grad[6] += 2 - (twoZij) / (c);
			}
		}
	}
}

void JPFITS::JPMath::alglib_Gauss_2d_LM_LS_CSTAT_grad_compound(array<double>^ p, double% f, array<double>^ grad, Object^ obj)
{
	f = 0;
	for (int i = 0; i < grad->Length; i++)
		grad[i] = 0;
	array<double, 2>^ Z = (array<double, 2>^)(((array<Object^>^)obj)[0]);
	array<double>^ x = (array<double>^)(((array<Object^>^)obj)[1]);
	array<double>^ y = (array<double>^)(((array<Object^>^)obj)[2]);
	int func = ((array<int>^)(((array<Object^>^)obj)[3]))[0];
	int count = ((array<int>^)(((array<Object^>^)obj)[3]))[1];
	int pindex;

	if (func == 5)
	{
		array<double>^ p3sq = gcnew array<double>(count);// = p[3] * p[3];
		array<double>^ twop3sq = gcnew array<double>(count);// = 2 * p3sq;
		array<double>^ p3cu = gcnew array<double>(count);// = p3sq * p[3];
		array<double>^ twop1 = gcnew array<double>(count); // = 2 * p[1];
		array<double>^ twop2 = gcnew array<double>(count);// = 2 * p[2];
		for (int k = 0; k < count; k++)
		{
			pindex = k * (func - 1);
			p3sq[k] = p[pindex + 3] * p[pindex + 3];
			twop3sq[k] = 2 * p3sq[k];
			p3cu[k] = p3sq[k] * p[pindex + 3];
			twop1[k] = 2 * p[pindex + 1];
			twop2[k] = 2 * p[pindex + 2];
		}

		double twoxi, exparg, p0expres, a, sumexpress;
		array<double>^ dxsq = gcnew array<double>(count);
		array<double>^ dysq = gcnew array<double>(count);
		array<double>^ dxsqdysq = gcnew array<double>(count);
		array<double>^ expres = gcnew array<double>(count);

		for (int i = 0; i < x->Length; i++)
		{
			for (int k = 0; k < count; k++)
			{
				pindex = k * (func - 1);

				dxsq[k] = p[pindex + 1] - x[i];
				dxsq[k] *= dxsq[k];
			}
			twoxi = 2 * x[i];

			for (int j = 0; j < y->Length; j++)
			{
				p0expres = 0;
				sumexpress = 0;
				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					dysq[k] = p[pindex + 2] - y[j];
					dysq[k] *= dysq[k];
					dxsqdysq[k] = dxsq[k] + dysq[k];
					exparg = -dxsqdysq[k] / twop3sq[k];
					expres[k] = Math::Exp(exparg);
					
					sumexpress += expres[k];
					p0expres += p[pindex + 0] * expres[k];
				}
				a = p[p->Length - 1] + p0expres;
				f += 2 * a - 2 * Z[i, j] * Math::Log(a);
				grad[grad->Length - 1] += 2 - (2 * Z[i, j]) / (a);

				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					grad[pindex + 0] += 2 * expres[k] - (2 * Z[i, j] * expres[k]) / (a);
					grad[pindex + 1] += (Z[i, j] * p[pindex + 0] * expres[k] * (twop1[k] - twoxi)) / (p3sq[k] * (a)) - (p[pindex + 0] * expres[k] * (twop1[k] - twoxi)) / p3sq[k];
					grad[pindex + 2] += (Z[i, j] * p[pindex + 0] * expres[k] * (twop2[k] - 2 * y[j])) / (p3sq[k] * (a)) - (p[pindex + 0] * expres[k] * (twop2[k] - 2 * y[j])) / p3sq[k];
					grad[pindex + 3] += (2 * p[pindex + 0] * expres[k] * dxsqdysq[k]) / p3cu[k] - (2 * Z[i, j] * p[pindex + 0] * expres[k] * dxsqdysq[k]) / (p3cu[k] * (a));
				}
			}
		}
		return;
	}

	if (func == 7)
	{
		array<double>^ p4sq = gcnew array<double>(count);// p[4] * p[4];
		array<double>^ p5sq = gcnew array<double>(count);// p[5] * p[5];
		array<double>^ cosp3 = gcnew array<double>(count);//  Math::Cos(p[3]);
		array<double>^ sinp3 = gcnew array<double>(count);// Math::Sin(p[3]);
		array<double>^ twop4sq = gcnew array<double>(count);// 2 * p4sq;
		array<double>^ twop5sq = gcnew array<double>(count);// 2 * p5sq;
		array<double>^ p4cu = gcnew array<double>(count);//  p4sq * p[4];
		array<double>^ p5cu = gcnew array<double>(count);// p5sq * p[5];
		for (int k = 0; k < count; k++)
		{
			pindex = k * (func - 1);

			p4sq[k] = p[pindex + 4] * p[pindex + 4];
			p5sq[k] = p[pindex + 5] * p[pindex + 5];
			cosp3[k] = Math::Cos(p[pindex + 3]);
			sinp3[k] = Math::Sin(p[pindex + 3]);
			twop4sq[k] = 2 * p4sq[k];
			twop5sq[k] = 2 * p5sq[k];
			p4cu[k] = p4sq[k] * p[pindex + 4];
			p5cu[k] = p5sq[k] * p[pindex + 5];
		}

		double dy, sinp3dy, a, b, c, ab, cosp3dy, exparg, p0expres, twop0expres, twoZijp0expres, twoZij, cosp3a, cosp3b, sinp3b, sinp3a;
		array<double>^ dx = gcnew array<double>(count);
		array<double>^ cosp3dx = gcnew array<double>(count);
		array<double>^ sinp3dx = gcnew array<double>(count);
		array<double>^ expres = gcnew array<double>(count);
		array<double>^ d = gcnew array<double>(count);
		array<double>^ g = gcnew array<double>(count);
		array<double>^ h = gcnew array<double>(count);
		array<double>^ asq = gcnew array<double>(count);
		array<double>^ bsq = gcnew array<double>(count);

		for (int i = 0; i < x->Length; i++)
		{
			for (int k = 0; k < count; k++)
			{
				pindex = k * (func - 1);

				dx[k] = p[pindex + 1] - x[i];
				cosp3dx[k] = cosp3[k] * dx[k];
				sinp3dx[k] = sinp3[k] * dx[k];
			}

			for (int j = 0; j < y->Length; j++)
			{
				twoZij = 2 * Z[i, j];
				p0expres = 0;

				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					dy = p[pindex + 2] - y[j];
					sinp3dy = sinp3[k] * dy;
					a = cosp3dx[k] + sinp3dy;
					cosp3dy = cosp3[k] * dy;
					b = cosp3dy - sinp3dx[k];
					asq[k] = a * a;
					bsq[k] = b * b;
					exparg = -asq[k] / twop4sq[k] - bsq[k] / twop5sq[k];
					expres[k] = Math::Exp(exparg);
					p0expres += p[pindex + 0] * expres[k];
					cosp3a = cosp3[k] * a;
					cosp3b = cosp3[k] * b;
					sinp3b = sinp3[k] * b;
					sinp3a = sinp3[k] * a;
					d[k] = cosp3a / p4sq[k] - sinp3b / p5sq[k];
					g[k] = cosp3b / p5sq[k] + sinp3a / p4sq[k];
					ab = a * b;
					h[k] = ab / p4sq[k] - ab / p5sq[k];
				}
				twop0expres = 2 * p0expres;
				twoZijp0expres = twoZij * p0expres;
				c = p[p->Length - 1] + p0expres;

				f += 2 * c - twoZij * Math::Log(c);
				grad[grad->Length - 1] += 2 - (twoZij) / (c);

				for (int k = 0; k < count; k++)
				{
					pindex = k * (func - 1);

					grad[pindex + 0] += 2 * expres[k] - (twoZij * expres[k]) / (c);
					grad[pindex + 1] += (twoZij * p[pindex + 0] * expres[k] * d[k]) / (c)-2 * p[pindex + 0] * expres[k] * d[k];
					grad[pindex + 2] += (twoZij * p[pindex + 0] * expres[k] * g[k]) / (c)-2 * p[pindex + 0] * expres[k] * g[k];
					grad[pindex + 3] += (twoZij * p[pindex + 0] * expres[k] * h[k]) / (c)-2 * p[pindex + 0] * expres[k] * h[k];
					grad[pindex + 4] += (2 * p[pindex + 0] * expres[k] * asq[k]) / p4cu[k] - (twoZij * p[pindex + 0] * expres[k] * asq[k]) / (p4cu[k] * (c));
					grad[pindex + 5] += (2 * p[pindex + 0] * expres[k] * bsq[k]) / p5cu[k] - (twoZij * p[pindex + 0] * expres[k] * bsq[k]) / (p5cu[k] * (c));
				}
			}
		}
	}
}

void JPFITS::JPMath::alglib_Moffat_2d_LM_LS_grad(array<double>^ p, double %f, array<double>^ grad, Object^ obj)
{
	f = 0;
	for (int i = 0; i < grad->Length; i++)
		grad[i] = 0;
	array<double, 2>^ Z = (array<double, 2>^)(((array<Object^>^)obj)[0]);
	array<double>^ x = (array<double>^)(((array<Object^>^)obj)[1]);
	array<double>^ y = (array<double>^)(((array<Object^>^)obj)[2]);

	if (p->Length == 6)
	{
		double p3sq = p[3] * p[3];
		double p41 = p[4] + 1;
		double twop1 = 2 * p[1];
		double twop2 = 2 * p[2];
		double p0p4 = p[0] * p[4];
		double p3cu = p3sq * p[3];
		double dx, dxsq, twoxi, twop1mtwoxi, dy, dysq, dxsqdysq, a, loga, apowp4, apowp41, b;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			dxsq = dx * dx;
			twoxi = 2 * x[i];
			twop1mtwoxi = twop1 - twoxi;

			for (int j = 0; j < y->Length; j++)
			{
				dy = p[2] - y[j];
				dysq = dy * dy;
				dxsqdysq = dxsq + dysq;
				a = dxsqdysq / p3sq + 1;
				loga = Math::Log(a);
				apowp4 = Math::Pow(a, p[4]);
				apowp41 = Math::Pow(a, p41);
				b = p[5] - Z[i, j] + p[0] / apowp4;

				f += b * b;

				grad[0] += (2 * b) / apowp4;
				grad[1] += -(2 * p0p4*twop1mtwoxi*b) / (p3sq * apowp41);
				grad[2] += -(2 * p0p4*(twop2 - 2 * y[j])*b) / (p3sq * apowp41);
				grad[3] += (4 * p0p4*dxsqdysq*b) / (p3cu * apowp41);
				grad[4] += -(2 * p[0]*loga*b) / apowp4;
				grad[5] += 2 * p[5] - 2 * Z[i, j] + (2 * p[0]) / apowp4;
			}
		}
		return;
	}

	if (p->Length == 8)
	{
		double cosp3 = Math::Cos(p[3]);
		double sinp3 = Math::Sin(p[3]);
		double p4sq = p[4] * p[4];
		double p5sq = p[5] * p[5];
		double p6p1 = p[6] + 1;
		double p0p6 = p[0] * p[6];
		double p4cu = p[4] * p4sq;
		double p5cu = p[5] * p5sq;
		double twop0p6 = 2 * p0p6;
		double fourp0p6 = 4 * p0p6;
		double twocosp3 = 2 * cosp3;
		double twosinp3 = 2 * sinp3;
		double twop7 = 2 * p[7];
		double twop0 = 2 * p[0];
		double dx, cosp3dx, sinp3dx, dy, sinp3dy, a, b, c, d, g, h, asq, cosp3dy, bsq, logc, p0logc;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			cosp3dx = cosp3 * dx;
			sinp3dx = sinp3 * dx;

			for (int j = 0; j < y->Length; j++)
			{
				dy = p[2] - y[j];
				sinp3dy = sinp3 * dy;
				a = cosp3dx + sinp3dy;
				asq = a * a;
				cosp3dy = cosp3 * dy;
				b = cosp3dy - sinp3dx;
				bsq = b * b;
				c = asq / p4sq + bsq / p5sq + 1;
				logc = Math::Log(c);
				d = Math::Pow(c, p6p1);
				g = Math::Pow(c, p[6]);
				h = p[7] - Z[i, j] + p[0] / g;
				p0logc = p[0] * logc;

				f += h * h;

				grad[0] += (2 * (h)) / g;
				grad[1] += -(twop0p6*((twocosp3*a) / p4sq - (twosinp3*b) / p5sq)*(h)) / d;
				grad[2] += -(twop0p6*((twocosp3*b) / p5sq + (twosinp3*a) / p4sq)*(h)) / d;
				grad[3] += -(twop0p6*((2 * a*b) / p4sq - (2 * a*b) / p5sq)*(h)) / d;
				grad[4] += (fourp0p6*asq * (h)) / (p4cu * d);
				grad[5] += (fourp0p6*bsq * (h)) / (p5cu * d);
				grad[6] += -(2 * p0logc*(h)) / g;
				grad[7] += twop7 - 2 * Z[i, j] + twop0 / g;
			}
		}
	}
}

void JPFITS::JPMath::alglib_Moffat_2d_LM_LS_CHISQ_grad(array<double>^ p, double %f, array<double>^ grad, Object^ obj)
{
	f = 0;
	for (int i = 0; i < grad->Length; i++)
		grad[i] = 0;
	array<double, 2>^ Z = (array<double, 2>^)(((array<Object^>^)obj)[0]);
	array<double>^ x = (array<double>^)(((array<Object^>^)obj)[1]);
	array<double>^ y = (array<double>^)(((array<Object^>^)obj)[2]);
	double den;

	if (p->Length == 6)
	{
		double p3sq = p[3] * p[3];
		double p41 = p[4] + 1;
		double twop1 = 2 * p[1];
		double twop2 = 2 * p[2];
		double p0p4 = p[0] * p[4];
		double p3cu = p3sq * p[3];
		double dx, dxsq, twoxi, twop1mtwoxi, dy, dysq, dxsqdysq, a, loga, apowp4, apowp41, b;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			dxsq = dx * dx;
			twoxi = 2 * x[i];
			twop1mtwoxi = twop1 - twoxi;

			for (int j = 0; j < y->Length; j++)
			{
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Z[i, j];

				dy = p[2] - y[j];
				dysq = dy * dy;
				dxsqdysq = dxsq + dysq;
				a = dxsqdysq / p3sq + 1;
				loga = Math::Log(a);
				apowp4 = Math::Pow(a, p[4]);
				apowp41 = Math::Pow(a, p41);
				b = p[5] - Z[i, j] + p[0] / apowp4;

				f += b * b / den;

				grad[0] += (2 * b) / (den*apowp4);
				grad[1] += -(2 * p0p4*twop1mtwoxi*b) / (den*p3sq * apowp41);
				grad[2] += -(2 * p0p4*(twop2 - 2 * y[j])*b) / (den*p3sq * apowp41);
				grad[3] += (4 * p0p4*dxsqdysq*b) / (den*p3cu * apowp41);
				grad[4] += -(2 * p[0]*loga*b) / (den*apowp4);
				grad[5] += (2 * p[5] - 2 * Z[i, j] + (2 * p[0]) / apowp4) / den;
			}
		}
		return;
	}

	if (p->Length == 8)
	{
		double cosp3 = Math::Cos(p[3]);
		double sinp3 = Math::Sin(p[3]);
		double p4sq = p[4] * p[4];
		double p5sq = p[5] * p[5];
		double p6p1 = p[6] + 1;
		double p0p6 = p[0] * p[6];
		double p4cu = p[4] * p4sq;
		double p5cu = p[5] * p5sq;
		double twop0p6 = 2 * p0p6;
		double fourp0p6 = 4 * p0p6;
		double twocosp3 = 2 * cosp3;
		double twosinp3 = 2 * sinp3;
		double twop7 = 2 * p[7];
		double twop0 = 2 * p[0];
		double dx, cosp3dx, sinp3dx, dy, sinp3dy, a, b, c, d, g, h, asq, cosp3dy, bsq, logc, p0logc;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			cosp3dx = cosp3 * dx;
			sinp3dx = sinp3 * dx;

			for (int j = 0; j < y->Length; j++)
			{
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Z[i, j];

				dy = p[2] - y[j];
				sinp3dy = sinp3 * dy;
				a = cosp3dx + sinp3dy;
				asq = a * a;
				cosp3dy = cosp3 * dy;
				b = cosp3dy - sinp3dx;
				bsq = b * b;
				c = asq / p4sq + bsq / p5sq + 1;
				logc = Math::Log(c);
				d = Math::Pow(c, p6p1);
				g = Math::Pow(c, p[6]);
				h = p[7] - Z[i, j] + p[0] / g;
				p0logc = p[0] * logc;

				f += h * h / den;

				grad[0] += (2 * (h)) / (den*g);
				grad[1] += -(twop0p6*((twocosp3*a) / p4sq - (twosinp3*b) / p5sq)*(h)) / (den*d);
				grad[2] += -(twop0p6*((twocosp3*b) / p5sq + (twosinp3*a) / p4sq)*(h)) / (den*d);
				grad[3] += -(twop0p6*((2 * a*b) / p4sq - (2 * a*b) / p5sq)*(h)) / (den*d);
				grad[4] += (fourp0p6*asq * (h)) / (den*p4cu * d);
				grad[5] += (fourp0p6*bsq * (h)) / (den*p5cu * d);
				grad[6] += -(2 * p0logc*(h)) / (den*g);
				grad[7] += (twop7 - 2 * Z[i, j] + twop0 / g) / den;
			}
		}
	}
}

void JPFITS::JPMath::alglib_Moffat_2d_LM_LS_ROBUST_grad(array<double>^ p, double %f, array<double>^ grad, Object^ obj)
{
	f = 0;
	for (int i = 0; i < grad->Length; i++)
		grad[i] = 0;
	array<double, 2>^ Z = (array<double, 2>^)(((array<Object^>^)obj)[0]);
	array<double>^ x = (array<double>^)(((array<Object^>^)obj)[1]);
	array<double>^ y = (array<double>^)(((array<Object^>^)obj)[2]);
	double den;

	if (p->Length == 6)
	{
		double p3sq = p[3] * p[3];
		double p41 = p[4] + 1;
		double twop1 = 2 * p[1];
		double twop2 = 2 * p[2];
		double p0p4 = p[0] * p[4];
		double p3cu = p3sq * p[3];
		double dx, dxsq, twoxi, twop1mtwoxi, dy, dysq, dxsqdysq, a, loga, apowp4, apowp41, b;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			dxsq = dx * dx;
			twoxi = 2 * x[i];
			twop1mtwoxi = twop1 - twoxi;

			for (int j = 0; j < y->Length; j++)
			{
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Math::Sqrt(Z[i, j]);

				dy = p[2] - y[j];
				dysq = dy * dy;
				dxsqdysq = dxsq + dysq;
				a = dxsqdysq / p3sq + 1;
				loga = Math::Log(a);
				apowp4 = Math::Pow(a, p[4]);
				apowp41 = Math::Pow(a, p41);
				b = p[5] - Z[i, j] + p[0] / apowp4;

				f += b * b / den;

				grad[0] += (2 * b) / (den*apowp4);
				grad[1] += -(2 * p0p4*twop1mtwoxi*b) / (den*p3sq * apowp41);
				grad[2] += -(2 * p0p4*(twop2 - 2 * y[j])*b) / (den*p3sq * apowp41);
				grad[3] += (4 * p0p4*dxsqdysq*b) / (den*p3cu * apowp41);
				grad[4] += -(2 * p[0]*loga*b) / (den*apowp4);
				grad[5] += (2 * p[5] - 2 * Z[i, j] + (2 * p[0]) / apowp4) / den;
			}
		}
		return;
	}

	if (p->Length == 8)
	{
		double cosp3 = Math::Cos(p[3]);
		double sinp3 = Math::Sin(p[3]);
		double p4sq = p[4] * p[4];
		double p5sq = p[5] * p[5];
		double p6p1 = p[6] + 1;
		double p0p6 = p[0] * p[6];
		double p4cu = p[4] * p4sq;
		double p5cu = p[5] * p5sq;
		double twop0p6 = 2 * p0p6;
		double fourp0p6 = 4 * p0p6;
		double twocosp3 = 2 * cosp3;
		double twosinp3 = 2 * sinp3;
		double twop7 = 2 * p[7];
		double twop0 = 2 * p[0];
		double dx, cosp3dx, sinp3dx, dy, sinp3dy, a, b, c, d, g, h, asq, cosp3dy, bsq, logc, p0logc;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			cosp3dx = cosp3 * dx;
			sinp3dx = sinp3 * dx;

			for (int j = 0; j < y->Length; j++)
			{
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Math::Sqrt(Z[i, j]);

				dy = p[2] - y[j];
				sinp3dy = sinp3 * dy;
				a = cosp3dx + sinp3dy;
				asq = a * a;
				cosp3dy = cosp3 * dy;
				b = cosp3dy - sinp3dx;
				bsq = b * b;
				c = asq / p4sq + bsq / p5sq + 1;
				logc = Math::Log(c);
				d = Math::Pow(c, p6p1);
				g = Math::Pow(c, p[6]);
				h = p[7] - Z[i, j] + p[0] / g;
				p0logc = p[0] * logc;

				f += h * h / den;

				grad[0] += (2 * (h)) / (den*g);
				grad[1] += -(twop0p6*((twocosp3*a) / p4sq - (twosinp3*b) / p5sq)*(h)) / (den*d);
				grad[2] += -(twop0p6*((twocosp3*b) / p5sq + (twosinp3*a) / p4sq)*(h)) / (den*d);
				grad[3] += -(twop0p6*((2 * a*b) / p4sq - (2 * a*b) / p5sq)*(h)) / (den*d);
				grad[4] += (fourp0p6*asq * (h)) / (den*p4cu * d);
				grad[5] += (fourp0p6*bsq * (h)) / (den*p5cu * d);
				grad[6] += -(2 * p0logc*(h)) / (den*g);
				grad[7] += (twop7 - 2 * Z[i, j] + twop0 / g) / den;
			}
		}
	}
}

void JPFITS::JPMath::alglib_Moffat_2d_LM_LS_CSTAT_grad(array<double>^ p, double %f, array<double>^ grad, Object^ obj)
{
	f = 0;
	for (int i = 0; i < grad->Length; i++)
		grad[i] = 0;
	array<double, 2>^ Z = (array<double, 2>^)(((array<Object^>^)obj)[0]);
	array<double>^ x = (array<double>^)(((array<Object^>^)obj)[1]);
	array<double>^ y = (array<double>^)(((array<Object^>^)obj)[2]);

	if (p->Length == 6)
	{
		double p3sq = p[3] * p[3];
		double p41 = p[4] + 1;
		double twop1 = 2 * p[1];
		double twop2 = 2 * p[2];
		double p0p4 = p[0] * p[4];
		double p3cu = p3sq * p[3];
		double dx, dxsq, twoxi, twop1mtwoxi, dy, dysq, dxsqdysq, a, loga, apowp4, apowp41, b, d;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			dxsq = dx * dx;
			twoxi = 2 * x[i];
			twop1mtwoxi = twop1 - twoxi;

			for (int j = 0; j < y->Length; j++)
			{
				dy = p[2] - y[j];
				dysq = dy * dy;
				dxsqdysq = dxsq + dysq;
				a = dxsqdysq / p3sq + 1;
				loga = Math::Log(a);
				apowp4 = Math::Pow(a, p[4]);
				apowp41 = Math::Pow(a, p41);
				b = p[5] - Z[i, j] + p[0] / apowp4;
				d = p[5] + p[0] / apowp4;

				f += 2 * p[5] - 2 * Z[i, j]*Math::Log(d) + (2 * p[0]) / apowp4;

				grad[0] += 2 / apowp4 - (2 * Z[i, j]) / (d*apowp4);
				grad[1] += (2 * Z[i, j] *p0p4*twop1mtwoxi) / (p3sq * d*apowp41) - (2 * p0p4*twop1mtwoxi) / (p3sq * apowp41);
				grad[2] += (2 * Z[i, j] *p0p4*(twop2 - 2 * y[j])) / (p3sq * d*apowp41) - (2 * p0p4*(twop2 - 2 * y[j])) / (p3sq * apowp41);
				grad[3] += (4 * p0p4*dxsqdysq) / (p3cu * apowp41) - (4 * Z[i, j] *p0p4*dxsqdysq) / (p3cu * d*apowp41);
				grad[4] += (2 * Z[i, j] *p[0]*loga) / (d*apowp4) - (2 * p[0]*loga) / apowp4;
				grad[5] += 2 - (2 * Z[i, j]) / d;
			}
		}
		return;
	}

	if (p->Length == 8)
	{
		double cosp3 = Math::Cos(p[3]);
		double sinp3 = Math::Sin(p[3]);
		double p4sq = p[4] * p[4];
		double p5sq = p[5] * p[5];
		double p61 = p[6] + 1;
		double p4cu = p[4] * p4sq;
		double p5cu = p[5] * p5sq;
		double p0p6 = p[0] * p[6];
		double dx, cosp3dx, sinp3dx, dy, sinp3dy, cosp3dy, cosp3dxsinp3dy, cosp3dysinp3dx, cosp3dxsinp3dysq, cosp3dysinp3dxsq, a, loga, apowp6, apowp61, b, twoZij;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			cosp3dx = cosp3 * dx;
			sinp3dx = sinp3 * dx;

			for (int j = 0; j < y->Length; j++)
			{
				dy = p[2] - y[j];
				sinp3dy = sinp3 * dy;
				cosp3dy = cosp3 * dy;
				cosp3dxsinp3dy = cosp3dx + sinp3dy;
				cosp3dysinp3dx = cosp3dy - sinp3dx;
				cosp3dxsinp3dysq = cosp3dxsinp3dy * cosp3dxsinp3dy;
				cosp3dysinp3dxsq = cosp3dysinp3dx * cosp3dysinp3dx;
				a = cosp3dxsinp3dysq / p4sq + cosp3dysinp3dxsq / p5sq + 1;
				loga = Math::Log(a);
				apowp6 = Math::Pow(a, p[6]);
				apowp61 = apowp6 * a;
				b = p[7] + p[0] / apowp6;
				twoZij = 2 * Z[i, j];

				f += 2 * p[7] - twoZij * log(b) + (2 * p[0]) / apowp6;

				grad[0] += 2 / apowp6 - (twoZij) / ((b)*apowp6);
				grad[1] += (twoZij*p0p6*((2 * cosp3*cosp3dxsinp3dy) / p4sq - (2 * sinp3*cosp3dysinp3dx) / p5sq)) / ((b)*apowp61) - (2 * p0p6*((2 * cosp3*cosp3dxsinp3dy) / p4sq - (2 * sinp3*cosp3dysinp3dx) / p5sq)) / apowp61;
				grad[2] += (twoZij*p0p6*((2 * cosp3*cosp3dysinp3dx) / p5sq + (2 * sinp3*cosp3dxsinp3dy) / p4sq)) / ((b)*apowp61) - (2 * p0p6*((2 * cosp3*cosp3dysinp3dx) / p5sq + (2 * sinp3*cosp3dxsinp3dy) / p4sq)) / apowp61;
				grad[3] += (twoZij*p0p6*((2 * cosp3dxsinp3dy*cosp3dysinp3dx) / p4sq - (2 * cosp3dxsinp3dy*cosp3dysinp3dx) / p5sq)) / ((b)*apowp61) - (2 * p0p6*((2 * cosp3dxsinp3dy*cosp3dysinp3dx) / p4sq - (2 * cosp3dxsinp3dy*cosp3dysinp3dx) / p5sq)) / apowp61;
				grad[4] += (4 * p0p6*cosp3dxsinp3dysq) / (p4cu * apowp61) - (4 * Z[i, j] *p0p6*cosp3dxsinp3dysq) / (p4cu * (b)*apowp61);
				grad[5] += (4 * p0p6*cosp3dysinp3dxsq) / (p5cu * apowp61) - (4 * Z[i, j] *p0p6*cosp3dysinp3dxsq) / (p5cu * (b)*apowp61);
				grad[6] += (twoZij*p[0]*loga) / ((b)*apowp6) - (2 * p[0]*loga) / apowp6;
				grad[7] += 2 - (twoZij) / (b);					
			}
		}
	}
}

array<double>^ JPFITS::JPMath::Gauss_2D_param_err(array<double>^ p, array<double>^ x, array<double>^ y, array<double, 2>^ Z)
{
	array<double, 2>^ hess = gcnew array<double, 2>(p->Length, p->Length);
	double den;

	if (p->Length == 5)
	{
		double p3sq = p[3] * p[3];
		double p3cu = p[3] * p3sq;
		double p3qu = p3sq * p3sq;
		double p3he = p3cu * p3cu;
		double p0sq = p[0] * p[0];
		double twop3sq = 2 * p3sq;
		double p3pe = p3sq * p3cu;
		double twop1 = 2 * p[1];
		double twop2 = 2 * p[2];
		double dx, dxsq, twoxi, twop1mtwoxi, twop1mtwoxisq, dy, dxsqpdysq, exparg, expres, exparg2, expres2, p0expres2, a, dxsqpdysqsq, twoyj, twop2mtwoyj, twop2mtwoyjsq;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			dxsq = dx * dx;
			twoxi = 2 * x[i];
			twop1mtwoxi = twop1 - twoxi;
			twop1mtwoxisq = twop1mtwoxi * twop1mtwoxi;

			for (int j = 0; j < y->Length; j++)
			{
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Z[i, j];

				dy = p[2] - y[j];
				dxsqpdysq = dxsq + dy * dy;
				exparg = -dxsqpdysq / p3sq;
				expres = Math::Exp(exparg);
				exparg2 = -dxsqpdysq / twop3sq;
				expres2 = Math::Exp(exparg2);
				p0expres2 = p[0] * expres2;
				a = p[4] - Z[i, j] + p0expres2;
				dxsqpdysqsq = dxsqpdysq * dxsqpdysq;
				twoyj = 2 * y[j];
				twop2mtwoyj = twop2 - twoyj;
				twop2mtwoyjsq = twop2mtwoyj * twop2mtwoyj;

				hess[0, 0] += (2 * expres) / den;
				hess[1, 1] += (p0sq * expres*twop1mtwoxisq) / (2 * den*p3qu) - (2 * p0expres2*a) / (den*p3sq) + (p0expres2*twop1mtwoxisq * a) / (2 * den*p3qu);
				hess[2, 2] += (p0sq * expres*twop2mtwoyjsq) / (2 * den*p3qu) - (2 * p0expres2*a) / (den*p3sq) + (p0expres2*twop2mtwoyjsq * a) / (2 * den*p3qu);
				hess[3, 3] += (2 * p0sq * expres*dxsqpdysqsq) / (den*p3he) - (6 * p0expres2*(dxsqpdysq)*a) / (den*p3qu) + (2 * p0expres2*dxsqpdysqsq * a) / (den*p3he);
				hess[4, 4] += 2 / den;
				hess[0, 1] += -(expres2*twop1mtwoxi*a) / (den*p3sq) - (p[0]*expres*twop1mtwoxi) / (den*p3sq);
				hess[0, 2] += -(expres2*twop2mtwoyj*a) / (den*p3sq) - (p[0]*expres*twop2mtwoyj) / (den*p3sq);
				hess[0, 3] += (2 * expres2*(dxsqpdysq)*a) / (den*p3cu) + (2 * p[0]*expres*(dxsqpdysq)) / (den*p3cu);
				hess[0, 4] += (2 * expres2) / den;
				hess[1, 2] += (p0sq * expres*twop1mtwoxi*twop2mtwoyj) / (2 * den*p3qu) + (p0expres2*twop1mtwoxi*twop2mtwoyj*a) / (2 * den*p3qu);
				hess[1, 3] += (2 * p0expres2*twop1mtwoxi*a) / (den*p3cu) - (p0sq * expres*twop1mtwoxi*(dxsqpdysq)) / (den*p3pe) - (p0expres2*twop1mtwoxi*(dxsqpdysq)*a) / (den*p3pe);
				hess[1, 4] += -(p0expres2*twop1mtwoxi) / (den*p3sq);
				hess[2, 3] += (2 * p0expres2*twop2mtwoyj*a) / (den*p3cu) - (p0sq * expres*twop2mtwoyj*(dxsqpdysq)) / (den*p3pe) - (p0expres2*twop2mtwoyj*(dxsqpdysq)*a) / (den*p3pe);
				hess[2, 4] += -(p0expres2*twop2mtwoyj) / (den*p3sq);
				hess[3, 4] += (2 * p0expres2*(dxsqpdysq)) / (den*p3cu);
			}
		}
		hess[1, 0] = hess[0, 1];
		hess[2, 0] = hess[0, 2];
		hess[3, 0] = hess[0, 3];
		hess[4, 0] = hess[0, 4];
		hess[2, 1] = hess[1, 2];
		hess[3, 1] = hess[1, 3];
		hess[4, 1] = hess[1, 4];
		hess[3, 2] = hess[2, 3];
		hess[4, 2] = hess[2, 4];
		hess[4, 3] = hess[3, 4];
	}
	else if (p->Length == 7)
	{
		double cosp3 = Math::Cos(p[3]);
		double sinp3 = Math::Sin(p[3]);
		double p4sq = p[4] * p[4];
		double p5sq = p[5] * p[5];
		double p0sq = p[0] * p[0];
		double cosp3sq = cosp3 * cosp3;
		double sinp3sq = sinp3 * sinp3;
		double twop0sq = 2 * p0sq;
		double twop0 = 2 * p[0];
		double p4cu = p4sq * p[4];
		double p5cu = p5sq * p[5];
		double p4he = p4cu * p4cu;
		double p5he = p5cu * p5cu;
		double p4qu = p4sq * p4sq;
		double p5qu = p5sq * p5sq;
		double dx, cosp3dx, dy, sinp3dy, a, b, asq, bsq, exarg1, expres1, c, exparg2, expres2, p0expres2, g, twop0sqexpres1, h, k;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			cosp3dx = cosp3 * dx;

			for (int j = 0; j < y->Length; j++)
			{
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Z[i, j];

				dy = p[2] - y[j];
				sinp3dy = sinp3 * dy;
				a = cosp3dx + sinp3dy;
				b = cosp3 * dy - sinp3 * dx;
				asq = a * a;
				bsq = b * b;
				exarg1 = -asq / p4sq - bsq / p5sq;
				expres1 = Math::Exp(exarg1);
				c = (cosp3*a) / p4sq - (sinp3*b) / p5sq;
				exparg2 = -asq / (2 * p4sq) - bsq / (2 * p5sq);
				expres2 = Math::Exp(exparg2);
				p0expres2 = p[0] * expres2;
				g = p[6] - Z[i, j] + p0expres2;
				twop0sqexpres1 = twop0sq * expres1;
				h = (cosp3*b) / p5sq + (sinp3*a) / p4sq;
				k = (a*b) / p4sq - (a*b) / p5sq;

				hess[0, 0] += (2 * expres1) / den;
				hess[1, 1] += (twop0sqexpres1*c * c) / den + (2 * p0expres2*c * c * g) / den - (2 * p0expres2*(cosp3sq / p4sq + sinp3sq / p5sq)*g) / den;
				hess[2, 2] += (twop0sqexpres1*h * h) / den + (2 * p0expres2*h * h * g) / den - (2 * p0expres2*(cosp3sq / p5sq + sinp3sq / p4sq)*g) / den;
				hess[3, 3] += (twop0sqexpres1*k * k) / den + (2 * p0expres2*g*(asq / p4sq - bsq / p4sq - asq / p5sq + bsq / p5sq)) / den + (2 * p0expres2*k * k * g) / den;
				hess[4, 4] += (twop0sqexpres1*a*a*a*a) / (den*p4he) - (6 * p0expres2*asq * g) / (den*p4qu) + (2 * p0expres2*a*a*a*a * g) / (den*p4he);
				hess[5, 5] += (twop0sqexpres1*b*b*b*b) / (den*p5he) - (6 * p0expres2*bsq * g) / (den*p5qu) + (2 * p0expres2*b*b*b*b * g) / (den*p5he);
				hess[6, 6] += 2 / den;
				hess[0, 1] += -(twop0*expres1*c) / den - (2 * expres2*c*g) / den;
				hess[0, 2] += -(twop0*expres1*h) / den - (2 * expres2*h*g) / den;
				hess[0, 3] += -(twop0*expres1*k) / den - (2 * expres2*k*g) / den;
				hess[0, 4] += (twop0*expres1*asq) / (den*p4cu) + (2 * expres2*asq * g) / (den*p4cu);
				hess[0, 5] += (twop0*expres1*bsq) / (den*p5cu) + (2 * expres2*bsq * g) / (den*p5cu);
				hess[0, 6] += (2 * expres2) / den;
				hess[1, 2] += (twop0sqexpres1*h*c) / den - (2 * p0expres2*((cosp3*sinp3) / p4sq - (cosp3*sinp3) / p5sq)*g) / den + (2 * p0expres2*h*c*g) / den;
				hess[1, 3] += (twop0sqexpres1*k*c) / den - (2 * p0expres2*g*((cosp3*b) / p4sq - (cosp3*b) / p5sq - (sinp3*a) / p4sq + (sinp3*a) / p5sq)) / den + (2 * p0expres2*k*c*g) / den;
				hess[1, 4] += (4 * p0expres2*cosp3*a*g) / (den*p4cu) - (2 * p0expres2*asq * c*g) / (den*p4cu) - (twop0sqexpres1*asq * c) / (den*p4cu);
				hess[1, 5] += -(twop0sqexpres1*bsq * c) / (den*p5cu) - (4 * p0expres2*sinp3*b*g) / (den*p5cu) - (2 * p0expres2*bsq * c*g) / (den*p5cu);
				hess[1, 6] += -(2 * p0expres2*c) / den;
				hess[2, 3] += (twop0sqexpres1*k*h) / den - (2 * p0expres2*g*((cosp3*a) / p4sq - (cosp3*a) / p5sq + (sinp3*b) / p4sq - (sinp3*b) / p5sq)) / den + (2 * p0expres2*k*h*g) / den;
				hess[2, 4] += (4 * p0expres2*sinp3*a*g) / (den*p4cu) - (twop0sqexpres1*asq * h) / (den*p4cu) - (2 * p0expres2*asq * h*g) / (den*p4cu);
				hess[2, 5] += (4 * p0expres2*cosp3*b*g) / (den*p5cu) - (2 * p0expres2*bsq * h*g) / (den*p5cu) - (twop0sqexpres1*bsq * h) / (den*p5cu);
				hess[2, 6] += -(2 * p0expres2*h) / den;
				hess[3, 4] += (4 * p0expres2*a*b*g) / (den*p4cu) - (twop0sqexpres1*asq * k) / (den*p4cu) - (2 * p0expres2*asq * k*g) / (den*p4cu);
				hess[3, 5] += -(twop0sqexpres1*bsq * k) / (den*p5cu) - (4 * p0expres2*a*b*g) / (den*p5cu) - (2 * p0expres2*bsq * k*g) / (den*p5cu);
				hess[3, 6] += -(2 * p0expres2*k) / den;
				hess[4, 5] += (twop0sqexpres1*asq * bsq) / (den*p4cu * p5cu) + (2 * p0expres2*asq * bsq * g) / (den*p4cu * p5cu);
				hess[4, 6] += (2 * p0expres2*asq) / (den*p4cu);
				hess[5, 6] += (2 * p0expres2*bsq) / (den*p5cu);
			}
		}
		hess[1, 0] = hess[0, 1];
		hess[2, 0] = hess[0, 2];
		hess[3, 0] = hess[0, 3];
		hess[4, 0] = hess[0, 4];
		hess[5, 0] = hess[0, 5];
		hess[6, 0] = hess[0, 6];
		hess[2, 1] = hess[1, 2];
		hess[3, 1] = hess[1, 3];
		hess[4, 1] = hess[1, 4];
		hess[5, 1] = hess[1, 5];
		hess[6, 1] = hess[1, 6];
		hess[3, 2] = hess[2, 3];
		hess[4, 2] = hess[2, 4];
		hess[5, 2] = hess[2, 5];
		hess[6, 2] = hess[2, 6];
		hess[4, 3] = hess[3, 4];
		hess[5, 3] = hess[3, 5];
		hess[6, 3] = hess[3, 6];
		hess[5, 4] = hess[4, 5];
		hess[6, 4] = hess[4, 6];
		hess[6, 5] = hess[5, 6];
	}

	int info;
	alglib::matinvreport^ rep;
	alglib::rmatrixinverse(hess, info, rep);
	array<double>^ errs = gcnew array<double>(p->Length);
	for (int i = 0; i < p->Length; i++)
		errs[i] = Math::Sqrt(hess[i, i]);
		
	return errs;
}

array<double>^ JPFITS::JPMath::Moffat_2D_param_err(array<double>^ p, array<double>^ x, array<double>^ y, array<double, 2>^ Z)
{
	array<double, 2>^ hess = gcnew array<double, 2>(p->Length, p->Length);
	double den;

	if (p->Length == 6)
	{
		double p3sq = p[3] * p[3];
		double p41 = p[4] + 1;
		double p42 = p[4] + 2;
		double p3qu = p3sq * p3sq;
		double p3pe = p3qu * p[3];
		double p3cu = p3sq * p[3];
		double p3he = p3cu * p3cu;
		double p0sq = p[0] * p[0];
		double p4sq = p[4] * p[4];
		double p0p4 = p[0] * p[4];
		double p0sqp4sq = p0sq * p4sq;
		double twop1 = 2 * p[1];
		double twop2 = 2 * p[2];

		for (int i = 0; i < x->Length; i++)
		{
			double dx = p[1] - x[i];
			double dxsq = dx * dx;
			double twoxi = 2 * x[i];

			for (int j = 0; j < y->Length; j++)
			{
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Z[i, j];

				double dy = p[2] - y[j];
				double dysq = dy * dy;
				double dxsqpdysq = dxsq + dysq;
				double logarg1 = dxsqpdysq / p3sq + 1;
				double logres1 = Math::Log(logarg1);
				double logarg1powp4 = Math::Pow(logarg1, p[4]);
				double logarg1powp41 = logarg1powp4 * logarg1;
				double logarg1powp42 = logarg1powp41 * logarg1;
				double a = p[5] - Z[i, j] + p[0] / logarg1powp4;
				double twoyj = 2 * y[j];
				double twop1mtwoxi = twop1 - twoxi;
				double twop2mtwoyj = twop2 - twoyj;
				double logarg1pow2p4 = logarg1powp4 * logarg1powp4;
				double logarg1pow2p42 = Math::Pow(logarg1, 2 * p[4] + 2);

				hess[0, 0] += (2 / logarg1pow2p4) / den;
				hess[1, 1] += (2 * p0sqp4sq * twop1mtwoxi * twop1mtwoxi) / (den*p3qu * logarg1pow2p42) - (4 * p0p4*a) / (den*p3sq * logarg1powp41) + (2 * p0p4*twop1mtwoxi * twop1mtwoxi * p41*a) / (den*p3qu * logarg1powp42);
				hess[2, 2] += (2 * p0sqp4sq * twop2mtwoyj * twop2mtwoyj) / (den*p3qu * logarg1pow2p42) - (4 * p0p4*a) / (den*p3sq * logarg1powp41) + (2 * p0p4*twop2mtwoyj * twop2mtwoyj * p41*a) / (den*p3qu * logarg1powp42);
				hess[3, 3] += (8 * p0sqp4sq * dxsqpdysq * dxsqpdysq) / (den*p3he * logarg1pow2p42) - (12 * p0p4*dxsqpdysq*a) / (den*p3qu * logarg1powp41) + (8 * p0p4*dxsqpdysq * dxsqpdysq * p41*a) / (den*p3he * logarg1powp42);
				hess[4, 4] += (2 * p0sq * logres1 * logres1 / logarg1pow2p4) / den + (2 * p[0]*logres1 * logres1 * a) / (den*logarg1powp4);
				hess[5, 5] += 2 / den;
				hess[0, 1] += -(2 * p[4]*twop1mtwoxi*a) / (den*p3sq * logarg1powp41) - (2 * p0p4*twop1mtwoxi) / (den*p3sq * logarg1powp4*logarg1powp41);
				hess[0, 2] += -(2 * p[4]*twop2mtwoyj*a) / (den*p3sq * logarg1powp41) - (2 * p0p4*twop2mtwoyj) / (den*p3sq * logarg1powp4*logarg1powp41);
				hess[0, 3] += (4 * p[4]*dxsqpdysq*a) / (den*p3cu * logarg1powp41) + (4 * p0p4*dxsqpdysq) / (den*p3cu * logarg1powp4*logarg1powp41);
				hess[0, 4] += -(2 * p[0]*logres1 / logarg1pow2p4) / den - (2 * logres1*a) / (den*logarg1powp4);
				hess[0, 5] += 2 / (den*logarg1powp4);
				hess[1, 2] += (2 * p0sqp4sq * twop1mtwoxi*twop2mtwoyj) / (den*p3qu * logarg1pow2p42) + (2 * p0p4*twop1mtwoxi*twop2mtwoyj*p41*a) / (den*p3qu * logarg1powp42);
				hess[1, 3] += (4 * p0p4*twop1mtwoxi*a) / (den*p3cu * logarg1powp41) - (4 * p0sqp4sq * twop1mtwoxi*dxsqpdysq) / (den*p3pe * logarg1pow2p42) - (4 * p0p4*twop1mtwoxi*dxsqpdysq*p41*a) / (den*p3pe * logarg1powp42);
				hess[1, 4] += (2 * p0p4*logres1*twop1mtwoxi*a) / (den*p3sq * logarg1powp41) - (2 * p[0]*twop1mtwoxi*a) / (den*p3sq * logarg1powp41) + (2 * p0sq * p[4]*logres1*twop1mtwoxi) / (den*p3sq * logarg1powp4*logarg1powp41);
				hess[1, 5] += -(2 * p0p4*twop1mtwoxi) / (den*p3sq * logarg1powp41);
				hess[2, 3] += (4 * p0p4*twop2mtwoyj*a) / (den*p3cu * logarg1powp41) - (4 * p0sqp4sq * twop2mtwoyj*dxsqpdysq) / (den*p3pe * logarg1pow2p42) - (4 * p0p4*twop2mtwoyj*dxsqpdysq*p41*a) / (den*p3pe * logarg1powp42);
				hess[2, 4] += (2 * p0p4*logres1*twop2mtwoyj*a) / (den*p3sq * logarg1powp41) - (2 * p[0]*twop2mtwoyj*a) / (den*p3sq * logarg1powp41) + (2 * p0sq * p[4]*logres1*twop2mtwoyj) / (den*p3sq * logarg1powp4*logarg1powp41);
				hess[2, 5] += -(2 * p0p4*twop2mtwoyj) / (den*p3sq * logarg1powp41);
				hess[3, 4] += (4 * p[0]*dxsqpdysq*a) / (den*p3cu * logarg1powp41) - (4 * p0p4*logres1*dxsqpdysq*a) / (den*p3cu * logarg1powp41) - (4 * p0sq * p[4]*logres1*dxsqpdysq) / (den*p3cu * logarg1powp4*logarg1powp41);
				hess[3, 5] += (4 * p0p4*dxsqpdysq) / (den*p3cu * logarg1powp41);
				hess[4, 5] += -(2 * p[0]*logres1) / (den*logarg1powp4);
			}
		}
		hess[1, 0] = hess[0, 1];
		hess[2, 0] = hess[0, 2];
		hess[3, 0] = hess[0, 3];
		hess[4, 0] = hess[0, 4];
		hess[5, 0] = hess[0, 5];
		hess[2, 1] = hess[1, 2];
		hess[3, 1] = hess[1, 3];
		hess[4, 1] = hess[1, 4];
		hess[5, 1] = hess[1, 5];
		hess[3, 2] = hess[2, 3];
		hess[4, 2] = hess[2, 4];
		hess[5, 2] = hess[2, 5];
		hess[4, 3] = hess[3, 4];
		hess[5, 3] = hess[3, 5];
		hess[5, 4] = hess[4, 5];
	}
	else if (p->Length == 8)
	{
		double cosp3 = Math::Cos(p[3]);
		double sinp3 = Math::Sin(p[3]);
		double p4sq = p[4] * p[4];
		double p5sq = p[5] * p[5];
		double twocosp3 = 2 * cosp3;
		double twosinp3 = 2 * sinp3;
		double p0sq = p[0] * p[0];
		double p6sq = p[6] * p[6];
		double p0p6 = p[0] * p[6];
		double p0sqp6sq = p0sq * p6sq;
		double p61 = p[6] + 1;
		double twocosp3sq = twocosp3 * twocosp3;
		double twosinp3sq = twosinp3 * twosinp3;
		double p4cu = p4sq * p[4];
		double p4he = p4cu * p4cu;
		double p5cu = p5sq * p[5];
		double p5he = p5cu * p5cu;
		double p4qu = p4sq * p4sq;
		double p5qu = p5sq * p5sq;
		double dx, cosp3dx, sinp3dx, dy, cosp3dxpsinp3dy, cosp3dymsinp3dx, logarg, logres, logargpowp6, logargpow2p6, logargpowp61, logargpowp62, logargpow2p62, a, b, c, d, dsq, f, fsq, g, gsq, h, k;

		for (int i = 0; i < x->Length; i++)
		{
			dx = p[1] - x[i];
			cosp3dx = cosp3 * dx;
			sinp3dx = sinp3 * dx;

			for (int j = 0; j < y->Length; j++)
			{
				if (Z[i, j] < 1)
					den = 1;
				else
					den = Z[i, j];

				dy = p[2] - y[j];
				cosp3dxpsinp3dy = cosp3dx + sinp3 * dy;
				cosp3dymsinp3dx = cosp3 * dy - sinp3dx;
				logarg = cosp3dxpsinp3dy * cosp3dxpsinp3dy / p4sq + cosp3dymsinp3dx * cosp3dymsinp3dx / p5sq + 1;
				logres = Math::Log(logarg);
				logargpowp6 = Math::Pow(logarg, p[6]);
				logargpow2p6 = logargpowp6 * logargpowp6;
				logargpowp61 = logargpowp6 * logarg;
				logargpowp62 = logargpowp61 * logarg;
				logargpow2p62 = logargpowp61 * logargpowp61;
				a = p[7] - Z[i, j] + p[0] / logargpowp6;
				b = cosp3dxpsinp3dy * cosp3dxpsinp3dy;
				c = cosp3dymsinp3dx * cosp3dymsinp3dx;
				d = twocosp3 * cosp3dxpsinp3dy / p4sq - twosinp3 * cosp3dymsinp3dx / p5sq;
				dsq = d * d;
				f = twocosp3 * cosp3dymsinp3dx / p5sq + twosinp3 * cosp3dxpsinp3dy / p4sq;
				fsq = f * f;
				g = 2 * cosp3dxpsinp3dy*cosp3dymsinp3dx / p4sq - 2 * cosp3dxpsinp3dy*cosp3dymsinp3dx / p5sq;
				gsq = g * g;
				h = cosp3dxpsinp3dy * cosp3dxpsinp3dy * cosp3dxpsinp3dy * cosp3dxpsinp3dy;
				k = cosp3dymsinp3dx * cosp3dymsinp3dx * cosp3dymsinp3dx * cosp3dymsinp3dx;

				hess[0, 0] += (2 / logargpow2p6) / den;
				hess[1, 1] += (2 * p0sqp6sq * dsq) / (den*logargpow2p62) - (2 * p0p6*((twocosp3sq) / p4sq + (twosinp3sq) / p5sq)*(a)) / (den*logargpowp61) + (2 * p0p6*dsq * p61*(a)) / (den*logargpowp62);
				hess[2, 2] += (2 * p0sqp6sq * fsq) / (den*logargpow2p62) - (2 * p0p6*((twocosp3sq) / p5sq + (twosinp3sq) / p4sq)*(a)) / (den*logargpowp61) + (2 * p0p6*fsq * p61*(a)) / (den*logargpowp62);
				hess[3, 3] += (2 * p0sqp6sq * gsq) / (den*logargpow2p62) + (2 * p0p6*(a)*((2 * b) / p4sq - (2 * c) / p4sq - (2 * b) / p5sq + (2 * c) / p5sq)) / (den*logargpowp61) + (2 * p0p6*gsq * p61*(a)) / (den*logargpowp62);
				hess[4, 4] += (8 * p0sqp6sq * h) / (den*p4he * logargpow2p62) - (12 * p0p6*b * (a)) / (den*p4qu * logargpowp61) + (8 * p0p6*h * p61*(a)) / (den*p4he * logargpowp62);
				hess[5, 5] += (8 * p0sqp6sq * k) / (den*p5he * logargpow2p62) - (12 * p0p6*c * (a)) / (den*p5qu * logargpowp61) + (8 * p0p6*k * p61*(a)) / (den*p5he * logargpowp62);
				hess[6, 6] += (2 * p0sq * logres * logres / logargpow2p6) / den + (2 * p[0] * logres * logres * (a)) / (den*logargpowp6);
				hess[7, 7] += 2 / den;
				hess[0, 1] += -(2 * p[6] * d*(a)) / (den*logargpowp61) - (2 * p0p6*d) / (den*logargpowp61*logargpowp6);
				hess[0, 2] += -(2 * p[6] * (f)*(a)) / (den*logargpowp61) - (2 * p0p6*(f)) / (den*logargpowp61*logargpowp6);
				hess[0, 3] += -(2 * p[6] * (g)*(a)) / (den*logargpowp61) - (2 * p0p6*(g)) / (den*logargpowp61*logargpowp6);
				hess[0, 4] += (4 * p[6] * b * (a)) / (den*p4cu * logargpowp61) + (4 * p0p6*b) / (den*p4cu * logargpowp61*logargpowp6);
				hess[0, 5] += (4 * p[6] * c * (a)) / (den*p5cu * logargpowp61) + (4 * p0p6*c) / (den*p5cu * logargpowp61*logargpowp6);
				hess[0, 6] += -(2 * logres*(a)) / (den*logargpowp6) - (2 * p[0] * logres / logargpow2p6) / den;
				hess[0, 7] += 2 / (den*logargpowp6);
				hess[1, 2] += (2 * p0sqp6sq * d*(f)) / (den*logargpow2p62) - (2 * p0p6*((twocosp3*sinp3) / p4sq - (twocosp3*sinp3) / p5sq)*(a)) / (den*logargpowp61) + (2 * p0p6*d*(f)*p61*(a)) / (den*logargpowp62);
				hess[1, 3] += (2 * p0sqp6sq * (g)*d) / (den*logargpow2p62) - (2 * p0p6*(a)*((twocosp3*cosp3dymsinp3dx) / p4sq - (twocosp3*cosp3dymsinp3dx) / p5sq - (twosinp3*cosp3dxpsinp3dy) / p4sq + (twosinp3*cosp3dxpsinp3dy) / p5sq)) / (den*logargpowp61) + (2 * p0p6*(g)*d*p61*(a)) / (den*logargpowp62);
				hess[1, 4] += (8 * p0p6*cosp3*cosp3dxpsinp3dy*(a)) / (den*p4cu * logargpowp61) - (4 * p0sqp6sq * b * d) / (den*p4cu * logargpow2p62) - (4 * p0p6*b * d*p61*(a)) / (den*p4cu * logargpowp62);
				hess[1, 5] += -(4 * p0sqp6sq * c * d) / (den*p5cu * logargpow2p62) - (8 * p0p6*sinp3*cosp3dymsinp3dx*(a)) / (den*p5cu * logargpowp61) - (4 * p0p6*c * d*p61*(a)) / (den*p5cu * logargpowp62);
				hess[1, 6] += (2 * p0p6*logres*d*(a)) / (den*logargpowp61) - (2 * p[0] * d*(a)) / (den*logargpowp61) + (2 * p0sq * p[6] * logres*d) / (den*logargpowp61*logargpowp6);
				hess[1, 7] += -(2 * p0p6*d) / (den*logargpowp61);
				hess[2, 3] += (2 * p0sqp6sq * (g)*(f)) / (den*logargpow2p62) - (2 * p0p6*(a)*((twocosp3*cosp3dxpsinp3dy) / p4sq - (twocosp3*cosp3dxpsinp3dy) / p5sq + (twosinp3*cosp3dymsinp3dx) / p4sq - (twosinp3*cosp3dymsinp3dx) / p5sq)) / (den*logargpowp61) + (2 * p0p6*(g)*(f)*p61*(a)) / (den*logargpowp62);
				hess[2, 4] += (8 * p0p6*sinp3*cosp3dxpsinp3dy*(a)) / (den*p4cu * logargpowp61) - (4 * p0sqp6sq * b * (f)) / (den*p4cu * logargpow2p62) - (4 * p0p6*b * (f)*p61*(a)) / (den*p4cu * logargpowp62);
				hess[2, 5] += (8 * p0p6*cosp3*cosp3dymsinp3dx*(a)) / (den*p5cu * logargpowp61) - (4 * p0sqp6sq * c * (f)) / (den*p5cu * logargpow2p62) - (4 * p0p6*c * (f)*p61*(a)) / (den*p5cu * logargpowp62);
				hess[2, 6] += (2 * p0p6*logres*(f)*(a)) / (den*logargpowp61) - (2 * p[0] * (f)*(a)) / (den*logargpowp61) + (2 * p0sq * p[6] * logres*(f)) / (den*logargpowp61*logargpowp6);
				hess[2, 7] += -(2 * p0p6*(f)) / (den*logargpowp61);
				hess[3, 4] += (8 * p0p6*cosp3dxpsinp3dy*cosp3dymsinp3dx*(a)) / (den*p4cu * logargpowp61) - (4 * p0sqp6sq * b * (g)) / (den*p4cu * logargpow2p62) - (4 * p0p6*b * (g)*p61*(a)) / (den*p4cu * logargpowp62);
				hess[3, 5] += -(4 * p0sqp6sq * c * (g)) / (den*p5cu * logargpow2p62) - (8 * p0p6*cosp3dxpsinp3dy*cosp3dymsinp3dx*(a)) / (den*p5cu * logargpowp61) - (4 * p0p6*c * (g)*p61*(a)) / (den*p5cu * logargpowp62);
				hess[3, 6] += (2 * p0p6*logres*(g)*(a)) / (den*logargpowp61) - (2 * p[0] * (g)*(a)) / (den*logargpowp61) + (2 * p0sq * p[6] * logres*(g)) / (den*logargpowp61*logargpowp6);
				hess[3, 7] += -(2 * p0p6*(g)) / (den*logargpowp61);
				hess[4, 5] += (8 * p0sqp6sq * b * c) / (den*p4cu * p5cu * logargpow2p62) + (8 * p0p6*b * c * p61*(a)) / (den*p4cu * p5cu * logargpowp62);
				hess[4, 6] += (4 * p[0] * b * (a)) / (den*p4cu * logargpowp61) - (4 * p0p6*logres*b * (a)) / (den*p4cu * logargpowp61) - (4 * p0sq * p[6] * logres*b) / (den*p4cu * logargpowp61*logargpowp6);
				hess[4, 7] += (4 * p0p6*b) / (den*p4cu * logargpowp61);
				hess[5, 6] += (4 * p[0] * c * (a)) / (den*p5cu * logargpowp61) - (4 * p0p6*logres*c * (a)) / (den*p5cu * logargpowp61) - (4 * p0sq * p[6] * logres*c) / (den*p5cu * logargpowp61*logargpowp6);
				hess[5, 7] += (4 * p0p6*c) / (den*p5cu * logargpowp61);
				hess[6, 7] += -(2 * p[0] * logres) / (den*logargpowp6);
			}
		}
		hess[1, 0] = hess[0, 1];
		hess[2, 0] = hess[0, 2];
		hess[3, 0] = hess[0, 3];
		hess[4, 0] = hess[0, 4];
		hess[5, 0] = hess[0, 5];
		hess[6, 0] = hess[0, 6];
		hess[7, 0] = hess[0, 7];
		hess[2, 1] = hess[1, 2];
		hess[3, 1] = hess[1, 3];
		hess[4, 1] = hess[1, 4];
		hess[5, 1] = hess[1, 5];
		hess[6, 1] = hess[1, 6];
		hess[7, 1] = hess[1, 7];
		hess[3, 2] = hess[2, 3];
		hess[4, 2] = hess[2, 4];
		hess[5, 2] = hess[2, 5];
		hess[6, 2] = hess[2, 6];
		hess[7, 2] = hess[2, 7];
		hess[4, 3] = hess[3, 4];
		hess[5, 3] = hess[3, 5];
		hess[6, 3] = hess[3, 6];
		hess[7, 3] = hess[3, 7];
		hess[5, 4] = hess[4, 5];
		hess[6, 4] = hess[4, 6];
		hess[7, 4] = hess[4, 7];
		hess[6, 5] = hess[5, 6];
		hess[7, 5] = hess[5, 7];
		hess[7, 6] = hess[6, 7];
	}

	int info;
	alglib::matinvreport^ rep;
	alglib::rmatrixinverse(hess, info, rep);
	array<double>^ errs = gcnew array<double>(p->Length);
	for (int i = 0; i < p->Length; i++)
		errs[i] = Math::Sqrt(hess[i, i]);

	return errs;
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

void JPFITS::JPMath::Fit_WCSTransform2d(array<double>^ x_intrmdt, array<double>^ y_intrmdt, array<double>^ x_pix, array<double>^ y_pix, array<double>^ &p, array<double>^ p_lbnd, array<double>^ p_ubnd, array<double>^ p_scale)
{
	double epsx = 0.000000001;
	int maxits = 1000;
	alglib::minlmstate^ state;
	alglib::minlmreport^ report;
	alglib::ndimensional_fvec^ fvec = gcnew alglib::ndimensional_fvec(alglib_WCSTransform2d_fvec);
	alglib::ndimensional_jac^ jac = gcnew alglib::ndimensional_jac(alglib_WCSTransform2d_jac);
	array<Object^>^ objj = gcnew array<Object^>(4);
	objj[0] = (Object^)x_pix;
	objj[1] = (Object^)y_pix;
	objj[2] = (Object^)x_intrmdt;
	objj[3] = (Object^)y_intrmdt;
	Object^ obj = (Object^)objj;

	alglib::minlmcreatevj(y_pix->Length, p, state);
	alglib::minlmsetcond(state, epsx, maxits);
	alglib::minlmsetscale(state, p_scale);
	alglib::minlmsetbc(state, p_lbnd, p_ubnd);
	alglib::minlmoptimize(state, fvec, jac, nullptr, obj);
	alglib::minlmresults(state, p, report);
}

void JPFITS::JPMath::alglib_WCSTransform2d_fvec(array<double>^ p, array<double>^ f, Object^ obj)
{
	array<Object^>^ objj = (array<Object^>^)obj;
	array<double>^ x_pix = (array<double>^)objj[0];
	array<double>^ y_pix = (array<double>^)objj[1];
	array<double>^ x_intrmdt = (array<double>^)objj[2];
	array<double>^ y_intrmdt = (array<double>^)objj[3];

	double xres, yres;

	if (p->Length == 4)
	{
		double cosp1 = Math::Cos(p[1]);
		double sinp1 = Math::Sin(p[1]);
		double xpiximp2, ypiximp3;

		for (int i = 0; i < x_intrmdt->Length; i++)
		{
			xpiximp2 = x_pix[i] - p[2];
			ypiximp3 = y_pix[i] - p[3];

			xres = p[0] * (cosp1 * xpiximp2 - sinp1 * ypiximp3) - x_intrmdt[i];
			yres = p[0] * (sinp1 * xpiximp2 + cosp1 * ypiximp3) - y_intrmdt[i];

			f[i] = Math::Sqrt(xres * xres + yres * yres);
		}
		return;
	}

	if (p->Length == 6)
	{
		double xpiximp4, ypiximp5;

		for (int i = 0; i < x_pix->Length; i++)
		{
			xpiximp4 = x_pix[i] - p[4];
			ypiximp5 = y_pix[i] - p[5];

			xres = p[0] * xpiximp4 + p[1] * ypiximp5 - x_intrmdt[i];
			yres = p[2] * xpiximp4 + p[3] * ypiximp5 - y_intrmdt[i];

			f[i] = Math::Sqrt(xres * xres + yres * yres);
		}
	}
}

void JPFITS::JPMath::alglib_WCSTransform2d_jac(array<double>^ p, array<double>^ f, array<double, 2>^ jac, Object^ obj)
{
	array<Object^>^ objj = (array<Object^>^)obj;
	array<double>^ x_pix = (array<double>^)objj[0];
	array<double>^ y_pix = (array<double>^)objj[1];
	array<double>^ x_intrmdt = (array<double>^)objj[2];
	array<double>^ y_intrmdt = (array<double>^)objj[3];
	
	alglib_WCSTransform2d_fvec(p, f, obj);

	if (p->Length == 4)
	{
		double cosp1 = Math::Cos(p[1]);
		double sinp1 = Math::Sin(p[1]);
		double p2mxpixi, p3mypixi, cosp1p2mxpixisinp1p3mypixi, twox_intrmdtip0cosp1p2mxpixisinp1p3mypixisqy_intrmdtip0cosp1p3mypixisinp1p2mxpixisqroot, x_intrmdtip0cosp1p2mxpixisinp1p3mypixisqy_intrmdtip0cosp1p3mypixisinp1p2mxpixisq, cosp1p3mypixisinp1p2mxpixi, x_intrmdtip0cosp1p2mxpixisinp1p3mypixi, y_intrmdtip0cosp1p3mypixisinp1p2mxpixi, x_intrmdtip0cosp1p2mxpixisinp1p3mypixisq, y_intrmdtip0cosp1p3mypixisinp1p2mxpixisq;

		for (int i = 0; i < x_intrmdt->Length; i++)
		{
			p2mxpixi = p[2] - x_pix[i];
			p3mypixi = p[3] - y_pix[i];
			cosp1p2mxpixisinp1p3mypixi = cosp1 * p2mxpixi - sinp1 * p3mypixi;
			cosp1p3mypixisinp1p2mxpixi = cosp1 * p3mypixi + sinp1 * p2mxpixi;
			x_intrmdtip0cosp1p2mxpixisinp1p3mypixi = x_intrmdt[i] + p[0] * cosp1p2mxpixisinp1p3mypixi;
			y_intrmdtip0cosp1p3mypixisinp1p2mxpixi = y_intrmdt[i] + p[0] * cosp1p3mypixisinp1p2mxpixi;
			x_intrmdtip0cosp1p2mxpixisinp1p3mypixisq = x_intrmdtip0cosp1p2mxpixisinp1p3mypixi * x_intrmdtip0cosp1p2mxpixisinp1p3mypixi;
			y_intrmdtip0cosp1p3mypixisinp1p2mxpixisq = y_intrmdtip0cosp1p3mypixisinp1p2mxpixi * y_intrmdtip0cosp1p3mypixisinp1p2mxpixi;
			x_intrmdtip0cosp1p2mxpixisinp1p3mypixisqy_intrmdtip0cosp1p3mypixisinp1p2mxpixisq = x_intrmdtip0cosp1p2mxpixisinp1p3mypixisq + y_intrmdtip0cosp1p3mypixisinp1p2mxpixisq;
			twox_intrmdtip0cosp1p2mxpixisinp1p3mypixisqy_intrmdtip0cosp1p3mypixisinp1p2mxpixisqroot = 2 * Math::Sqrt(x_intrmdtip0cosp1p2mxpixisinp1p3mypixisqy_intrmdtip0cosp1p3mypixisinp1p2mxpixisq);

			jac[i, 0] = (2 * cosp1p2mxpixisinp1p3mypixi*x_intrmdtip0cosp1p2mxpixisinp1p3mypixi + 2 * cosp1p3mypixisinp1p2mxpixi*y_intrmdtip0cosp1p3mypixisinp1p2mxpixi) / twox_intrmdtip0cosp1p2mxpixisinp1p3mypixisqy_intrmdtip0cosp1p3mypixisinp1p2mxpixisqroot;
			jac[i, 1] = -(2 * p[0] * cosp1p3mypixisinp1p2mxpixi*x_intrmdtip0cosp1p2mxpixisinp1p3mypixi - 2 * p[0] * cosp1p2mxpixisinp1p3mypixi*y_intrmdtip0cosp1p3mypixisinp1p2mxpixi) / twox_intrmdtip0cosp1p2mxpixisinp1p3mypixisqy_intrmdtip0cosp1p3mypixisinp1p2mxpixisqroot;
			jac[i, 2] = (2 * p[0] * cosp1*x_intrmdtip0cosp1p2mxpixisinp1p3mypixi + 2 * p[0] * sinp1*y_intrmdtip0cosp1p3mypixisinp1p2mxpixi) / twox_intrmdtip0cosp1p2mxpixisinp1p3mypixisqy_intrmdtip0cosp1p3mypixisinp1p2mxpixisqroot;
			jac[i, 3] = (2 * p[0] * cosp1*y_intrmdtip0cosp1p3mypixisinp1p2mxpixi - 2 * p[0] * sinp1*x_intrmdtip0cosp1p2mxpixisinp1p3mypixi) / twox_intrmdtip0cosp1p2mxpixisinp1p3mypixisqy_intrmdtip0cosp1p3mypixisinp1p2mxpixisqroot;
		}
		return;
	}

	if (p->Length == 6)
	{
		double p4mx_pixi, p5my_pixi, p0p4mx_pixi, rootx_intrmdtip0p4mx_pixip1p5my_pixisqy_intrmdtip2p4mx_pixip3p5my_pixisq, x_intrmdtip0p4mx_pixip1p5my_pixisqy_intrmdtip2p4mx_pixip3p5my_pixisq, p1p5my_pixi, p2p4mx_pixi, p3p5my_pixi, x_intrmdtip0p4mx_pixip1p5my_pixi, y_intrmdtip2p4mx_pixip3p5my_pixi, x_intrmdtip0p4mx_pixip1p5my_pixisq, y_intrmdtip2p4mx_pixip3p5my_pixisq;

		for (int i = 0; i < x_pix->Length; i++)
		{
			p4mx_pixi = p[4] - x_pix[i];
			p5my_pixi = p[5] - y_pix[i];
			p0p4mx_pixi = p[0] * p4mx_pixi;
			p1p5my_pixi = p[1] * p5my_pixi;
			p2p4mx_pixi = p[2] * p4mx_pixi;
			p3p5my_pixi = p[3] * p5my_pixi;
			x_intrmdtip0p4mx_pixip1p5my_pixi = x_intrmdt[i] + p0p4mx_pixi + p1p5my_pixi;
			y_intrmdtip2p4mx_pixip3p5my_pixi = y_intrmdt[i] + p2p4mx_pixi + p3p5my_pixi;
			x_intrmdtip0p4mx_pixip1p5my_pixisq = x_intrmdtip0p4mx_pixip1p5my_pixi * x_intrmdtip0p4mx_pixip1p5my_pixi;
			y_intrmdtip2p4mx_pixip3p5my_pixisq = y_intrmdtip2p4mx_pixip3p5my_pixi * y_intrmdtip2p4mx_pixip3p5my_pixi;
			x_intrmdtip0p4mx_pixip1p5my_pixisqy_intrmdtip2p4mx_pixip3p5my_pixisq = x_intrmdtip0p4mx_pixip1p5my_pixisq + y_intrmdtip2p4mx_pixip3p5my_pixisq;
			rootx_intrmdtip0p4mx_pixip1p5my_pixisqy_intrmdtip2p4mx_pixip3p5my_pixisq = Math::Sqrt(x_intrmdtip0p4mx_pixip1p5my_pixisqy_intrmdtip2p4mx_pixip3p5my_pixisq);

			jac[i, 0] = (p4mx_pixi*x_intrmdtip0p4mx_pixip1p5my_pixi) / rootx_intrmdtip0p4mx_pixip1p5my_pixisqy_intrmdtip2p4mx_pixip3p5my_pixisq;
			jac[i, 1] = (p5my_pixi*x_intrmdtip0p4mx_pixip1p5my_pixi) / rootx_intrmdtip0p4mx_pixip1p5my_pixisqy_intrmdtip2p4mx_pixip3p5my_pixisq;
			jac[i, 2] = (p4mx_pixi*y_intrmdtip2p4mx_pixip3p5my_pixi) / rootx_intrmdtip0p4mx_pixip1p5my_pixisqy_intrmdtip2p4mx_pixip3p5my_pixisq;
			jac[i, 3] = (p5my_pixi*y_intrmdtip2p4mx_pixip3p5my_pixi) / rootx_intrmdtip0p4mx_pixip1p5my_pixisqy_intrmdtip2p4mx_pixip3p5my_pixisq;
			jac[i, 4] = (2 * p[0]*x_intrmdtip0p4mx_pixip1p5my_pixi + 2 * p[2]*y_intrmdtip2p4mx_pixip3p5my_pixi) / (2 * rootx_intrmdtip0p4mx_pixip1p5my_pixisqy_intrmdtip2p4mx_pixip3p5my_pixisq);
			jac[i, 5] = (2 * p[1]*x_intrmdtip0p4mx_pixip1p5my_pixi + 2 * p[3]*y_intrmdtip2p4mx_pixip3p5my_pixi) / (2 * rootx_intrmdtip0p4mx_pixip1p5my_pixisqy_intrmdtip2p4mx_pixip3p5my_pixisq);
		}
	}
}

void JPFITS::JPMath::Fit_GeneralTransform2d(array<double>^ x_ref, array<double>^ y_ref, array<double>^ x_tran, array<double>^ y_tran, array<double>^ &p, array<double>^ p_lbnd, array<double>^ p_ubnd, array<double>^ p_scale)
{
	double epsx = 0.00000000001;
	int maxits = 0;
	alglib::minlmstate^ state;
	alglib::minlmreport^ report;
	alglib::ndimensional_fvec^ fvec = gcnew alglib::ndimensional_fvec(alglib_GeneralTransform2d_fvec);
	array<Object^>^ objj = gcnew array<Object^>(4);
	objj[0] = (Object^)x_tran;
	objj[1] = (Object^)y_tran;
	objj[2] = (Object^)x_ref;
	objj[3] = (Object^)y_ref;
	Object^ obj = (Object^)objj;

	alglib::minlmcreatev(x_ref->Length, p, 0.000001, state);
	alglib::minlmsetcond(state, epsx, maxits);
	alglib::minlmsetscale(state, p_scale);
	alglib::minlmsetbc(state, p_lbnd, p_ubnd);
	alglib::minlmoptimize(state, fvec, nullptr, obj);
	alglib::minlmresults(state, p, report);
}

void JPFITS::JPMath::alglib_GeneralTransform2d_fvec(array<double>^ p, array<double>^ f, Object^ obj)
{
	array<Object^>^ objj = (array<Object^>^)obj;
	array<double>^ x_tran = (array<double>^)objj[0];
	array<double>^ y_tran = (array<double>^)objj[1];
	array<double>^ x_ref = (array<double>^)objj[2];
	array<double>^ y_ref = (array<double>^)objj[3];

	if (p->Length == 6)
	{
		double xres, yres;

		for (int i = 0; i < x_ref->Length; i++)
		{
			xres = p[0] * (Math::Cos(p[1]) * (x_tran[i] - p[2]) - Math::Sin(p[1]) * (y_tran[i] - p[3])) + p[2] + p[4] - x_ref[i];
			yres = p[0] * (Math::Sin(p[1]) * (x_tran[i] - p[2]) + Math::Cos(p[1]) * (y_tran[i] - p[3])) + p[3] + p[5] - y_ref[i];

			f[i] = Math::Sqrt(xres * xres + yres * yres);
		}
		return;
	}

	if (p->Length == 8)
	{
		double xres, yres;

		for (int i = 0; i < x_tran->Length; i++)
		{
			xres = p[0] * (x_tran[i] - p[4]) + p[1] * (y_tran[i] - p[5]) + p[4] + p[6] - x_ref[i];
			yres = p[2] * (x_tran[i] - p[4]) + p[3] * (y_tran[i] - p[5]) + p[5] + p[7] - y_ref[i];

			f[i] = Math::Sqrt(xres * xres + yres * yres);
		}
	}
}

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

