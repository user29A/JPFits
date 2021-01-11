/*Copyright 2017 Joseph Edwin Postma

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

You should have received a copy of the GNU General Public License
along with JPFITS. If not, see http://www.gnu.org/licenses/. */


#include "stdafx.h"
#include "JPFITS.h"

using namespace JPFITS;


array<double,2>^ JPFITS::FITSImage::operator +(FITSImage^ lhs_img, FITSImage^ rhs_img)
{
	if (lhs_img->Width != rhs_img->Width || lhs_img->Height != rhs_img->Height)
	{
		throw gcnew System::ArrayTypeMismatchException("Image Data Matrices not the Same Size...Discontinuing.");
	}
	int W = lhs_img->Width;
	int H = lhs_img->Height;
	array<double,2>^ result = gcnew array<double,2>(W,H);
	#pragma omp parallel for
	for (int i = 0; i < W; i++)
		for (int j = 0; j < H; j++)
			result[i,j] = lhs_img[i,j] + rhs_img[i,j];

	return result;
}

array<double,2>^ JPFITS::FITSImage::operator +(FITSImage^ lhs_img, double scalar)
{
	int W = lhs_img->Width;
	int H = lhs_img->Height;
	array<double,2>^ result = gcnew array<double,2>(W,H);
	#pragma omp parallel for
	for (int i=0; i<W; i++)
		for (int j=0; j<H; j++)
			result[i,j] = lhs_img[i,j] + scalar;

	return result;
}

array<double,2>^ JPFITS::FITSImage::operator -(FITSImage^ lhs_img, FITSImage^ rhs_img)
{
	if (lhs_img->Width != rhs_img->Width || lhs_img->Height != rhs_img->Height)
	{
		throw gcnew System::ArrayTypeMismatchException("Image Data Matrices not the Same Size...Discontinuing.");
	}
	int W = lhs_img->Width;
	int H = lhs_img->Height;
	array<double,2>^ result = gcnew array<double,2>(W,H);
	#pragma omp parallel for
	for (int i=0; i<W; i++)
		for (int j=0; j<H; j++)
			result[i,j] = lhs_img[i,j] - rhs_img[i,j];

	return result;
}

array<double,2>^ JPFITS::FITSImage::operator -(FITSImage^ lhs_img,double scalar)
{
	int W = lhs_img->Width;
	int H = lhs_img->Height;
	array<double,2>^ result = gcnew array<double,2>(W,H);
	#pragma omp parallel for
	for (int i=0; i<W; i++)
		for (int j=0; j<H; j++)
			result[i,j] = lhs_img[i,j]-scalar;
	return result;
}

array<double,2>^ JPFITS::FITSImage::operator /(FITSImage^ lhs_img,FITSImage^ rhs_img)
{
	if (lhs_img->Width != rhs_img->Width || lhs_img->Height != rhs_img->Height)
	{
		throw gcnew System::ArrayTypeMismatchException("Image Data Matrices not the Same Size...Discontinuing.");
	}
	int W = lhs_img->Width;
	int H = lhs_img->Height;
	array<double,2>^ result = gcnew array<double,2>(W,H);
	#pragma omp parallel for
	for (int i=0; i<W; i++)
		for (int j=0; j<H; j++)
			result[i,j] = lhs_img[i,j]/rhs_img[i,j];

	return result;
}

array<double,2>^ JPFITS::FITSImage::operator /(FITSImage^ lhs_img,double scalar)
{
	int W = lhs_img->Width;
	int H = lhs_img->Height;
	array<double,2>^ result = gcnew array<double,2>(W,H);
	scalar = 1/scalar;
	#pragma omp parallel for
	for (int i=0; i<W; i++)
		for (int j=0; j<H; j++)
			result[i,j] = lhs_img[i,j]*scalar;

	return result;
}

array<double,2>^ JPFITS::FITSImage::operator *(FITSImage^ lhs_img,FITSImage^ rhs_img)
{
	if (lhs_img->Width != rhs_img->Width || lhs_img->Height != rhs_img->Height)
	{
		throw gcnew System::ArrayTypeMismatchException("Image Data Matrices not the Same Size...Discontinuing.");
	}
	int W = lhs_img->Width;
	int H = lhs_img->Height;
	array<double,2>^ result = gcnew array<double,2>(W,H);
	#pragma omp parallel for
	for (int i=0; i<W; i++)
		for (int j=0; j<H; j++)
			result[i,j] = lhs_img[i,j]*rhs_img[i,j];

	return result;
}

array<double,2>^ JPFITS::FITSImage::operator *(FITSImage^ lhs_img,double scalar)
{
	int W = lhs_img->Width;
	int H = lhs_img->Height;
	array<double,2>^ result = gcnew array<double,2>(W,H);
	#pragma omp parallel for
	for (int i=0; i<W; i++)
		for (int j=0; j<H; j++)
			result[i,j] = lhs_img[i,j]*scalar;

	return result;
}

array<double,2>^ JPFITS::FITSImage::operator ^(FITSImage^ lhs_img,FITSImage^ rhs_img)
{
	if (lhs_img->Width != rhs_img->Width || lhs_img->Height != rhs_img->Height)
	{
		throw gcnew System::ArrayTypeMismatchException("Image Data Matrices not the Same Size...Discontinuing.");
	}
	int W = lhs_img->Width;
	int H = lhs_img->Height;
	array<double,2>^ result = gcnew array<double,2>(W,H);
	#pragma omp parallel for
	for (int i=0; i<W; i++)
		for (int j=0; j<H; j++)
			result[i,j] = Math::Pow(lhs_img[i,j],rhs_img[i,j]);

	return result;
}

array<double,2>^ JPFITS::FITSImage::operator ^(FITSImage^ lhs_img,double scalar)
{
	int W = lhs_img->Width;
	int H = lhs_img->Height;
	array<double,2>^ result = gcnew array<double,2>(W,H);
	#pragma omp parallel for
	for (int i=0; i<W; i++)
		for (int j=0; j<H; j++)
			result[i,j] = Math::Pow(lhs_img[i,j],scalar);

	return result;
}

void JPFITS::FITSImage::StatsUpD(bool do_parallel)
{
	int L = NAXIS1*NAXIS2;
	double N = double(L);
	MIN = ::Double::MaxValue;
	MAX = ::Double::MinValue;
	MEDIAN = 0;
	MEAN = 0;
	STD = 0;
	SUM = 0;

	#pragma omp parallel sections if(do_parallel)
	{
		#pragma omp section
		{
			double sum = 0;
			//#pragma omp parallel for reduction(+:sum)
			for (int i=0; i < NAXIS1; i++)
				for (int j=0; j < NAXIS2; j++)
					{
						sum = sum + DIMAGE[i,j];
						if (MIN > DIMAGE[i,j])
							MIN = DIMAGE[i,j];
						if (MAX < DIMAGE[i,j])
							MAX = DIMAGE[i,j];
					}
			SUM = sum;
			MEAN = SUM/N;
		}

		#pragma omp section
		{
			MEDIAN = JPMath::Median(DIMAGE);
		}
	}

	double std = 0.0;
	#pragma omp parallel for if (do_parallel) reduction(+:std)
		for (int i=0; i < NAXIS1; i++)
			for (int j=0; j < NAXIS2; j++)
				std += (DIMAGE[i,j]-MEAN)*(DIMAGE[i,j]-MEAN);
	STD = Math::Sqrt(std/(N-1.0));

	/*double sum = 0;
	#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < NAXIS1; i++)
	{
		for (int j = 0; j < NAXIS2; j++)
		{
			sum += DIMAGE[i, j];

			if (DIMAGE[i, j] < MIN)
			{
				#pragma omp critical  
				{
					if (DIMAGE[i, j] < MIN)
						MIN = DIMAGE[i, j];
				}
			}
			if (DIMAGE[i, j] > MAX)
			{
				#pragma omp critical  
				{
					if (DIMAGE[i, j] > MAX)
						MAX = DIMAGE[i, j];
				}
			}
		}
	}
	SUM = sum;
	MEAN = SUM / N;		

	#pragma omp parallel sections
	{
		#pragma omp section
		{
			MEDIAN = JPMath::Median(DIMAGE);
		}

		#pragma omp section
		{
			double std = 0;
			#pragma omp parallel for reduction(+:std)
			for (int i = 0; i < NAXIS1; i++)
				for (int j = 0; j < NAXIS2; j++)
					std += (DIMAGE[i, j] - MEAN)*(DIMAGE[i, j] - MEAN);
			STD = Math::Sqrt(std / (N - 1.0));
		}
	}*/
	
}

void JPFITS::FITSImage::SETBITPIX(System::TypeCode Precision)
{
	if (DIMAGE == nullptr || DIMAGE->Length == 0)
	{
		BITPIX = 8;
		NAXIS = 0;
		NAXIS1 = 0;
		NAXIS2 = 0;
		SetKey("BITPIX", BITPIX.ToString(), false, 1);
		SetKey("NAXIS", NAXIS.ToString(), false, 2);
		return;
	}

	switch (Precision)
	{
		case TypeCode::SByte:
		{
			BITPIX = 8;
			BZERO = 0;
			BSCALE = 1;
			break;
		}

		case TypeCode::Byte:
		{
			BITPIX = 8;
			BZERO = 128;
			BSCALE = 1;
			break;
		}

		case TypeCode::Int16:
		{
			BITPIX = 16;
			BZERO = 0;
			BSCALE = 1;
			break;
		}

		case TypeCode::UInt16:
		{
			BITPIX = 16;
			BZERO = 32768;
			BSCALE = 1;
			break;
		}

		case TypeCode::Int32:
		{
			BITPIX = 32;
			BZERO = 0;
			BSCALE = 1;
			break;
		}

		case TypeCode::UInt32:
		{
			BITPIX = 32;
			BZERO = 2147483648;
			BSCALE = 1;
			break;
		}

		case TypeCode::Int64:
		{
			BITPIX = 64;
			BZERO = 0;
			BSCALE = 1;
			break;
		}

		case TypeCode::Double:
		{
			BITPIX = -64;
			BZERO = 0;
			BSCALE = 1;
			break;
		}

		default:
		{
			BITPIX = 0;
			BZERO = 0;
			BSCALE = 0;
			break;
		}
	}

	SetKey("BITPIX", BITPIX.ToString(), false, 0);
	SetKey("BZERO", BZERO.ToString(), "Data Offset; pixel = pixel*BSCALE+BZERO", true, 4);
	SetKey("BSCALE", BSCALE.ToString(), "Data Scaling; pixel = pixel*BSCALE+BZERO", true, 5);
}

void JPFITS::FITSImage::RotateCW(bool CW)
{
	array<double,2>^ rotimg = gcnew array<double,2>(NAXIS2,NAXIS1);
	if (CW)
	{
		#pragma omp parallel for
		for (int i = 0; i < int((NAXIS1+1)/2); i++)
		{
			double dum1, dum2, dum3, dum4;
			for (int j = 0; j < int((NAXIS2+1)/2); j++)
			{
				dum1 = DIMAGE[i,j];
				dum2 = DIMAGE[i,NAXIS2-1-j];
				dum3 = DIMAGE[NAXIS1-1-i,NAXIS2-1-j];
				dum4 = DIMAGE[NAXIS1-1-i,j];

				rotimg[j,i] = dum2;
				rotimg[NAXIS2-1-j,i] = dum1;
				rotimg[NAXIS2-1-j,NAXIS1-1-i] = dum4;
				rotimg[j,NAXIS1-1-i] = dum3;
			}
		}
		this->AddKey("ROTATN","-90","Image Rotated CW 90",-1);
	}
	else
	{
		#pragma omp parallel for
		for (int i = 0; i < int((NAXIS1+1)/2); i++)
		{
			double dum1, dum2, dum3, dum4;
			for (int j = 0; j < int((NAXIS2+1)/2); j++)
			{
				dum1 = DIMAGE[i,j];
				dum2 = DIMAGE[i,NAXIS2-1-j];
				dum3 = DIMAGE[NAXIS1-1-i,NAXIS2-1-j];
				dum4 = DIMAGE[NAXIS1-1-i,j];

				rotimg[j,i] = dum4;
				rotimg[NAXIS2-1-j,i] = dum3;
				rotimg[NAXIS2-1-j,NAXIS1-1-i] = dum2;
				rotimg[j,NAXIS1-1-i] = dum1;
			}
		}
		this->AddKey("ROTATN","90","Image Rotated CCW 90",-1);
	}

	DIMAGE = rotimg;
	delete rotimg;
	this->SetKey("NAXIS1",NAXIS2.ToString(),0,1);
	this->SetKey("NAXIS2",NAXIS1.ToString(),0,2);
	int dumn = NAXIS1;
	NAXIS1 = NAXIS2;
	NAXIS2 = dumn;
}

void JPFITS::FITSImage::FlipVertical()
{
	#pragma omp parallel for
	for (int i = 0; i < NAXIS1; i++)
	{
		double dum1, dum2;
		for (int j = 0; j < int(NAXIS2/2); j++)
		{
			dum1 = DIMAGE[i,j];
			dum2 = DIMAGE[i,NAXIS2-1-j];
			DIMAGE[i,j] = dum2;
			DIMAGE[i,NAXIS2-1-j] = dum1;
		}
	}
	this->AddKey("VFLIP","true","Image Vertically Flipped",-1);
}

void JPFITS::FITSImage::FlipHorizontal()
{
	#pragma omp parallel for
	for (int j = 0; j < NAXIS2; j++)
	{
		double dum1, dum2;
		for (int i = 0; i < int(NAXIS1/2); i++)
		{
			dum1 = DIMAGE[i,j];
			dum2 = DIMAGE[NAXIS1-1-i,j];
			DIMAGE[i,j] = dum2;
			DIMAGE[NAXIS1-1-i,j] = dum1;
		}
	}
	this->AddKey("HFLIP","true","Image Horizontally Flipped",-1);
}

