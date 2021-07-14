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

using namespace JPFITS;


JPFITS::WorldCoordinateSolution::~WorldCoordinateSolution()
{

}

JPFITS::WorldCoordinateSolution::WorldCoordinateSolution()
{

}

JPFITS::WorldCoordinateSolution::WorldCoordinateSolution(JPFITS::FITSImageHeader^ header)
{
	EATHEADERFORWCS(header);
}

void JPFITS::WorldCoordinateSolution::EATHEADERFORWCS(JPFITS::FITSImageHeader^ header)
{
	CD1_1 = (double)header->GetKeyIndex("CD1_1", false);
	if (CD1_1 < 0)
	{
		WCSEXISTS = false;
		return;
	}
	CD1_1 = Convert::ToDouble(header->HeaderKeyValues[(int)Math::Round(CD1_1)]);

	CTYPEN = gcnew array<String^>(2);
	CRPIXN = gcnew array<double>(2);
	CRVALN = gcnew array<double>(2);
	
	#pragma omp parallel sections
	{
		#pragma omp section
		CTYPEN[0] = header->GetKeyValue("CTYPE1");
		#pragma omp section
		CTYPEN[1] = header->GetKeyValue("CTYPE2");
		#pragma omp section
		CD1_2 = Convert::ToDouble(header->GetKeyValue("CD1_2"));
		#pragma omp section
		CD2_1 = Convert::ToDouble(header->GetKeyValue("CD2_1"));
		#pragma omp section
		CD2_2 = Convert::ToDouble(header->GetKeyValue("CD2_2"));
		#pragma omp section
		CRPIXN[0] = Convert::ToDouble(header->GetKeyValue("CRPIX1"));
		#pragma omp section
		CRPIXN[1] = Convert::ToDouble(header->GetKeyValue("CRPIX2"));
		#pragma omp section
		CRVALN[0] = Convert::ToDouble(header->GetKeyValue("CRVAL1"));
		#pragma omp section
		CRVALN[1] = Convert::ToDouble(header->GetKeyValue("CRVAL2"));
	}

	CDMATRIX = gcnew array<double, 2>(2, 2);
	CDMATRIXINV = gcnew array<double, 2>(2, 2);
	CDMATRIX[0, 0] = CD1_1;
	CDMATRIX[1, 0] = CD1_2;
	CDMATRIX[0, 1] = CD2_1;
	CDMATRIX[1, 1] = CD2_2;
	SET_CDMATRIXINV();

	CDELTN = gcnew array<double>(2);
	CDELTN[0] = Math::Sqrt(CD1_1 * CD1_1 + CD1_2 * CD1_2) * 3600;
	CDELTN[1] = Math::Sqrt(CD2_1 * CD2_1 + CD2_2 * CD2_2) * 3600;

	CROTAN = gcnew array<double>(2);
	CROTAN[0] = Math::Atan2(CD1_2, -CD1_1) * 180 / Math::PI;
	CROTAN[1] = Math::Atan2(-CD2_1, -CD2_2) * 180 / Math::PI;

	WCSEXISTS = true;

	//optionally populate this?
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CCVALD1", false);
			if (ind != -1)
				CCVALD1 = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CCVALD2", false);
			if (ind != -1)
				CCVALD2 = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CCVALS1", false);
			if (ind != -1)
				CCVALS1 = header->GetKeyValue("CCVALS1");
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CCVALS2", false);
			if (ind != -1)
				CCVALS2 = header->GetKeyValue("CCVALS2");
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CPIX1RM", false);
			if (ind != -1)
				CPIX1RM = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CPIX1RS", false);
			if (ind != -1)
				CPIX1RS = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CVAL1RM", false);
			if (ind != -1)
				CVAL1RM = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CVAL1RS", false);
			if (ind != -1)
				CVAL1RS = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CPIX2RM", false);
			if (ind != -1)
				CPIX2RM = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CPIX2RS", false);
			if (ind != -1)
				CPIX2RS = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CVAL2RM", false);
			if (ind != -1)
				CVAL2RM = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CVAL2RS", false);
			if (ind != -1)
				CVAL2RS = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CPIXRM", false);
			if (ind != -1)
				CPIXRM = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CPIXRS", false);
			if (ind != -1)
				CPIXRS = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CVALRM", false);
			if (ind != -1)
				CVALRM = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
		#pragma omp section
		{
			int ind = header->GetKeyIndex("CVALRS", false);
			if (ind != -1)
				CVALRS = Convert::ToDouble(header->HeaderKeyValues[ind]);
		}
	}

	int num = 0, key = 0;
	while (key != -1)
	{
		num++;
		key = header->GetKeyIndex("WCP1_" + num.ToString("000"), false);
	}
	num--;
	CPIX1 = gcnew array<double>(num);
	CPIX2 = gcnew array<double>(num);
	CVAL1 = gcnew array<double>(num);
	CVAL2 = gcnew array<double>(num);
	DVAL1 = gcnew array<double>(num);
	DVAL2 = gcnew array<double>(num);

	#pragma omp parallel for
	for (int i = 1; i <= CPIX1->Length; i++)
	{
		int ind = header->GetKeyIndex("WCP1_" + i.ToString("000"), false);
		if (ind != -1)
			CPIX1[i - 1] = Convert::ToDouble(header->HeaderKeyValues[ind]);

		ind = header->GetKeyIndex("WCP2_" + i.ToString("000"), false);
		if (ind != -1)
			CPIX2[i - 1] = Convert::ToDouble(header->HeaderKeyValues[ind]);

		ind = header->GetKeyIndex("WCV1_" + i.ToString("000"), false);
		if (ind != -1)
			CVAL1[i - 1] = Convert::ToDouble(header->HeaderKeyValues[ind]);

		ind = header->GetKeyIndex("WCV2_" + i.ToString("000"), false);
		if (ind != -1)
			CVAL2[i - 1] = Convert::ToDouble(header->HeaderKeyValues[ind]);
		
		ind = header->GetKeyIndex("WCD1_" + i.ToString("000"), false);
		if (ind != -1)
			DVAL1[i - 1] = Convert::ToDouble(header->HeaderKeyValues[ind]);

		ind = header->GetKeyIndex("WCD2_" + i.ToString("000"), false);
		if (ind != -1)
			DVAL2[i - 1] = Convert::ToDouble(header->HeaderKeyValues[ind]);
		
		/*CPIX1[i - 1] = Convert::ToDouble(header->GetKeyValue("WCP1_" + i.ToString("000")));
		CPIX2[i - 1] = Convert::ToDouble(header->GetKeyValue("WCP2_" + i.ToString("000")));
		CVAL1[i - 1] = Convert::ToDouble(header->GetKeyValue("WCV1_" + i.ToString("000")));
		CVAL2[i - 1] = Convert::ToDouble(header->GetKeyValue("WCV2_" + i.ToString("000")));
		DVAL1[i - 1] = Convert::ToDouble(header->GetKeyValue("WCD1_" + i.ToString("000")));
		DVAL2[i - 1] = Convert::ToDouble(header->GetKeyValue("WCD2_" + i.ToString("000")));*/
	}
}

void JPFITS::WorldCoordinateSolution::Solve_WCS(String^ WCS_Type, array<double>^ X_pix, array<double>^ Y_pix, bool zero_based_pixels, array<double>^ cval1, array<double>^ cval2, JPFITS::FITSImageHeader^ header)
{
	//should first do a check of WCS type to make sure it is valid
	CTYPEN = gcnew array<String^>(2);
	CTYPEN[0] = "RA---" + WCS_Type;
	CTYPEN[1] = "DEC--" + WCS_Type;

	CPIX1 = gcnew array<double>(X_pix->Length);
	CPIX2 = gcnew array<double>(X_pix->Length);
	CVAL1 = gcnew array<double>(X_pix->Length);
	CVAL2 = gcnew array<double>(X_pix->Length);
	DVAL1 = gcnew array<double>(X_pix->Length);
	DVAL2 = gcnew array<double>(X_pix->Length);
	for (int i = 0; i < X_pix->Length; i++)
	{
		CPIX1[i] = X_pix[i];
		CPIX2[i] = Y_pix[i];
		if (zero_based_pixels)
		{
			CPIX1[i]++;
			CPIX2[i]++;
		}
		CVAL1[i] = cval1[i];
		CVAL2[i] = cval2[i];
	}
	CRPIXN = gcnew array<double>(2);
	CRVALN = gcnew array<double>(2);
	CRPIXN[0] = JPMath::Mean(CPIX1, true);//fix this in the fit boundaries? - NO...let it get the best one...but they should be close
	CRPIXN[1] = JPMath::Mean(CPIX2, true);//fix this in the fit boundaries? - NO...let it get the best one...but they should be close
	CRVALN[0] = JPMath::Mean(CVAL1, true);//these are fixed as the coordinate reference value
	CRVALN[1] = JPMath::Mean(CVAL2, true);//these are fixed as the coordinate reference value

	array<double>^ X_intrmdt = gcnew array<double>(CPIX1->Length);//intermediate coords (degrees)
	array<double>^ Y_intrmdt = gcnew array<double>(CPIX1->Length);//intermediate coords (degrees)
	double a0 = CRVALN[0] * Math::PI / 180, d0 = CRVALN[1] * Math::PI / 180, a, d;
	for (int i = 0; i < CPIX1->Length; i++)
	{
		a = CVAL1[i] * Math::PI / 180;//radians
		d = CVAL2[i] * Math::PI / 180;//radians

		//for tangent plane Gnomic
		if (WCS_Type == "TAN")
		{
			X_intrmdt[i] = Math::Cos(d)*Math::Sin(a - a0) / (Math::Cos(d0)*Math::Cos(d)*Math::Cos(a - a0) + Math::Sin(d0)*Math::Sin(d));
			Y_intrmdt[i] = (Math::Cos(d0)*Math::Sin(d) - Math::Cos(d)*Math::Sin(d0)*Math::Cos(a - a0)) / (Math::Cos(d0)*Math::Cos(d)*Math::Cos(a - a0) + Math::Sin(d0)*Math::Sin(d));
		}
	}

	array<double>^ P0 = gcnew array<double>(6) { 0, 0, 0, 0, CRPIXN[0], CRPIXN[1] };
	array<double>^ plb = gcnew array<double>(6) { -0.1, -0.1, -0.1, -0.1, JPMath::Min(CPIX1, false), JPMath::Min(CPIX2, false) };
	array<double>^ pub = gcnew array<double>(6) { 0.1, 0.1, 0.1, 0.1, JPMath::Max(CPIX1, false), JPMath::Max(CPIX2, false) };
	array<double>^ scale = gcnew array<double>(6) { 1e-6, 1e-6, 1e-6, 1e-6, JPMath::Mean(CPIX1, true), JPMath::Mean(CPIX2, true) };
	JPMath::Fit_WCSTransform2d(X_intrmdt, Y_intrmdt, CPIX1, CPIX2, P0, plb, pub, scale);

	CDMATRIX = gcnew array<double, 2>(2, 2);
	CDMATRIX[0, 0] = P0[0] * 180 / Math::PI;
	CDMATRIX[1, 0] = P0[1] * 180 / Math::PI;
	CDMATRIX[0, 1] = P0[2] * 180 / Math::PI;
	CDMATRIX[1, 1] = P0[3] * 180 / Math::PI;
	CRPIXN[0] = P0[4];
	CRPIXN[1] = P0[5];
	CD1_1 = CDMATRIX[0, 0];
	CD1_2 = CDMATRIX[1, 0];
	CD2_1 = CDMATRIX[0, 1];
	CD2_2 = CDMATRIX[1, 1];

	CDELTN = gcnew array<double>(2);
	CDELTN[0] = Math::Sqrt(CD1_1 * CD1_1 + CD1_2 * CD1_2) * 3600;
	CDELTN[1] = Math::Sqrt(CD2_1 * CD2_1 + CD2_2 * CD2_2) * 3600;

	CROTAN = gcnew array<double>(2);
	CROTAN[0] = Math::Atan2(CD1_2, -CD1_1) * 180 / Math::PI;
	CROTAN[1] = Math::Atan2(-CD2_1, -CD2_2) * 180 / Math::PI;

	SET_CDMATRIXINV();

	array<double>^ dxpix = gcnew array<double>(CPIX1->Length);
	array<double>^ dypix = gcnew array<double>(CPIX1->Length);
	double xpix, ypix;
	for (int i = 0; i < CPIX1->Length; i++)
	{
		this->Get_Pixel(CVAL1[i], CVAL2[i], "TAN", xpix, ypix, false);
		dxpix[i] = xpix - CPIX1[i];
		dypix[i] = ypix - CPIX2[i];

		DVAL1[i] = dxpix[i] * CDELTN[0];
		DVAL2[i] = dypix[i] * CDELTN[1];
	}

	CPIX1RM = JPMath::Mean(dxpix, true);
	CPIX1RS = JPMath::Stdv(dxpix, true);
	CVAL1RM = CPIX1RM * CDELTN[0];
	CVAL1RS = CPIX1RS * CDELTN[0];
	CPIX2RM = JPMath::Mean(dypix, true);
	CPIX2RS = JPMath::Stdv(dypix, true);
	CVAL2RM = CPIX2RM * CDELTN[1];
	CVAL2RS = CPIX2RS * CDELTN[1];

	CPIXRM = Math::Sqrt(CPIX1RM * CPIX1RM + CPIX2RM * CPIX2RM);
	CPIXRS = Math::Sqrt(CPIX1RS * CPIX1RS + CPIX2RS * CPIX2RS);
	CVALRM = Math::Sqrt(CVAL1RM * CVAL1RM + CVAL2RM * CVAL2RM);
	CVALRS = Math::Sqrt(CVAL1RS * CVAL1RS + CVAL2RS * CVAL2RS);

	WCSEXISTS = true;

	if (header == nullptr)
		return;

	double ccvald1, ccvald2;
	String^ ccvals1;
	String^ ccvals2;
	double width = Convert::ToDouble(header->GetKeyValue("NAXIS1"));
	double height = Convert::ToDouble(header->GetKeyValue("NAXIS2"));
	Get_Coordinate(width / 2, height / 2, false, "TAN", ccvald1, ccvald2, ccvals1, ccvals2);
	CCVALD1 = ccvald1;
	CCVALD2 = ccvald2;
	CCVALS1 = ccvals1;
	CCVALS2 = ccvals2;

	Clear(header);
	this->CopyTo(header);	
}

void JPFITS::WorldCoordinateSolution::Get_Pixel(double cval1, double cval2, String^ WCS_Type, double &X_pix, double &Y_pix, bool return_zero_based_pixels)
{
	double a0 = CRVALN[0] * Math::PI / 180, d0 = CRVALN[1] * Math::PI / 180;
	double a = cval1 * Math::PI / 180, d = cval2 * Math::PI / 180;//radians
	double X_intrmdt = Math::Cos(d)*Math::Sin(a - a0) / (Math::Cos(d0)*Math::Cos(d)*Math::Cos(a - a0) + Math::Sin(d0)*Math::Sin(d));
	double Y_intrmdt = (Math::Cos(d0)*Math::Sin(d) - Math::Cos(d)*Math::Sin(d0)*Math::Cos(a - a0)) / (Math::Sin(d0)*Math::Sin(d) + Math::Cos(d0)*Math::Cos(d)*Math::Cos(a - a0));
	X_pix = CDMATRIXINV[0, 0] * X_intrmdt + CDMATRIXINV[1, 0] * Y_intrmdt + CRPIXN[0];
	Y_pix = CDMATRIXINV[0, 1] * X_intrmdt + CDMATRIXINV[1, 1] * Y_intrmdt + CRPIXN[1];
	if (return_zero_based_pixels)
	{
		X_pix--;
		Y_pix--;
	}
}

void JPFITS::WorldCoordinateSolution::Get_Pixels(array<double>^ cval1, array<double>^ cval2, String^ WCS_Type, array<double>^ &X_pix, array<double>^ &Y_pix, bool return_zero_based_pixels)
{
	double xpix, ypix;
	for (int i = 0; i < cval1->Length; i++)
	{
		this->Get_Pixel(cval1[i], cval2[i], WCS_Type, xpix, ypix, return_zero_based_pixels);
		X_pix[i] = xpix;
		Y_pix[i] = ypix;
	}
}

void JPFITS::WorldCoordinateSolution::Get_Coordinate(double X_pix, double Y_pix, bool zero_based_pixels, String^ WCS_Type, double &cval1, double &cval2)
{
	String^ sx1;
	String^ sx2;

	Get_Coordinate(X_pix, Y_pix, zero_based_pixels, WCS_Type, cval1, cval2, sx1, sx2);
}

void JPFITS::WorldCoordinateSolution::Get_Coordinate(double X_pix, double Y_pix, bool zero_based_pixels, String^ WCS_Type, String^ &cval1_sxgsml, String^ &cval2_sxgsml)
{
	double cv1, cv2;

	Get_Coordinate(X_pix, Y_pix, zero_based_pixels, WCS_Type, cv1, cv2, cval1_sxgsml, cval2_sxgsml);
}

void JPFITS::WorldCoordinateSolution::Get_Coordinate(double X_pix, double Y_pix, bool zero_based_pixels, String^ WCS_Type, double &cval1, double &cval2, String^ &cval1_sxgsml, String^ &cval2_sxgsml)
{
	if (zero_based_pixels)
	{
		X_pix++;
		Y_pix++;
	}
	double X_intrmdt = CDMATRIX[0, 0] * (X_pix - CRPIXN[0]) * Math::PI / 180 + CDMATRIX[1, 0] * (Y_pix - CRPIXN[1]) * Math::PI / 180;
	double Y_intrmdt = CDMATRIX[0, 1] * (X_pix - CRPIXN[0]) * Math::PI / 180 + CDMATRIX[1, 1] * (Y_pix - CRPIXN[1]) * Math::PI / 180;
	double a = CRVALN[0] * Math::PI / 180 + Math::Atan(X_intrmdt / (Math::Cos(CRVALN[1] * Math::PI / 180) - Y_intrmdt * Math::Sin(CRVALN[1] * Math::PI / 180)));
	double d = Math::Asin((Math::Sin(CRVALN[1] * Math::PI / 180) + Y_intrmdt * Math::Cos(CRVALN[1] * Math::PI / 180)) / Math::Sqrt(1 + X_intrmdt * X_intrmdt + Y_intrmdt * Y_intrmdt));
	a = a * 180 / Math::PI;
	d = d * 180 / Math::PI;

	if (a < 0)
		a += 360;

	cval1 = a;
	cval2 = d;

	double h = Math::Floor(a / 360 * 24);
	double m = Math::Floor((a / 360 * 24 - h) * 60);
	double s = Math::Round((a / 360 * 24 - h - m / 60) * 3600, 2);

	double decdeg = Math::Abs(d);
	double deg = Math::Floor(decdeg);
	double am = Math::Floor((decdeg - deg) * 60);
	double as = Math::Round((decdeg - deg - am / 60) * 3600, 2);

	String^ sign = "+";
	if (d < 0)
		sign = "-";

	cval1_sxgsml = h.ToString("00") + ":" + m.ToString("00") + ":" + s.ToString("00.00");
	cval2_sxgsml = sign + deg.ToString("00") + ":" + am.ToString("00") + ":" + as.ToString("00.00");
}

void JPFITS::WorldCoordinateSolution::Get_Coordinates(array<double>^ X_pix, array<double>^ Y_pix, bool zero_based_pixels, String^ WCS_Type, array<double>^ &cval1, array<double>^ &cval2)
{
	array<String^>^ sx1 = gcnew array<String^>(X_pix->Length);
	array<String^>^ sx2 = gcnew array<String^>(X_pix->Length);

	Get_Coordinates(X_pix, Y_pix, zero_based_pixels, WCS_Type, cval1, cval2, sx1, sx2);
}

void JPFITS::WorldCoordinateSolution::Get_Coordinates(array<double>^ X_pix, array<double>^ Y_pix, bool zero_based_pixels, String^ WCS_Type, array<String^>^ &cval1_sxgsml, array<String^>^ &cval2_sxgsml)
{
	array<double>^ cv1 = gcnew array<double>(X_pix->Length);
	array<double>^ cv2 = gcnew array<double>(X_pix->Length);

	Get_Coordinates(X_pix, Y_pix, zero_based_pixels, WCS_Type, cv1, cv2, cval1_sxgsml, cval2_sxgsml);
}

void JPFITS::WorldCoordinateSolution::Get_Coordinates(array<double>^ X_pix, array<double>^ Y_pix, bool zero_based_pixels, String^ WCS_Type, array<double>^ &cval1, array<double>^ &cval2, array<String^>^ &cval1_sxgsml, array<String^>^ &cval2_sxgsml)
{
	double radeg, decdeg;
	String^ rasx;
	String^ decsx;
	for (int i = 0; i < cval1->Length; i++)
	{
		this->Get_Coordinate(X_pix[i], Y_pix[i], zero_based_pixels, WCS_Type, radeg, decdeg, rasx, decsx);
		cval1[i] = radeg;
		cval2[i] = decdeg;
		cval1_sxgsml[i] = rasx;
		cval2_sxgsml[i] = decsx;
	}
}

void JPFITS::WorldCoordinateSolution::SET_CDMATRIXINV()
{
	CDMATRIXINV = gcnew array<double, 2>(2, 2);
	double det = 1 / ((CDMATRIX[0, 0] * CDMATRIX[1, 1] - CDMATRIX[1, 0] * CDMATRIX[0, 1])  * Math::PI / 180);
	CDMATRIXINV[0, 0] = det * CDMATRIX[1, 1];
	CDMATRIXINV[1, 0] = -det * CDMATRIX[1, 0];
	CDMATRIXINV[0, 1] = -det * CDMATRIX[0, 1];
	CDMATRIXINV[1, 1] = det * CDMATRIX[0, 0];
}

void JPFITS::WorldCoordinateSolution::CopyFrom(JPFITS::WorldCoordinateSolution^ wcs_source)
{
	try
	{
		WCSEXISTS = true;

		this->CD_Matrix = wcs_source->CD_Matrix;

		this->CTYPEN = gcnew array<String^>(2);
		this->CTYPEN[0] = wcs_source->CTYPEn[1];
		this->CTYPEN[1] = wcs_source->CTYPEn[2];

		this->CRPIXN = gcnew array<double>(2);
		this->CRPIXN[0] = wcs_source->CRPIXn[1];
		this->CRPIXN[1] = wcs_source->CRPIXn[2];

		this->CRVALN = gcnew array<double>(2);
		this->CRVALN[0] = wcs_source->CRVALn[1];
		this->CRVALN[1] = wcs_source->CRVALn[2];

		this->CDELTN = gcnew array<double>(2);
		this->CDELTN[0] = wcs_source->CDELTn[1];
		this->CDELTN[1] = wcs_source->CDELTn[2];

		this->CROTAN = gcnew array<double>(2);
		this->CROTAN[0] = wcs_source->CROTAn[1];
		this->CROTAN[1] = wcs_source->CROTAn[2];

		this->CPIX1 = wcs_source->Coordinate_Pixels[1];
		this->CPIX2 = wcs_source->Coordinate_Pixels[2];
		this->CVAL1 = wcs_source->Coordinate_Values[1];
		this->CVAL2 = wcs_source->Coordinate_Values[2];

		this->CCVALD1 = wcs_source->CCVALD1;
		this->CCVALD2 = wcs_source->CCVALD2;
		this->CCVALS1 = wcs_source->CCVALS1;
		this->CCVALS2 = wcs_source->CCVALS2;

		this->CPIX1RM = wcs_source->CPIX1RM;
		this->CPIX1RS = wcs_source->CPIX1RS;
		this->CVAL1RM = wcs_source->CVAL1RM;
		this->CVAL1RS = wcs_source->CVAL1RS;
		this->CPIX2RM = wcs_source->CPIX2RM;
		this->CPIX2RS = wcs_source->CPIX2RS;
		this->CVAL2RM = wcs_source->CVAL2RM;
		this->CVAL2RS = wcs_source->CVAL2RS;
		this->CPIXRM = wcs_source->CPIXRM;
		this->CPIXRS = wcs_source->CPIXRS;
		this->CVALRM = wcs_source->CVALRM;
		this->CVALRS = wcs_source->CVALRS;
	}
	catch (Exception^ e)
	{
		MessageBox::Show(e->Data + "	" + e->InnerException + "	" + e->Message + "	" + e->Source + "	" + e->StackTrace + "	" + e->TargetSite);
	}
}

void JPFITS::WorldCoordinateSolution::CopyTo(JPFITS::FITSImageHeader^ header)
{
	try
	{
		header->SetKey("CTYPE1", CTYPEN[0], "WCS type of horizontal coordinate transformation", true, -1);
		header->SetKey("CTYPE2", CTYPEN[1], "WCS type of vertical coordinate transformation", true, -1);
		header->SetKey("CRPIX1", CRPIXN[0].ToString("F5"), "WCS coordinate reference pixel on axis 1", true, -1);
		header->SetKey("CRPIX2", CRPIXN[1].ToString("F5"), "WCS coordinate reference pixel on axis 2", true, -1);
		header->SetKey("CRVAL1", CRVALN[0].ToString("F8"), "WCS coordinate reference value on axis 1 (deg)", true, -1);
		header->SetKey("CRVAL2", CRVALN[1].ToString("F8"), "WCS coordinate reference value on axis 2 (deg)", true, -1);
		header->SetKey("CD1_1", CDMATRIX[0, 0].ToString("0.0#########e+00"), "WCS rotation and scaling matrix", true, -1);
		header->SetKey("CD1_2", CDMATRIX[1, 0].ToString("0.0#########e+00"), "WCS rotation and scaling matrix", true, -1);
		header->SetKey("CD2_1", CDMATRIX[0, 1].ToString("0.0#########e+00"), "WCS rotation and scaling matrix", true, -1);
		header->SetKey("CD2_2", CDMATRIX[1, 1].ToString("0.0#########e+00"), "WCS rotation and scaling matrix", true, -1);
		header->SetKey("CDELT1", CDELTN[0].ToString("F8"), "WCS plate scale on axis 1 (arcsec per pixel)", true, -1);
		header->SetKey("CDELT2", CDELTN[1].ToString("F8"), "WCS plate Scale on axis 2 (arcsec per pixel)", true, -1);
		header->SetKey("CROTA1", CROTAN[0].ToString("F8"), "WCS field rotation angle on axis 1 (degrees)", true, -1);
		header->SetKey("CROTA2", CROTAN[1].ToString("F8"), "WCS field rotation angle on axis 2 (degrees)", true, -1);
		header->SetKey("CCVALD1", CCVALD1.ToString("F8"), "WCS field center on axis 1 (degrees)", true, -1);
		header->SetKey("CCVALD2", CCVALD2.ToString("F8"), "WCS field center on axis 2 (degrees)", true, -1);
		header->SetKey("CCVALS1", CCVALS1, "WCS field center on axis 1 (sexigesimal h m s)", true, -1);
		header->SetKey("CCVALS2", CCVALS2, "WCS field center on axis 2 (sexigesimal d am as)", true, -1);
		header->SetKey("CPIX1RM", CPIX1RM.ToString("G"), "Mean of WCS residuals on axis 1 (pixels)", true, -1);
		header->SetKey("CPIX1RS", CPIX1RS.ToString("G"), "Standard dev of WCS residuals on axis 1 (pixels)", true, -1);
		header->SetKey("CVAL1RM", CVAL1RM.ToString("G"), "Mean of WCS residuals on axis 1 (arcsec)", true, -1);
		header->SetKey("CVAL1RS", CVAL1RS.ToString("G"), "Standard dev of WCS residuals on axis 1 (arcsec)", true, -1);
		header->SetKey("CPIX2RM", CPIX2RM.ToString("G"), "Mean of WCS residuals on axis 2 (pixels)", true, -1);
		header->SetKey("CPIX2RS", CPIX2RS.ToString("G"), "Standard dev of WCS residuals on axis 2 (pixels)", true, -1);
		header->SetKey("CVAL2RM", CVAL2RM.ToString("G"), "Mean of WCS residuals on axis 2 (arcsec)", true, -1);
		header->SetKey("CVAL2RS", CVAL2RS.ToString("G"), "Standard dev of WCS residuals on axis 2 (arcsec)", true, -1);
		header->SetKey("CPIXRM", CPIXRM.ToString("G"), "Mean of WCS residuals (pixels)", true, -1);
		header->SetKey("CPIXRS", CPIXRS.ToString("G"), "Standard dev of WCS residuals (pixels)", true, -1);
		header->SetKey("CVALRM", CVALRM.ToString("G"), "Mean of WCS residuals (arcsec)", true, -1);
		header->SetKey("CVALRS", CVALRS.ToString("G"), "Standard dev of WCS residuals (arcsec)", true, -1);

		int key = 0, num = 1;
		while (key != -1)
		{
			key = header->GetKeyIndex("WCD1_" + num.ToString("000"), false);
			header->RemoveKey(key);
			key = header->GetKeyIndex("WCD2_" + num.ToString("000"), false);
			header->RemoveKey(key);
			key = header->GetKeyIndex("WCP1_" + num.ToString("000"), false);
			header->RemoveKey(key);
			key = header->GetKeyIndex("WCP2_" + num.ToString("000"), false);
			header->RemoveKey(key);
			key = header->GetKeyIndex("WCV1_" + num.ToString("000"), false);
			header->RemoveKey(key);
			key = header->GetKeyIndex("WCV2_" + num.ToString("000"), false);
			header->RemoveKey(key);
			num++;
		}

		num = 1;
		for (int i = 0; i < CPIX1->Length; i++)
		{
			header->SetKey("WCP1_" + num.ToString("000"), CPIX1[i].ToString("F5"), "WCS coordinate pixel on axis 1", true, -1);
			header->SetKey("WCP2_" + num.ToString("000"), CPIX2[i].ToString("F5"), "WCS coordinate pixel on axis 2", true, -1);
			header->SetKey("WCV1_" + num.ToString("000"), CVAL1[i].ToString("F8"), "WCS coordinate value on axis 1 (degrees)", true, -1);
			header->SetKey("WCV2_" + num.ToString("000"), CVAL2[i].ToString("F8"), "WCS coordinate value on axis 2 (degrees)", true, -1);
			header->SetKey("WCD1_" + num.ToString("000"), DVAL1[i].ToString("F8"), "WCS coordinate delta on axis 1 (arcsec)", true, -1);
			header->SetKey("WCD2_" + num.ToString("000"), DVAL2[i].ToString("F8"), "WCS coordinate delta on axis 2 (arcsec)", true, -1);
			num++;
		}
	}
	catch (Exception^ e)
	{
		MessageBox::Show(e->Data + "	" + e->InnerException + "	" + e->Message + "	" + e->Source + "	" + e->StackTrace + "	" + e->TargetSite);
	}
}

void JPFITS::WorldCoordinateSolution::Clear()
{
	CDMATRIX = nullptr;
	CDMATRIXINV = nullptr;

	CD1_1 = 0;
	CD1_2 = 0;
	CD2_1 = 0;
	CD2_2 = 0;

	CTYPEN = nullptr;
	CRPIXN = nullptr;
	CRVALN = nullptr;
	CDELTN = nullptr;
	CROTAN = nullptr;

	CPIX1RM = 0; 
	CPIX1RS = 0; 
	CVAL1RM = 0; 
	CVAL1RS = 0; 
	CPIX2RM = 0; 
	CPIX2RS = 0; 
	CVAL2RM = 0; 
	CVAL2RS = 0; 
	CPIXRM = 0; 
	CPIXRS = 0; 
	CVALRM = 0; 
	CVALRS = 0; 
	CCVALD1 = 0; 
	CCVALD2 = 0;
	CCVALS1 = "";
	CCVALS2 = "";

	WCSEXISTS = false;
}

void JPFITS::WorldCoordinateSolution::Clear(JPFITS::FITSImageHeader^ header)
{
	header->RemoveKey("CTYPE1");
	header->RemoveKey("CTYPE2");
	header->RemoveKey("CRPIX1");
	header->RemoveKey("CRPIX2");
	header->RemoveKey("CRVAL1");
	header->RemoveKey("CRVAL2");
	header->RemoveKey("CD1_1");
	header->RemoveKey("CD1_2");
	header->RemoveKey("CD2_1");
	header->RemoveKey("CD2_2");
	header->RemoveKey("CDELT1");
	header->RemoveKey("CDELT2");
	header->RemoveKey("CROTA1");
	header->RemoveKey("CROTA2");
	header->RemoveKey("CCVALD1");
	header->RemoveKey("CCVALD2");
	header->RemoveKey("CCVALS1");
	header->RemoveKey("CCVALS2");
	header->RemoveKey("CPIX1RM");
	header->RemoveKey("CPIX1RS");
	header->RemoveKey("CVAL1RM");
	header->RemoveKey("CVAL1RS");
	header->RemoveKey("CPIX2RM");
	header->RemoveKey("CPIX2RS");
	header->RemoveKey("CVAL2RM");
	header->RemoveKey("CVAL2RS");
	header->RemoveKey("CPIXRM");
	header->RemoveKey("CPIXRS");
	header->RemoveKey("CVALRM");
	header->RemoveKey("CVALRS");

	int key = 0, num = 1;
	while (key != -1)
	{
		key = header->GetKeyIndex("WCD1_" + num.ToString("000"), false);
		header->RemoveKey(key);
		key = header->GetKeyIndex("WCD2_" + num.ToString("000"), false);
		header->RemoveKey(key);
		key = header->GetKeyIndex("WCP1_" + num.ToString("000"), false);
		header->RemoveKey(key);
		key = header->GetKeyIndex("WCP2_" + num.ToString("000"), false);
		header->RemoveKey(key);
		key = header->GetKeyIndex("WCV1_" + num.ToString("000"), false);
		header->RemoveKey(key);
		key = header->GetKeyIndex("WCV2_" + num.ToString("000"), false);
		header->RemoveKey(key);
		num++;
	}
}

bool JPFITS::WorldCoordinateSolution::Exists(FITSImageHeader^ header, array<String^>^ WCS_CTYPEN)
{
	for (int i = 1; i <= WCS_CTYPEN->Length; i++)
		if (!header->GetKeyValue("CTYPE" + i.ToString())->Contains(WCS_CTYPEN[i - 1]))
			return false;
	
	return true;
}

