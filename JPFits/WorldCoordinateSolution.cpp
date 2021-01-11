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

void JPFITS::WorldCoordinateSolution::Solve_WCS(String^ WCS_Type, array<double>^ X_pix, array<double>^ Y_pix, bool zero_based_pixels, array<double>^ cval1, array<double>^ cval2, FITSImage^ FITS)
{
	//should first do a check of WCS type to make sure it is valid
	CTYPEN = gcnew array<String^>(2);
	CTYPEN[0] = "RA---" + WCS_Type;
	CTYPEN[1] = "DEC--" + WCS_Type;

	CPIX1 = gcnew array<double>(X_pix->Length);
	CPIX2 = gcnew array<double>(X_pix->Length);
	CVAL1 = gcnew array<double>(X_pix->Length);
	CVAL2 = gcnew array<double>(X_pix->Length);
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
	}

	CPIX1RM = JPMath::Mean(dxpix, true);
	CPIX1RS = JPMath::Stdv(dxpix, true);
	CVAL1RM = JPMath::Mean(dxpix, true) * CDELTN[0];
	CVAL1RS = JPMath::Stdv(dxpix, true) * CDELTN[0];
	CPIX2RM = JPMath::Mean(dypix, true);
	CPIX2RS = JPMath::Stdv(dypix, true);
	CVAL2RM = JPMath::Mean(dypix, true) * CDELTN[1];
	CVAL2RS = JPMath::Stdv(dypix, true) * CDELTN[1];

	CPIXRM = Math::Sqrt(CPIX1RM * CPIX1RM + CPIX2RM * CPIX2RM);
	CPIXRS = Math::Sqrt(CPIX1RS * CPIX1RS + CPIX2RS * CPIX2RS);
	CVALRM = Math::Sqrt(CVAL1RM * CVAL1RM + CVAL2RM * CVAL2RM);
	CVALRS = Math::Sqrt(CVAL1RS * CVAL1RS + CVAL2RS * CVAL2RS);

	if (FITS == nullptr)
		return;

	Clear(FITS);
	this->CopyTo(FITS);
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

	//linear
	//int xpix = int(WCS_CDMATRIX_INV[0, 0] * (RA - WCS_CRVAL1) + WCS_CDMATRIX_INV[1, 0] * (Dec - WCS_CRVAL2) + WCS_CRPIX1);
	//int ypix = int(WCS_CDMATRIX_INV[0, 1] * (RA - WCS_CRVAL1) + WCS_CDMATRIX_INV[1, 1] * (Dec - WCS_CRVAL2) + WCS_CRPIX2);
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
	//double d = Math::Atan((Math::Cos(a - WCS_CRVAL1 * Math::PI / 180) * (Y + Math::Tan(WCS_CRVAL2 * Math::PI / 180))) / (1 - Y / Math::Tan(WCS_CRVAL2 * Math::PI / 180)));
	a = a * 180 / Math::PI;
	d = d * 180 / Math::PI;

	//linear
	//double a = WCS_CDMATRIX[0, 0]*(XPOS_CURSOR - WCS_CRPIX1) + WCS_CDMATRIX[1, 0]*(YPOS_CURSOR - WCS_CRPIX2) + WCS_CRVAL1;
	//double d = WCS_CDMATRIX[0, 1]*(XPOS_CURSOR - WCS_CRPIX1) + WCS_CDMATRIX[1, 1]*(YPOS_CURSOR - WCS_CRPIX2) + WCS_CRVAL2;

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

	cval1_sxgsml = h.ToString("00") + "h " + m.ToString("00") + "m " + s.ToString("00.##") + "s";
	cval2_sxgsml = sign + deg.ToString("00") + "d " + am.ToString("00") + "' " + as.ToString("00.##") + "''";
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

void JPFITS::WorldCoordinateSolution::Get_WCS(FITSImage^ FITS)
{
	if (FITS->GetKeyValue("CD1_1") == "")
		throw gcnew Exception("Error: CD matrix not found...");

	CDMATRIX = gcnew array<double, 2>(2, 2);
	CDMATRIXINV = gcnew array<double, 2>(2, 2);

	CTYPEN = gcnew array<String^>(2);
	CTYPEN[0] = FITS->GetKeyValue("CTYPE1");
	CTYPEN[1] = FITS->GetKeyValue("CTYPE2");

	CRPIXN = gcnew array<double>(2);
	CRPIXN[0] = Convert::ToDouble(FITS->GetKeyValue("CRPIX1"));
	CRPIXN[1] = Convert::ToDouble(FITS->GetKeyValue("CRPIX2"));

	CRVALN = gcnew array<double>(2);
	CRVALN[0] = Convert::ToDouble(FITS->GetKeyValue("CRVAL1"));
	CRVALN[1] = Convert::ToDouble(FITS->GetKeyValue("CRVAL2"));

	CDMATRIX = gcnew array<double, 2>(2, 2);
	CDMATRIX[0, 0] = Convert::ToDouble(FITS->GetKeyValue("CD1_1"));
	CDMATRIX[1, 0] = Convert::ToDouble(FITS->GetKeyValue("CD1_2"));
	CDMATRIX[0, 1] = Convert::ToDouble(FITS->GetKeyValue("CD2_1"));
	CDMATRIX[1, 1] = Convert::ToDouble(FITS->GetKeyValue("CD2_2"));

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

	try
	{
		CDELTN[0] = Convert::ToDouble(FITS->GetKeyValue("CDELT1"));
		CDELTN[1] = Convert::ToDouble(FITS->GetKeyValue("CDELT2"));
		CROTAN[0] = Convert::ToDouble(FITS->GetKeyValue("CROTA1"));
		CROTAN[1] = Convert::ToDouble(FITS->GetKeyValue("CROTA2"));
		CPIX1RM = Convert::ToDouble(FITS->GetKeyValue("CPIX1RM"));
		CPIX1RS = Convert::ToDouble(FITS->GetKeyValue("CPIX1RS"));
		CVAL1RM = Convert::ToDouble(FITS->GetKeyValue("CVAL1RM"));
		CVAL1RS = Convert::ToDouble(FITS->GetKeyValue("CVAL1RS"));
		CPIX2RM = Convert::ToDouble(FITS->GetKeyValue("CPIX2RM"));
		CPIX2RS = Convert::ToDouble(FITS->GetKeyValue("CPIX2RS"));
		CVAL2RM = Convert::ToDouble(FITS->GetKeyValue("CVAL2RM"));
		CVAL2RS = Convert::ToDouble(FITS->GetKeyValue("CVAL2RS"));
		CPIXRM = Convert::ToDouble(FITS->GetKeyValue("CPIXRM"));
		CPIXRS = Convert::ToDouble(FITS->GetKeyValue("CPIXRS"));
		CVALRM = Convert::ToDouble(FITS->GetKeyValue("CVALRM"));
		CVALRS = Convert::ToDouble(FITS->GetKeyValue("CVALRS"));

		int num = 0, key = 0;
		while (key != -1)
		{
			num++;
			key = FITS->GetKeyIndex("WCP1_" + num.ToString("000"));
		}
		num--;
		CPIX1 = gcnew array<double>(num);
		CPIX2 = gcnew array<double>(num);
		CVAL1 = gcnew array<double>(num);
		CVAL2 = gcnew array<double>(num);
		num = 1;
		for (int i = 0; i < CPIX1->Length; i++)
		{
			CPIX1[i] = Convert::ToDouble(FITS->GetKeyValue("WCP1_" + num.ToString("000")));
			CPIX2[i] = Convert::ToDouble(FITS->GetKeyValue("WCP2_" + num.ToString("000")));
			CVAL1[i] = Convert::ToDouble(FITS->GetKeyValue("WCV1_" + num.ToString("000")));
			CVAL2[i] = Convert::ToDouble(FITS->GetKeyValue("WCV2_" + num.ToString("000")));
			num++;
		}
	}
	catch (... /*Exception^ e*/)
	{
		//MessageBox::Show(e->Data + "	" + e->InnerException + "	" + e->Message + "	" + e->Source + "	" + e->StackTrace + "	" + e->TargetSite);
	}
}

void JPFITS::WorldCoordinateSolution::Get_WCS(String^ FITS_filename)
{
	FITSImage^ FITS = gcnew FITSImage(FITS_filename, nullptr, true, false, false, false);
	this->Get_WCS(FITS);
}

void JPFITS::WorldCoordinateSolution::CopyTo(FITSImage^ FITS)
{
	FITS->SetKey("CTYPE1", CTYPEN[0], "WCS type of horizontal coordinate transformation", true, -1);
	FITS->SetKey("CTYPE2", CTYPEN[1], "WCS type of vertical coordinate transformation", true, -1);
	FITS->SetKey("CRPIX1", CRPIXN[0].ToString("F5"), "WCS coordinate reference pixel on axis 1", true, -1);
	FITS->SetKey("CRPIX2", CRPIXN[1].ToString("F5"), "WCS coordinate reference pixel on axis 2", true, -1);
	FITS->SetKey("CRVAL1", CRVALN[0].ToString("F8"), "WCS coordinate reference value on axis 1 (deg)", true, -1);
	FITS->SetKey("CRVAL2", CRVALN[1].ToString("F8"), "WCS coordinate reference value on axis 2 (deg)", true, -1);
	FITS->SetKey("CD1_1", CDMATRIX[0, 0].ToString("0.0#######e+00"), "WCS rotation and scaling matrix", true, -1);
	FITS->SetKey("CD1_2", CDMATRIX[1, 0].ToString("0.0#######e+00"), "WCS rotation and scaling matrix", true, -1);
	FITS->SetKey("CD2_1", CDMATRIX[0, 1].ToString("0.0#######e+00"), "WCS rotation and scaling matrix", true, -1);
	FITS->SetKey("CD2_2", CDMATRIX[1, 1].ToString("0.0#######e+00"), "WCS rotation and scaling matrix", true, -1);
	FITS->SetKey("CDELT1", CDELTN[0].ToString("F8"), "WCS plate scale on axis 1 (arcsec per pixel)", true, -1);
	FITS->SetKey("CDELT2", CDELTN[1].ToString("F8"), "WCS plate Scale on axis 2 (arcsec per pixel)", true, -1);
	FITS->SetKey("CROTA1", CROTAN[0].ToString("F8"), "WCS field rotation angle on axis 1 (degrees)", true, -1);
	FITS->SetKey("CROTA2", CROTAN[1].ToString("F8"), "WCS field rotation angle on axis 2 (degrees)", true, -1);

	double cval1, cval2;
	String^ scval1;
	String^ scval2;
	Get_Coordinate(FITS->Width / 2, FITS->Height / 2, false, "TAN", cval1, cval2, scval1, scval2);
	FITS->SetKey("CCVALD1", cval1.ToString("F8"), "WCS field center on axis 1 (degrees)", true, -1);
	FITS->SetKey("CCVALD2", cval2.ToString("F8"), "WCS field center on axis 2 (degrees)", true, -1);
	FITS->SetKey("CCVALS1", scval1, "WCS field center on axis 1 (sexigesimal h m s)", true, -1);
	FITS->SetKey("CCVALS2", scval2, "WCS field center on axis 2 (sexigesimal d am as)", true, -1);
	FITS->SetKey("CPIX1RM", CPIX1RM.ToString("G"), "Mean of WCS residuals on axis 1 (pixels)", true, -1);
	FITS->SetKey("CPIX1RS", CPIX1RS.ToString("G"), "Standard dev of WCS residuals on axis 1 (pixels)", true, -1);
	FITS->SetKey("CVAL1RM", CVAL1RM.ToString("G"), "Mean of WCS residuals on axis 1 (arcsec)", true, -1);
	FITS->SetKey("CVAL1RS", CVAL1RS.ToString("G"), "Standard dev of WCS residuals on axis 1 (arcsec)", true, -1);
	FITS->SetKey("CPIX2RM", CPIX2RM.ToString("G"), "Mean of WCS residuals on axis 2 (pixels)", true, -1);
	FITS->SetKey("CPIX2RS", CPIX2RS.ToString("G"), "Standard dev of WCS residuals on axis 2 (pixels)", true, -1);
	FITS->SetKey("CVAL2RM", CVAL2RM.ToString("G"), "Mean of WCS residuals on axis 2 (arcsec)", true, -1);
	FITS->SetKey("CVAL2RS", CVAL2RS.ToString("G"), "Standard dev of WCS residuals on axis 2 (arcsec)", true, -1);
	FITS->SetKey("CPIXRM", CPIXRM.ToString("G"), "Mean of WCS residuals (pixels)", true, -1);
	FITS->SetKey("CPIXRS", CPIXRS.ToString("G"), "Standard dev of WCS residuals (pixels)", true, -1);
	FITS->SetKey("CVALRM", CVALRM.ToString("G"), "Mean of WCS residuals (arcsec)", true, -1);
	FITS->SetKey("CVALRS", CVALRS.ToString("G"), "Standard dev of WCS residuals (arcsec)", true, -1);

	int key = 0, num = 1;
	while (key != -1)
	{
		key = FITS->GetKeyIndex("WCP1_" + num.ToString("000"));
		FITS->RemoveKey(key);
		key = FITS->GetKeyIndex("WCP2_" + num.ToString("000"));
		FITS->RemoveKey(key);
		key = FITS->GetKeyIndex("WCV1_" + num.ToString("000"));
		FITS->RemoveKey(key);
		key = FITS->GetKeyIndex("WCV2_" + num.ToString("000"));
		FITS->RemoveKey(key);
		num++;
	}

	num = 1;
	for (int i = 0; i < CPIX1->Length; i++)
	{
		FITS->SetKey("WCP1_" + num.ToString("000"), CPIX1[i].ToString("F5"), "WCS coordinate pixel on axis 1", true, -1);
		FITS->SetKey("WCP2_" + num.ToString("000"), CPIX2[i].ToString("F5"), "WCS coordinate pixel on axis 2", true, -1);
		FITS->SetKey("WCV1_" + num.ToString("000"), CVAL1[i].ToString("F8"), "WCS coordinate value on axis 1 (degrees)", true, -1);
		FITS->SetKey("WCV2_" + num.ToString("000"), CVAL2[i].ToString("F8"), "WCS coordinate value on axis 2 (degrees)", true, -1);
		num++;
	}
}

void JPFITS::WorldCoordinateSolution::Clear(FITSImage^ FITS)
{
	FITS->RemoveKey("CTYPE1");
	FITS->RemoveKey("CTYPE2");
	FITS->RemoveKey("CRPIX1");
	FITS->RemoveKey("CRPIX2");
	FITS->RemoveKey("CRVAL1");
	FITS->RemoveKey("CRVAL2");
	FITS->RemoveKey("CD1_1");
	FITS->RemoveKey("CD1_2");
	FITS->RemoveKey("CD2_1");
	FITS->RemoveKey("CD2_2");
	FITS->RemoveKey("CDELT1");
	FITS->RemoveKey("CDELT2");
	FITS->RemoveKey("CROTA1");
	FITS->RemoveKey("CROTA2");
	FITS->RemoveKey("CCVALD1");
	FITS->RemoveKey("CCVALD2");
	FITS->RemoveKey("CCVALS1");
	FITS->RemoveKey("CCVALS2");
	FITS->RemoveKey("CPIX1RM");
	FITS->RemoveKey("CPIX1RS");
	FITS->RemoveKey("CVAL1RM");
	FITS->RemoveKey("CVAL1RS");
	FITS->RemoveKey("CPIX2RM");
	FITS->RemoveKey("CPIX2RS");
	FITS->RemoveKey("CVAL2RM");
	FITS->RemoveKey("CVAL2RS");
	FITS->RemoveKey("CPIXRM");
	FITS->RemoveKey("CPIXRS");
	FITS->RemoveKey("CVALRM");
	FITS->RemoveKey("CVALRS");

	int key = 0, num = 1;
	while (key != -1)
	{
		key = FITS->GetKeyIndex("WCP1_" + num.ToString("000"));
		FITS->RemoveKey(key);
		key = FITS->GetKeyIndex("WCP2_" + num.ToString("000"));
		FITS->RemoveKey(key);
		key = FITS->GetKeyIndex("WCV1_" + num.ToString("000"));
		FITS->RemoveKey(key);
		key = FITS->GetKeyIndex("WCV2_" + num.ToString("000"));
		FITS->RemoveKey(key);
		num++;
	}
}

bool JPFITS::WorldCoordinateSolution::Exists(FITSImage^ FITS, array<String^>^ WCS_CTYPEN)
{
	for (int i = 1; i <= WCS_CTYPEN->Length; i++)
		if (!FITS->GetKeyValue("CTYPE" + i.ToString())->Contains(WCS_CTYPEN[i - 1]))
			return false;
	
	return true;
}

