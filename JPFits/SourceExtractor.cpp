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

JPFITS::SourceExtractor::~SourceExtractor()
{

}

void JPFITS::SourceExtractor::INITTHIS()
{
	this->BGWRKR = gcnew BackgroundWorker();
	this->BGWRKR->WorkerReportsProgress = true;
	this->BGWRKR->WorkerSupportsCancellation = true;
	this->BGWRKR->DoWork += gcnew System::ComponentModel::DoWorkEventHandler(this, &SourceExtractor::BGWRKR_DoWork);
	this->BGWRKR->ProgressChanged += gcnew System::ComponentModel::ProgressChangedEventHandler(this, &SourceExtractor::BGWRKR_ProgressChanged);
	this->BGWRKR->RunWorkerCompleted += gcnew System::ComponentModel::RunWorkerCompletedEventHandler(this, &SourceExtractor::BGWRKR_RunWorkerCompleted);
}

JPFITS::SourceExtractor::SourceExtractor()
{
	INITTHIS();
}

JPFITS::SourceExtractor::SourceExtractor(array<double>^ XCoords, array<double>^ YCoords)
{
	INITTHIS();
	this->N_SRC = XCoords->Length;
	this->INITARRAYS();
	this->Centroids_X = XCoords;
	this->Centroids_Y = YCoords;
}

inline array<double, 2>^ JPFITS::SourceExtractor::GetKernel(array<double, 2>^ image, int x0, int y0, int radius)
{
	int width = radius * 2 + 1, kx = -1, ky = 0, x, y, xmin = x0 - radius, ymin = y0 - radius;
	array<double, 2>^ kernel = gcnew array<double, 2>(width, width);

	for (x = 0; x < width; x++)
	{
		kx = xmin + x;
		for (y = 0; y < width; y++)
			kernel[x, y] = image[kx, ymin + y];
	}
	
	return kernel;
}

inline void JPFITS::SourceExtractor::Centroid(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ kernel, double &x_centroid, double& y_centroid)
{
	int xw = kernel->GetLength(0);
	int yh = kernel->GetLength(1);

	if (xdata == nullptr)
	{
		xdata = gcnew array<int>(xw);
		for (int i = 0; i < xw; i++)
			xdata[i] = i;
	}
	if (ydata == nullptr)
	{
		ydata = gcnew array<int>(yh);
		for (int i = 0; i < yh; i++)
			ydata[i] = i;
	}

	double xweighted = 0, yweighted = 0, kernel_sum = 0;
	for (int x = 0; x < xw; x++)
		for (int y = 0; y < yh; y++)
		{
			xweighted += kernel[x, y] * double(xdata[x]);
			yweighted += kernel[x, y] * double(ydata[y]);
			kernel_sum += kernel[x, y];
		}

	x_centroid = xweighted / kernel_sum;
	y_centroid = yweighted / kernel_sum;
}

void JPFITS::SourceExtractor::Extract_Sources(array<double, 2>^ image, array<double>^ XCoords, array<double>^ YCoords, int kernel_radius, bool auto_background, String^ kernel_filename_template)
{
	KERNEL_RADIUS = kernel_radius;
	KERNEL_WIDTH = KERNEL_RADIUS * 2 + 1;
	N_SRC = XCoords->Length;
	INITARRAYS();
	AUTO_BG = auto_background;
	SAVE_PS = kernel_filename_template != nullptr && kernel_filename_template != "";
	SAVE_PS_FILENAME = kernel_filename_template;
	IMAGE = image;
	//IMAGE_KERNEL_BOOL_SOURCE = gcnew array<bool, 2>(IMAGE->GetLength(0), IMAGE->GetLength(1));
	FITTED = false;

	int HW = KERNEL_RADIUS;
	int intprog = 0;

	#pragma omp parallel for
	for (int i = 0; i < N_SRC; i++)
	{
		array<double, 2>^ kernel = GetKernel(image, (int)XCoords[i], (int)YCoords[i], KERNEL_RADIUS);
		/*int xmax, ymax;
		JPMath::Max(kernel, xmax, ymax);
		if (xmax != KERNEL_RADIUS || ymax != KERNEL_RADIUS)
		{
			XCoords[i] += double(xmax - KERNEL_RADIUS);
			YCoords[i] += double(ymax - KERNEL_RADIUS);
			kernel = GetKernel(image, (int)XCoords[i], (int)YCoords[i], KERNEL_RADIUS);
		}*/

		double bg_est = 0;//default
		if (AUTO_BG)
		{
			bg_est = ESTIMATELOCALBACKGROUND((int)XCoords[i], (int)YCoords[i], HW);
			kernel = JPMath::MatrixAddScalar(kernel, -bg_est, false);
		}

		array<double>^ xcoords = gcnew array<double>(KERNEL_WIDTH);// x coords
		array<double>^ ycoords = gcnew array<double>(KERNEL_WIDTH);// y coords
		for (int j = 0; j < KERNEL_WIDTH; j++)
		{
			xcoords[j] = double((int)XCoords[i] - HW + j);
			ycoords[j] = double((int)YCoords[i] - HW + j);
		}

		double xweighted = 0, yweighted = 0;
		for (int x = 0; x < KERNEL_WIDTH; x++)
			for (int y = 0; y < KERNEL_WIDTH; y++)
			{
				xweighted += kernel[x, y] * xcoords[x];
				yweighted += kernel[x, y] * ycoords[y];
			}

		//centroids
		double kernel_sum = JPMath::Sum(kernel, false);
		CENTROIDS_X[i] = xweighted / kernel_sum;
		CENTROIDS_Y[i] = yweighted / kernel_sum;
		CENTROIDS_VOLUME[i] = kernel_sum;
		CENTROIDS_AMPLITUDE[i] = kernel[KERNEL_RADIUS, KERNEL_RADIUS];
		CENTROIDS_BGESTIMATE[i] = bg_est;

		if (SAVE_PS)
		{
			String^ file = SAVE_PS_FILENAME;
			int ind = file->LastIndexOf(".");//for saving PS
			file = String::Concat(file->Substring(0, ind), "_", (i + 1).ToString("00000000"), ".fits");

			JPFITS::FITSImage^ f = gcnew ::JPFITS::FITSImage(file, kernel, false, false);
			f->WriteImage(TypeCode::Double, false);
		}
	}
}

void JPFITS::SourceExtractor::Extract_Sources(array<double, 2>^ image, double pix_saturation, double pix_min, double pix_max, double kernel_min, double kernel_max, bool threshholds_as_SN, int kernel_radius, int source_separation, bool auto_background, String^ kernel_filename_template, array<bool, 2>^ ROI_region, bool show_waitbar)
{
	IMAGE = image;
	PIX_SAT = pix_saturation;
	IMAGEWIDTH = IMAGE->GetLength(0);
	IMAGEHEIGHT = IMAGE->GetLength(1);
	KERNEL_RADIUS = kernel_radius;
	KERNEL_WIDTH = KERNEL_RADIUS * 2 + 1;
	SOURCE_SEPARATION = source_separation;
	PIX_MAX = pix_max;
	PIX_MIN = pix_min;
	N_SRC = 0;
	KERNEL_MIN = kernel_min;
	KERNEL_MAX = kernel_max;
	AUTO_BG = auto_background;
	SAVE_PS = kernel_filename_template != nullptr && kernel_filename_template != "";
	SAVE_PS_FILENAME = kernel_filename_template;
	SOURCE_BOOLEAN_MAP = gcnew array<bool, 2>(IMAGEWIDTH, IMAGEHEIGHT);
	SOURCE_INDEX_MAP = gcnew array<int, 2>(IMAGEWIDTH, IMAGEHEIGHT);
	#pragma omp parallel for
	for (int i = 0; i < IMAGE->GetLength(0); i++)
		for (int j = 0; j < IMAGE->GetLength(1); j++)
			SOURCE_INDEX_MAP[i, j] = -1;
	FITTED = false;
	THRESHHOLDS_AS_SN = threshholds_as_SN;
	SEARCH_ROI = ROI_region != nullptr;
	ROI_REGION = ROI_region;
	SHOWWAITBAR = show_waitbar;

	if (SHOWWAITBAR)
	{
		WAITBAR = gcnew JPWaitBar::WaitBar();
		WAITBAR->ProgressBar->Maximum = 100;
		WAITBAR->Text = "Scanning Image...";
		BGWRKR->RunWorkerAsync(1);
		WAITBAR->ShowDialog();

		if (WAITBAR->DialogResult == DialogResult::Cancel)
			N_SRC = -1;
	}
	else
	{
		Object^ sender = gcnew Object();
		System::ComponentModel::DoWorkEventArgs^ e = gcnew System::ComponentModel::DoWorkEventArgs(1);
		BGWRKR_DoWork(sender, e);
	}
}

void JPFITS::SourceExtractor::BGWRKR_DoWork(System::Object^  sender, System::ComponentModel::DoWorkEventArgs^  e)
{
	if (::Convert::ToInt32(e->Argument) == 1)//Extract
	{
		ArrayList^ Xs = gcnew ArrayList();// x positions
		ArrayList^ Ys = gcnew ArrayList();// y positions
		ArrayList^ Ks = gcnew ArrayList();// kernel sums
		ArrayList^ Ps = gcnew ArrayList();// pixel values
		ArrayList^ Bs = gcnew ArrayList();// background estimates
		double SSEP2 = SOURCE_SEPARATION * SOURCE_SEPARATION, KRAD2 = KERNEL_RADIUS * KERNEL_RADIUS;
		int src_index = 0;
		int intprog = 0;

		if (PIX_SAT > 0)//check for saturation islands
		{
			#pragma omp parallel for
			for (int x = SOURCE_SEPARATION; x < IMAGEWIDTH - SOURCE_SEPARATION; x++)
			{
				for (int y = SOURCE_SEPARATION; y < IMAGEHEIGHT - SOURCE_SEPARATION; y++)
				{
					if (SEARCH_ROI)
						if (!ROI_REGION[x, y])
							continue;

					if (IMAGE[x, y] < PIX_SAT)
						continue;

					if (SOURCE_BOOLEAN_MAP[x, y])
						continue;

					//map saturation island
					#pragma omp critical
					{
						if (!SOURCE_BOOLEAN_MAP[x, y] && IMAGE[x, y] >= PIX_SAT)
						{
							int Xmin = x, Xmax = x, Ymin = y, Ymax = y;
							MAPSATURATIONISLAND(x, y, src_index, Xmin, Xmax, Ymin, Ymax);

							if (Xmax - Xmin == 0)//single pixel, expand to 3 pixels
							{
								Xmin--;
								Xmax++;
							}
							else if (!JPMath::IsEven(Xmax - Xmin))//544-543 = 1 = 2 pixels, so make odd number of pixels...with Xmax I guess
								Xmax++;
							if (Ymax - Ymin == 0)//single pixel, expand to 3 pixels
							{
								Ymin--;
								Ymax++;
							}
							else if (!JPMath::IsEven(Ymax - Ymin))//544-543 = 1 = 2 pixels, so make odd number of pixels...with Xmax I guess
								Ymax++;

							array<double, 2>^ kernel = gcnew array<double, 2>(Xmax - Xmin + 1, Ymax - Ymin + 1);
							double x_centroid = 0, y_centroid = 0, kernel_sum = 0;
							for (int i = Xmin; i <= Xmax; i++)
								for (int j = Ymin; j <= Ymax; j++)
								{
									kernel[i - Xmin, j - Ymin] = IMAGE[i, j];
									x_centroid += (IMAGE[i, j] * double(i));
									y_centroid += (IMAGE[i, j] * double(j));
									kernel_sum += (IMAGE[i, j] /*- bg_est*/);
									/*IMAGE_KERNEL_BOOL_SOURCE[i, j] = true;
									IMAGE_KERNEL_INDEX_SOURCE[i, j] = src_index;*/
								}
							x_centroid /= kernel_sum;
							y_centroid /= kernel_sum;

							//MAPSATURATIONISLAND((int)x_centroid, (int)y_centroid, src_index, Xmin, Xmax, Ymin, Ymax);
							SOURCE_BOOLEAN_MAP[(int)x_centroid, (int)y_centroid] = true;
							SOURCE_INDEX_MAP[(int)x_centroid, (int)y_centroid] = src_index;

							src_index++;
							N_SATURATED++;
							Xs->Add(x_centroid);
							Ys->Add(y_centroid);
							Ks->Add(kernel_sum);
							Ps->Add(IMAGE[x, y]/*pixel*/);
							Bs->Add(0.0/*bg_est*/);

							if (SAVE_PS)
							{
								String^ file = SAVE_PS_FILENAME;
								int ind = file->LastIndexOf(".");//for saving PS
								file = String::Concat(file->Substring(0, ind), "_", Xs->Count.ToString("00000000"), ".fits");

								JPFITS::FITSImage^ f = gcnew JPFITS::FITSImage(file, kernel, false, false);
								f->WriteImage(TypeCode::Double, false);
							}
						}
					}
				}
			}

			/*JPFITS::FITSImage^ ff = gcnew FITSImage("C:\\Users\\Joseph E Postma\\Desktop\\atest.fits", IMAGE_KERNEL_INDEX_SOURCE, false);
			ff->WriteFile(TypeCode::Int32);*/
		}

		#pragma omp parallel for
		for (int x = SOURCE_SEPARATION; x < IMAGEWIDTH - SOURCE_SEPARATION; x++)
		{
			if (SHOWWAITBAR)
			{
				if (WAITBAR->DialogResult == DialogResult::Cancel)
					break;

				int nthreads = omp_get_num_threads();
				if (x < IMAGE->GetLength(0) / nthreads)
					if (x * 100 * nthreads / IMAGE->GetLength(0) > intprog)
					{
						intprog = x * 100 * nthreads / IMAGE->GetLength(0);
						BGWRKR->ReportProgress(intprog + 1);
					}
			}

			for (int y = SOURCE_SEPARATION; y < IMAGEHEIGHT - SOURCE_SEPARATION; y++)
			{
				if (SEARCH_ROI)
					if (!ROI_REGION[x, y])
						continue;

				if (SOURCE_BOOLEAN_MAP[x, y])
					continue;

				double bg_est = 0;
				if (AUTO_BG)
					bg_est = ESTIMATELOCALBACKGROUND(x, y, SOURCE_SEPARATION);

				double pixel = IMAGE[x, y] - bg_est;

				if (!THRESHHOLDS_AS_SN)
				{
					if (pixel < PIX_MIN || pixel > PIX_MAX)
						continue;
				}
				else//check as S/N
				{
					double Nbg = 0;
					if (bg_est < 1)
						Nbg = 1;
					else
						Nbg = Math::Sqrt(bg_est);
					double SNval = pixel / Nbg;
					if (SNval < PIX_MIN || SNval > PIX_MAX)
						continue;
				}

				bool brek = false;
				for (int i = x - SOURCE_SEPARATION; i <= x + SOURCE_SEPARATION; i++)
				{
					double sx2 = double(x - i);
					sx2 *= sx2;
					for (int j = y - SOURCE_SEPARATION; j <= y + SOURCE_SEPARATION; j++)
					{
						double sy = double(y - j);
						double r2 =  (sx2 +  sy * sy) / SSEP2;

						if (r2 > 1) //outside the source separation circle
							continue;

						if (IMAGE[i, j] - bg_est > pixel) // max failure, the pixel isn't the maximum in the source separation circle
						{
							brek = true;
							break;
						}

						if (r2 < 0.25 && SOURCE_BOOLEAN_MAP[i, j]) //a source was already found within the source separation
						{
							brek = true;
							break;
						}
					}
					if (brek)
						break;
				}
				if (brek)
					continue;

				//if got to here then x,y is a possible source depending on total sum of kernel
				array<double, 2>^ kernel = GetKernel(IMAGE, x, y, KERNEL_RADIUS);
				double kernel_sum = JPMath::Sum(kernel, false) - double(KERNEL_WIDTH*KERNEL_WIDTH)*bg_est;///square kernel sum
				
				/////do PSF kernel sum????????????????????????????????????????????????????????
				/*double kernel_psf_sum = 0, n_psf_pixels = 0;
				for (int i = x - KERNEL_RADIUS; i <= x + KERNEL_RADIUS; i++)
					for (int j = y - KERNEL_RADIUS; j <= y + KERNEL_RADIUS; j++)
					{
						double r2 = double((i - x) * (i - x) + (j - y) * (j - y));
						if (r2 > KRAD2)
							continue;

						kernel_psf_sum += IMAGE[i, j];
						n_psf_pixels++;
					}
				kernel_psf_sum -= (bg_est * n_psf_pixels);*/

				if (!THRESHHOLDS_AS_SN)
				{
					if (kernel_sum < KERNEL_MIN || kernel_sum > KERNEL_MAX)
						continue;
				}
				else//check as S/N
				{
					double Nbg = 0;
					if (bg_est < 1)
						Nbg = Math::Sqrt(double(KERNEL_WIDTH*KERNEL_WIDTH));
					else
						Nbg = Math::Sqrt(bg_est*double(KERNEL_WIDTH*KERNEL_WIDTH));
					double SNenergy = kernel_sum / Nbg;
					if (kernel_sum < KERNEL_MIN || kernel_sum > KERNEL_MAX)
						continue;
				}

				//if got to here then must centroid at this pixel
				double x_centroid, y_centroid;
				array<int>^ xdata = gcnew array<int>(KERNEL_WIDTH);
				array<int>^ ydata = gcnew array<int>(KERNEL_WIDTH);
				for (int i = -KERNEL_RADIUS; i <= KERNEL_RADIUS; i++)
				{
					xdata[i + KERNEL_RADIUS] = x + i;
					ydata[i + KERNEL_RADIUS] = y + i;
				}
				kernel = JPMath::MatrixSubScalar(kernel, bg_est, false);
				Centroid(xdata, ydata, kernel, x_centroid, y_centroid);

				int r_x_cent = (int)Math::Round(x_centroid);
				int r_y_cent = (int)Math::Round(y_centroid);

				#pragma omp critical
				{
					for (int ii = r_x_cent - KERNEL_RADIUS; ii <= r_x_cent + KERNEL_RADIUS; ii++)
					{
						double sx2 = double(r_x_cent - ii);
						sx2 *= sx2;
						for (int jj = r_y_cent - KERNEL_RADIUS; jj <= r_y_cent + KERNEL_RADIUS; jj++)
						{
							double sy = double(r_y_cent - jj);
							double r2 = (sx2 +  sy * sy);
							if (r2 > KRAD2)
								continue;

							if (ii > 0 && jj > 0 && ii < SOURCE_BOOLEAN_MAP->GetLength(0) && jj < SOURCE_BOOLEAN_MAP->GetLength(1))
							{
								SOURCE_BOOLEAN_MAP[ii, jj] = true;//this kernel radius position has a source detected
								SOURCE_INDEX_MAP[ii, jj] = src_index;//this is the source index of the given detection kernel radius position
							}
						}
					}

					src_index++;
					Xs->Add(x_centroid);
					Ys->Add(y_centroid);
					Ks->Add(kernel_sum);
					Ps->Add(pixel);
					Bs->Add(bg_est);

					if (SAVE_PS)
					{
						String^ file = SAVE_PS_FILENAME;
						int ind = file->LastIndexOf(".");//for saving PS
						file = String::Concat(file->Substring(0, ind), "_", Xs->Count.ToString("00000000"), ".fits");

						JPFITS::FITSImage^ f = gcnew JPFITS::FITSImage(file, kernel, false, false);
						f->WriteImage(TypeCode::Double, false);
					}
				}
			}
		}

		if (SHOWWAITBAR)
			if (WAITBAR->DialogResult == DialogResult::Cancel)
				return;

		this->N_SRC = Xs->Count;
		this->INITARRAYS();

		for (int i = 0; i < N_SRC; i++)
		{
			CENTROIDS_X[i] = Convert::ToDouble(Xs[i]);
			CENTROIDS_Y[i] = Convert::ToDouble(Ys[i]);
			CENTROIDS_AMPLITUDE[i] = Convert::ToDouble(Ps[i]);
			CENTROIDS_VOLUME[i] = Convert::ToDouble(Ks[i]);
			CENTROIDS_BGESTIMATE[i] = Convert::ToDouble(Bs[i]);
		}
		return;
	}
	//returned if after Source Extraction

	int intprog = 0;
	array<double>^ empty = gcnew array<double>(0);
	if (LBND == nullptr)
		LBND = gcnew array<double>(0);
	if (UBND == nullptr)
		UBND = gcnew array<double>(0);

	#pragma omp parallel for
	for (int k = 0; k < N_SRC; k++)
	{
		if (WAITBAR->DialogResult == DialogResult::Cancel)
		{
				BGWRKR->CancelAsync();
				break;
		}

		int nthreads = omp_get_num_threads();
		if (k < N_SRC / nthreads)
			if (k * 100 * nthreads / N_SRC > intprog)
			{
				intprog = k * 100 * nthreads / N_SRC;
				BGWRKR->ReportProgress(intprog + 100);
			}

		array<double, 2>^ kernel = GetKernel(IMAGE, int(CENTROIDS_X[k] + .5), int(CENTROIDS_Y[k] + .5), KERNEL_RADIUS);
		array<int>^ xcoords = gcnew array<int>(KERNEL_WIDTH);
		array<int>^ ycoords = gcnew array<int>(KERNEL_WIDTH);
		for (int i = 0; i < KERNEL_WIDTH; i++)
		{
			xcoords[i] = int(CENTROIDS_X[k] + .5) - KERNEL_RADIUS + i;
			ycoords[i] = int(CENTROIDS_Y[k] + .5) - KERNEL_RADIUS + i;
		}
		
		array<double, 2>^ fit_resid = gcnew array<double, 2>(KERNEL_WIDTH, KERNEL_WIDTH);
		array<double>^ P0;
		array<double>^ lb = gcnew array<double>(LBND->Length);
		array<double>^ ub = gcnew array<double>(UBND->Length);
		if (LBND->Length > 0)//set bounds to make sense
		{
			lb[0] = 0;
			lb[1] = CENTROIDS_X[k] - KERNEL_RADIUS;
			lb[2] = CENTROIDS_Y[k] - KERNEL_RADIUS;
			for (int i = 3; i < LBND->Length; i++)
				lb[i] = LBND[i];
		}
		if (UBND->Length > 0)//set bounds to make sense
		{
			ub[0] = Math::Abs(CENTROIDS_AMPLITUDE[k] * 2);
			ub[1] = CENTROIDS_X[k] + KERNEL_RADIUS;
			ub[2] = CENTROIDS_Y[k] + KERNEL_RADIUS;
			for (int i = 3; i < UBND->Length; i++)
				ub[i] = UBND[i];
		}

		if (::Convert::ToInt32(e->Argument) == 2)//Fit circular Gaussian
		{
			if (PINI != nullptr)
				P0 = gcnew array<double, 1>(5) { Math::Abs(CENTROIDS_AMPLITUDE[k]), CENTROIDS_X[k], CENTROIDS_Y[k], PINI[3], PINI[4] };
			else
				P0 = gcnew array<double, 1>(5) { Math::Abs(CENTROIDS_AMPLITUDE[k]), CENTROIDS_X[k], CENTROIDS_Y[k], 2, CENTROIDS_BGESTIMATE[k] };

			JPMath::Fit_Gaussian2d(xcoords, ycoords, kernel, P0, lb, ub, empty, fit_resid);
			FITS_VOLUME[k] = 2 * Math::PI * P0[0] * P0[3] * P0[3];
			FITS_FWHM_X[k] = 2.355 * P0[3];
			FITS_FWHM_Y[k] = FITS_FWHM_X[k];
			for (int i = 0; i < P0->Length; i++)
				FITS_PARAMS[i, k] = P0[i];
		}

		if (::Convert::ToInt32(e->Argument) == 3)//Fit elliptical Gaussian
		{
			if (PINI != nullptr)
				P0 = gcnew array<double, 1>(7) { Math::Abs(CENTROIDS_AMPLITUDE[k]), CENTROIDS_X[k], CENTROIDS_Y[k], PINI[3], PINI[4], PINI[5], PINI[6] };
			else
				P0 = gcnew array<double, 1>(7) { Math::Abs(CENTROIDS_AMPLITUDE[k]), CENTROIDS_X[k], CENTROIDS_Y[k], 0, 2, 2, CENTROIDS_BGESTIMATE[k] };

			JPMath::Fit_Gaussian2d(xcoords, ycoords, kernel, P0, lb, ub, empty, fit_resid);
			FITS_VOLUME[k] = 2 * Math::PI * P0[0] * P0[4] * P0[5];
			FITS_FWHM_X[k] = 2.355 * P0[4];
			FITS_FWHM_Y[k] = 2.355 * P0[5];
			FITS_PHI[k] = P0[3];
			for (int i = 0; i < P0->Length; i++)
				FITS_PARAMS[i, k] = P0[i];
		}

		if (::Convert::ToInt32(e->Argument) == 4)//Fit circular Moffat
		{
			if (PINI != nullptr)
				P0 = gcnew array<double, 1>(6) { Math::Abs(CENTROIDS_AMPLITUDE[k]), CENTROIDS_X[k], CENTROIDS_Y[k], PINI[3], PINI[4], PINI[5] };
			else
				P0 = gcnew array<double, 1>(6) { Math::Abs(CENTROIDS_AMPLITUDE[k]), CENTROIDS_X[k], CENTROIDS_Y[k], 2, 2, CENTROIDS_BGESTIMATE[k] };

			JPMath::Fit_Moffat2d(xcoords, ycoords, kernel, P0, lb, ub, empty, fit_resid);
			FITS_VOLUME[k] = Math::PI * P0[3] * P0[3] * P0[0] / (P0[4] - 1);
			FITS_FWHM_X[k] = 2 * P0[3] * Math::Sqrt(Math::Pow(2, 1 / (P0[4])) - 1);
			FITS_FWHM_Y[k] = FITS_FWHM_X[k];
			for (int i = 0; i < P0->Length; i++)
				FITS_PARAMS[i, k] = P0[i];
		}

		if (::Convert::ToInt32(e->Argument) == 5)//Fit elliptical Moffat
		{
			if (PINI != nullptr)
				P0 = gcnew array<double, 1>(8) { Math::Abs(CENTROIDS_AMPLITUDE[k]), CENTROIDS_X[k], CENTROIDS_Y[k], PINI[3], PINI[4], PINI[5], PINI[6], PINI[7] };
			else
				P0 = gcnew array<double, 1>(8) { Math::Abs(CENTROIDS_AMPLITUDE[k]), CENTROIDS_X[k], CENTROIDS_Y[k], 0, 2, 2, 2, CENTROIDS_BGESTIMATE[k] };

			JPMath::Fit_Moffat2d(xcoords, ycoords, kernel, P0, lb, ub, empty, fit_resid);
			FITS_VOLUME[k] = Math::PI * P0[4] * P0[5] * P0[0] / (P0[6] - 1);
			FITS_FWHM_X[k] = 2 * P0[4] * Math::Sqrt(Math::Pow(2, 1 / (P0[6])) - 1);
			FITS_FWHM_Y[k] = 2 * P0[5] * Math::Sqrt(Math::Pow(2, 1 / (P0[6])) - 1);
			FITS_PHI[k] = P0[3];
			for (int i = 0; i < P0->Length; i++)
				FITS_PARAMS[i, k] = P0[i];
		}

		FITS_AMPLITUDE[k] = P0[0];
		FITS_X[k] = P0[1];
		FITS_Y[k] = P0[2];
		FITS_BGESTIMATE[k] = P0[P0->Length - 1];

		double chisq_norm = 0;
		for (int i = 0; i < KERNEL_WIDTH; i++)
			for (int j = 0; j < KERNEL_WIDTH; j++)
			{
				if (kernel[i, j] - P0[P0->Length - 1] == 0)
					chisq_norm += fit_resid[i, j] * fit_resid[i, j];
				else
					chisq_norm += fit_resid[i, j] * fit_resid[i, j] / Math::Abs(kernel[i, j] - P0[P0->Length - 1]);
			}
		chisq_norm /= (kernel->Length - P0->Length);
		FITS_CHISQNORM[k] = chisq_norm;
	}
}

void JPFITS::SourceExtractor::BGWRKR_ProgressChanged(System::Object^  sender, System::ComponentModel::ProgressChangedEventArgs^  e)
{
	if (e->ProgressPercentage >= 0 && e->ProgressPercentage <= 100)
	{
		WAITBAR->ProgressBar->Value = e->ProgressPercentage;
		WAITBAR->TextMsg->Text = "Scanned " + e->ProgressPercentage.ToString() + "% of the image...";
	}

	/*if (e->ProgressPercentage < 0)
	{
		WAITBAR->ProgressBar->Value = -e->ProgressPercentage;
		WAITBAR->TextMsg->Text = "Please wait while I remove adjacent sources..." + (-e->ProgressPercentage).ToString() + "%";
	}*/

	if (e->ProgressPercentage > 100)
	{
		WAITBAR->ProgressBar->Value = e->ProgressPercentage - 100;
		WAITBAR->TextMsg->Text = "Fitted " + (e->ProgressPercentage - 100).ToString() + "% of the sources...";
	}

	WAITBAR->Refresh();
}

void JPFITS::SourceExtractor::BGWRKR_RunWorkerCompleted(System::Object^  sender, System::ComponentModel::RunWorkerCompletedEventArgs^  e)
{
	if (WAITBAR != nullptr)
		WAITBAR->Close();
}

void JPFITS::SourceExtractor::Extract_Attempt_N_Sources(int N, array<double, 2>^ image, double pix_saturation, double pix_min, double pix_max, double kernel_min, double kernel_max, bool threshholds_as_SN, int kernel_radius, int source_separation, bool auto_background, String^ kernel_filename_template, array<bool, 2>^ ROI_region, bool show_waitbar)
{
	JPFITS::FITSImage^ FITS = gcnew FITSImage("", image, true, true);

	double immax = JPMath::Max(image, true);
	double pixthresh = immax / 16;
	double div = 2;
	double amp = pixthresh;
	int PSEiters = 0;
	int maxPSEiters = 20;
	int nPSEpts_min = N;
	int nPSEpts_max = N + 1;
	while (this->N_Sources < nPSEpts_min || this->N_Sources > nPSEpts_max)
	{
		PSEiters++;
		if (PSEiters > maxPSEiters)
			break;

		if (this->N_SaturatedSources >= nPSEpts_min)
			break;

		if (this->N_Sources >= nPSEpts_min)
			break;

		Extract_Sources(image, pix_saturation, pixthresh, pix_max, kernel_min, kernel_max, threshholds_as_SN, kernel_radius, source_separation, auto_background, kernel_filename_template, ROI_region, show_waitbar);

		if (this->N_Sources < nPSEpts_min)
			pixthresh -= amp / div;
		if (this->N_Sources > nPSEpts_max)
			pixthresh += amp / div;
		div *= 2;

		if (pixthresh < pix_min)
		{
			Extract_Sources(image, pix_saturation, pix_min, pix_max, kernel_min, kernel_max, threshholds_as_SN, kernel_radius, source_separation, auto_background, kernel_filename_template, ROI_region, show_waitbar);
			break;
		}
	}
	if (this->N_Sources > nPSEpts_min)
		this->ClipToNBrightest(nPSEpts_min);
}

void JPFITS::SourceExtractor::MAPSATURATIONISLAND(int X, int Y, int sourceindex, int &xmin, int &xmax, int &ymin, int &ymax)
{
	SOURCE_BOOLEAN_MAP[X, Y] = true;
	SOURCE_INDEX_MAP[X, Y] = sourceindex;

	if (X < xmin)
		xmin = X;
	if (X > xmax)
		xmax = X;
	if (Y < ymin)
		ymin = Y;
	if (Y > ymax)
		ymax = Y;

	for (int x = X - 1; x <= X + 1; x++)
		for (int y = Y - 1; y <= Y + 1; y++)
			if (SAFETOMAPSATURATION(x, y))
				MAPSATURATIONISLAND(x, y, sourceindex, xmin, xmax, ymin, ymax);
}

bool JPFITS::SourceExtractor::SAFETOMAPSATURATION(int x, int y)
{
	return (x >= 0) && (x < IMAGEWIDTH) && (y >= 0) && (y < IMAGEHEIGHT) && (IMAGE[x, y] > PIX_SAT) && /*!SOURCE_BOOLEAN_MAP[x, y]*/ (SOURCE_INDEX_MAP[x, y] == -1);
}

void JPFITS::SourceExtractor::DEMAP(int X, int Y, int sourceindex)
{
	SOURCE_BOOLEAN_MAP[X, Y] = false;
	SOURCE_INDEX_MAP[X, Y] = -1;

	for (int x = X - 1; x <= X + 1; x++)
		for (int y = Y - 1; y <= Y + 1; y++)
			if (SAFETODEMAP(x, y, sourceindex))
				DEMAP(x, y, sourceindex);
}

bool JPFITS::SourceExtractor::SAFETODEMAP(int x, int y, int sourceindex)
{
	return (x >= 0) && (x < IMAGEWIDTH) && (y >= 0) && (y < IMAGEHEIGHT) && SOURCE_BOOLEAN_MAP[x, y] && (sourceindex == SOURCE_INDEX_MAP[x, y]);
}

void JPFITS::SourceExtractor::REMAP(int X, int Y, int oldindex, int newindex)
{
	if (oldindex == newindex)
		return;

	SOURCE_INDEX_MAP[X, Y] = newindex;

	for (int x = X - 1; x <= X + 1; x++)
		for (int y = Y - 1; y <= Y + 1; y++)
			if (SAFETOREMAP(x, y, oldindex))
				REMAP(x, y, oldindex, newindex);
}

bool JPFITS::SourceExtractor::SAFETOREMAP(int x, int y, int oldindex)
{
	return (x >= 0) && (x < IMAGEWIDTH) && (y >= 0) && (y < IMAGEHEIGHT) && SOURCE_BOOLEAN_MAP[x, y] && (oldindex == SOURCE_INDEX_MAP[x, y]);
}

void JPFITS::SourceExtractor::Extract_Source_LSFits_Gaussian_Circular(array<double>^ Pinit, array<double>^ LBnds, array<double>^ UBnds)
{
	//G = P(0) * exp( -((X-P(1)).^2 + (Y-P(2)).^2 ) / (2*P(3)^2)) + P(4);

	FITS_PARAMS = gcnew array<double, 2>(5, N_SRC);
	FITTED = true;
	FIT_EQUATION = "Gaussian (Circular): P(0) * exp( -((X-P(1)).^2 + (Y-P(2)).^2 ) / (2*P(3)^2)) + P(4)";
	//VIEWFITS = view;
	//P_INIT = P0;
	LBND = LBnds;
	UBND = UBnds;
	PINI = Pinit;

	WAITBAR = gcnew JPWaitBar::WaitBar();
	WAITBAR->ProgressBar->Maximum = 100;
	WAITBAR->Text = "Fitting Sources...";
	BGWRKR->RunWorkerAsync(2);
	WAITBAR->ShowDialog();

	if (WAITBAR->DialogResult == DialogResult::Cancel)
		FITTED = false;
}

void JPFITS::SourceExtractor::Extract_Source_LSFits_Gaussian_Elliptical(array<double>^ Pinit, array<double>^ LBnds, array<double>^ UBnds)// 2-D Elliptical Gaussian
{
	//G = P(0) * exp( -((x-P(1))*cosd(P(3)) + (y-P(2))*sind(P(3))).^2 / (2*P(4)^2) - ( -(x-P(1))*sind(P(3)) + (y-P(2))*cosd(P(3))).^2 / (2*P(5)^2) ) + P(6);

	FITS_PARAMS = gcnew array<double, 2>(7, N_SRC);
	FITTED = true;
	FIT_EQUATION = "Gaussian (Elliptical): P(0) * exp( -((x-P(1))*cos(P(3)) + (y-P(2))*sin(P(3))).^2 / (2*P(4)^2) - ( -(x-P(1))*sin(P(3)) + (y-P(2))*cos(P(3))).^2 / (2*P(5)^2) ) + P(6)";
	/*VIEWFITS = view;
	P_INIT = P0;*/
	LBND = LBnds;
	UBND = UBnds;
	PINI = Pinit;

	WAITBAR = gcnew JPWaitBar::WaitBar();
	WAITBAR->ProgressBar->Maximum = 100;
	WAITBAR->Text = "Fitting Sources...";
	BGWRKR->RunWorkerAsync(3);
	WAITBAR->ShowDialog();

	if (WAITBAR->DialogResult == DialogResult::Cancel)
		FITTED = false;
}

void JPFITS::SourceExtractor::Extract_Source_LSFits_Moffat_Circular(array<double>^ Pinit, array<double>^ LBnds, array<double>^ UBnds)// 2-D Circular Moffat
{
	// M = P(0) * ( 1 + { (X-P(1))^2 + (Y-P(2))^2 } / P(3)^2 ) ^ (-P(4)) + P(5)

	FITS_PARAMS = gcnew array<double, 2>(6, N_SRC);
	FITTED = true;
	FIT_EQUATION = "Moffat (Circular): P(0) * ( 1 + { (X-P(1))^2 + (Y-P(2))^2 } / P(3)^2 ) ^ (-P(4)) + P(5)";
	/*VIEWFITS = view;
	P_INIT = P0;*/
	LBND = LBnds;
	UBND = UBnds;
	PINI = Pinit;

	WAITBAR = gcnew JPWaitBar::WaitBar();
	WAITBAR->ProgressBar->Maximum = 100;
	WAITBAR->Text = "Fitting Sources...";
	BGWRKR->RunWorkerAsync(4);
	WAITBAR->ShowDialog();

	if (WAITBAR->DialogResult == DialogResult::Cancel)
		FITTED = false;
}

void JPFITS::SourceExtractor::Extract_Source_LSFits_Moffat_Elliptical(array<double>^ Pinit, array<double>^ LBnds, array<double>^ UBnds)// 2-D Elliptical Moffat
{
	//M = P(0) * ( 1 + { ((X-P(1))*cosd(P(3)) + (Y-P(2))*sind(P(3)))^2 } / P(4)^2 + { (-(X-P(1))*sind(P(3)) + (Y-P(2))*cosd(P(3)))^2 } / P(5)^2 ) ^ (-P(6)) + P(7);

	FITS_PARAMS = gcnew array<double, 2>(8, N_SRC);
	FITTED = true;
	FIT_EQUATION = "Moffat (Elliptical): P(0) * ( 1 + { ((X-P(1))*cos(P(3)) + (Y-P(2))*sin(P(3)))^2 } / P(4)^2 + { (-(X-P(1))*sin(P(3)) + (Y-P(2))*cos(P(3)))^2 } / P(5)^2 ) ^ (-P(6)) + P(7)";
	/*VIEWFITS = view;
	P_INIT = P0;*/
	LBND = LBnds;
	UBND = UBnds;
	PINI = Pinit;

	WAITBAR = gcnew JPWaitBar::WaitBar();
	WAITBAR->ProgressBar->Maximum = 100;
	WAITBAR->Text = "Fitting Sources...";
	BGWRKR->RunWorkerAsync(5);
	WAITBAR->ShowDialog();

	if (WAITBAR->DialogResult == DialogResult::Cancel)
		FITTED = false;
}

inline double JPFITS::SourceExtractor::ESTIMATELOCALBACKGROUND(int x, int y, int HW)
{
	double fcmin = IMAGE[x - HW, y - HW];
	if (IMAGE[x - HW, y + HW] < fcmin)
		fcmin = IMAGE[x - HW, y + HW];
	if (IMAGE[x + HW, y - HW] < fcmin)
		fcmin = IMAGE[x + HW, y - HW];
	if (IMAGE[x + HW, y + HW] < fcmin)
		fcmin = IMAGE[x + HW, y + HW];

	double fcmin2 = IMAGE[x - HW, y - HW];
	if (fcmin2 != fcmin && IMAGE[x - HW, y + HW] < fcmin2)
		fcmin2 = IMAGE[x - HW, y + HW];
	if (fcmin2 != fcmin && IMAGE[x + HW, y - HW] < fcmin2)
		fcmin2 = IMAGE[x + HW, y - HW];
	if (fcmin2 != fcmin && IMAGE[x + HW, y + HW] < fcmin2)
		fcmin2 = IMAGE[x + HW, y + HW];

	return fcmin2;

	/*int xmin = x - HW, xmax = x + HW, ymin = y - HW, ymax = y + HW;

	array<double>^ cmins = gcnew array<double>(4) { IMAGE[xmin, ymin], IMAGE[xmin, ymax], IMAGE[xmax, ymin], IMAGE[xmax, ymax] };
	Array::Sort(cmins);
	return cmins[1];*/
}

void JPFITS::SourceExtractor::Generate_Source_RADec_Coords(JPFITS::WorldCoordinateSolution^ wcs)
{
	for (int i = 0; i < N_SRC; i++)
	{
		double a, d;
		String^ RAhms;
		String^ DECdamas;

		wcs->Get_Coordinate(CENTROIDS_X[i], CENTROIDS_Y[i], true, "TAN", a, d, RAhms, DECdamas);

		CENTROIDS_RA_DEG[i] = a;
		CENTROIDS_RA_HMS[i] = RAhms;
		CENTROIDS_DEC_DEG[i] = d;
		CENTROIDS_DEC_DMS[i] = DECdamas;

		if (FITTED)
		{
			wcs->Get_Coordinate(FITS_X[i], FITS_Y[i], true, "TAN", a, d, RAhms, DECdamas);

			FITS_RA_DEG[i] = a;
			FITS_RA_HMS[i] = RAhms;
			FITS_DEC_DEG[i] = d;
			FITS_DEC_DMS[i] = DECdamas;
		}
	}

	WCS_GENERATED = true;
}

void JPFITS::SourceExtractor::Save_Source_Table(String^ delimit)
{
	SaveFileDialog^ sfd = gcnew ::SaveFileDialog();
	sfd->Filter = "Delimited Text file (*.txt)|*.txt";
	if (sfd->ShowDialog() == ::DialogResult::Cancel)
		return;

	this->GENERATEPSETABLE();

	if (delimit == "tab")
		delimit = "\t";

	String^ line;
	StreamWriter^ sw = gcnew StreamWriter(sfd->FileName);
	for (int j = 0; j < PSE_TABLE->GetLength(1); j++)
	{
		line = "";
		for (int i = 0; i < PSE_TABLE->GetLength(0); i++)
			line += PSE_TABLE[i, j] + delimit;
		line = line->Substring(0, line->LastIndexOf(delimit) + 1);
		sw->WriteLine(line);
	}
	sw->Close();
}

void JPFITS::SourceExtractor::GENERATEPSETABLE()
{
	if (!FITTED)
		PSE_TABLE = gcnew array<String^, 2>(9, N_SRC + 1);
	else
		PSE_TABLE = gcnew array<String^, 2>(22 + FITS_PARAMS->GetLength(0), N_SRC + 1);

	int c = 0;

	PSE_TABLE[c++, 0] = "Pixel X-Centroid";
	PSE_TABLE[c++, 0] = "Pixel Y-Centroid";
	PSE_TABLE[c++, 0] = "Pixel Amplitude";
	PSE_TABLE[c++, 0] = "Pixel Volume";
	PSE_TABLE[c++, 0] = "Pixel Background";
	PSE_TABLE[c++, 0] = "Pixel Right Ascension (deg)";
	PSE_TABLE[c++, 0] = "Pixel Declination (deg)";
	PSE_TABLE[c++, 0] = "Pixel Right Ascension (hms)";
	PSE_TABLE[c++, 0] = "Pixel Declination (d ' '')";

	if (FITTED)
	{
		PSE_TABLE[c++, 0] = "Fit X-Centroid";
		PSE_TABLE[c++, 0] = "Fit Y-Centroid";
		PSE_TABLE[c++, 0] = "Fit Amplitude";
		PSE_TABLE[c++, 0] = "Fit Volume";
		PSE_TABLE[c++, 0] = "Fit Background";
		PSE_TABLE[c++, 0] = "Fit Right Ascension (deg)";
		PSE_TABLE[c++, 0] = "Fit Declination (deg)";
		PSE_TABLE[c++, 0] = "Fit Right Ascension (hms)";
		PSE_TABLE[c++, 0] = "Fit Declination (d ' '')";
		PSE_TABLE[c++, 0] = "Fit FWHM_X";
		PSE_TABLE[c++, 0] = "Fit FWHM_Y";
		PSE_TABLE[c++, 0] = "Fit Phi";
		PSE_TABLE[c++, 0] = "Fit ChiSqNorm";

		for (int j = 0; j < FITS_PARAMS->GetLength(0); j++)
			PSE_TABLE[c++, 0] = "Fit P(" + j + ")";
	}

	for (int i = 0; i < N_SRC; i++)
	{
		c = 0;
		PSE_TABLE[c++, i + 1] = CENTROIDS_X[i].ToString();
		PSE_TABLE[c++, i + 1] = CENTROIDS_Y[i].ToString();
		PSE_TABLE[c++, i + 1] = CENTROIDS_AMPLITUDE[i].ToString();
		PSE_TABLE[c++, i + 1] = CENTROIDS_VOLUME[i].ToString();
		PSE_TABLE[c++, i + 1] = CENTROIDS_BGESTIMATE[i].ToString();
		PSE_TABLE[c++, i + 1] = CENTROIDS_RA_DEG[i].ToString();
		PSE_TABLE[c++, i + 1] = CENTROIDS_DEC_DEG[i].ToString();
		PSE_TABLE[c++, i + 1] = CENTROIDS_RA_HMS[i];
		PSE_TABLE[c++, i + 1] = CENTROIDS_DEC_DMS[i];

		if (FITTED)
		{
			PSE_TABLE[c++, i + 1] = FITS_X[i].ToString();
			PSE_TABLE[c++, i + 1] = FITS_Y[i].ToString();
			PSE_TABLE[c++, i + 1] = FITS_AMPLITUDE[i].ToString();
			PSE_TABLE[c++, i + 1] = FITS_VOLUME[i].ToString();
			PSE_TABLE[c++, i + 1] = FITS_BGESTIMATE[i].ToString();
			PSE_TABLE[c++, i + 1] = FITS_RA_DEG[i].ToString();
			PSE_TABLE[c++, i + 1] = FITS_DEC_DEG[i].ToString();
			PSE_TABLE[c++, i + 1] = FITS_RA_HMS[i];
			PSE_TABLE[c++, i + 1] = FITS_DEC_DMS[i];
			PSE_TABLE[c++, i + 1] = FITS_FWHM_X[i].ToString();
			PSE_TABLE[c++, i + 1] = FITS_FWHM_Y[i].ToString();
			PSE_TABLE[c++, i + 1] = FITS_PHI[i].ToString();
			PSE_TABLE[c++, i + 1] = FITS_CHISQNORM[i].ToString();

			for (int j = 0; j < FITS_PARAMS->GetLength(0); j++)
				PSE_TABLE[c++, i + 1] = FITS_PARAMS[j, i].ToString();
		}
	}
}

void JPFITS::SourceExtractor::INITARRAYS()
{
	CENTROIDS_X = gcnew array<double, 1>(N_SRC);
	CENTROIDS_Y = gcnew array<double, 1>(N_SRC);
	FITS_X = gcnew array<double, 1>(N_SRC);
	FITS_Y = gcnew array<double, 1>(N_SRC);
	FITS_FWHM_X = gcnew array<double, 1>(N_SRC);
	FITS_FWHM_Y = gcnew array<double, 1>(N_SRC);
	FITS_PHI = gcnew array<double, 1>(N_SRC);
	FITS_CHISQNORM = gcnew array<double, 1>(N_SRC);
	CENTROIDS_RA_DEG = gcnew array<double, 1>(N_SRC);
	CENTROIDS_RA_HMS = gcnew array<String^, 1>(N_SRC);
	CENTROIDS_DEC_DEG = gcnew array<double, 1>(N_SRC);
	CENTROIDS_DEC_DMS = gcnew array<String^, 1>(N_SRC);
	CENTROIDS_AMPLITUDE = gcnew array<double, 1>(N_SRC);
	CENTROIDS_VOLUME = gcnew array<double, 1>(N_SRC);
	CENTROIDS_BGESTIMATE = gcnew array<double, 1>(N_SRC);
	FITS_AMPLITUDE = gcnew array<double, 1>(N_SRC);
	FITS_VOLUME = gcnew array<double, 1>(N_SRC);
	FITS_BGESTIMATE = gcnew array<double, 1>(N_SRC);
	FITS_RA_DEG = gcnew array<double, 1>(N_SRC);
	FITS_RA_HMS = gcnew array<String^, 1>(N_SRC);
	FITS_DEC_DEG = gcnew array<double, 1>(N_SRC);
	FITS_DEC_DMS = gcnew array<String^, 1>(N_SRC);
}

void JPFITS::SourceExtractor::ClipToNBrightest(int NBright)
{
	if (NBright >= N_SRC)
		return;

	array<double>^ volkey = gcnew array<double>(N_SRC);
	Array::Copy(CENTROIDS_VOLUME, volkey, N_SRC);

	array<int>^ indices = gcnew array<int>(N_SRC);
	for (int i = 0; i < N_SRC; i++)
	{
		indices[i] = i; //IMAGE_KERNEL_INDEX_SOURCE[(int)CENTROIDS_X[i], (int)CENTROIDS_Y[i]];
		if (indices[i] == -1)
			MessageBox::Show("-1 a " + CENTROIDS_X[i] + " " + Math::Round(CENTROIDS_X[i]) + " " + CENTROIDS_Y[i] + " " + Math::Round(CENTROIDS_Y[i]));
	}

	Array::Sort(volkey, indices);//by increasing brightness
	Array::Reverse(indices);//by decreasing brightness; now all location at indices at index >= NBright are no longer wanted

	double dum;
	String^ dumstr;
	for (int i = 0; i < NBright; i++)
	{
		dum = CENTROIDS_X[i];
		CENTROIDS_X[i] = CENTROIDS_X[indices[i]];
		CENTROIDS_X[indices[i]] = dum;

		dum = CENTROIDS_Y[i];
		CENTROIDS_Y[i] = CENTROIDS_Y[indices[i]];
		CENTROIDS_Y[indices[i]] = dum;

		REMAP(int(CENTROIDS_X[i]), int(CENTROIDS_Y[i]), indices[i], i);
		REMAP(int(CENTROIDS_X[indices[i]]), int(CENTROIDS_Y[indices[i]]), i, indices[i]);

		dum = FITS_X[i];
		FITS_X[i] = FITS_X[indices[i]];
		FITS_X[indices[i]] = dum;

		dum = FITS_Y[i];
		FITS_Y[i] = FITS_Y[indices[i]];
		FITS_Y[indices[i]] = dum;

		dum = FITS_FWHM_X[i];
		FITS_FWHM_X[i] = FITS_FWHM_X[indices[i]];
		FITS_FWHM_X[indices[i]] = dum;

		dum = FITS_FWHM_Y[i];
		FITS_FWHM_Y[i] = FITS_FWHM_Y[indices[i]];
		FITS_FWHM_Y[indices[i]] = dum;

		dum = FITS_PHI[i];
		FITS_PHI[i] = FITS_PHI[indices[i]];
		FITS_PHI[indices[i]] = dum;

		dum = FITS_CHISQNORM[i];
		FITS_CHISQNORM[i] = FITS_CHISQNORM[indices[i]];
		FITS_CHISQNORM[indices[i]] = dum;

		dum = CENTROIDS_RA_DEG[i];
		CENTROIDS_RA_DEG[i] = CENTROIDS_RA_DEG[indices[i]];
		CENTROIDS_RA_DEG[indices[i]] = dum;

		dumstr = CENTROIDS_RA_HMS[i];
		CENTROIDS_RA_HMS[i] = CENTROIDS_RA_HMS[indices[i]];
		CENTROIDS_RA_HMS[indices[i]] = dumstr;

		dum = CENTROIDS_DEC_DEG[i];
		CENTROIDS_DEC_DEG[i] = CENTROIDS_DEC_DEG[indices[i]];
		CENTROIDS_DEC_DEG[indices[i]] = dum;

		dumstr = CENTROIDS_DEC_DMS[i];
		CENTROIDS_DEC_DMS[i] = CENTROIDS_DEC_DMS[indices[i]];
		CENTROIDS_DEC_DMS[indices[i]] = dumstr;

		dum = CENTROIDS_AMPLITUDE[i];
		CENTROIDS_AMPLITUDE[i] = CENTROIDS_AMPLITUDE[indices[i]];
		CENTROIDS_AMPLITUDE[indices[i]] = dum;

		dum = CENTROIDS_BGESTIMATE[i];
		CENTROIDS_BGESTIMATE[i] = CENTROIDS_BGESTIMATE[indices[i]];
		CENTROIDS_BGESTIMATE[indices[i]] = dum;

		dum = FITS_AMPLITUDE[i];
		FITS_AMPLITUDE[i] = FITS_AMPLITUDE[indices[i]];
		FITS_AMPLITUDE[indices[i]] = dum;

		dum = FITS_VOLUME[i];
		FITS_VOLUME[i] = FITS_VOLUME[indices[i]];
		FITS_VOLUME[indices[i]] = dum;

		dum = FITS_BGESTIMATE[i];
		FITS_BGESTIMATE[i] = FITS_BGESTIMATE[indices[i]];
		FITS_BGESTIMATE[indices[i]] = dum;

		dum = FITS_RA_DEG[i];
		FITS_RA_DEG[i] = FITS_RA_DEG[indices[i]];
		FITS_RA_DEG[indices[i]] = dum;

		dumstr = FITS_RA_HMS[i];
		FITS_RA_HMS[i] = FITS_RA_HMS[indices[i]];
		FITS_RA_HMS[indices[i]] = dumstr;

		dum = FITS_DEC_DEG[i];
		FITS_DEC_DEG[i] = FITS_DEC_DEG[indices[i]];
		FITS_DEC_DEG[indices[i]] = dum;

		dumstr = FITS_DEC_DMS[i];
		FITS_DEC_DMS[i] = FITS_DEC_DMS[indices[i]];
		FITS_DEC_DMS[indices[i]] = dumstr;

		dum = CENTROIDS_VOLUME[i];
		CENTROIDS_VOLUME[i] = CENTROIDS_VOLUME[indices[i]];
		CENTROIDS_VOLUME[indices[i]] = dum;
	}

	for (int i = NBright; i < N_SRC; i++)//all location at indices[i] where i >= NBright are no longer wanted
		DEMAP((int)CENTROIDS_X[i], (int)CENTROIDS_Y[i], SOURCE_INDEX_MAP[(int)CENTROIDS_X[i], (int)CENTROIDS_Y[i]]);

	Array::Resize(CENTROIDS_X, NBright);
	Array::Resize(CENTROIDS_Y, NBright);
	Array::Resize(FITS_X, NBright);
	Array::Resize(FITS_Y, NBright);
	Array::Resize(FITS_FWHM_X, NBright);
	Array::Resize(FITS_FWHM_Y, NBright);
	Array::Resize(FITS_PHI, NBright);
	Array::Resize(FITS_CHISQNORM, NBright);
	Array::Resize(CENTROIDS_RA_DEG, NBright);
	Array::Resize(CENTROIDS_RA_HMS, NBright);
	Array::Resize(CENTROIDS_DEC_DEG, NBright);
	Array::Resize(CENTROIDS_DEC_DMS, NBright);
	Array::Resize(CENTROIDS_AMPLITUDE, NBright);
	Array::Resize(CENTROIDS_VOLUME, NBright);
	Array::Resize(CENTROIDS_BGESTIMATE, NBright);
	Array::Resize(FITS_AMPLITUDE, NBright);
	Array::Resize(FITS_VOLUME, NBright);
	Array::Resize(FITS_BGESTIMATE, NBright);
	Array::Resize(FITS_RA_DEG, NBright);
	Array::Resize(FITS_RA_HMS, NBright);
	Array::Resize(FITS_DEC_DEG, NBright);
	Array::Resize(FITS_DEC_DMS, NBright);
	N_SRC = NBright;

	/*JPFITS::FITSImage^ ff = gcnew FITSImage("C:\\Users\\Joseph E Postma\\Desktop\\test.fits", IMAGE_KERNEL_INDEX_SOURCE, false);
	ff->WriteFile(TypeCode::Int32);*/
}

