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

#include "StdAfx.h"
#include "JPFITS.h"

array<double, 2>^ JPFITS::JPBitMap::Bin(cli::array<double, 2> ^data, int Nx, int Ny)
{
	int Lx = data->GetLength(0) / Nx;
	int Ly = data->GetLength(1) / Ny;
	array<double, 2>^ result = gcnew array<double, 2>(Lx, Ly);
	double inv_size = 1 / double(Nx * Ny), s;

	#pragma omp parallel for private(s)
	for (int i = 0; i < Lx; i++)
		for (int j = 0; j < Ly; j++)
		{
			s = 0;
			for (int k = i * Nx; k < i*Nx + Nx; k++)
				for (int l = j * Ny; l < j*Ny + Ny; l++)
					s += data[k, l];

			result[i, j] = s * inv_size;
		}

	return result;
}

Bitmap^ JPFITS::JPBitMap::ArrayToBmp(cli::array<double, 2>^ image, int scaling, int colour, bool invert, array<double, 1>^ DImCLim, int WinWidth, int WinHeight, bool invertYaxis)
{
	/*if (image->GetLength(0) > 8192 || image->GetLength(1) > 8192)
	{
		int Nx = 1;
		int Ny = 1;
		if (image->GetLength(0) > 8192)
			Nx = image->GetLength(0)/8192;
		if (image->GetLength(1) > 8192)
			Ny = image->GetLength(1)/8192;
		image = Bin(image, Nx, Ny);
	}*/

	if (image->GetLength(0) > WinWidth * 2 || image->GetLength(1) > WinHeight * 2)
	{
		int Nx = 1;
		int Ny = 1;
		if (image->GetLength(0) > WinWidth * 2)
			Nx = image->GetLength(0) / WinWidth;
		if (image->GetLength(1) > WinHeight * 2)
			Ny = image->GetLength(1) / WinHeight;
		if (Nx > 1 || Ny > 1)
			image = Bin(image, Nx, Ny);
	}

	Bitmap^ bmp = gcnew Bitmap(image->GetLength(0), image->GetLength(1), PixelFormat::Format24bppRgb);
	BitmapData^ data = bmp->LockBits(Drawing::Rectangle(0, 0, image->GetLength(0), image->GetLength(1)), ImageLockMode::WriteOnly, PixelFormat::Format24bppRgb);
	unsigned char *bits = (unsigned char *)data->Scan0.ToPointer();
	int bytesPerPixel = 3; // 3 bytes per pixel for 24 bpp rgb

	int height = data->Height;
	int width = data->Width;
	int stride = data->Stride;
	int bytesWidth = width * bytesPerPixel;
	double invDImCLimRange = 1 / (DImCLim[1] - DImCLim[0]);

	#pragma omp parallel for
	for (int i = 0; i < height; i++)
	{
		int istride = i * stride;
		int jcounter = -1;
		for (int j = 0; j < bytesWidth; j += bytesPerPixel)
		{
			jcounter++;
			double val = image[jcounter, i];
			if (val < DImCLim[0])
				val = DImCLim[0];
			else if (val > DImCLim[1])
				val = DImCLim[1];
			val = (val - DImCLim[0]) * invDImCLimRange;

			switch (scaling)
			{
				case (0):
				{
					val = val * 255;
					break;
				}
				case (1):
				{
					val = Math::Sqrt(val) * 255;
					break;
				}
				case (2):
				{
					val = val * val * 255;
					break;
				}
				case (3):
				{
					val = Math::Log(Math::Sqrt(val*val) + 1) * 255;//log
					break;
				}
			}

			if (invert)
				val = 255 - val;

			switch (colour)
			{
				case (0)://grayscale
				{
					bits[istride + j + 0] = int(val);   // blue
					bits[istride + j + 1] = int(val); // green
					bits[istride + j + 2] = int(val);   // red
					break;
				}
				case (1)://Jet
				{
					bits[istride + j + 0] = int(JetB(val));   // blue
					bits[istride + j + 1] = int(JetG(val)); // green
					bits[istride + j + 2] = int(JetR(val));   // red
					break;
				}
				case (2)://Winter
				{
					bits[istride + j + 0] = int(255 - 0.5*val);   // blue
					bits[istride + j + 1] = int(val); // green
					bits[istride + j + 2] = 0;   // red
					break;
				}
				case (3)://Lines
				{
					bits[istride + j + 0] = int(LinesB(val));   // blue
					bits[istride + j + 1] = int(LinesG(val)); // green
					bits[istride + j + 2] = int(LinesR(val));   // red
					break;
				}
			}
		}
	}

	bmp->UnlockBits(data);
	if (invertYaxis)
		bmp->RotateFlip(::Drawing::RotateFlipType::RotateNoneFlipY);

	return bmp;
}

Bitmap^ JPFITS::JPBitMap::RGBBitMap(cli::array<double, 2> ^R, cli::array<double, 2> ^G, cli::array<double, 2> ^B)
{
	bool codim = true;
	if (R->GetLength(0) != G->GetLength(0) || R->GetLength(0) != B->GetLength(0) || G->GetLength(0) != B->GetLength(0))
		codim = false;
	if (R->GetLength(1) != G->GetLength(1) || R->GetLength(1) != B->GetLength(1) || G->GetLength(1) != B->GetLength(1))
		codim = false;
	if (codim == false)
	{
		::MessageBox::Show("Error: RGB array set not co-dimensional...", "Error...");
		return nullptr;
	}

	Bitmap^ bmp = gcnew Bitmap(R->GetLength(0), R->GetLength(1), PixelFormat::Format24bppRgb);

	BitmapData^ data = bmp->LockBits(Drawing::Rectangle(0, 0, R->GetLength(0), R->GetLength(1)), ImageLockMode::WriteOnly, PixelFormat::Format24bppRgb);
	unsigned char *bits = (unsigned char *)data->Scan0.ToPointer();
	int bytesPerPixel = 3; // 3 bytes per pixel for 24 bpp rgb

	int height = data->Height;
	int width = data->Width;
	int stride = data->Stride;
	int bytesWidth = width * bytesPerPixel;

	#pragma omp parallel for
	for (int i = 0; i < height; i++)
	{
		int istride = i * stride;
		int jcounter = -1;
		for (int j = 0; j < bytesWidth; j += bytesPerPixel)
		{
			jcounter++;

			bits[istride + j + 0] = int(B[jcounter, i]);   // blue
			bits[istride + j + 1] = int(G[jcounter, i]); // green
			bits[istride + j + 2] = int(R[jcounter, i]);   // red
		}
	}
	bmp->UnlockBits(data);
	return bmp;
}

//Bitmap^ JPBitMap::ArrayToBmp(cli::array<int,2> ^image, int scaling, int colour, bool invert, cli::array<int,1> ^ImCLim)
//{
//	Bitmap^ bmp = gcnew Bitmap(image->GetLength(0),image->GetLength(1),PixelFormat::Format24bppRgb);
//	BitmapData^ data = bmp->LockBits(Drawing::Rectangle(0,0,image->GetLength(0),image->GetLength(1)),ImageLockMode::WriteOnly,PixelFormat::Format24bppRgb);
//	unsigned char *bits = (unsigned char *)data->Scan0.ToPointer();
//    int bytesPerPixel = 3; // 3 bytes per pixel for 24 bpp rgb
//
//	int height = data->Height;
//	int width = data->Width;
//	int stride = data->Stride;
//
//	#pragma omp parallel for
//	for (int i = 0; i < height; i++)
//	{
//		for (int j = 0; j < width * bytesPerPixel; j += bytesPerPixel)
//		{
//			int val = image[j/bytesPerPixel,i];
//			if (val<ImCLim[0]) val = ImCLim[0];
//			if (val>ImCLim[1]) val = ImCLim[1];
//			val = (val - ImCLim[0])/(ImCLim[1]-ImCLim[0]);
//
//			switch (scaling)
//			{
//			case (0):
//				{
//					val = val*255;
//					break;
//				}
//			case (1):
//				{
//					val = (int)Math::Sqrt(val)*255;
//					break;
//				}
//			case (2):
//				{
//					val = val*val*255;
//					break;
//				}
//			case (3):
//				{
//					val = (int)Math::Log(Math::Sqrt(val*val)+1)*255;//log
//					break;
//				}
//			}
//
//			if (invert)
//				val = 255-val;
//
//			switch (colour)
//			{
//			case (0)://grayscale
//				{
//					bits[i*stride + j + 0] = val;   // blue
//					bits[i*stride + j + 1] = val; // green
//					bits[i*stride + j + 2] = val;   // red
//					break;
//				}
//			case (1)://Jet
//				{
//					bits[i*stride + j + 0] = int(JetB(val));   // blue
//					bits[i*stride + j + 1] = int(JetG(val)); // green
//					bits[i*stride + j + 2] = int(JetR(val));   // red
//					break;
//				}
//			case (2)://Winter, see functional form Winter(B,G,R)
//				{
//					bits[i*stride + j + 0] = 255-val/2;   // blue
//					bits[i*stride + j + 1] = val; // green
//					bits[i*stride + j + 2] = 0;   // red
//					break;
//				}
//			case (3)://Lines
//				{
//					bits[i*stride + j + 0] = int(LinesB(val));   // blue
//					bits[i*stride + j + 1] = int(LinesG(val)); // green
//					bits[i*stride + j + 2] = int(LinesR(val));   // red
//					break;
//				}
//			}
//		}
//	}
//	bmp->UnlockBits(data);
//	return bmp;
//}

/*array<double, 2>^ JPFITS::JPBitMap::Sample(cli::array<double, 2> ^img, String^ Method, int WinWidth, int WinHeight)
{
	int Lx = img->GetLength(0);//length x
	int Ly = img->GetLength(1);//length y

	double qpx = Math::Max(double(Lx) / double(WinWidth), 1.0);
	double qpy = Math::Max(double(Ly) / double(WinHeight), 1.0);

	int sx_max = Math::Min(WinWidth, Lx);
	int sy_max = Math::Min(WinHeight, Ly);

	array<int>^ sx = gcnew array<int>(sx_max + 1);
	array<int>^ sy = gcnew array<int>(sy_max + 1);
	for (int i = 0; i <= sx_max; i++)
		sx[i] = (int)((double(i)*qpx) + 0.5);
	for (int j = 0; j <= sy_max; j++)
		sy[j] = (int)((double(j)*qpy) + 0.5);
	array<double, 2>^ s_img = gcnew array<double, 2>(sx_max, sy_max);

	#pragma omp parallel for
	for (int i = 0; i < sx_max; i++)
		for (int j = 0; j < sy_max; j++)
		{
			double samp = 0;
			double c = 0;
			for (int ii = sx[i]; ii < sx[i + 1]; ii++)
			{
				for (int jj = sy[j]; jj < sy[j + 1]; jj++)
				{
					c++;
					samp = samp + img[ii, jj];
				}
			}
			s_img[i, j] = samp / c;
		}

	return s_img;
}*/

/*array<double, 2>^ JPFITS::JPBitMap::Sample2(cli::array<double, 2> ^img, int WinWidth, int WinHeight)
{
	array<double, 2>^ sampled_img;

	int ImgPixX = img->GetLength(0);
	int ImgPixY = img->GetLength(1);

	double scaleX = Math::Max(double(ImgPixX) / double(WinWidth), 1.0);//if the scale is less than 1, then the pixels dont need to be binned - i.e. small image size
	double scaleY = Math::Max(double(ImgPixY) / double(WinWidth), 1.0);

	//if the scale is greater than 1, then interpolated binning required

	//.......................... ... hmmmmmm  :)







	return sampled_img;
}*/

double JPFITS::JPBitMap::WinterR(double val)
{
	//this channel dead
	return 0;
}

double JPFITS::JPBitMap::WinterG(double val)
{
	//this channel linear
	return val;
}

double JPFITS::JPBitMap::WinterB(double val)
{
	//this channel half linear backwards
	return 255 - 0.5*val;
}

inline double JPFITS::JPBitMap::JetR(double val)
{
	if (val < 96)
		return 0;
	if (val >= 96 && val < 160)
		return 4 * (val - 96) - 1;
	if (val >= 160 && val < 224)
		return 255;
	return 255 - (val - 224) * 4;
}

inline double JPFITS::JPBitMap::JetG(double val)
{
	if (val < 32)
		return 0;
	if (val >= 32 && val < 96)
		return 4 * (val - 32) - 1;
	if (val >= 96 && val < 160)
		return 255;
	if (val >= 160 && val < 224)
		return 255 - (val - 160) * 4;
	return 0;
}

inline double JPFITS::JPBitMap::JetB(double val)
{
	if (val < 28)
		return 144 + (val) * 4 - 1;
	if (val >= 28 && val < 96)
		return 255;
	if (val >= 96 && val < 160)
		return 255 - (val - 160) * 4;
	return 0;
}

double JPFITS::JPBitMap::LinesR(double val)
{
	double res = 0;
	int mod = (int(val - 2)) % 7;
	switch (mod)
	{
		case (0):
		{
			res = 255;
			break;
		}
		case(4):
		{
			res = 63;
			break;
		}
		case(2):
		case(3):
		{
			res = 191;
			break;
		}
	}
	return res;
}

double JPFITS::JPBitMap::LinesG(double val)
{
	double res = 0;
	int mod = (int(val - 1)) % 7;
	switch (mod)
	{
		case (0):
		{
			res = 127;
			break;
		}
		case(2):
		case(4):
		{
			res = 191;
			break;
		}
		case(5):
		{
			res = 63;
			break;
		}
	}
	return res;
}

double JPFITS::JPBitMap::LinesB(double val)
{
	double res = 0;
	int mod = (int(val)) % 7;
	switch (mod)
	{
		case (0):
		{
			res = 255;
			break;
		}
		case(3):
		case(4):
		{
			res = 191;
			break;
		}
		case(6):
		{
			res = 63;
			break;
		}
	}
	return res;
}

