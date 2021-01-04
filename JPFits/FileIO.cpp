/*Copyright 2017 Joseph Edwin Postma

joepostma@live.ca

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


void JPFITS::FITSImage::READHEADER(FileStream^ fs, bool pop)
{
	//read and format new header keys/vals/comms and advance basestream to end of header block.  Should also set basic image properties
	//here like BZERO and BSCALE etc, whatever's important
	ArrayList^ lines = gcnew ArrayList();
	ArrayList^ keys = gcnew ArrayList();
	ArrayList^ vals = gcnew ArrayList();
	ArrayList^ coms = gcnew ArrayList();

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (headerlines == 1)
		{
			if (strheaderline->Substring(0, 8) != "SIMPLE  ")
				throw gcnew Exception("Error: File is not a FITS. SIMPLE keyword non-existent");
			if (strheaderline->Substring(10, 20)->Trim() != "T")
				throw gcnew Exception("Error: File does not conform to FITS. SIMPLE keyvalue not 'T'");

			fs->Close();
			return;
		}

		if (strheaderline->Substring(0,8) == "BZERO   ")
			BZERO = (__int64)::Convert::ToDouble(strheaderline->Substring(10,20));
		if (strheaderline->Substring(0,8) == "BSCALE  ")
			BSCALE = (__int64)::Convert::ToDouble(strheaderline->Substring(10,20));
		if (strheaderline->Substring(0,8) == "NAXIS   ")
			NAXIS = ::Convert::ToInt32(strheaderline->Substring(10,20));
		if (strheaderline->Substring(0,8) == "NAXIS1  ")
			NAXIS1 = ::Convert::ToInt32(strheaderline->Substring(10,20));
		if (strheaderline->Substring(0,8) == "NAXIS2  ")
			NAXIS2 = ::Convert::ToInt32(strheaderline->Substring(10,20));
		if (strheaderline->Substring(0,8) == "BITPIX  ")
			BITPIX = ::Convert::ToInt32(strheaderline->Substring(10,20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;

		if (pop)
		{
			lines->Add(strheaderline);
			keys->Add(strheaderline->Substring(0,8));
			if (strheaderline->Substring(0,7)->Equals("COMMENT"))
			{
				vals->Add(strheaderline->Substring(8,18));
				coms->Add(strheaderline->Substring(26));
			}
			else
			{
				vals->Add(strheaderline->Substring(10,20));
				coms->Add(strheaderline->Substring(32));
			}
		}
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	int NKeys = keys->Count;

	if (!pop)
		return;//exit this without formatting header, to save some time if only data array wanted

	HEADERLINES = gcnew array<String^>(NKeys);
	HEADERKEYS = gcnew array<String^>(NKeys);
	HEADERKEYVALS = gcnew array<String^>(NKeys);
	HEADERKEYCOMS = gcnew array<String^>(NKeys);
	String^ s;
	String^ nock = "'";

	for (int i = 0; i < NKeys; i++)
	{
		HEADERLINES[i] = (String^)lines[i];

		s = (String^)keys[i];
		s = s->Trim();
		HEADERKEYS[i] = s;

		s = (String^)vals[i];
		if (JPMath::IsNumeric(s))//this has to work if it is supposed to be a numeric value here
			s = s->Trim();//get rid of leading and trailing white space
		else//is a string, possibly with ' and '.  If it was a comment...
		{
			s = s->Trim();//get rid of white space at beginning and end which doesnt need to be there, and so ' is visible
			s = s->Trim(nock->ToCharArray());
			s = s->Trim();//might be more white space now, so get rid of it...why are people so sloppy...just make a standrd and stick to it!
		}
		if (s->Trim()->Equals("T"))//happens for SIMPLE and EXTEND keywords..maybe others
			s = "T";//same thing as s = s->Trim(); 		
		HEADERKEYVALS[i] = s;

		s = (String^)coms[i];
		HEADERKEYCOMS[i] = s->Trim();//trim starting and trailing blanks from key comment cause they're not supposed to be there
	}
}

void JPFITS::FITSImage::READDATA(FileStream^ fs, array<int,1>^ Range, bool do_parallel)//assumed READHEADER has been completed, and bs is at beginning of data block
{
	/*if (NAXIS == 0)
	{
		NAXIS1 = 1;
		NAXIS2 = 1;
		if (HEADER_POP)
		{
			SetKey("NAXIS1","1",false,0);
			SetKey("NAXIS2","1",false,0);
		}
		DIMAGE = gcnew array<double,2>(1,1);
		return;
	}*/

	double bscale = (double)BSCALE;
	double bzero = (double)BZERO;
	int bpix = ::Math::Abs(BITPIX);
	int NBytes = NAXIS1*NAXIS2*(bpix/8);
	array<int>^ R = gcnew array<int>(4);
	
	if (Range != nullptr)
	{
		if (Range[1] > NAXIS1-1 || Range[1] == -1)
			Range[1] = NAXIS1-1;
		if (Range[3] > NAXIS2-1 || Range[3] == -1)
			Range[3] = NAXIS2-1;
	}

	if (Range == nullptr || Range[0] == -1)//then it is a full frame read
	{
		R[0] = 0;
		R[1] = NAXIS1-1;
		R[2] = 0;
		R[3] = NAXIS2-1;
	}
	else//else it is a sub-frame read
	{
		R[0] = Range[0];
		R[1] = Range[1];
		R[2] = Range[2];
		R[3] = Range[3];
	}

	int W = R[1] - R[0] + 1;
	int H = R[3] - R[2] + 1;
	DIMAGE = gcnew array<double,2>(W,H);
	
	array<unsigned char>^ arr = gcnew array<unsigned char>(NBytes);
	fs->Read(arr,0,NBytes);//seems to be fastest to just read the entire data even if only subimage will be used

	if (BITPIX == 8)
	{
		int cc;

		#pragma omp parallel for if (do_parallel) private(cc)
		for (int j = R[2]; j <= R[3]; j++)
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*NAXIS1 + i);
				DIMAGE[i-R[0],j-R[2]] = double(arr[cc])*bscale + bzero;
			}
	}

	if (BITPIX == 16)
	{
		int cc;
		__int16 val;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for (int j = R[2]; j <= R[3]; j++)
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*NAXIS1 + i)*2;
				val = (arr[cc] << 8) | arr[cc + 1];
				DIMAGE[i-R[0],j-R[2]] = double(val)*bscale + bzero;
			}
	}

	if (BITPIX == 32)
	{
		int cc;
		__int32 val;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for (int j = R[2]; j <= R[3]; j++)
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*NAXIS1 + i)*4;
				val = (arr[cc] << 24) | (arr[cc + 1] << 16) | (arr[cc + 2] << 8) | arr[cc + 3];
				DIMAGE[i-R[0],j-R[2]] = double(val)*bscale + bzero;
			}
	}

	if (BITPIX == -32)
	{
		float val = 0;
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(val, cc)
		for (int j = R[2]; j <= R[3]; j++)
		{
			array<unsigned char>^ flt = gcnew array<unsigned char>(4);
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*NAXIS1 + i) * 4;
				flt[3] = arr[cc];
				flt[2] = arr[cc + 1];
				flt[1] = arr[cc + 2];
				flt[0] = arr[cc + 3];
				val = BitConverter::ToSingle(flt, 0);
				DIMAGE[i - R[0], j - R[2]] = double(val)*bscale + bzero;
			}
		}
	}

	if (BITPIX == 64)
	{
		__int64 val = 0;
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(val, cc)
		for (int j = R[2]; j <= R[3]; j++)
		{
			array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*NAXIS1 + i) * 8;
				dbl[7] = arr[cc];
				dbl[6] = arr[cc + 1];
				dbl[5] = arr[cc + 2];
				dbl[4] = arr[cc + 3];
				dbl[3] = arr[cc + 4];
				dbl[2] = arr[cc + 5];
				dbl[1] = arr[cc + 6];
				dbl[0] = arr[cc + 7];
				val = BitConverter::ToInt64(dbl, 0);
				DIMAGE[i - R[0], j - R[2]] = double(val)*bscale + bzero;
			}
		}
	}

	if (BITPIX == -64)
	{
		double val = 0;
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(val, cc)
		for (int j = R[2]; j <= R[3]; j++)
		{
			array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*NAXIS1 + i) * 8;
				dbl[7] = arr[cc];
				dbl[6] = arr[cc + 1];
				dbl[5] = arr[cc + 2];
				dbl[4] = arr[cc + 3];
				dbl[3] = arr[cc + 4];
				dbl[2] = arr[cc + 5];
				dbl[1] = arr[cc + 6];
				dbl[0] = arr[cc + 7];
				val = BitConverter::ToDouble(dbl, 0);
				DIMAGE[i - R[0], j - R[2]] = (val)*bscale + bzero;
			}
		}
	}

	/*
	//subimage reading only...seems much slower than just reading entire image
	if (BITPIX == -64)
	{
		array<unsigned char>^ arr = gcnew array<unsigned char>(W*8);
		double val;
		int cc;
		array<unsigned char>^ dbl = gcnew array<unsigned char>(8);

		bs->Seek((R[2]*NAXIS1 + R[0])*8,::IO::SeekOrigin::Current);

		for (int j = R[2]; j <= R[3]; j++)
		{
			bs->Read(arr,0,W*8);

			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (i - R[0])*8;
				dbl[7] = arr[cc]; 
				dbl[6] = arr[cc + 1]; 
				dbl[5] = arr[cc + 2]; 
				dbl[4] = arr[cc + 3]; 
				dbl[3] = arr[cc + 4]; 
				dbl[2] = arr[cc + 5]; 
				dbl[1] = arr[cc + 6]; 
				dbl[0] = arr[cc + 7];
				val = BitConverter::ToDouble(dbl, 0);
				DIMAGE[i-R[0],j-R[2]] = (val)*bscale + bzero;
			}

			bs->Seek((NAXIS1 - R[1] - 1 + R[0])*8,::IO::SeekOrigin::Current);
		}
	}*/

	NAXIS1 = W;
	NAXIS2 = H;
	if (HEADER_POP)
	{
		SetKey("NAXIS1", W.ToString(), false, 0);
		SetKey("NAXIS2", H.ToString(), false, 0);
	}
}

void JPFITS::FITSImage::WRITEFILE(TypeCode Precision, bool do_parallel)
{
	//first need to set header BZERO and BCSALE key values depending on prec type.  BSCALE always assumed to be one I guess
	//also need to set HEADERKEYVALS
	SETBITPIX(Precision);
	FORMATHEADER();//make header ready to write...even adds padding to 2880 blocks, so just go through

	int NHead = HEADER->Length*80;//number of bytes in fitsheader...should always be multiple of 2880.
	__int64 NIm = __int64(NAXIS1)*__int64(NAXIS2)*__int64(Math::Abs(BITPIX/8));//number of bytes in image
	__int64 NImNHead = __int64(NHead) + NIm;//this is the number of bytes in the file + header...but need to write file so multiple of 2880 bytes
	int NBlocks = int(Math::Ceiling(double(NImNHead)/2880.0));
	int NBytesTot = NBlocks*2880;

	double bscale = (double)BSCALE;
	double bzero = (double)BZERO;

	FileStream^ fs = gcnew FileStream(FULLFILENAME,IO::FileMode::Create);

	array<unsigned char>^ data = gcnew array<unsigned char>(NBytesTot);
	
	for (int i = 0; i < HEADER->Length; i++)
		for (int j = 0; j < 80; j++)
			data[i * 80 + j] = (unsigned char)HEADER[i][j];

	if (Precision == TypeCode::Char)
	{
		char val = 0;
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for(int j = 0; j < NAXIS2; j++)
			for(int i = 0; i < NAXIS1; i++)
			{
				cc = NHead + (j*NAXIS1 + i)*2;
				val = char((DIMAGE[i,j] - bzero)/bscale);
				data[cc] = val;
			}
	}

	if (Precision == TypeCode::UInt16 || Precision == TypeCode::Int16)//use signed Int16 in FITS standard...bzero is used to scale to UInt.
	{
		__int16 val = 0;
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for(int j = 0; j < NAXIS2; j++)
			for(int i = 0; i < NAXIS1; i++)
			{
				cc = NHead + (j*NAXIS1 + i)*2;
				val = __int16((DIMAGE[i,j] - bzero)/bscale);
				data[cc] = ((val>>8)&0xff);
				data[cc + 1] = (val&0xff);
			}
	}

	if (Precision == TypeCode::UInt32 || Precision == TypeCode::Int32)
	{
		__int32 val = 0;
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for(int j = 0; j < NAXIS2; j++)
			for(int i = 0; i < NAXIS1; i++)
			{
				cc = NHead + (j*NAXIS1 + i)*4;
				val = __int32((DIMAGE[i,j] - bzero)/bscale);
				data[cc] = ((val>>24)&0xff);
				data[cc + 1] = ((val>>16)&0xff);
				data[cc + 2] = ((val>>8)&0xff);
				data[cc + 3] = (val&0xff);
			}
	}

	if (Precision == TypeCode::Double)
	{
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(cc)
		for (int j = 0; j < NAXIS2; j++)
		{
			array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
			for (int i = 0; i < NAXIS1; i++)
			{
				cc = NHead + (j*NAXIS1 + i) * 8;
				dbl = BitConverter::GetBytes((DIMAGE[i, j] - bzero) / bscale);
				data[cc] = dbl[7];
				data[cc + 1] = dbl[6];
				data[cc + 2] = dbl[5];
				data[cc + 3] = dbl[4];
				data[cc + 4] = dbl[3];
				data[cc + 5] = dbl[2];
				data[cc + 6] = dbl[1];
				data[cc + 7] = dbl[0];
			}
		}
	}

	fs->Write(data, 0, NBytesTot);
	fs->Close();
}

array<double,2>^ JPFITS::FITSImage::ReadImageArrayOnly(System::String ^file, cli::array<int,1> ^Range, bool do_parallel)
{
	FileStream^ fs = gcnew FileStream(file,IO::FileMode::Open);

	double bzero = 0;
	double bscale = 1;
	int naxis1, naxis2, bitpix;

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (strheaderline->Substring(0, 8) == "BZERO   ")
			bzero = ::Convert::ToDouble(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "BSCALE  ")
			bscale = ::Convert::ToDouble(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "NAXIS1  ")
			naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "NAXIS2  ")
			naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	int bpix = ::Math::Abs(bitpix);
	int NBytes = naxis1*naxis2*(bpix/8);
	array<int>^ R = gcnew array<int>(4);

	if (Range == nullptr || Range[0] == -1)//then it is a full frame read
	{
		R[0] = 0;
		R[1] = naxis1-1;
		R[2] = 0;
		R[3] = naxis2-1;
	}
	else//else it is a sub-frame read
	{
		R[0] = Range[0];
		R[1] = Range[1];
		R[2] = Range[2];
		R[3] = Range[3];
	}

	int W = R[1] - R[0] + 1;
	int H = R[3] - R[2] + 1;
	array<double,2>^ result = gcnew array<double,2>(W,H);
	array<unsigned char>^ arr = gcnew array<unsigned char>(NBytes);
	fs->Read(arr,0,NBytes);
	fs->Close();

	if (bitpix == 8)
	{
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(cc)
		for (int j = R[2]; j <= R[3]; j++)
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*naxis1 + i);
				result[i-R[0],j-R[2]] = double(arr[cc])*bscale + bzero;
			}
	}

	if (bitpix == 16)
	{
		int cc = 0;
		__int16 val = 0;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for (int j = R[2]; j <= R[3]; j++)
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*naxis1 + i)*2;
				val = (((arr[cc])<<8) | arr[cc + 1]);
				result[i-R[0],j-R[2]] = double(val)*bscale + bzero;
			}
	}

	if (bitpix == 32)
	{
		int cc = 0;
		__int32 val = 0;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for (int j = R[2]; j <= R[3]; j++)
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*naxis1 + i)*4;
				val = ((arr[cc])<<24) | ((arr[cc + 1])<<16) | (arr[cc + 2]<<8) | arr[cc + 3];
				result[i-R[0],j-R[2]] = double(val)*bscale + bzero;
			}
	}

	if (bitpix == -32)
	{
		float val = 0;
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(val, cc)
		for (int j = R[2]; j <= R[3]; j++)
		{
			array<unsigned char>^ flt = gcnew array<unsigned char>(4);
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*naxis1 + i) * 4;
				flt[3] = arr[cc];
				flt[2] = arr[cc + 1];
				flt[1] = arr[cc + 2];
				flt[0] = arr[cc + 3];
				val = BitConverter::ToSingle(flt, 0);
				result[i - R[0], j - R[2]] = double(val)*bscale + bzero;
			}
		}
	}

	if (bitpix == -64)
	{
		double val = 0;
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(val, cc)
		for (int j = R[2]; j <= R[3]; j++)
		{
			array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*naxis1 + i) * 8;
				dbl[7] = arr[cc];
				dbl[6] = arr[cc + 1];
				dbl[5] = arr[cc + 2];
				dbl[4] = arr[cc + 3];
				dbl[3] = arr[cc + 4];
				dbl[2] = arr[cc + 5];
				dbl[1] = arr[cc + 6];
				dbl[0] = arr[cc + 7];
				val = BitConverter::ToDouble(dbl, 0);
				result[i - R[0], j - R[2]] = (val)*bscale + bzero;
			}
		}
	}

	return result;
}

array<double>^ JPFITS::FITSImage::ReadImageVectorOnly(System::String ^file, cli::array<int,1> ^Range, bool do_parallel)
{
	if (Range != nullptr && (Range[1] - Range[0]) > 0 && (Range[3] - Range[2]) > 0 )
	{
		::MessageBox::Show("Requested Vector Output but Specified Range is 2D","Error");
		return nullptr;
	}

	FileStream^ fs = gcnew FileStream(file,IO::FileMode::Open);
	array<unsigned char>^ c = gcnew array<unsigned char>(2880);

	double bzero = 0;
	double bscale = 1;
	int naxis1, naxis2, bitpix;

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (strheaderline->Substring(0, 8) == "BZERO   ")
			bzero = ::Convert::ToDouble(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "BSCALE  ")
			bscale = ::Convert::ToDouble(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "NAXIS1  ")
			naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "NAXIS2  ")
			naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	int bpix = ::Math::Abs(bitpix);
	int NBytes = naxis1*naxis2*(bpix/8);
	array<int>^ R = gcnew array<int>(4);

	if (Range == nullptr || Range[0] == -1)//then it is a full frame read
	{
		R[0] = 0;
		R[1] = naxis1-1;
		R[2] = 0;
		R[3] = naxis2-1;
	}
	else//else it is a sub-frame read
	{
		R[0] = Range[0];
		R[1] = Range[1];
		R[2] = Range[2];
		R[3] = Range[3];
	}

	if ((R[1] - R[0]) > 0 && (R[3] - R[2]) > 0)
	{
		::MessageBox::Show("Requested default Vector Output but Image is 2D", "Error");
		return nullptr;
	}

	int W = R[1] - R[0] + 1;
	int H = R[3] - R[2] + 1;
	array<double>^ result = gcnew array<double>(W*H);
	array<unsigned char>^ arr = gcnew array<unsigned char>(NBytes);
	fs->Read(arr,0,NBytes);
	fs->Close();

	if (bitpix == 8)
	{
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(cc)
		for (int j = R[2]; j <= R[3]; j++)
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*naxis1 + i);
				result[(i-R[0])*(j-R[2]) + j-R[2]] = double(arr[cc])*bscale + bzero;
			}
	}

	if (bitpix == 16)
	{
		int cc = 0;
		__int16 val = 0;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for (int j = R[2]; j <= R[3]; j++)
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*naxis1 + i)*2;
				val = (((arr[cc])<<8) | arr[cc + 1]);
				result[(i-R[0])*(j-R[2]) + j-R[2]] = double(val)*bscale + bzero;
			}
	}

	if (bitpix == 32)
	{
		int cc = 0;
		__int32 val = 0;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for (int j = R[2]; j <= R[3]; j++)
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*naxis1 + i)*4;
				val = ((arr[cc])<<24) | ((arr[cc + 1])<<16) | (arr[cc + 2]<<8) | arr[cc + 3];
				result[(i-R[0])*(j-R[2]) + j-R[2]] = double(val)*bscale + bzero;
			}
	}

	if (bitpix == -32)
	{
		float val = 0;
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for (int j = R[2]; j <= R[3]; j++)
		{
			array<unsigned char>^ flt = gcnew array<unsigned char>(4);
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*naxis1 + i) * 4;
				flt[3] = arr[cc];
				flt[2] = arr[cc + 1];
				flt[1] = arr[cc + 2];
				flt[0] = arr[cc + 3];
				val = BitConverter::ToSingle(flt, 0);
				result[(i - R[0])*(j - R[2]) + j - R[2]] = double(val)*bscale + bzero;
			}
		}
	}

	if (bitpix == -64)
	{
		double val = 0;
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for (int j = R[2]; j <= R[3]; j++)
		{
			array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*naxis1 + i) * 8;
				dbl[7] = arr[cc];
				dbl[6] = arr[cc + 1];
				dbl[5] = arr[cc + 2];
				dbl[4] = arr[cc + 3];
				dbl[3] = arr[cc + 4];
				dbl[2] = arr[cc + 5];
				dbl[1] = arr[cc + 6];
				dbl[0] = arr[cc + 7];
				val = BitConverter::ToDouble(dbl, 0);
				result[(i - R[0])*(j - R[2]) + j - R[2]] = (val)*bscale + bzero;
			}
		}
	}

	return result;
}

