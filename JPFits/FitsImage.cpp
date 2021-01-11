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

JPFITS::FITSImage::~FITSImage()
{
	delete DIMAGE, HEADERKEYS, HEADERKEYVALS, HEADERKEYCOMS;
}

JPFITS::FITSImage::FITSImage(String^ FullFileName, String^ DiskUCharBufferName, TypeCode Precision, int NAxis1, int NAxis2)
{
	ISEXTENSION = false;
	EXTNAME = nullptr;

	FULLFILENAME = FullFileName;
	int index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	DISKBUFFERFULLNAME = DiskUCharBufferName;
	FROMDISKBUFFER = true;

	DIMAGE = gcnew array<double,2>(1,1);//default fill

	HEADERKEYS = gcnew array<String^>(8);
	HEADERKEYVALS = gcnew array<String^>(8);
	HEADERKEYCOMS = gcnew array<String^>(8);
	NAXIS = 2;
	NAXIS1 = NAxis1; 
	NAXIS2 = NAxis2;
	HEADERKEYS[0] = "SIMPLE";
	HEADERKEYS[1] = "BITPIX";
	HEADERKEYS[2] = "NAXIS";
	HEADERKEYS[3] = "NAXIS1";
	HEADERKEYS[4] = "NAXIS2";
	HEADERKEYS[5] = "BZERO";
	HEADERKEYS[6] = "BSCALE";
	HEADERKEYS[7] = "END";
	HEADERKEYVALS[0] = "T";
	HEADERKEYVALS[1] = "-64";
	HEADERKEYVALS[2] = "2";
	HEADERKEYVALS[3] = NAXIS1.ToString();
	HEADERKEYVALS[4] = NAXIS2.ToString();
	HEADERKEYVALS[5] = "0";
	HEADERKEYVALS[6] = "1";
	HEADERKEYVALS[7] = "";
	HEADERKEYCOMS[0] = "Valid Fits File";
	HEADERKEYCOMS[1] = "Bits per Pixel";
	HEADERKEYCOMS[2] = "Number of Axes";
	HEADERKEYCOMS[3] = "Width (No. of Columns)";
	HEADERKEYCOMS[4] = "Height (No. of Rows)";
	HEADERKEYCOMS[5] = "Data Offset";
	HEADERKEYCOMS[6] = "Data Scaling";
	HEADERKEYCOMS[7] = "";
	SETBITPIX(Precision);
}

JPFITS::FITSImage::FITSImage(System::String ^FullFileName, array<int,1>^ Range, bool Populate_Header, bool Populate_Data, bool Do_Stats, bool do_parallel)
{
	HEADER_POP = Populate_Header;
	DATA_POP = Populate_Data;
	STATS_POP = Do_Stats;
	ISEXTENSION = false;

	if (DATA_POP == false)
	{
		Do_Stats = false;//can't do stats on data that isn't read
		STATS_POP = false;
	}

	ISEXTENSION = false;
	EXTNAME = nullptr;

	FULLFILENAME = FullFileName;
	int index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	
	FileStream^ fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Open);
	ArrayList^ header = gcnew ArrayList();
	bool hasext;
	if (!FITSFILEOPS::SCANPRIMARYUNIT(fs, false, header, hasext))
	{
		fs->Close();
		throw gcnew Exception("File not formatted as FITS file.");
		return;
	}
	EATRAWIMAGEHEADER(header, HEADER_POP);

	if (DATA_POP == true)
		READIMAGE(fs, Range, do_parallel);

	fs->Close();

	if (STATS_POP)
		StatsUpD(do_parallel);
}

JPFITS::FITSImage::FITSImage(String^ FullFileName, String^ extensionName, array<int, 1>^ Range, bool Populate_Header, bool Populate_Data, bool Do_Stats, bool do_parallel)
{
	FileStream^ fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Open);
	bool hasext;
	if (!FITSFILEOPS::SCANPRIMARYUNIT(fs, true, nullptr, hasext) || !hasext)
	{
		throw gcnew Exception("File not formatted as FITS file, or indicates no extensions present.");
	}

	ArrayList^ header = gcnew ArrayList();
	__int64 extensionstartposition, extensionendposition;
	if (!FITSFILEOPS::SEEKEXTENSION(fs, "IMAGE", extensionName, header, extensionstartposition, extensionendposition))
	{
		fs->Close();
		throw gcnew Exception("Could not find IMAGE extension with name '" + extensionName + "'");
		return;
	}
	EATRAWIMAGEHEADER(header, Populate_Header);

	HEADER_POP = Populate_Header;
	DATA_POP = Populate_Data;
	STATS_POP = Do_Stats;
	ISEXTENSION = true;
	EXTNAME = extensionName;

	if (DATA_POP == false)
	{
		Do_Stats = false;//can't do stats on data that isn't read
		STATS_POP = false;
	}

	FULLFILENAME = FullFileName;
	int index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index + 1);
	FILEPATH = FULLFILENAME->Substring(0, index + 1);

	if (DATA_POP == true)
		READIMAGE(fs, Range, do_parallel);

	fs->Close();

	if (STATS_POP)
		StatsUpD(do_parallel);
}

JPFITS::FITSImage::FITSImage(String^ FullFileName, Object^ ImageData, bool Do_Stats, bool do_parallel)
{
	if (!(ImageData->GetType()->IsArray))
	{
		throw gcnew Exception("Error: Object 'ImageData' is not an array.");
		return;
	}
	
	int rank = ((Array^)ImageData)->Rank;

	if (rank > 2)
	{
		throw gcnew Exception("Error: Rank of Object 'ImageData' is not a 1D or 2D array. Use other implementations for creating larger dimensional FITS primary data.");
		return;
	}

	NAXIS = 2;
	if (rank == 1)
	{
		NAXIS1 = 1;
		NAXIS2 = ((Array^)ImageData)->Length;
	}
	else
	{
		NAXIS1 = ((Array^)ImageData)->GetLength(0);
		NAXIS2 = ((Array^)ImageData)->GetLength(1);
	}

	DIMAGE = gcnew array<double, 2>(NAXIS1, NAXIS2);

	TypeCode type = Type::GetTypeCode((((Array^)ImageData)->GetType())->GetElementType());

	switch (type)
	{
		case TypeCode::Double:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS2; x++)
					DIMAGE[0, x] = ((array<double>^)ImageData)[x];
			}
			else
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS1; x++)
					for (int y = 0; y < NAXIS2; y++)
						DIMAGE[x, y] = ((array<double, 2>^)ImageData)[x, y];
			}

			break;
		}

		case TypeCode::UInt16:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS2; x++)
					DIMAGE[0, x] = (double)((array<unsigned __int16>^)ImageData)[x];
			}
			else
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS1; x++)
					for (int y = 0; y < NAXIS2; y++)
						DIMAGE[x, y] = (double)((array<unsigned __int16, 2>^)ImageData)[x, y];
			}

			break;
		}

		case TypeCode::Int16:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS2; x++)
					DIMAGE[0, x] = (double)((array<__int16>^)ImageData)[x];
			}
			else
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS1; x++)
					for (int y = 0; y < NAXIS2; y++)
						DIMAGE[x, y] = (double)((array<__int16, 2>^)ImageData)[x, y];
			}

			break;
		}

		case TypeCode::UInt32:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS2; x++)
					DIMAGE[0, x] = (double)((array<unsigned __int32>^)ImageData)[x];
			}
			else
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS1; x++)
					for (int y = 0; y < NAXIS2; y++)
						DIMAGE[x, y] = (double)((array<unsigned __int32, 2>^)ImageData)[x, y];
			}

			break;
		}

		case TypeCode::Int32:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS2; x++)
					DIMAGE[0, x] = (double)((array<__int32>^)ImageData)[x];
			}
			else
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS1; x++)
					for (int y = 0; y < NAXIS2; y++)
						DIMAGE[x, y] = (double)((array<__int32, 2>^)ImageData)[x, y];
			}

			break;
		}

		default:
		{
			throw gcnew Exception("Data of TypeCode:" + type.ToString() + " not suported.");
			return;
		}

	}

	ISEXTENSION = false;
	EXTNAME = nullptr;

	FULLFILENAME = FullFileName;
	int index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index + 1);
	FILEPATH = FULLFILENAME->Substring(0, index + 1);

	MIN = 0.0;
	MAX = 0.0;
	MEDIAN = 0.0;
	MEAN = 0.0;
	STD = 0.0;
	SUM = 0.0;

	HEADER_POP = true;
	DATA_POP = true;
	STATS_POP = Do_Stats;

	//create default header
	MAKEDEFHEADER();

	if (STATS_POP)
		StatsUpD(do_parallel);
}

void JPFITS::FITSImage::SetImage(array<double,2>^ ImageData, bool Do_Stats, bool do_parallel)//sometimes want to not bother with stats, so cant use property setter for images
{
	DIMAGE = ImageData;
	NAXIS1 = ImageData->GetLength(0);
	NAXIS2 = ImageData->GetLength(1);
	SetKey("NAXIS1",NAXIS1.ToString(), false,-1);
	SetKey("NAXIS2",NAXIS2.ToString(), false,-1);
	STATS_POP = Do_Stats;
	if (STATS_POP)
		StatsUpD(do_parallel);
}

void JPFITS::FITSImage::WriteImage(System::TypeCode Precision, bool do_parallel)
{
	ISEXTENSION = false;
	EXTNAME = nullptr;

	WRITEIMAGE(Precision, do_parallel);
}

void JPFITS::FITSImage::WriteImage(System::String ^FullFileName, System::TypeCode Precision, bool do_parallel)
{
	ISEXTENSION = false;
	EXTNAME = nullptr;

	FULLFILENAME = FullFileName;
	int index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index + 1);
	FILEPATH = FULLFILENAME->Substring(0, index + 1);
	
	WRITEIMAGE(Precision, do_parallel);
}

void JPFITS::FITSImage::WriteImage(String^ FullFileName, String^ extensionName, bool overwriteIfExists, TypeCode Precision, bool do_parallel)
{
	ISEXTENSION = true;
	EXTNAME = extensionName;
	EXTNAMEOVERWRITE = overwriteIfExists;

	FULLFILENAME = FullFileName;
	int index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index + 1);
	FILEPATH = FULLFILENAME->Substring(0, index + 1);
	
	WRITEIMAGE(Precision, do_parallel);
}

void JPFITS::FITSImage::WriteFileFromDiskBuffer(bool DeleteOrigDiskBuffer)
{
	if (FROMDISKBUFFER != true)
	{
		::MessageBox::Show("This FITS not created from ""JPFITS::FITSImage(String^ FullFileName, String^ DiskUCharBufferName, TypeCode Precision, int NAxis1, int NAxis2)"" constructor.","Error");
		return;
	}

	array<String^>^ HEADER = FITSFILEOPS::GETFORMATTEDIMAGEHEADER(HEADERKEYS, HEADERKEYVALS, HEADERKEYCOMS, false);
	//FORMATHEADER();
	FileStream^ fits_fs = gcnew FileStream(FULLFILENAME,IO::FileMode::Create);
	
	array<unsigned char>^ head = gcnew array<unsigned char>(HEADER->Length * 80);
	for (int i = 0; i < HEADER->Length; i++)
		for (int j = 0; j < 80; j++)
			head[i * 80 + j] = (unsigned char)HEADER[i][j];

	fits_fs->Write(head, 0, head->Length);//header is written

	FileStream^ buff_fs = gcnew FileStream(DISKBUFFERFULLNAME,FileMode::Open,::FileAccess::Read);

	int NBytes = (int)buff_fs->Length;//size of disk data buffer
	int buffsize = 1024*1024*64;//size of memory array buffer
	array<unsigned char>^ buff_arr = gcnew array<unsigned char>(buffsize);

	int NBuffArrs = int(::Math::Ceiling(double(NBytes)/double(buffsize) - ::Double::Epsilon*100));
	for (int i = 0; i < NBuffArrs; i++)
	{
		int bytestoread = buffsize;
		if (i == NBuffArrs-1 && NBuffArrs != 1 || NBuffArrs == 1)
			bytestoread = NBytes - i*buffsize;
		buff_fs->Read(buff_arr, 0, bytestoread);
		fits_fs->Write(buff_arr, 0, bytestoread);
	}
	buff_fs->Close();
	int resid = int(Math::Ceiling(double(fits_fs->Position)/2880.0))*2880 - int(fits_fs->Position);
	array<unsigned char>^ resds = gcnew array<unsigned char>(resid);
	fits_fs->Write(resds, 0, resid);
	fits_fs->Close();

	if (DeleteOrigDiskBuffer)
		::File::Delete(DISKBUFFERFULLNAME);
}

JPFITS::FITSImage::FITSImage(String^ FullFileName, bool mayContainExtensions)
{
	ISEXTENSION = false;
	EXTNAME = nullptr;

	int index;
	FULLFILENAME = FullFileName;
	index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);

	HEADERKEYS = gcnew array<String^>(5);
	HEADERKEYVALS = gcnew array<String^>(5);
	HEADERKEYCOMS = gcnew array<String^>(5);

	HEADERKEYS[0] = "SIMPLE";
	HEADERKEYS[1] = "BITPIX";
	HEADERKEYS[2] = "NAXIS";
	HEADERKEYS[3] = "EXTEND";
	HEADERKEYS[4] = "END";

	HEADERKEYVALS[0] = "T";
	HEADERKEYVALS[1] = "-64";
	HEADERKEYVALS[2] = "0";
	if (mayContainExtensions)
		HEADERKEYVALS[3] = "T";
	else
		HEADERKEYVALS[3] = "F";
	HEADERKEYVALS[4] = "";

	HEADERKEYCOMS[0] = "Valid Fits File";
	HEADERKEYCOMS[1] = "Bits per Pixel";
	HEADERKEYCOMS[2] = "Number of Axes";
	if (mayContainExtensions)
		HEADERKEYCOMS[3] = "File may contain extensions";
	else
		HEADERKEYCOMS[3] = "File does not contain extensions";
	HEADERKEYCOMS[4] = "";
}

array<double,2>^ JPFITS::FITSImage::ConvertTxtToDblArray(String^ FullFileName, bool IsPoorlyFormatted)
{
	if (IsPoorlyFormatted)
	{
		StreamReader^ sr = gcnew StreamReader(FullFileName);
		String^ line;

		int j = 0;
		int cols = 0;
		line = sr->ReadLine();
		line = line->Trim();
		bool foundcol = false;
		String^ chstr;

		for (int i = 0; i < line->Length; i++)
		{
			chstr = line[i].ToString();
			if (chstr == "+" || chstr == "-" || chstr == "E" || chstr == "e" || chstr == ".")
				continue;

			if (JPMath::IsNumeric(chstr))
				foundcol = false;

			if (!JPMath::IsNumeric(chstr) && !foundcol)
			{
				cols++;
				foundcol = true;
				continue;
			}
		}
		cols++;//because there will be one more data point after the last delimit

		int rows = 1;//because one row was already read
		while (!sr->EndOfStream)
		{
			line = sr->ReadLine();
			rows++;
		}
		sr->Close();

		sr = gcnew StreamReader(FullFileName);
		array<double, 2>^ arr = gcnew array<double, 2>(cols, rows);

		for (int i = 0; i < rows; i++)
		{
			line = sr->ReadLine();
			line = line->Trim();
			int numstartind, col = 0;
			bool startindfnd = false;

			for (int j = 0; j < line->Length; j++)
			{
				chstr = line[j].ToString();

				if (chstr == "+" || chstr == "-" || chstr == "E" || chstr == "e" || chstr == ".")
					continue;

				if (JPMath::IsNumeric(chstr) && !startindfnd)
				{
					startindfnd = true;
					numstartind = j;
					continue;
				}

				if (!JPMath::IsNumeric(chstr) && startindfnd)
				{
					arr[col, i] = Convert::ToDouble(line->Substring(numstartind, j - numstartind));
					col++;
					startindfnd = false;
				}
			}
			arr[col, i] = Convert::ToDouble(line->Substring(numstartind));
		}
		sr->Close();
		return arr;
	}
	else
	{
		StreamReader^ sr = gcnew StreamReader(FullFileName);
		String^ line;

		int j = 0;
		int cols = 0;
		line = sr->ReadLine();
		while (line->IndexOf("	", j + 1) != -1)//"	" is the tab delimit; could also use comma "," if required etc
		{
			cols++;
			j = line->IndexOf("	", j + 1);
		}
		cols++; //because there will be one more data point after the last delimit

		int rows = 1;//because one row was already read
		while (!sr->EndOfStream)
		{
			line = sr->ReadLine();
			rows++;
		}
		sr->Close();

		sr = gcnew StreamReader(FullFileName);
		array<double, 2>^ arr = gcnew array<double, 2>(cols, rows);
		array<int>^ data_inds = gcnew array<int>(cols);
		double datum = 0;

		for (int i = 0; i < rows; i++)
		{
			line = sr->ReadLine();
			int j = 0;
			int c = 0;
			data_inds[0] = 0;
			while (line->IndexOf("	", j + 1) != -1)//"	" is the tab delimit
			{
				j = line->IndexOf("	", j + 1);
				c++;
				data_inds[c] = j;
			}
			for (int j = 0; j < cols; j++)
			{
				if (j < cols - 1)
					datum = ::Convert::ToDouble(line->Substring(data_inds[j], data_inds[j + 1] - data_inds[j]));
				else
					datum = ::Convert::ToDouble(line->Substring(data_inds[j]));

				arr[j, i] = datum;
			}
		}
		sr->Close();
		return arr;
	}
}

array<double,2>^ JPFITS::FITSImage::GetSubImage(int X_Center, int Y_Center, int X_HalfWidth, int Y_HalfWidth)
{
	array<double,2>^ result = gcnew array<double,2>(X_HalfWidth*2+1,Y_HalfWidth*2+1);

	/*if (X_Center - X_HalfWidth < 0 || Y_Center - Y_HalfWidth < 0 || X_Center - X_HalfWidth + result->GetLength(0) > NAXIS1 || Y_Center - Y_HalfWidth + result->GetLength(1) > NAXIS2)
	{

	}*/

	for (int x = 0; x < result->GetLength(0); x++)
		for (int y = 0; y < result->GetLength(1); y++)
			result[x,y] = DIMAGE[X_Center - X_HalfWidth + x, Y_Center - Y_HalfWidth + y];

	return result;
}

array<double, 2>^ JPFITS::FITSImage::GetSubImage(int X_Center, int Y_Center, int X_HalfWidth, int Y_HalfWidth, array<int>^& xdata, array<int>^& ydata)
{
	array<double, 2>^ result = gcnew array<double, 2>(X_HalfWidth * 2 + 1, Y_HalfWidth * 2 + 1);

	/*if (X_Center - X_HalfWidth < 0 || Y_Center - Y_HalfWidth < 0 || X_Center - X_HalfWidth + result->GetLength(0) > NAXIS1 || Y_Center - Y_HalfWidth + result->GetLength(1) > NAXIS2)
	{

	}*/

	for (int x = 0; x < result->GetLength(0); x++)
		for (int y = 0; y < result->GetLength(1); y++)
			result[x, y] = DIMAGE[X_Center - X_HalfWidth + x, Y_Center - Y_HalfWidth + y];

	for (int x = 0; x < result->GetLength(0); x++)
		xdata[x] = X_Center - X_HalfWidth + x;
	for (int y = 0; y < result->GetLength(1); y++)
		ydata[y] = Y_Center - Y_HalfWidth + y;

	return result;
}

array<double, 2>^ JPFITS::FITSImage::GetSubImage(array<int, 1>^ Range)
{
	array<double, 2>^ result = gcnew array<double, 2>(Range[1] - Range[0] + 1, Range[3] - Range[2] + 1);

	for (int x = 0; x < result->GetLength(0); x++)
		for (int y = 0; y < result->GetLength(1); y++)
			result[x, y] = DIMAGE[Range[0] + x, Range[2] + y];

	return result;
}

void JPFITS::FITSImage::EATRAWIMAGEHEADER(ArrayList^ header, bool populate_nonessential)
{
	//HEADERLINES = gcnew array<String^>(header->Count);
	HEADERKEYS = gcnew array<String^>(header->Count);
	HEADERKEYVALS = gcnew array<String^>(header->Count);
	HEADERKEYCOMS = gcnew array<String^>(header->Count);

	int Nnaxisn = 1;
	array<int>^ naxisn;

	for (int i = 0; i < header->Count; i++)
	{
		String^ line = (String^)header[i];

		if (BITPIX == -1)
			if (line->Substring(0, 8)->Trim() == "BITPIX")
				BITPIX = ::Convert::ToInt32(line->Substring(10, 20));

		if (NAXIS == -1)
			if (line->Substring(0, 8)->Trim() == "NAXIS")
			{
				NAXIS = ::Convert::ToInt32(line->Substring(10, 20));
				naxisn = gcnew array<int>(NAXIS);
			}

		if (Nnaxisn <= NAXIS)
			if (line->Substring(0, 8)->Trim() == ("NAXIS" + Nnaxisn.ToString()))
			{
				naxisn[Nnaxisn - 1] = ::Convert::ToInt32(line->Substring(10, 20));
				Nnaxisn++;
			}

		if (BZERO == -1)
			if (line->Substring(0, 8)->Trim() == "BZERO")
				BZERO = (__int64)::Convert::ToDouble(line->Substring(10, 20));

		if (BSCALE == -1)
			if (line->Substring(0, 8)->Trim() == "BSCALE")
				BSCALE = (__int64)::Convert::ToDouble(line->Substring(10, 20));

		if (populate_nonessential)
		{
			//HEADERLINES[i] = line;

			HEADERKEYS[i] = line->Substring(0, 8)->Trim();

			if (HEADERKEYS[i] == "COMMENT")
			{
				HEADERKEYVALS[i] = line->Substring(8, 18);
				HEADERKEYCOMS[i] = line->Substring(26);
			}
			else
			{
				if (JPMath::IsNumeric(line->Substring(10, 20)))//this has to work if it is supposed to be a numeric value here
					HEADERKEYVALS[i] = line->Substring(10, 20)->Trim();//get rid of leading and trailing white space
				else
				{
					String^ nock = "'";
					HEADERKEYVALS[i] = line->Substring(10, 20)->Trim();
					HEADERKEYVALS[i] = HEADERKEYVALS[i]->Trim(nock->ToCharArray());
					HEADERKEYVALS[i] = HEADERKEYVALS[i]->Trim();
				}

				HEADERKEYCOMS[i] = line->Substring(32)->Trim();
			}
		}
	}

	if (NAXIS == 0)
	{
		NAXIS1 = 0;
		NAXIS2 = 0;
	}
	else if (NAXIS == 1)
	{
		NAXIS1 = 1;
		NAXIS2 = naxisn[0];
	}
	else if (NAXIS > 1)
	{
		NAXIS1 = naxisn[0];
		NAXIS2 = naxisn[1];
	}

	if (BZERO == -1)
		BZERO = 0;
	if (BSCALE == -1)
		BSCALE = 1;
}

void JPFITS::FITSImage::READIMAGE(FileStream^ fs, array<int, 1>^ Range, bool do_parallel)//assumed READHEADER has been completed, and bs is at beginning of data block
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
	int NBytes = NAXIS1 * NAXIS2*(bpix / 8);
	array<int>^ R = gcnew array<int>(4);

	if (Range != nullptr)
	{
		if (Range[1] > NAXIS1 - 1 || Range[1] == -1)
			Range[1] = NAXIS1 - 1;
		if (Range[3] > NAXIS2 - 1 || Range[3] == -1)
			Range[3] = NAXIS2 - 1;
	}

	if (Range == nullptr || Range[0] == -1)//then it is a full frame read
	{
		R[0] = 0;
		R[1] = NAXIS1 - 1;
		R[2] = 0;
		R[3] = NAXIS2 - 1;
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
	DIMAGE = gcnew array<double, 2>(W, H);

	array<unsigned char>^ arr = gcnew array<unsigned char>(NBytes);
	fs->Read(arr, 0, NBytes);//seems to be fastest to just read the entire data even if only subimage will be used

	if (BITPIX == 8)
	{
		int cc;

		#pragma omp parallel for if (do_parallel) private(cc)
		for (int j = R[2]; j <= R[3]; j++)
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*NAXIS1 + i);
				DIMAGE[i - R[0], j - R[2]] = double(arr[cc])*bscale + bzero;
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
				cc = (j*NAXIS1 + i) * 2;
				val = (arr[cc] << 8) | arr[cc + 1];
				DIMAGE[i - R[0], j - R[2]] = double(val)*bscale + bzero;
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
				cc = (j*NAXIS1 + i) * 4;
				val = (arr[cc] << 24) | (arr[cc + 1] << 16) | (arr[cc + 2] << 8) | arr[cc + 3];
				DIMAGE[i - R[0], j - R[2]] = double(val)*bscale + bzero;
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

void JPFITS::FITSImage::WRITEIMAGE(TypeCode Precision, bool do_parallel)
{
	//FileStream^ fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Create);
	
	
	
	
	
	
	
	
	FileStream^ fs;
	bool filexists = File::Exists(FULLFILENAME);
	array<unsigned char>^ prependdata;
	array<unsigned char>^ appenddata;
	if (!filexists && !ISEXTENSION)//then just write to a new file
		fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Create);
	if (filexists && !ISEXTENSION)//then write the primary unit, and append any extensions if they already exist on the existing file
	{
		fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Open);
		//check for extensions
		bool hasext = false;
		FITSFILEOPS::SCANPRIMARYUNIT(fs, true, nullptr, hasext);
		if (hasext)
		{
			appenddata = gcnew array<unsigned char>(int(fs->Length - fs->Position));
			fs->Read(appenddata, 0, appenddata->Length);
		}
		fs->Close();
		fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Create);//write the primary unit, don't forget to append the extensions after if it isn't null
	}
	if (filexists && ISEXTENSION)//then get the primary data unit and also check for other extensions, etc
	{
		fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Open);
		bool hasext = false;
		if (!FITSFILEOPS::SCANPRIMARYUNIT(fs, true, nullptr, hasext))
		{
			fs->Close();
			throw gcnew Exception("Primary data unit of file is not a FITS file.");
			return;
		}

		__int64 extstartpos, extendpos;
		bool extexists = FITSFILEOPS::SEEKEXTENSION(fs, "IMAGE", EXTNAME, nullptr, extstartpos, extendpos);
		if (extexists && !EXTNAMEOVERWRITE)
		{
			fs->Close();
			throw gcnew Exception("IMAGE extension '" + EXTNAME + "' exists but told to not overwrite it.");
			return;
		}
		if (!extexists)//then the fs stream will be at the end of the file, get all prependdata
		{
			prependdata = gcnew array<unsigned char>((int)fs->Position);
			fs->Position = 0;
			fs->Read(prependdata, 0, prependdata->Length);
			fs->Close();
			fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Create);
			fs->Write(prependdata, 0, prependdata->Length);
		}
		else//then get the prepend units and any append data
		{
			prependdata = gcnew array<unsigned char>((int)extstartpos);
			fs->Position = 0;
			fs->Read(prependdata, 0, prependdata->Length);
			appenddata = gcnew array<unsigned char>(int(fs->Length - extendpos));
			fs->Position = extendpos;
			fs->Read(appenddata, 0, appenddata->Length);
			fs->Close();
			fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Create);
			fs->Write(prependdata, 0, prependdata->Length);
		}
	}
	if (!filexists && ISEXTENSION)//then write the extension to a new file with an empty primary unit
	{
		array<String^>^ pheader = FITSFILEOPS::MAKEPRIMARYDEFFORMATTEDHEADER(true);
		int NHead = pheader->Length * 80;
		prependdata = gcnew array<unsigned char>(NHead);
		for (int i = 0; i < pheader->Length; i++)
			for (int j = 0; j < 80; j++)
				prependdata[i * 80 + j] = (unsigned char)pheader[i][j];
		fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Create);
		fs->Write(prependdata, 0, prependdata->Length);
	}









	//set header BZERO and BCSALE key values depending on prec type.
	SETBITPIX(Precision);

	//get formatted header block
	array<String^>^ HEADER = FITSFILEOPS::GETFORMATTEDIMAGEHEADER(HEADERKEYS, HEADERKEYVALS, HEADERKEYCOMS, ISEXTENSION);
	int NHead = HEADER->Length * 80;//number of bytes in fitsheader...should always be multiple of 2880.
	__int64 NIm = __int64(NAXIS1)*__int64(NAXIS2)*__int64(Math::Abs(BITPIX / 8));//number of bytes in image
	__int64 NImNHead = __int64(NHead) + NIm;//this is the number of bytes in the file + header...but need to write file so multiple of 2880 bytes
	int NBlocks = int(Math::Ceiling(double(NImNHead) / 2880.0));
	int NBytesTot = NBlocks * 2880;

	double bscale = (double)BSCALE;
	double bzero = (double)BZERO;

	array<unsigned char>^ data = gcnew array<unsigned char>(NBytesTot);

	for (int i = 0; i < HEADER->Length; i++)
		for (int j = 0; j < 80; j++)
			data[i * 80 + j] = (unsigned char)HEADER[i][j];

	if (Precision == TypeCode::Char)
	{
		char val = 0;
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for (int j = 0; j < NAXIS2; j++)
			for (int i = 0; i < NAXIS1; i++)
			{
				cc = NHead + (j*NAXIS1 + i) * 2;
				val = char((DIMAGE[i, j] - bzero) / bscale);
				data[cc] = val;
			}
	}

	if (Precision == TypeCode::UInt16 || Precision == TypeCode::Int16)//use signed Int16 in FITS standard...bzero is used to scale to UInt.
	{
		__int16 val = 0;
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for (int j = 0; j < NAXIS2; j++)
			for (int i = 0; i < NAXIS1; i++)
			{
				cc = NHead + (j*NAXIS1 + i) * 2;
				val = __int16((DIMAGE[i, j] - bzero) / bscale);
				data[cc] = ((val >> 8) & 0xff);
				data[cc + 1] = (val & 0xff);
			}
	}

	if (Precision == TypeCode::UInt32 || Precision == TypeCode::Int32)
	{
		__int32 val = 0;
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(cc, val)
		for (int j = 0; j < NAXIS2; j++)
			for (int i = 0; i < NAXIS1; i++)
			{
				cc = NHead + (j*NAXIS1 + i) * 4;
				val = __int32((DIMAGE[i, j] - bzero) / bscale);
				data[cc] = ((val >> 24) & 0xff);
				data[cc + 1] = ((val >> 16) & 0xff);
				data[cc + 2] = ((val >> 8) & 0xff);
				data[cc + 3] = (val & 0xff);
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
	if (appenddata != nullptr)
		fs->Write(appenddata, 0, appenddata->Length);
	fs->Close();
}

array<double, 2>^ JPFITS::FITSImage::ReadImageArrayOnly(System::String ^file, cli::array<int, 1> ^Range, bool do_parallel)
{
	FileStream^ fs = gcnew FileStream(file, IO::FileMode::Open);

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
	int NBytes = naxis1 * naxis2*(bpix / 8);
	array<int>^ R = gcnew array<int>(4);

	if (Range == nullptr || Range[0] == -1)//then it is a full frame read
	{
		R[0] = 0;
		R[1] = naxis1 - 1;
		R[2] = 0;
		R[3] = naxis2 - 1;
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
	array<double, 2>^ result = gcnew array<double, 2>(W, H);
	array<unsigned char>^ arr = gcnew array<unsigned char>(NBytes);
	fs->Read(arr, 0, NBytes);
	fs->Close();

	if (bitpix == 8)
	{
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(cc)
		for (int j = R[2]; j <= R[3]; j++)
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*naxis1 + i);
				result[i - R[0], j - R[2]] = double(arr[cc])*bscale + bzero;
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
				cc = (j*naxis1 + i) * 2;
				val = (((arr[cc]) << 8) | arr[cc + 1]);
				result[i - R[0], j - R[2]] = double(val)*bscale + bzero;
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
				cc = (j*naxis1 + i) * 4;
				val = ((arr[cc]) << 24) | ((arr[cc + 1]) << 16) | (arr[cc + 2] << 8) | arr[cc + 3];
				result[i - R[0], j - R[2]] = double(val)*bscale + bzero;
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

array<double>^ JPFITS::FITSImage::ReadImageVectorOnly(System::String ^file, cli::array<int, 1> ^Range, bool do_parallel)
{
	if (Range != nullptr && (Range[1] - Range[0]) > 0 && (Range[3] - Range[2]) > 0)
	{
		::MessageBox::Show("Requested Vector Output but Specified Range is 2D", "Error");
		return nullptr;
	}

	FileStream^ fs = gcnew FileStream(file, IO::FileMode::Open);
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
	int NBytes = naxis1 * naxis2*(bpix / 8);
	array<int>^ R = gcnew array<int>(4);

	if (Range == nullptr || Range[0] == -1)//then it is a full frame read
	{
		R[0] = 0;
		R[1] = naxis1 - 1;
		R[2] = 0;
		R[3] = naxis2 - 1;
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
		::MessageBox::Show("Requested defualt Vector Output but Image is 2D", "Error");
		return nullptr;
	}

	int W = R[1] - R[0] + 1;
	int H = R[3] - R[2] + 1;
	array<double>^ result = gcnew array<double>(W*H);
	array<unsigned char>^ arr = gcnew array<unsigned char>(NBytes);
	fs->Read(arr, 0, NBytes);
	fs->Close();

	if (bitpix == 8)
	{
		int cc = 0;

		#pragma omp parallel for if (do_parallel) private(cc)
		for (int j = R[2]; j <= R[3]; j++)
			for (int i = R[0]; i <= R[1]; i++)
			{
				cc = (j*naxis1 + i);
				result[(i - R[0])*(j - R[2]) + j - R[2]] = double(arr[cc])*bscale + bzero;
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
				cc = (j*naxis1 + i) * 2;
				val = (((arr[cc]) << 8) | arr[cc + 1]);
				result[(i - R[0])*(j - R[2]) + j - R[2]] = double(val)*bscale + bzero;
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
				cc = (j*naxis1 + i) * 4;
				val = ((arr[cc]) << 24) | ((arr[cc + 1]) << 16) | (arr[cc + 2] << 8) | arr[cc + 3];
				result[(i - R[0])*(j - R[2]) + j - R[2]] = double(val)*bscale + bzero;
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

/*Object^ JPFITS::FITSImage::ReadPrimaryNDimensionalData(String^ fullFileName, TypeCode^ convertToTypeCode, int &nAaxis, array<int>^ &nAxisN)
{
	FileStream^ fs = gcnew FileStream(fullFileName, FileMode::Open);
	ArrayList^ header = gcnew ArrayList();
	bool hasext;
	if (!FITSFILEOPS::SCANPRIMARYUNIT(fs, false, header, hasext))
	{
		fs->Close();
		throw gcnew Exception("File not formatted as FITS file.");
		return nullptr;
	}

	int BITPIX = -1, nAxis = -1 , nAxisn = 1;
	for (int i = 0; i < header->Count; i++)
	{
		String^ line = (String^)header[i];

		if (BITPIX == -1)
			if (line->Substring(0, 8)->Trim() == "BITPIX")
				BITPIX = ::Convert::ToInt32(line->Substring(10, 20));

		if (nAxis == -1)
			if (line->Substring(0, 8)->Trim() == "NAXIS")
			{
				nAxis = ::Convert::ToInt32(line->Substring(10, 20));
				nAxisN = gcnew array<int>(nAxis);
			}

		if (nAxisn <= nAxis)
			if (line->Substring(0, 8)->Trim() == ("NAXIS" + nAxisn.ToString()))
			{
				nAxisN[nAxisn - 1] = ::Convert::ToInt32(line->Substring(10, 20));
				nAxisn++;
			}
	}






	return gcnew Object();
}*/

void JPFITS::FITSImage::MAKEDEFHEADER()
{
	HEADERKEYS = gcnew array<String^>(9);
	HEADERKEYVALS = gcnew array<String^>(9);
	HEADERKEYCOMS = gcnew array<String^>(9);
	BZERO = 0;
	BSCALE = 1;
	BITPIX = -64;
	if (DIMAGE != nullptr)
	{
		NAXIS = 2;
		NAXIS1 = DIMAGE->GetLength(0);
		NAXIS2 = DIMAGE->GetLength(1);
	}
	else
	{
		NAXIS = 0;
		NAXIS1 = 0;
		NAXIS2 = 0;
	}

	HEADERKEYS[0] = "SIMPLE";
	HEADERKEYS[1] = "BITPIX";
	HEADERKEYS[2] = "NAXIS";
	HEADERKEYS[3] = "NAXIS1";
	HEADERKEYS[4] = "NAXIS2";
	HEADERKEYS[5] = "BZERO";
	HEADERKEYS[6] = "BSCALE";
	HEADERKEYS[7] = "EXTEND";
	HEADERKEYS[8] = "END";

	HEADERKEYVALS[0] = "T";
	HEADERKEYVALS[1] = "-64";
	HEADERKEYVALS[2] = NAXIS.ToString();
	HEADERKEYVALS[3] = NAXIS1.ToString();
	HEADERKEYVALS[4] = NAXIS2.ToString();
	HEADERKEYVALS[5] = "0";
	HEADERKEYVALS[6] = "1";
	HEADERKEYVALS[7] = "T";
	HEADERKEYVALS[8] = "";

	HEADERKEYCOMS[0] = "Valid Fits File";
	HEADERKEYCOMS[1] = "Bits per Pixel";
	HEADERKEYCOMS[2] = "Number of Axes";
	HEADERKEYCOMS[3] = "Width (No. of Columns)";
	HEADERKEYCOMS[4] = "Height (No. of Rows)";
	HEADERKEYCOMS[5] = "Data Offset";
	HEADERKEYCOMS[6] = "Data Scaling";
	HEADERKEYCOMS[7] = "File may contain extensions";
	HEADERKEYCOMS[8] = "";
}

void JPFITS::FITSImage::SetKey(String^ Key, String^ Value, bool AddIfNotFound, int AddAtIndex)
{
	Key = Key->ToUpper();
	for (int i = 0; i < HEADERKEYS->Length; i++)
	{
		if (Key == HEADERKEYS[i])
		{
			HEADERKEYVALS[i] = Value;
			return;
		}
	}
	if (AddIfNotFound)
		AddKey(Key, Value, "", AddAtIndex);
}

void JPFITS::FITSImage::SetKey(String^ Key, String^ Value, String^ Comment, bool AddIfNotFound, int AddAtIndex)
{
	Key = Key->ToUpper();
	for (int i = 0; i < HEADERKEYS->Length; i++)
	{
		if (Key == HEADERKEYS[i])
		{
			HEADERKEYVALS[i] = Value;
			HEADERKEYCOMS[i] = Comment;
			return;
		}
	}
	if (AddIfNotFound)
		AddKey(Key, Value, Comment, AddAtIndex);
}

void JPFITS::FITSImage::SetKey(int index, String^ Key, String^ Value, String^ Comment)
{
	Key = Key->ToUpper();
	if (index < 0 || index > HEADERKEYS->Length - 1)
		return;

	HEADERKEYS[index] = Key;
	HEADERKEYVALS[index] = Value;
	HEADERKEYCOMS[index] = Comment;
}

void JPFITS::FITSImage::AddKey(String^ NewKey, String^ NewValue, String^ NewComment, int KeyIndex)
{
	int L = HEADERKEYS->Length;
	if (KeyIndex < 0 || KeyIndex >= L)
		KeyIndex = L - 1;//add to end of header (before END)
	array<String^>^ headerkeys = gcnew array<String^>(L + 1);
	array<String^>^ headerkeyvals = gcnew array<String^>(L + 1);
	array<String^>^ headerkeycoms = gcnew array<String^>(L + 1);
	int c = 0;
	for (int i = 0; i < L + 1; i++)
	{
		if (i == KeyIndex)
		{
			headerkeys[i] = NewKey;
			headerkeyvals[i] = NewValue;
			headerkeycoms[i] = NewComment;
			continue;
		}
		headerkeys[i] = HEADERKEYS[c];
		headerkeyvals[i] = HEADERKEYVALS[c];
		headerkeycoms[i] = HEADERKEYCOMS[c];
		c++;
	}
	HEADERKEYS = gcnew array<String^>(L + 1);
	HEADERKEYVALS = gcnew array<String^>(L + 1);
	HEADERKEYCOMS = gcnew array<String^>(L + 1);
	HEADERKEYS = headerkeys;
	HEADERKEYVALS = headerkeyvals;
	HEADERKEYCOMS = headerkeycoms;
}

String^ JPFITS::FITSImage::GetKeyValue(String^ key)
{
	String^ result = "";
	for (int i = 0; i < HEADERKEYS->Length; i++)
	{
		if (key->CompareTo(HEADERKEYS[i]) == 0)
		{
			result = HEADERKEYVALS[i];
			break;
		}
	}
	return result;
}

String^ JPFITS::FITSImage::GetKeyComment(String^ key)
{
	String^ result = "";
	for (int i = 0; i < HEADERKEYS->Length; i++)
	{
		if (key->CompareTo(HEADERKEYS[i]) == 0)
		{
			result = HEADERKEYCOMS[i];
			break;
		}
	}
	return result;
}

int JPFITS::FITSImage::GetKeyIndex(String^ key, String^ keyvalue, String^ keycomment)
{
	int result = -1;
	for (int i = 0; i < HEADERKEYS->Length; i++)
		if (key->CompareTo(HEADERKEYS[i]) == 0)
			if (keyvalue->CompareTo(HEADERKEYVALS[i]) == 0)
				if (keycomment->CompareTo(HEADERKEYCOMS[i]) == 0)
				{
					result = i;
					break;
				}
	return result;
}

int JPFITS::FITSImage::GetKeyIndex(String^ key, String^ keyvalue)
{
	int result = -1;
	for (int i = 0; i < HEADERKEYS->Length; i++)
		if (key->CompareTo(HEADERKEYS[i]) == 0)
			if (keyvalue->CompareTo(HEADERKEYVALS[i]) == 0)
			{
				result = i;
				break;
			}
	return result;
}

int JPFITS::FITSImage::GetKeyIndex(String^ key)
{
	int result = -1;
	for (int i = 0; i < HEADERKEYS->Length; i++)
		if (key == HEADERKEYS[i])
			return i;

	return result;
}

String^ JPFITS::FITSImage::GetKeyName(int index)
{
	if (index < HEADERKEYS->Length)
		return HEADERKEYS[index];
	else
		return "";
}
String^ JPFITS::FITSImage::GetKeyValue(int index)
{
	if (index < HEADERKEYVALS->Length)
		return HEADERKEYVALS[index];
	else
		return "";
}
String^ JPFITS::FITSImage::GetKeyComment(int index)
{
	if (index < HEADERKEYCOMS->Length)
		return HEADERKEYCOMS[index];
	else
		return "";
}

void JPFITS::FITSImage::RemoveKey(int KeyIndex)
{
	if (KeyIndex > HEADERKEYS->Length - 1 || KeyIndex < 0)
		return;

	int c = -1;
	array<String^>^ keys = gcnew array<String^>(HEADERKEYS->Length - 1);
	array<String^>^ vals = gcnew array<String^>(HEADERKEYS->Length - 1);
	array<String^>^ coms = gcnew array<String^>(HEADERKEYS->Length - 1);
	for (int i = 0; i < HEADERKEYS->Length; i++)
	{
		if (i == KeyIndex)
			continue;
		c++;
		keys[c] = HEADERKEYS[i];
		vals[c] = HEADERKEYVALS[i];
		coms[c] = HEADERKEYCOMS[i];
	}
	HEADERKEYS = keys;
	HEADERKEYVALS = vals;
	HEADERKEYCOMS = coms;
}

void JPFITS::FITSImage::RemoveKey(String^ Key)
{
	RemoveKey(GetKeyIndex(Key));
}

void JPFITS::FITSImage::RemoveKey(String^ Key, String^ Value)
{
	RemoveKey(GetKeyIndex(Key, Value));
}

void JPFITS::FITSImage::RemoveAllKeys()
{
	MAKEDEFHEADER();
}

bool JPFITS::FITSImage::ValidKeyEdit(String^ EditingKey)
{
	return VALIDKEYEDIT(EditingKey);
}
bool JPFITS::FITSImage::VALIDKEYEDIT(String^ checkkey)
{
	bool result = true;
	array<String^>^ checks = { "SIMPLE","EXTEND","BITPIX","NAXIS","NAXIS1","NAXIS2","BZERO","BSCALE","END" };
	for (int i = 0; i < checks->Length; i++)
		if (checkkey->CompareTo(checks[i]) == 0)
			result = false;

	return result;
}

void JPFITS::FITSImage::CopyHeader(JPFITS::FITSImage ^source)
{
	for (int i = 0; i < source->HeaderKeys->Length; i++)
	{
		if (!VALIDKEYEDIT(source->HeaderKeys[i]))
			continue;

		this->AddKey(source->HeaderKeys[i], source->HeaderKeyValues[i], source->HeaderKeyComments[i], -1);
	}
}

