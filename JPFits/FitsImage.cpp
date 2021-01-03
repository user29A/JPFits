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

JPFITS::FITSImage::~FITSImage()
{
	delete DIMAGE, HEADERKEYS, HEADERKEYVALS, HEADERKEYCOMS;
}

JPFITS::FITSImage::FITSImage(String^ FullFileName, String^ DiskUCharBufferName, TypeCode Precision, int NAxis1, int NAxis2)
{
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

	if (DATA_POP == false)
	{
		Do_Stats = false;//can't do stats on data that isn't read
		STATS_POP = false;
	}

	FULLFILENAME = FullFileName;
	int index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	MIN = 0.0;
	MAX = 0.0;
	MEDIAN = 0.0;
	MEAN = 0.0;
	STD = 0.0;
	SUM = 0.0;

	FileStream^ fs = gcnew FileStream(FULLFILENAME,IO::FileMode::Open);

	READHEADER(fs,HEADER_POP);//reads and sets and populates the header ans properties and sets the position of fs to the beginning of the image data

	if (DATA_POP == true)
	{
		READDATA(fs, Range, do_parallel);
	}

	fs->Close();

	if (STATS_POP)
		StatsUpD(do_parallel);
}

JPFITS::FITSImage::FITSImage(String^ FullFileName, array<double,2>^ ImageData, bool Do_Stats, bool do_parallel)//create new fits image from new data array, using path_file filename
{
	int index;
	DIMAGE = ImageData;
	FULLFILENAME = FullFileName;
	index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	NAXIS = 2;
	NAXIS1 = ImageData->GetLength(0);
	NAXIS2 = ImageData->GetLength(1);
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

JPFITS::FITSImage::FITSImage(String^ FullFileName, array<double>^ ImageData, bool Do_Stats, bool do_parallel)//create new fits image from new data array, using path_file filename
{
	int index;
	FULLFILENAME = FullFileName;
	index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	NAXIS = 2;
	NAXIS1 = 1;
	NAXIS2 = ImageData->Length;
	MIN = 0.0;
	MAX = 0.0;
	MEDIAN = 0.0;
	MEAN = 0.0;
	STD = 0.0;
	SUM = 0.0;

	//convert integer image to double image
	DIMAGE = gcnew array<double,2>(NAXIS1,NAXIS2);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < NAXIS2; i++)
			DIMAGE[0,i] = (double)ImageData[i];

	HEADER_POP = true;
	DATA_POP = true;
	STATS_POP = Do_Stats;

	//create default header
	MAKEDEFHEADER();

	if (STATS_POP)
		StatsUpD(do_parallel);
}

JPFITS::FITSImage::FITSImage(System::String ^FullFileName, cli::array<unsigned __int32,2> ^ImageData, bool Do_Stats, bool do_parallel)
{
	int index;
	FULLFILENAME = FullFileName;
	index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	NAXIS = 2;
	NAXIS1 = ImageData->GetLength(0);
	NAXIS2 = ImageData->GetLength(1);
	MIN = 0.0;
	MAX = 0.0;
	MEDIAN = 0.0;
	MEAN = 0.0;
	STD = 0.0;
	SUM = 0.0;

	//convert integer image to double image
	DIMAGE = gcnew array<double,2>(NAXIS1,NAXIS2);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < NAXIS1; i++)
		for (int j = 0; j < NAXIS2; j++)
			DIMAGE[i, j] = (double)ImageData[i, j];

	HEADER_POP = true;
	DATA_POP = true;
	STATS_POP = Do_Stats;

	MAKEDEFHEADER();
	BITPIX = 32;
	SetKey("BITPIX",BITPIX.ToString(),false,-1);
	BZERO = 2147483648;
	SetKey("BZERO",BZERO.ToString(),false,-1);
	BSCALE = 1;
	SetKey("BSCALE",BSCALE.ToString(),false,-1);
	if (STATS_POP)
		StatsUpD(do_parallel);
}

JPFITS::FITSImage::FITSImage(System::String ^FullFileName, cli::array<unsigned __int32> ^ImageData, bool Do_Stats, bool do_parallel)
{
	int index;
	FULLFILENAME = FullFileName;
	index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	NAXIS = 2;
	NAXIS1 = 1; 
	NAXIS2 = ImageData->Length;
	MIN = 0.0;
	MAX = 0.0;
	MEDIAN = 0.0;
	MEAN = 0.0;
	STD = 0.0;
	SUM = 0.0;

	//convert integer image to double image
	DIMAGE = gcnew array<double,2>(NAXIS1,NAXIS2);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < NAXIS2; i++)
			DIMAGE[0,i] = (double)ImageData[i];

	HEADER_POP = true;
	DATA_POP = true;
	STATS_POP = Do_Stats;

	MAKEDEFHEADER();
	BITPIX = 32;
	SetKey("BITPIX",BITPIX.ToString(),false,-1);
	BZERO = 2147483648;
	SetKey("BZERO",BZERO.ToString(),false,-1);
	BSCALE = 1;
	SetKey("BSCALE",BSCALE.ToString(),false,-1);
	if (STATS_POP)
		StatsUpD(do_parallel);
}

JPFITS::FITSImage::FITSImage(System::String ^FullFileName, cli::array<__int32,2> ^ImageData, bool Do_Stats, bool do_parallel)
{
	int index;
	FULLFILENAME = FullFileName;
	index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	NAXIS = 2;
	NAXIS1 = ImageData->GetLength(0);
	NAXIS2 = ImageData->GetLength(1);
	MIN = 0.0;
	MAX = 0.0;
	MEDIAN = 0.0;
	MEAN = 0.0;
	STD = 0.0;
	SUM = 0.0;

	//convert integer image to double image
	DIMAGE = gcnew array<double,2>(NAXIS1,NAXIS2);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < NAXIS1; i++)
		for (int j = 0; j < NAXIS2; j++)
			DIMAGE[i,j] = (double)ImageData[i,j];

	HEADER_POP = true;
	DATA_POP = true;
	STATS_POP = Do_Stats;

	MAKEDEFHEADER();
	BITPIX = 32;
	SetKey("BITPIX",BITPIX.ToString(),false,-1);
	BZERO = 0;
	SetKey("BZERO",BZERO.ToString(),false,-1);
	BSCALE = 1;
	SetKey("BSCALE",BSCALE.ToString(),false,-1);
	if (STATS_POP)
		StatsUpD(do_parallel);
}

JPFITS::FITSImage::FITSImage(System::String ^FullFileName, cli::array<unsigned __int16,2> ^ImageData, bool Do_Stats, bool do_parallel)
{
	int index;
	FULLFILENAME = FullFileName;
	index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	NAXIS = 2;
	NAXIS1 = ImageData->GetLength(0);  
	NAXIS2 = ImageData->GetLength(1);
	MIN = 0.0;
	MAX = 0.0;
	MEDIAN = 0.0;
	MEAN = 0.0;
	STD = 0.0;
	SUM = 0.0;

	//convert integer image to double image
	DIMAGE = gcnew array<double,2>(NAXIS1,NAXIS2);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < NAXIS1; i++)
		for (int j = 0; j < NAXIS2; j++)
			DIMAGE[i,j] = (double)ImageData[i,j];

	HEADER_POP = true;
	DATA_POP = true;
	STATS_POP = Do_Stats;

	//to do...create dummy header
	MAKEDEFHEADER();
	BITPIX = 16;
	SetKey("BITPIX",BITPIX.ToString(),false,-1);
	BZERO = 32768;
	SetKey("BZERO",BZERO.ToString(),false,-1);
	BSCALE = 1;
	SetKey("BSCALE",BSCALE.ToString(),false,-1);
	if (STATS_POP)
		StatsUpD(do_parallel);
}

JPFITS::FITSImage::FITSImage(System::String ^FullFileName, cli::array<__int16,2> ^ImageData, bool Do_Stats, bool do_parallel)
{
	int index;
	FULLFILENAME = FullFileName;
	index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	NAXIS = 2;
	NAXIS1 = ImageData->GetLength(0);
	NAXIS2 = ImageData->GetLength(1);
	MIN = 0.0;
	MAX = 0.0;
	MEDIAN = 0.0;
	MEAN = 0.0;
	STD = 0.0;
	SUM = 0.0;

	//convert integer image to double image
	DIMAGE = gcnew array<double,2>(NAXIS1,NAXIS2);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < NAXIS1; i++)
		for (int j = 0; j < NAXIS2; j++)
			DIMAGE[i,j] = (double)ImageData[i,j];

	HEADER_POP = true;
	DATA_POP = true;
	STATS_POP = Do_Stats;

	//to do...create dummy header
	MAKEDEFHEADER();
	BITPIX = 16;
	SetKey("BITPIX",BITPIX.ToString(),false,-1);
	BZERO = 0;
	SetKey("BZERO",BZERO.ToString(),false,-1);
	BSCALE = 1;
	SetKey("BSCALE",BSCALE.ToString(),false,-1);
	if (STATS_POP)
		StatsUpD(do_parallel);
}

JPFITS::FITSImage::FITSImage(String^ FullFileName, array<unsigned __int16>^ ImageData, bool Do_Stats, bool do_parallel)
{
	int index;
	FULLFILENAME = FullFileName;
	index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	NAXIS = 2;
	NAXIS1 = 1;
	NAXIS2 = ImageData->Length;
	MIN = 0.0;
	MAX = 0.0;
	MEDIAN = 0.0;
	MEAN = 0.0;
	STD = 0.0;
	SUM = 0.0;

	//convert integer image to double image
	DIMAGE = gcnew array<double,2>(NAXIS1,NAXIS2);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < NAXIS2; i++)
			DIMAGE[0,i] = (double)ImageData[i];

	HEADER_POP = true;
	DATA_POP = true;
	STATS_POP = Do_Stats;

	MAKEDEFHEADER();
	BITPIX = 16;
	SetKey("BITPIX",BITPIX.ToString(),false,-1);
	BZERO = 32768;
	SetKey("BZERO",BZERO.ToString(),false,-1);
	BSCALE = 1;
	SetKey("BSCALE",BSCALE.ToString(),false,-1);
	if (STATS_POP)
		StatsUpD(do_parallel);
}

JPFITS::FITSImage::FITSImage(String^ FullFileName, array<__int16>^ ImageData, bool Do_Stats, bool do_parallel)
{
	int index;
	FULLFILENAME = FullFileName;
	index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	NAXIS = 2;
	NAXIS1 = 1;
	NAXIS2 = ImageData->Length;
	MIN = 0.0;
	MAX = 0.0;
	MEDIAN = 0.0;
	MEAN = 0.0;
	STD = 0.0;
	SUM = 0.0;

	//convert integer image to double image
	DIMAGE = gcnew array<double,2>(NAXIS1,NAXIS2);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < NAXIS2; i++)
			DIMAGE[0,i] = (double)ImageData[i];

	HEADER_POP = true;
	DATA_POP = true;
	STATS_POP = Do_Stats;

	MAKEDEFHEADER();
	BITPIX = 16;
	SetKey("BITPIX",BITPIX.ToString(),false,-1);
	BZERO = 0;
	SetKey("BZERO",BZERO.ToString(),false,-1);
	BSCALE = 1;
	SetKey("BSCALE",BSCALE.ToString(),false,-1);
	if (STATS_POP)
		StatsUpD(do_parallel);
}

JPFITS::FITSImage::FITSImage(System::String ^FullFileName, cli::array<__int8,2> ^ImageData, bool Do_Stats, bool do_parallel)
{
	int index;
	FULLFILENAME = FullFileName;
	index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	NAXIS = 2;
	NAXIS1 = ImageData->GetLength(0);
	NAXIS2 = ImageData->GetLength(1);
	MIN = 0.0;
	MAX = 0.0;
	MEDIAN = 0.0;
	MEAN = 0.0;
	STD = 0.0;
	SUM = 0.0;

	//convert integer image to double image
	DIMAGE = gcnew array<double,2>(NAXIS1,NAXIS2);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < NAXIS1; i++)
		for (int j = 0; j < NAXIS2; j++)
			DIMAGE[i,j] = (double)ImageData[i,j];

	HEADER_POP = true;
	DATA_POP = true;
	STATS_POP = Do_Stats;

	//to do...create dummy header
	MAKEDEFHEADER();
	BITPIX = 8;
	SetKey("BITPIX",BITPIX.ToString(),false,-1);
	BZERO = 0;
	SetKey("BZERO",BZERO.ToString(),false,-1);
	BSCALE = 1;
	SetKey("BSCALE",BSCALE.ToString(),false,-1);
	if (STATS_POP)
		StatsUpD(do_parallel);
}

JPFITS::FITSImage::FITSImage(System::String ^FullFileName, cli::array<unsigned __int8,2> ^ImageData, bool Do_Stats, bool do_parallel)
{
	int index;
	FULLFILENAME = FullFileName;
	index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	NAXIS = 2;
	NAXIS1 = ImageData->GetLength(0);
	NAXIS2 = ImageData->GetLength(1);
	MIN = 0.0;
	MAX = 0.0;
	MEDIAN = 0.0;
	MEAN = 0.0;
	STD = 0.0;
	SUM = 0.0;

	//convert integer image to double image
	DIMAGE = gcnew array<double,2>(NAXIS1,NAXIS2);
	#pragma omp parallel for if (do_parallel)
	for (int i = 0; i < NAXIS1; i++)
		for (int j = 0; j < NAXIS2; j++)
			DIMAGE[i,j] = (double)ImageData[i,j];

	HEADER_POP = true;
	DATA_POP = true;
	STATS_POP = Do_Stats;

	//to do...create dummy header
	MAKEDEFHEADER();
	BITPIX = 8;
	SetKey("BITPIX",BITPIX.ToString(),false,-1);
	BZERO = 128;
	SetKey("BZERO",BZERO.ToString(),false,-1);
	BSCALE = 1;
	SetKey("BSCALE",BSCALE.ToString(),false,-1);
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

void JPFITS::FITSImage::WriteFile(System::String ^FullFileName, System::TypeCode Precision, bool do_parallel)
{
	FULLFILENAME = FullFileName;
	int index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);
	WRITEFILE(Precision, do_parallel);
}

void JPFITS::FITSImage::WriteFile(System::TypeCode Precision, bool do_parallel)
{
	WRITEFILE(Precision, do_parallel);
}

void JPFITS::FITSImage::WriteFileFromDiskBuffer(bool DeleteOrigDiskBuffer)
{
	if (FROMDISKBUFFER != true)
	{
		::MessageBox::Show("This FITS not created from ""JPFITS::FITSImage(String^ FullFileName, String^ DiskUCharBufferName, TypeCode Precision, int NAxis1, int NAxis2)"" constructor.","Error");
		return;
	}
	FORMATHEADER();
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

JPFITS::FITSImage::FITSImage(String^ FullFileName)
{
	int index;
	FULLFILENAME = FullFileName;
	index = FULLFILENAME->LastIndexOf("\\");
	FILENAME = FULLFILENAME->Substring(index+1);
	FILEPATH = FULLFILENAME->Substring(0,index+1);

	HEADERKEYS = gcnew array<String^>(8);
	HEADERKEYVALS = gcnew array<String^>(8);
	HEADERKEYCOMS = gcnew array<String^>(8);

	DIMAGE = gcnew array<double, 2>(0, 0);
	NAXIS = 0;
	BZERO = 0;
	BSCALE = 1;
	BITPIX = -64;
	NAXIS1 = 0;
	NAXIS2 = 0;

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
	HEADERKEYVALS[2] = NAXIS.ToString();
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

void JPFITS::FITSImage::ConvertToImage(String^ Source_FullFileName, String^ Destination_FullFileName, String^ contrast_scaling, bool invert_colormap, bool do_parallel)
{
	FITSImage^ fi = gcnew FITSImage(Source_FullFileName, nullptr, true, true, true, do_parallel);
	//JPBitMap::ArrayToBmp(fi->Image, 0, 0, invert_colormap, )
}

