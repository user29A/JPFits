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
	delete DIMAGE;
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

	DIMAGE = gcnew array<double, 2>(NAxis1, NAxis2);//default fill

	HEADER = gcnew FITSImageHeader(true, DIMAGE);
	this->Header->SetBITPIXNAXISBSCZ(Precision, DIMAGE);
	EATIMAGEHEADER();

	WORLDCOORDINATESOLUTION = gcnew JPFITS::WorldCoordinateSolution();
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
	HEADER = gcnew JPFITS::FITSImageHeader(header, HEADER_POP);
	EATIMAGEHEADER();

	if (DATA_POP == true)
		READIMAGE(fs, Range, do_parallel);

	fs->Close();

	if (STATS_POP)
		StatsUpD(do_parallel);

	WORLDCOORDINATESOLUTION = gcnew JPFITS::WorldCoordinateSolution();
}

JPFITS::FITSImage::FITSImage(String^ FullFileName, String^ extensionName, array<int, 1>^ Range, bool Populate_Header, bool Populate_Data, bool Do_Stats, bool do_parallel)
{
	FileStream^ fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Open);
	bool hasext;
	if (!FITSFILEOPS::SCANPRIMARYUNIT(fs, true, nullptr, hasext) || !hasext)
	{
		fs->Close();
		throw gcnew Exception("File not formatted as FITS file, or indicates no extensions present.");
	}

	ArrayList^ header = gcnew ArrayList();
	__int64 extensionstartposition, extensionendposition, tableendposition, pcount, theap;
	if (!FITSFILEOPS::SEEKEXTENSION(fs, "IMAGE", extensionName, header, extensionstartposition, extensionendposition, tableendposition, pcount, theap))
	{
		fs->Close();
		throw gcnew Exception("Could not find IMAGE extension with name '" + extensionName + "'");
		return;
	}
	HEADER = gcnew JPFITS::FITSImageHeader(header, HEADER_POP);
	EATIMAGEHEADER();

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

	WORLDCOORDINATESOLUTION = gcnew JPFITS::WorldCoordinateSolution();
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
		throw gcnew Exception("Error: Rank of Object 'ImageData' is not a 1D or 2D array. Use other interface functions for creating larger dimensional FITS primary data.");
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

		case TypeCode::Single:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS2; x++)
					DIMAGE[0, x] = (double)((array<float>^)ImageData)[x];
			}
			else
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS1; x++)
					for (int y = 0; y < NAXIS2; y++)
						DIMAGE[x, y] = (double)((array<float, 2>^)ImageData)[x, y];
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

		case TypeCode::UInt64:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS2; x++)
					DIMAGE[0, x] = (double)((array<unsigned __int64>^)ImageData)[x];
			}
			else
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS1; x++)
					for (int y = 0; y < NAXIS2; y++)
						DIMAGE[x, y] = (double)((array<unsigned __int64, 2>^)ImageData)[x, y];
			}

			break;
		}

		case TypeCode::Int64:
		{
			if (rank == 1)
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS2; x++)
					DIMAGE[0, x] = (double)((array<__int64>^)ImageData)[x];
			}
			else
			{
				#pragma omp parallel for if (do_parallel)
				for (int x = 0; x < NAXIS1; x++)
					for (int y = 0; y < NAXIS2; y++)
						DIMAGE[x, y] = (double)((array<__int64, 2>^)ImageData)[x, y];
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
	HEADER = gcnew FITSImageHeader(true, DIMAGE);
	EATIMAGEHEADER();

	if (STATS_POP)
		StatsUpD(do_parallel);

	WORLDCOORDINATESOLUTION = gcnew JPFITS::WorldCoordinateSolution();
}

void JPFITS::FITSImage::SetImage(array<double,2>^ ImageData, bool Do_Stats, bool do_parallel)//sometimes want to not bother with stats, so cant use property setter for images
{
	DIMAGE = ImageData;
	NAXIS1 = ImageData->GetLength(0);
	NAXIS2 = ImageData->GetLength(1);
	this->Header->SetKey("NAXIS1", NAXIS1.ToString(), false, -1);
	this->Header->SetKey("NAXIS2", NAXIS2.ToString(), false, -1);
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

void JPFITS::FITSImage::WriteImage(String^ FullFileName, String^ extensionName, bool overwriteExtensionIfExists, TypeCode Precision, bool do_parallel)
{
	ISEXTENSION = true;
	EXTNAME = extensionName;
	EXTNAMEOVERWRITE = overwriteExtensionIfExists;

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

	array<String^>^ HEADER = this->Header->GetFormattedHeaderBlock(false, false);
	FileStream^ fits_fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Create);
	
	array<unsigned char>^ head = gcnew array<unsigned char>(HEADER->Length * 80);
	for (int i = 0; i < HEADER->Length; i++)
		for (int j = 0; j < 80; j++)
			head[i * 80 + j] = (unsigned char)HEADER[i][j];

	fits_fs->Write(head, 0, head->Length);//header is written

	FileStream^ buff_fs = gcnew FileStream(DISKBUFFERFULLNAME, FileMode::Open, ::FileAccess::Read);

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

	HEADER = gcnew FITSImageHeader(true, nullptr);
	EATIMAGEHEADER();
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
	array<double, 2>^ result = gcnew array<double, 2>(X_HalfWidth * 2 + 1, Y_HalfWidth * 2 + 1);

	/*if (X_Center - X_HalfWidth < 0 || Y_Center - Y_HalfWidth < 0 || X_Center - X_HalfWidth + result->GetLength(0) > NAXIS1 || Y_Center - Y_HalfWidth + result->GetLength(1) > NAXIS2)
	{

	}*/

	for (int x = 0; x < result->GetLength(0); x++)
		for (int y = 0; y < result->GetLength(1); y++)
			result[x, y] = DIMAGE[X_Center - X_HalfWidth + x, Y_Center - Y_HalfWidth + y];

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

void JPFITS::FITSImage::EATIMAGEHEADER()
{
	BITPIX = Convert::ToInt32(HEADER->GetKeyValue("BITPIX"));
	NAXIS = Convert::ToInt32(HEADER->GetKeyValue("NAXIS"));
	if (HEADER->GetKeyValue("NAXIS1") != "")
		NAXIS1 = Convert::ToInt32(HEADER->GetKeyValue("NAXIS1"));
	if (HEADER->GetKeyValue("NAXIS2") != "")
		NAXIS2 = Convert::ToInt32(HEADER->GetKeyValue("NAXIS2"));
	if (HEADER->GetKeyValue("BZERO") != "")
		BZERO = Convert::ToDouble(HEADER->GetKeyValue("BZERO"));
	if (HEADER->GetKeyValue("BSCALE") != "")
		BSCALE = Convert::ToDouble(HEADER->GetKeyValue("BSCALE"));

	if (NAXIS == 0)
	{
		NAXIS1 = 0;
		NAXIS2 = 0;
	}
	else if (NAXIS == 1)
	{
		NAXIS2 = NAXIS1;
		NAXIS1 = 1;		
	}
	if (BZERO == -1)
		BZERO = 0;
	if (BSCALE == -1)
		BSCALE = 1;
}

void JPFITS::FITSImage::WRITEIMAGE(TypeCode Precision, bool do_parallel)
{	
	//try
	{
		FileStream^ fs;
		bool filexists = File::Exists(FULLFILENAME);
		array<unsigned char>^ prependdata;
		array<unsigned char>^ appenddata;
		if (!filexists && !ISEXTENSION)//then just write to a new file
			fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Create);
		else if (filexists && !ISEXTENSION)//then write the primary unit, and append any extensions if they already exist on the existing file
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
		else if (filexists && ISEXTENSION)//then get the primary data unit and also check for other extensions, etc
		{
			fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Open);
			bool hasext = false;
			if (!FITSFILEOPS::SCANPRIMARYUNIT(fs, true, nullptr, hasext))
			{
				fs->Close();
				throw gcnew Exception("Primary data unit of file is not a FITS file.");
				return;
			}

			__int64 extensionstartposition, extensionendposition, tableendposition, pcount, theap;
			bool extexists = FITSFILEOPS::SEEKEXTENSION(fs, "IMAGE", EXTNAME, nullptr, extensionstartposition, extensionendposition, tableendposition, pcount, theap);
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
				prependdata = gcnew array<unsigned char>((int)extensionstartposition);
				fs->Position = 0;
				fs->Read(prependdata, 0, prependdata->Length);
				appenddata = gcnew array<unsigned char>(int(fs->Length - extensionendposition));
				fs->Position = extensionendposition;
				fs->Read(appenddata, 0, appenddata->Length);
				fs->Close();
				fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Create);
				fs->Write(prependdata, 0, prependdata->Length);
			}
		}
		else if (!filexists && ISEXTENSION)//then write the extension to a new file with an empty primary unit
		{
			FITSImageHeader^ hed = gcnew FITSImageHeader(true, nullptr);
			array<String^>^ pheader = hed->GetFormattedHeaderBlock(true, false);//   MAKEPRIMARYDEFFORMATTEDHEADER(true);// FITSHEADER::MAKEPRIMARYDEFFORMATTEDHEADER(true);
			int NHead = pheader->Length * 80;
			prependdata = gcnew array<unsigned char>(NHead);
			for (int i = 0; i < pheader->Length; i++)
				for (int j = 0; j < 80; j++)
					prependdata[i * 80 + j] = (unsigned char)pheader[i][j];
			fs = gcnew FileStream(FULLFILENAME, IO::FileMode::Create);
			fs->Write(prependdata, 0, prependdata->Length);
		}

		//set header BZERO and BCSALE key values depending on prec type.
		this->Header->SetBITPIXNAXISBSCZ(Precision, DIMAGE);
		EATIMAGEHEADER();

		//get formatted header block
		array<String^>^ HEADER = this->Header->GetFormattedHeaderBlock(ISEXTENSION, false);// FITSHEADER::MAKEFORMATTEDIMAGEHEADER(this->Header->HeaderKeys, this->Header->HeaderKeyValues, this->Header->HeaderKeyComments, ISEXTENSION);
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

		int cc = 0;

		switch (Precision)
		{
			case TypeCode::Byte:
			case TypeCode::SByte:
			{
				__int8 val = 0;

				#pragma omp parallel for if (do_parallel) private(cc, val)
				for (int j = 0; j < NAXIS2; j++)
					for (int i = 0; i < NAXIS1; i++)
					{
						cc = NHead + (j*NAXIS1 + i) * 2;
						val = __int8((DIMAGE[i, j] - bzero) / bscale);
						data[cc] = val;
					}
				break;
			}

			case TypeCode::UInt16:
			case TypeCode::Int16:
			{
				__int16 val = 0;

				#pragma omp parallel for if (do_parallel) private(cc, val)
				for (int j = 0; j < NAXIS2; j++)
					for (int i = 0; i < NAXIS1; i++)
					{
						cc = NHead + (j*NAXIS1 + i) * 2;
						val = __int16((DIMAGE[i, j] - bzero) / bscale);
						data[cc] = ((val >> 8) & 0xff);
						data[cc + 1] = (val & 0xff);
					}
				break;
			}

			case TypeCode::UInt32:
			case TypeCode::Int32:
			{
				__int32 val = 0;

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
				break;
			}

			case TypeCode::UInt64:
			case TypeCode::Int64:
			{
				__int64 val = 0;

				#pragma omp parallel for if (do_parallel) private(cc, val)
				for (int j = 0; j < NAXIS2; j++)
					for (int i = 0; i < NAXIS1; i++)
					{
						cc = NHead + (j*NAXIS1 + i) * 8;
						val = __int64((DIMAGE[i, j] - bzero) / bscale);
						data[cc] = ((val >> 56) & 0xff);
						data[cc + 1] = ((val >> 48) & 0xff);
						data[cc + 2] = ((val >> 40) & 0xff);
						data[cc + 3] = ((val >> 32) & 0xff);
						data[cc + 4] = ((val >> 24) & 0xff);
						data[cc + 5] = ((val >> 16) & 0xff);
						data[cc + 6] = ((val >> 8) & 0xff);
						data[cc + 7] = (val & 0xff);
					}
				break;
			}

			case TypeCode::Single:
			{
				#pragma omp parallel for if (do_parallel) private(cc)
				for (int j = 0; j < NAXIS2; j++)
				{
					array<unsigned char>^ sng = gcnew array<unsigned char>(4);
					for (int i = 0; i < NAXIS1; i++)
					{
						cc = NHead + (j*NAXIS1 + i) * 4;
						sng = BitConverter::GetBytes(float((DIMAGE[i, j] - bzero) / bscale));	
						data[cc] = sng[3];
						data[cc + 1] = sng[2];
						data[cc + 2] = sng[1];
						data[cc + 3] = sng[0];
					}
				}
				break;
			}

			case TypeCode::Double:
			{
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
				break;
			}
		}

		fs->Write(data, 0, NBytesTot);
		if (appenddata != nullptr)
			fs->Write(appenddata, 0, appenddata->Length);
		fs->Close();
	}
	/*catch (Exception^ e)
	{
		MessageBox::Show(e->Data + "	" + e->InnerException + "	" + e->Message + "	" + e->Source + "	" + e->StackTrace + "	" + e->TargetSite);
	}*/
}

/*array<double, 2>^ JPFITS::FITSImage::READIMAGE(FileStream^ fs, array<int, 1>^ Range, int bitpix, int naxis1, int naxis2, double bzero, double bscale, bool do_parallel)
{

}*/

void JPFITS::FITSImage::READIMAGE(FileStream^ fs, array<int, 1>^ Range, bool do_parallel)//assumed READHEADER has been completed, and bs is at beginning of data block
{
	if (NAXIS <= 0)
	{
		NAXIS1 = 0;
		NAXIS2 = 0;
		if (HEADER_POP)
		{
			this->Header->SetKey("NAXIS1", "0", false, 0);
			this->Header->SetKey("NAXIS2", "0", false, 0);
		}
		DIMAGE = gcnew array<double, 2>(0, 0);
		return;
	}

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

	DIMAGE = gcnew array<double, 2>(R[1] - R[0] + 1, R[3] - R[2] + 1);
	array<unsigned char>^ arr = gcnew array<unsigned char>(NBytes);
	fs->Read(arr, 0, NBytes);//seems to be fastest to just read the entire data even if only subimage will be used
	int cc;

	switch (BITPIX)
	{
		case 8:
		{
			#pragma omp parallel for if (do_parallel) private(cc)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = j * NAXIS1 + R[0];
				for (int i = R[0]; i <= R[1]; i++)
					DIMAGE[i - R[0], j - R[2]] = double(arr[cc + i])*BSCALE + BZERO;
			}
			break;
		}

		case 16:
		{
			__int16 val;

			#pragma omp parallel for if (do_parallel) private(cc, val)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * NAXIS1 + R[0]) * 2;
				for (int i = R[0]; i <= R[1]; i++)
				{
					val = (arr[cc] << 8) | arr[cc + 1];
					DIMAGE[i - R[0], j - R[2]] = double(val)*BSCALE + BZERO;
					cc += 2;
				}
			}
			break;
		}

		case 32:
		{
			__int32 val;

			#pragma omp parallel for if (do_parallel) private(cc, val)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * NAXIS1 + R[0]) * 4;
				for (int i = R[0]; i <= R[1]; i++)
				{					
					val = (arr[cc] << 24) | (arr[cc + 1] << 16) | (arr[cc + 2] << 8) | arr[cc + 3];
					DIMAGE[i - R[0], j - R[2]] = double(val)*BSCALE + BZERO;
					cc += 4;
				}
			}
			break;
		}

		case 64:
		{
			#pragma omp parallel for if (do_parallel) private(cc)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * NAXIS1 + R[0]) * 8;
				array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
				for (int i = R[0]; i <= R[1]; i++)
				{					
					dbl[7] = arr[cc];
					dbl[6] = arr[cc + 1];
					dbl[5] = arr[cc + 2];
					dbl[4] = arr[cc + 3];
					dbl[3] = arr[cc + 4];
					dbl[2] = arr[cc + 5];
					dbl[1] = arr[cc + 6];
					dbl[0] = arr[cc + 7];
					DIMAGE[i - R[0], j - R[2]] = double(BitConverter::ToInt64(dbl, 0))*BSCALE + BZERO;
					cc += 8;
				}
			}
			break;
		}

		case -32:
		{
			#pragma omp parallel for if (do_parallel) private(cc)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * NAXIS1 + R[0]) * 4;
				array<unsigned char>^ flt = gcnew array<unsigned char>(4);
				for (int i = R[0]; i <= R[1]; i++)
				{					
					flt[3] = arr[cc];
					flt[2] = arr[cc + 1];
					flt[1] = arr[cc + 2];
					flt[0] = arr[cc + 3];
					DIMAGE[i - R[0], j - R[2]] = double(BitConverter::ToSingle(flt, 0))*BSCALE + BZERO;
					cc += 4;
				}
			}
			break;
		}

		case -64:
		{
			#pragma omp parallel for if (do_parallel) private(cc)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * NAXIS1 + R[0]) * 8;
				array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
				for (int i = R[0]; i <= R[1]; i++)
				{					
					dbl[7] = arr[cc];
					dbl[6] = arr[cc + 1];
					dbl[5] = arr[cc + 2];
					dbl[4] = arr[cc + 3];
					dbl[3] = arr[cc + 4];
					dbl[2] = arr[cc + 5];
					dbl[1] = arr[cc + 6];
					dbl[0] = arr[cc + 7];
					DIMAGE[i - R[0], j - R[2]] = (BitConverter::ToDouble(dbl, 0))*BSCALE + BZERO;
					cc += 8;
				}
			}
			break;
		}

		/*//subimage reading only...seems much slower than just reading entire image
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
	}

	if (Range != nullptr && Range[0] != -1)//then it is a full frame read
	{
		NAXIS1 = DIMAGE->GetLength(0);
		NAXIS2 = DIMAGE->GetLength(1);
		if (NAXIS2 == 1)
		{
			NAXIS2 = NAXIS1;
			NAXIS1 = 1;
		}
		this->Header->SetKey("NAXIS1", NAXIS1.ToString(), false, 0);
		this->Header->SetKey("NAXIS2", NAXIS2.ToString(), false, 0);
	}
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

	array<double, 2>^ result = gcnew array<double, 2>(R[1] - R[0] + 1, R[3] - R[2] + 1);
	array<unsigned char>^ arr = gcnew array<unsigned char>(NBytes);
	fs->Read(arr, 0, NBytes);//seems to be fastest to just read the entire data even if only subimage will be used
	fs->Close();
	int cc;

	switch (bitpix)
	{
		case 8:
		{
			#pragma omp parallel for if (do_parallel) private(cc)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = j * naxis1 + R[0];
				for (int i = R[0]; i <= R[1]; i++)
					result[i - R[0], j - R[2]] = double(arr[cc + i])*bscale + bzero;
			}
			break;
		}

		case 16:
		{
			__int16 val;

			#pragma omp parallel for if (do_parallel) private(cc, val)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * naxis1 + R[0]) * 2;
				for (int i = R[0]; i <= R[1]; i++)
				{
					val = (arr[cc] << 8) | arr[cc + 1];
					result[i - R[0], j - R[2]] = double(val)*bscale + bzero;
					cc += 2;
				}
			}
			break;
		}

		case 32:
		{
			__int32 val;

			#pragma omp parallel for if (do_parallel) private(cc, val)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * naxis1 + R[0]) * 4;
				for (int i = R[0]; i <= R[1]; i++)
				{
					val = (arr[cc] << 24) | (arr[cc + 1] << 16) | (arr[cc + 2] << 8) | arr[cc + 3];
					result[i - R[0], j - R[2]] = double(val)*bscale + bzero;
					cc += 4;
				}
			}
			break;
		}

		case 64:
		{
			#pragma omp parallel for if (do_parallel) private(cc)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * naxis1 + R[0]) * 8;
				array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
				for (int i = R[0]; i <= R[1]; i++)
				{
					dbl[7] = arr[cc];
					dbl[6] = arr[cc + 1];
					dbl[5] = arr[cc + 2];
					dbl[4] = arr[cc + 3];
					dbl[3] = arr[cc + 4];
					dbl[2] = arr[cc + 5];
					dbl[1] = arr[cc + 6];
					dbl[0] = arr[cc + 7];
					result[i - R[0], j - R[2]] = double(BitConverter::ToInt64(dbl, 0))*bscale + bzero;
					cc += 8;
				}
			}
			break;
		}

		case -32:
		{
			#pragma omp parallel for if (do_parallel) private(cc)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * naxis1 + R[0]) * 4;
				array<unsigned char>^ flt = gcnew array<unsigned char>(4);
				for (int i = R[0]; i <= R[1]; i++)
				{
					flt[3] = arr[cc];
					flt[2] = arr[cc + 1];
					flt[1] = arr[cc + 2];
					flt[0] = arr[cc + 3];
					result[i - R[0], j - R[2]] = double(BitConverter::ToSingle(flt, 0))*bscale + bzero;
					cc += 4;
				}
			}
			break;
		}

		case -64:
		{
			#pragma omp parallel for if (do_parallel) private(cc)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * naxis1 + R[0]) * 8;
				array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
				for (int i = R[0]; i <= R[1]; i++)
				{
					dbl[7] = arr[cc];
					dbl[6] = arr[cc + 1];
					dbl[5] = arr[cc + 2];
					dbl[4] = arr[cc + 3];
					dbl[3] = arr[cc + 4];
					dbl[2] = arr[cc + 5];
					dbl[1] = arr[cc + 6];
					dbl[0] = arr[cc + 7];
					result[i - R[0], j - R[2]] = (BitConverter::ToDouble(dbl, 0))*bscale + bzero;
					cc += 8;
				}
			}
			break;
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

	array<double>^ result = gcnew array<double>((R[1] - R[0] + 1) * (R[3] - R[2] + 1));
	array<unsigned char>^ arr = gcnew array<unsigned char>(NBytes);
	fs->Read(arr, 0, NBytes);//seems to be fastest to just read the entire data even if only subimage will be used
	fs->Close();
	int cc, jmr2;

	switch (bitpix)
	{
		case 8:
		{
			#pragma omp parallel for if (do_parallel) private(cc, jmr2)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = j * naxis1 + R[0];
				jmr2 = j - R[2];
				for (int i = R[0]; i <= R[1]; i++)
					result[(i - R[0]) * jmr2 + jmr2] = double(arr[cc + i])*bscale + bzero;
			}
			break;
		}

		case 16:
		{
			__int16 val;

			#pragma omp parallel for if (do_parallel) private(cc, val, jmr2)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * naxis1 + R[0]) * 2;
				jmr2 = j - R[2];
				for (int i = R[0]; i <= R[1]; i++)
				{
					val = (arr[cc] << 8) | arr[cc + 1];
					result[(i - R[0]) * jmr2 + jmr2] = double(val)*bscale + bzero;
					cc += 2;
				}
			}
			break;
		}

		case 32:
		{
			__int32 val;

			#pragma omp parallel for if (do_parallel) private(cc, val, jmr2)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * naxis1 + R[0]) * 4;
				jmr2 = j - R[2];
				for (int i = R[0]; i <= R[1]; i++)
				{
					val = (arr[cc] << 24) | (arr[cc + 1] << 16) | (arr[cc + 2] << 8) | arr[cc + 3];
					result[(i - R[0]) * jmr2 + jmr2] = double(val)*bscale + bzero;
					cc += 4;
				}
			}
			break;
		}

		case 64:
		{
			#pragma omp parallel for if (do_parallel) private(cc, jmr2)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * naxis1 + R[0]) * 8;
				jmr2 = j - R[2];
				array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
				for (int i = R[0]; i <= R[1]; i++)
				{
					dbl[7] = arr[cc];
					dbl[6] = arr[cc + 1];
					dbl[5] = arr[cc + 2];
					dbl[4] = arr[cc + 3];
					dbl[3] = arr[cc + 4];
					dbl[2] = arr[cc + 5];
					dbl[1] = arr[cc + 6];
					dbl[0] = arr[cc + 7];
					result[(i - R[0]) * jmr2 + jmr2] = double(BitConverter::ToInt64(dbl, 0))*bscale + bzero;
					cc += 8;
				}
			}
			break;
		}

		case -32:
		{
			#pragma omp parallel for if (do_parallel) private(cc, jmr2)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * naxis1 + R[0]) * 4;
				jmr2 = j - R[2];
				array<unsigned char>^ flt = gcnew array<unsigned char>(4);
				for (int i = R[0]; i <= R[1]; i++)
				{
					flt[3] = arr[cc];
					flt[2] = arr[cc + 1];
					flt[1] = arr[cc + 2];
					flt[0] = arr[cc + 3];
					result[(i - R[0]) * jmr2 + jmr2] = double(BitConverter::ToSingle(flt, 0))*bscale + bzero;
					cc += 4;
				}
			}
			break;
		}

		case -64:
		{
			#pragma omp parallel for if (do_parallel) private(cc, jmr2)
			for (int j = R[2]; j <= R[3]; j++)
			{
				cc = (j * naxis1 + R[0]) * 8;
				jmr2 = j - R[2];
				array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
				for (int i = R[0]; i <= R[1]; i++)
				{
					dbl[7] = arr[cc];
					dbl[6] = arr[cc + 1];
					dbl[5] = arr[cc + 2];
					dbl[4] = arr[cc + 3];
					dbl[3] = arr[cc + 4];
					dbl[2] = arr[cc + 5];
					dbl[1] = arr[cc + 6];
					dbl[0] = arr[cc + 7];
					result[(i - R[0]) * jmr2 + jmr2] = (BitConverter::ToDouble(dbl, 0))*bscale + bzero;
					cc += 8;
				}
			}
			break;
		}
	}

	return result;
}

array<double>^ JPFITS::FITSImage::ReadPrimaryNDimensionalData(String^ fullFileName, array<int>^ &nAxisN)
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

	int bitpix = -1, nAxis = -1, nAxisn = 1;
	double bscale = 1, bzero = 0;
	for (int i = 0; i < header->Count; i++)
	{
		String^ line = (String^)header[i];

		if (bitpix == -1)
			if (line->Substring(0, 8)->Trim() == "BITPIX")
				bitpix = ::Convert::ToInt32(line->Substring(10, 20));

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

		if (line->Substring(0, 8) == "BZERO   ")
			bzero = ::Convert::ToDouble(line->Substring(10, 20));

		if (line->Substring(0, 8) == "BSCALE  ")
			bscale = ::Convert::ToDouble(line->Substring(10, 20));
	}

	int NBytes = Math::Abs(bitpix) / 8;
	for (int i = 0; i < nAxisN->Length; i++)
		NBytes *= nAxisN[i];
	
	array<unsigned char>^ arr = gcnew array<unsigned char>(NBytes);
	fs->Read(arr, 0, NBytes);
	fs->Close();

	array<double>^ result = gcnew array<double>(NBytes / Math::Abs(bitpix));

	if (bitpix == 8)
	{
		#pragma omp parallel for
		for (int i = 0; i < result->Length; i++)
			result[i] = (double)arr[i] * bscale + bzero;
	}

	if (bitpix == 16)
	{
		int cc = 0;
		__int16 val = 0;

		#pragma omp parallel for private(cc, val)
		for (int i = 0; i < result->Length; i++)
		{
			cc = i * 2;
			val = (arr[cc] << 8) | arr[cc + 1];
			result[i] = (double)val * bscale + bzero;
		}
	}

	if (bitpix == 32)
	{
		int cc = 0;
		__int32 val = 0;

		#pragma omp parallel for private(cc, val)
		for (int i = 0; i < result->Length; i++)
		{
			cc = i * 4;
			val = (arr[cc] << 24) | (arr[cc + 1] << 16) | (arr[cc + 2] << 8) | arr[cc + 3];
			result[i] = double(val) * bscale + bzero;
		}
	}

	if (bitpix == -32)
	{
		float val = 0;
		int cc = 0;

		#pragma omp parallel for private(cc, val)
		for (int i = 0; i < result->Length; i++)
		{
			array<unsigned char>^ flt = gcnew array<unsigned char>(4);
			cc = i * 4;
			flt[3] = arr[cc];
			flt[2] = arr[cc + 1];
			flt[1] = arr[cc + 2];
			flt[0] = arr[cc + 3];
			val = BitConverter::ToSingle(flt, 0);
			result[i] = double(val) * bscale + bzero;
		}
	}

	if (bitpix == -64)
	{
		double val = 0;
		int cc = 0;

		#pragma omp parallel for private(cc, val)
		for (int i = 0; i < result->Length; i++)
		{
			cc = i * 8;
			array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
			dbl[7] = arr[cc];
			dbl[6] = arr[cc + 1];
			dbl[5] = arr[cc + 2];
			dbl[4] = arr[cc + 3];
			dbl[3] = arr[cc + 4];
			dbl[2] = arr[cc + 5];
			dbl[1] = arr[cc + 6];
			dbl[0] = arr[cc + 7];
			val = BitConverter::ToDouble(dbl, 0);
			result[i] = val * bscale + bzero;
		}
	}

	return result;
}

void JPFITS::FITSImage::ExtendizePrimaryImageLayerCube(String^ sourceFullFileName, String^ destFullFileName, array<String^>^ layerExtensionNames)
{
	for (int i = 0; i < layerExtensionNames->Length - 1; i++)
		for (int j = i + 1; j < layerExtensionNames->Length; j++)
			if (layerExtensionNames[i] == layerExtensionNames[j])
			{
				throw gcnew Exception("layerExtensionNames are not all unique");
				return;
			}
	for (int i = 0; i < layerExtensionNames->Length; i++)
		if (layerExtensionNames[i] == "")
		{
			throw gcnew Exception("layerExtensionNames cannot contain a nameless extension (empty string)");
			return;
		}

	array<int>^ axesN;
	array<double>^ cube = FITSImage::ReadPrimaryNDimensionalData(sourceFullFileName, axesN);

	if (layerExtensionNames->Length != axesN->Length)
	{
		throw gcnew Exception("layerExtensionNames array not equal in length to the number of layers");
		return;
	}

	if (destFullFileName == sourceFullFileName)
		File::Delete(destFullFileName);
	JPFITS::FITSImage^ fi = gcnew JPFITS::FITSImage(destFullFileName, true);
	fi->WriteImage(TypeCode::Double, false);

	for (int z = 0; z < axesN[2]; z++)//z is each layer of the cube
	{
		array<double, 2>^ layer = gcnew array<double, 2>(axesN[0], axesN[1]);

		#pragma omp parallel for
		for (int y = 0; y < axesN[1]; y++)
			for (int x = 0; x < axesN[0]; x++)
				layer[x, y] = cube[z * axesN[1] * axesN[0] + y * axesN[0] + x];

		fi = gcnew FITSImage(destFullFileName, layer, false, true);
		fi->WriteImage(destFullFileName, layerExtensionNames[z], false, TypeCode::Double, true);
	}
}

