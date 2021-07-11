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

JPFITS::FITSBinTable::FITSBinTable()
{

}

JPFITS::FITSBinTable::FITSBinTable(String^ fileName, String^ extensionName)
{
	FileStream^ fs = gcnew FileStream(fileName, IO::FileMode::Open);
	bool hasext = false;
	if (!FITSFILEOPS::SCANPRIMARYUNIT(fs, true, nullptr, hasext) || !hasext)
	{
		fs->Close();
		if (!hasext)
			throw gcnew Exception("File indicates no extensions present.");
		else
			throw gcnew Exception("File not formatted as FITS file.");
		return;
	}

	ArrayList^ header = gcnew ArrayList();
	__int64 extensionstartposition, extensionendposition, tableendposition, pcount, theap;
	if (!FITSFILEOPS::SEEKEXTENSION(fs, "BINTABLE", extensionName, header, extensionstartposition, extensionendposition, tableendposition, pcount, theap))
	{
		fs->Close();
		throw gcnew Exception("Could not find BINTABLE with name '" + extensionName + "'");
		return;
	}

	FILENAME = fileName;
	EXTENSIONNAME = extensionName;

	BINTABLE = gcnew array<unsigned char>(int(tableendposition - fs->Position));

	fs->Read(BINTABLE, 0, BINTABLE->Length);

	if (pcount != 0)
	{
		fs->Position = fs->Position + theap - BINTABLE->Length;
		HEAPDATA = gcnew array<unsigned char>(int(tableendposition + pcount - fs->Position));
		fs->Read(HEAPDATA, 0, HEAPDATA->Length);
	}

	fs->Close();

	EATRAWBINTABLEHEADER(header);
}

bool JPFITS::FITSBinTable::TTYPEEntryExists(String^ ttypeEntry)
{
	for (int i = 0; i < TTYPES->Length; i++)
		if (TTYPES[i] == ttypeEntry)
			return true;
	return false;
}

array<String^>^ JPFITS::FITSBinTable::GetAllExtensionNames(String^ FileName)
{
	return FITSFILEOPS::GETALLEXTENSIONNAMES(FileName, "BINTABLE");
}

array<double>^ JPFITS::FITSBinTable::GetTTYPEEntry(String^ ttypeEntry)
{
	int ttypeindex = -1;
	for (int i = 0; i < TTYPES->Length; i++)
		if (TTYPES[i] == ttypeEntry)
		{
			ttypeindex = i;
			break;
		}

	if (ttypeindex == -1)
	{
		throw gcnew Exception("Extension Entry TTYPE Label wasn't found: '" + ttypeEntry + "'");
		return nullptr;
	}

	if (TTYPEISCOMPLEX[ttypeindex] || TTYPEISHEAPARRAYDESC[ttypeindex])
	{
		String^ type;
		if (TTYPEISCOMPLEX[ttypeindex])
			type = "COMPLEX";
		else
			type = "HEAP ARRAY DESCRIPTOR";

		throw gcnew Exception("Cannot return entry '" + ttypeEntry + "' as a vector because it is " + type + ". Use an overload to get as an Object with dimensions returned to user.");
		return nullptr;
	}

	array<int>^ dimNElements;
	return GetTTYPEEntry(ttypeEntry, dimNElements);
}

array<double>^ JPFITS::FITSBinTable::GetTTYPEEntry(String^ ttypeEntry, array<int>^ &dimNElements)
{
	int ttypeindex = -1;
	for (int i = 0; i < TTYPES->Length; i++)
		if (TTYPES[i] == ttypeEntry)
			if (TCODES[i] == TypeCode::Char)
			{
				throw gcnew Exception("Cannot return entry '" + ttypeEntry + "' as double array because it is a char (String) array. Use an overload to get as an Object and then cast the object as a (array<String^>^)Object.");
				return nullptr;
			}	

	TypeCode tcode;
	Object^ obj = GetTTYPEEntry(ttypeEntry, tcode, dimNElements);
	int rank = ((Array^)obj)->Rank;
	int width, height;
	if (rank == 1)
	{
		width = 1;
		height = ((Array^)obj)->Length;
	}
	else
	{
		width = ((Array^)obj)->GetLength(0);
		height = ((Array^)obj)->GetLength(1);
	}
	array<double>^ result = gcnew array<double>(width * height);

	switch (tcode)
	{
		case TypeCode::Double:
		{
			if (rank == 1)
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					result[y] = ((array<double>^)(obj))[y];
			}
			else
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					for (int x = 0; x < width; x++)
						result[y * width + x] = ((array<double, 2>^)(obj))[x, y];
			}
			break;
		}

		case (TypeCode::Int64):
		{
			if (rank == 1)
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					result[y] = (double)((array<__int64>^)(obj))[y];
			}
			else
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					for (int x = 0; x < width; x++)
						result[y * width + x] = (double)((array<__int64, 2>^)(obj))[x, y];
			}
			break;
		}

		case (TypeCode::UInt64):
		{
			if (rank == 1)
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					result[y] = (double)((array<unsigned __int64>^)(obj))[y];
			}
			else
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					for (int x = 0; x < width; x++)
						result[y * width + x] = (double)((array<unsigned __int64, 2>^)(obj))[x, y];
			}
			break;
		}

		case TypeCode::Single:
		{
			if (rank == 1)
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					result[y] = (double)((array<float>^)(obj))[y];
			}
			else
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					for (int x = 0; x < width; x++)
						result[y * width + x] = (double)((array<float, 2>^)(obj))[x, y];
			}
			break;
		}

		case TypeCode::UInt32:
		{
			if (rank == 1)
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					result[y] = (double)((array<unsigned __int32>^)(obj))[y];
			}
			else
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					for (int x = 0; x < width; x++)
						result[y * width + x] = (double)((array<unsigned __int32, 2>^)(obj))[x, y];
			}
			break;
		}

		case TypeCode::Int32:
		{
			if (rank == 1)
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					result[y] = (double)((array<__int32>^)(obj))[y];
			}
			else
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					for (int x = 0; x < width; x++)
						result[y * width + x] = (double)((array<__int32, 2>^)(obj))[x, y];
			}
			break;
		}

		case TypeCode::UInt16:
		{
			if (rank == 1)
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					result[y] = (double)((array<unsigned __int16>^)(obj))[y];
			}
			else
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					for (int x = 0; x < width; x++)
						result[y * width + x] = (double)((array<unsigned __int16, 2>^)(obj))[x, y];
			}
			break;
		}

		case TypeCode::Int16:
		{
			if (rank == 1)
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					result[y] = (double)((array<__int16>^)(obj))[y];
			}
			else
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					for (int x = 0; x < width; x++)
						result[y * width + x] = (double)((array<__int16, 2>^)(obj))[x, y];
			}
			break;
		}

		case TypeCode::Byte:
		{
			if (rank == 1)
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					result[y] = (double)((array<unsigned __int8>^)(obj))[y];
			}
			else
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					for (int x = 0; x < width; x++)
						result[y * width + x] = (double)((array<unsigned __int8, 2>^)(obj))[x, y];
			}
			break;
		}

		case TypeCode::SByte:
		{
			if (rank == 1)
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					result[y] = (double)((array<__int8>^)(obj))[y];
			}
			else
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					for (int x = 0; x < width; x++)
						result[y * width + x] = (double)((array<__int8, 2>^)(obj))[x, y];
			}
			break;
		}

		case TypeCode::Boolean:
		{
			if (rank == 1)
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					result[y] = (double)((array<bool>^)(obj))[y];
			}
			else
			{
				#pragma omp parallel for
				for (int y = 0; y < height; y++)
					for (int x = 0; x < width; x++)
						result[y * width + x] = (double)((array<bool, 2>^)(obj))[x, y];
			}
			break;
		}

		default:
			throw gcnew Exception("Unrecognized TypeCode: '" + tcode.ToString() + "'");
	}

	return result;
}

Object^ JPFITS::FITSBinTable::GetTTYPEEntry(String^ ttypeEntry, TypeCode &objectTypeCode, array<int>^ &dimNElements)
{
	int ttypeindex = -1;
	for (int i = 0; i < TTYPES->Length; i++)
		if (TTYPES[i] == ttypeEntry)
		{
			ttypeindex = i;
			break;
		}

	if (ttypeindex == -1)
	{
		throw gcnew Exception("Extension Entry TTYPE Label wasn't found: '" + ttypeEntry + "'");
		return nullptr;
	}

	if (TTYPEISHEAPARRAYDESC[ttypeindex])//get from heap
		return GETHEAPTTYPE(ttypeindex, objectTypeCode, dimNElements);

	objectTypeCode = TCODES[ttypeindex];

	if (TDIMS[ttypeindex] != nullptr)
		dimNElements = TDIMS[ttypeindex];
	else
		if (TREPEATS[ttypeindex] == 1 || TCODES[ttypeindex] == TypeCode::Char)
			dimNElements = gcnew array<int>(1) { NAXIS2 };
		else
			dimNElements = gcnew array<int>(2) { TREPEATS[ttypeindex], NAXIS2 };

	int byteoffset = 0;
	for (int i = 0; i < ttypeindex; i++)
		byteoffset += TBYTES[i];
	int currentbyte;

	switch (TCODES[ttypeindex])
	{
		case ::TypeCode::Double:
		{
			if (TREPEATS[ttypeindex] == 1)
			{
				array<double>^ vector = gcnew array<double>(NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
					currentbyte = byteoffset + i * NAXIS1;
					dbl[7] = BINTABLE[currentbyte];
					dbl[6] = BINTABLE[currentbyte + 1];
					dbl[5] = BINTABLE[currentbyte + 2];
					dbl[4] = BINTABLE[currentbyte + 3];
					dbl[3] = BINTABLE[currentbyte + 4];
					dbl[2] = BINTABLE[currentbyte + 5];
					dbl[1] = BINTABLE[currentbyte + 6];
					dbl[0] = BINTABLE[currentbyte + 7];
					vector[i] = BitConverter::ToDouble(dbl, 0);
				}
				return vector;
			}
			else
			{
				array<double, 2>^ arrya = gcnew array<double, 2>(TREPEATS[ttypeindex], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 8;
						dbl[7] = BINTABLE[currentbyte];
						dbl[6] = BINTABLE[currentbyte + 1];
						dbl[5] = BINTABLE[currentbyte + 2];
						dbl[4] = BINTABLE[currentbyte + 3];
						dbl[3] = BINTABLE[currentbyte + 4];
						dbl[2] = BINTABLE[currentbyte + 5];
						dbl[1] = BINTABLE[currentbyte + 6];
						dbl[0] = BINTABLE[currentbyte + 7];
						arrya[j, i] = BitConverter::ToDouble(dbl, 0);
					}
				}
				return arrya;
			}
			break;
		}

		case ::TypeCode::Single:
		{
			if (TREPEATS[ttypeindex] == 1)
			{
				array<float>^ vector = gcnew array<float>(NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ sng = gcnew array<unsigned char>(4);
					currentbyte = byteoffset + i * NAXIS1;
					sng[3] = BINTABLE[currentbyte];
					sng[2] = BINTABLE[currentbyte + 1];
					sng[1] = BINTABLE[currentbyte + 2];
					sng[0] = BINTABLE[currentbyte + 3];
					vector[i] = BitConverter::ToSingle(sng, 0);
				}
				return vector;
			}
			else
			{
				array<float, 2>^ arrya = gcnew array<float, 2>(TREPEATS[ttypeindex], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ sng = gcnew array<unsigned char>(4);
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 4;
						sng[3] = BINTABLE[currentbyte];
						sng[2] = BINTABLE[currentbyte + 1];
						sng[1] = BINTABLE[currentbyte + 2];
						sng[0] = BINTABLE[currentbyte + 3];
						arrya[j, i] = BitConverter::ToSingle(sng, 0);
					}
				}
				return arrya;
			}
			break;
		}

		case (::TypeCode::Int64):
		{
			if (TREPEATS[ttypeindex] == 1)
			{
				array<__int64>^ vector = gcnew array<__int64>(NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ i64 = gcnew array<unsigned char>(8);
					currentbyte = byteoffset + i * NAXIS1;
					i64[7] = BINTABLE[currentbyte];
					i64[6] = BINTABLE[currentbyte + 1];
					i64[5] = BINTABLE[currentbyte + 2];
					i64[4] = BINTABLE[currentbyte + 3];
					i64[3] = BINTABLE[currentbyte + 4];
					i64[2] = BINTABLE[currentbyte + 5];
					i64[1] = BINTABLE[currentbyte + 6];
					i64[0] = BINTABLE[currentbyte + 7];
					vector[i] = BitConverter::ToInt64(i64, 0);
				}
				return vector;
			}
			else
			{
				array<__int64, 2>^ arrya = gcnew array<__int64, 2>(TREPEATS[ttypeindex], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ i64 = gcnew array<unsigned char>(8);
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 8;
						i64[7] = BINTABLE[currentbyte];
						i64[6] = BINTABLE[currentbyte + 1];
						i64[5] = BINTABLE[currentbyte + 2];
						i64[4] = BINTABLE[currentbyte + 3];
						i64[3] = BINTABLE[currentbyte + 4];
						i64[2] = BINTABLE[currentbyte + 5];
						i64[1] = BINTABLE[currentbyte + 6];
						i64[0] = BINTABLE[currentbyte + 7];
						arrya[j, i] = BitConverter::ToInt64(i64, 0);
					}
				}
				return arrya;
			}
			break;
		}

		case (::TypeCode::UInt64):
		{
			unsigned __int64 bzero = 9223372036854775808;
			if (TREPEATS[ttypeindex] == 1)
			{
				array<unsigned __int64>^ vector = gcnew array<unsigned __int64>(NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ ui64 = gcnew array<unsigned char>(8);
					currentbyte = byteoffset + i * NAXIS1;
					ui64[7] = BINTABLE[currentbyte];
					ui64[6] = BINTABLE[currentbyte + 1];
					ui64[5] = BINTABLE[currentbyte + 2];
					ui64[4] = BINTABLE[currentbyte + 3];
					ui64[3] = BINTABLE[currentbyte + 4];
					ui64[2] = BINTABLE[currentbyte + 5];
					ui64[1] = BINTABLE[currentbyte + 6];
					ui64[0] = BINTABLE[currentbyte + 7];
					vector[i] = BitConverter::ToInt64(ui64, 0) + bzero;
				}
				return vector;
			}
			else
			{
				array<unsigned __int64, 2>^ arrya = gcnew array<unsigned __int64, 2>(TREPEATS[ttypeindex], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ ui64 = gcnew array<unsigned char>(8);
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 8;
						ui64[7] = BINTABLE[currentbyte];
						ui64[6] = BINTABLE[currentbyte + 1];
						ui64[5] = BINTABLE[currentbyte + 2];
						ui64[4] = BINTABLE[currentbyte + 3];
						ui64[3] = BINTABLE[currentbyte + 4];
						ui64[2] = BINTABLE[currentbyte + 5];
						ui64[1] = BINTABLE[currentbyte + 6];
						ui64[0] = BINTABLE[currentbyte + 7];
						arrya[j, i] = BitConverter::ToInt64(ui64, 0) + bzero;
					}
				}
				return arrya;
			}
			break;
		}

		case ::TypeCode::UInt32:
		{
			unsigned __int32 bzero = 2147483648;
			if (TREPEATS[ttypeindex] == 1)
			{
				array<unsigned __int32>^ vector = gcnew array<unsigned __int32>(NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ uint32 = gcnew array<unsigned char>(4);
					currentbyte = byteoffset + i * NAXIS1;
					uint32[3] = BINTABLE[currentbyte];
					uint32[2] = BINTABLE[currentbyte + 1];
					uint32[1] = BINTABLE[currentbyte + 2];
					uint32[0] = BINTABLE[currentbyte + 3];
					vector[i] = BitConverter::ToInt32(uint32, 0) + bzero;
				}
				return vector;
			}
			else
			{
				array<unsigned __int32, 2>^ arrya = gcnew array<unsigned __int32, 2>(TREPEATS[ttypeindex], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ uint32 = gcnew array<unsigned char>(4);
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 4;
						uint32[3] = BINTABLE[currentbyte];
						uint32[2] = BINTABLE[currentbyte + 1];
						uint32[1] = BINTABLE[currentbyte + 2];
						uint32[0] = BINTABLE[currentbyte + 3];
						arrya[j, i] = BitConverter::ToInt32(uint32, 0) + bzero;
					}
				}
				return arrya;
			}
			break;
		}

		case ::TypeCode::Int32:
		{
			if (TREPEATS[ttypeindex] == 1)
			{
				array<__int32>^ vector = gcnew array<__int32>(NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ int32 = gcnew array<unsigned char>(4);
					currentbyte = byteoffset + i * NAXIS1;
					int32[3] = BINTABLE[currentbyte];
					int32[2] = BINTABLE[currentbyte + 1];
					int32[1] = BINTABLE[currentbyte + 2];
					int32[0] = BINTABLE[currentbyte + 3];
					vector[i] = BitConverter::ToInt32(int32, 0);
				}
				return vector;
			}
			else			
			{
				array<__int32, 2>^ arrya = gcnew array<__int32, 2>(TREPEATS[ttypeindex], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ int32 = gcnew array<unsigned char>(4);
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 4;
						int32[3] = BINTABLE[currentbyte];
						int32[2] = BINTABLE[currentbyte + 1];
						int32[1] = BINTABLE[currentbyte + 2];
						int32[0] = BINTABLE[currentbyte + 3];
						arrya[j, i] = BitConverter::ToInt32(int32, 0);
					}
				}
				return arrya;
			}
			break;
		}

		case ::TypeCode::UInt16:
		{
			unsigned __int16 bzero = 32768;
			if (TREPEATS[ttypeindex] == 1)
			{
				array<unsigned __int16>^ vector = gcnew array<unsigned __int16>(NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ uint16 = gcnew array<unsigned char>(2);
					currentbyte = byteoffset + i * NAXIS1;
					uint16[1] = BINTABLE[currentbyte];
					uint16[0] = BINTABLE[currentbyte + 1];
					vector[i] = BitConverter::ToInt16(uint16, 0) + bzero;
				}
				return vector;
			}
			else			
			{
				array<unsigned __int16, 2>^ arrya = gcnew array<unsigned __int16, 2>(TREPEATS[ttypeindex], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ uint16 = gcnew array<unsigned char>(2);
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 2;
						uint16[1] = BINTABLE[currentbyte];
						uint16[0] = BINTABLE[currentbyte + 1];
						arrya[j, i] = BitConverter::ToInt16(uint16, 0) + bzero;
					}
				}
				return arrya;
			}
			break;
		}

		case ::TypeCode::Int16:
		{
			if (TREPEATS[ttypeindex] == 1)
			{
				array<__int16>^ vector = gcnew array<__int16>(NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ int16 = gcnew array<unsigned char>(2);
					currentbyte = byteoffset + i * NAXIS1;
					int16[1] = BINTABLE[currentbyte];
					int16[0] = BINTABLE[currentbyte + 1];
					vector[i] = BitConverter::ToInt16(int16, 0);
				}
				return vector;
			}
			else			
			{
				array<__int16, 2>^ arrya = gcnew array<__int16, 2>(TREPEATS[ttypeindex], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ int16 = gcnew array<unsigned char>(2);
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 2;
						int16[1] = BINTABLE[currentbyte];
						int16[0] = BINTABLE[currentbyte + 1];
						arrya[j, i] = BitConverter::ToInt16(int16, 0);
					}
				}
				return arrya;
			}
			break;
		}

		case ::TypeCode::Byte:
		{
			if (TREPEATS[ttypeindex] == 1)
			{
				array<unsigned __int8>^ vector = gcnew array<unsigned __int8>(NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					currentbyte = byteoffset + i * NAXIS1;
					vector[i] = (unsigned __int8)BINTABLE[currentbyte];
				}
				return vector;
			}
			else			
			{
				array<unsigned __int8, 2>^ arrya = gcnew array<unsigned __int8, 2>(TREPEATS[ttypeindex], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j;
						arrya[j, i] = (unsigned __int8)BINTABLE[currentbyte];
					}
				return arrya;
			}
			break;
		}

		case ::TypeCode::SByte:
		{
			if (TREPEATS[ttypeindex] == 1)
			{
				array<__int8>^ vector = gcnew array<__int8>(NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					currentbyte = byteoffset + i * NAXIS1;
					vector[i] = (__int8)(BINTABLE[currentbyte]);
				}
				return vector;
			}
			else			
			{
				array<__int8, 2>^ arrya = gcnew array<__int8, 2>(TREPEATS[ttypeindex], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j;
						arrya[j, i] = (__int8)(BINTABLE[currentbyte]);
					}
				return arrya;
			}
			break;
		}

		case ::TypeCode::Boolean:
		{
			if (TREPEATS[ttypeindex] == 1)
			{
				array<bool>^ vector = gcnew array<bool>(NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					currentbyte = byteoffset + i * NAXIS1;
					vector[i] = Convert::ToBoolean(BINTABLE[currentbyte]);
				}
				return vector;
			}
			else			
			{
				array<bool, 2>^ arrya = gcnew array<bool, 2>(TREPEATS[ttypeindex], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j;
						arrya[j, i] = Convert::ToBoolean(BINTABLE[currentbyte]);
					}
				return arrya;
			}
			break;
		}

		case TypeCode::Char:
		{
			array<String^>^ vector = gcnew array<String^>(NAXIS2);
			#pragma omp parallel for private(currentbyte)
			for (int i = 0; i < NAXIS2; i++)
			{
				array<unsigned char>^ charstr = gcnew array<unsigned char>(TREPEATS[ttypeindex]);
				for (int j = 0; j < TREPEATS[ttypeindex]; j++)
				{
					currentbyte = byteoffset + i * NAXIS1 + j;
					charstr[j] = BINTABLE[currentbyte];
				}
				vector[i] = System::Text::Encoding::ASCII->GetString(charstr);
			}
			return vector;
			break;
		}

		default:
		{
			throw gcnew Exception("Unrecognized TypeCode: '" + TCODES[ttypeindex].ToString() + "'");
			return nullptr;
		}
	}
}

array<int, 2>^ JPFITS::FITSBinTable::GETHEAPTTYPENELSPOS(int ttypeindex)
{
	int byteoffset = 0;
	for (int i = 0; i < ttypeindex; i++)
		byteoffset += TBYTES[i];

	array<int, 2>^ result = gcnew array<int, 2>(2, NAXIS2);

	if (TFORMS[ttypeindex]->Contains("P"))
	{
		for (int i = 0; i < NAXIS2; i++)
		{
			int cc = byteoffset + i * NAXIS1;
			array<unsigned char>^ int32 = gcnew array<unsigned char>(4);
			int32[3] = BINTABLE[cc];
			int32[2] = BINTABLE[cc + 1];
			int32[1] = BINTABLE[cc + 2];
			int32[0] = BINTABLE[cc + 3];
			result[0, i] = BitConverter::ToInt32(int32, 0);//nelements
			int32[3] = BINTABLE[cc + 4];
			int32[2] = BINTABLE[cc + 5];
			int32[1] = BINTABLE[cc + 6];
			int32[0] = BINTABLE[cc + 7];
			result[1, i] = BitConverter::ToInt32(int32, 0);//position
		}
	}
	else//"Q"
	{
		for (int i = 0; i < NAXIS2; i++)
		{
			int cc = byteoffset + i * NAXIS1;
			array<unsigned char>^ i64 = gcnew array<unsigned char>(8);
			i64[7] = BINTABLE[cc];
			i64[6] = BINTABLE[cc + 1];
			i64[5] = BINTABLE[cc + 2];
			i64[4] = BINTABLE[cc + 3];
			i64[3] = BINTABLE[cc + 4];
			i64[2] = BINTABLE[cc + 5];
			i64[1] = BINTABLE[cc + 6];
			i64[0] = BINTABLE[cc + 7];
			result[0, i] = (int)BitConverter::ToInt64(i64, 0);//nelements
			i64[7] = BINTABLE[cc + 8];
			i64[6] = BINTABLE[cc + 9];
			i64[5] = BINTABLE[cc + 10];
			i64[4] = BINTABLE[cc + 11];
			i64[3] = BINTABLE[cc + 12];
			i64[2] = BINTABLE[cc + 13];
			i64[1] = BINTABLE[cc + 14];
			i64[0] = BINTABLE[cc + 15];
			result[1, i] = (int)BitConverter::ToInt64(i64, 0);//position
		}
	}
	return result;
}

Object^ JPFITS::FITSBinTable::GETHEAPTTYPE(int ttypeindex, TypeCode &objectTypeCode, array<int>^ &dimNElements)
{
	objectTypeCode = HEAPTCODES[ttypeindex];

	if (TDIMS[ttypeindex] != nullptr)
		dimNElements = TDIMS[ttypeindex];
	else
	{
		dimNElements = gcnew array<int>(2);
		dimNElements[1] = NAXIS2;
		int max = 0;

		for (int i = 0; i < NAXIS2; i++)
			if (TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i] > max)
				max = TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i];
		dimNElements[0] = max;
		if (TTYPEISCOMPLEX[ttypeindex])
			dimNElements[0] /= 2;
	}

	switch (objectTypeCode)
	{
		case ::TypeCode::Double:
		{
			array<array<double>^>^ arrya = gcnew array<array<double>^>(NAXIS2);
			//#pragma omp parallel for
			for (int i = 0; i < NAXIS2; i++)
			{
				array<double>^ row = gcnew array<double>(TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i]);
				int pos = (int)TTYPEHEAPARRAYNELSPOS[ttypeindex][1, i];
				array<unsigned char>^ dbl = gcnew array<unsigned char>(8);

				for (int j = 0; j < row->Length; j++)
				{					
					dbl[7] = HEAPDATA[pos];
					dbl[6] = HEAPDATA[pos + 1];
					dbl[5] = HEAPDATA[pos + 2];
					dbl[4] = HEAPDATA[pos + 3];
					dbl[3] = HEAPDATA[pos + 4];
					dbl[2] = HEAPDATA[pos + 5];
					dbl[1] = HEAPDATA[pos + 6];
					dbl[0] = HEAPDATA[pos + 7];
					row[j] = BitConverter::ToDouble(dbl, 0);
					pos += 8;
				}
				arrya[i] = row;
			}
			return arrya;
		}

		case ::TypeCode::Single:
		{
			array<array<float>^>^ arrya = gcnew array<array<float>^>(NAXIS2);
			//#pragma omp parallel for
			for (int i = 0; i < NAXIS2; i++)
			{
				array<float>^ row = gcnew array<float>(TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i]);
				int pos = (int)TTYPEHEAPARRAYNELSPOS[ttypeindex][1, i];
				array<unsigned char>^ sng = gcnew array<unsigned char>(4);
					
				for (int j = 0; j < row->Length; j++)
				{					
					sng[3] = HEAPDATA[pos];
					sng[2] = HEAPDATA[pos + 1];
					sng[1] = HEAPDATA[pos + 2];
					sng[0] = HEAPDATA[pos + 3];
					row[j] = BitConverter::ToSingle(sng, 0);
					pos += 4;
				}
				arrya[i] = row;
			}
			return arrya;
		}

		case (::TypeCode::Int64):
		{
			array<array<__int64>^>^ arrya = gcnew array<array<__int64>^>(NAXIS2);
			//#pragma omp parallel for
			for (int i = 0; i < NAXIS2; i++)
			{
				array<__int64>^ row = gcnew array<__int64>(TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i]);
				int pos = (int)TTYPEHEAPARRAYNELSPOS[ttypeindex][1, i];
				array<unsigned char>^ i64 = gcnew array<unsigned char>(8);

				for (int j = 0; j < row->Length; j++)
				{
					i64[7] = HEAPDATA[pos];
					i64[6] = HEAPDATA[pos + 1];
					i64[5] = HEAPDATA[pos + 2];
					i64[4] = HEAPDATA[pos + 3];
					i64[3] = HEAPDATA[pos + 4];
					i64[2] = HEAPDATA[pos + 5];
					i64[1] = HEAPDATA[pos + 6];
					i64[0] = HEAPDATA[pos + 7];
					row[j] = BitConverter::ToInt64(i64, 0);
					pos += 8;
				}
				arrya[i] = row;
			}
			return arrya;
		}

		case (::TypeCode::UInt64):
		{
			unsigned __int64 bzero = 9223372036854775808;
			array<array<unsigned __int64>^>^ arrya = gcnew array<array<unsigned __int64>^>(NAXIS2);
			//#pragma omp parallel for
			for (int i = 0; i < NAXIS2; i++)
			{
				array<unsigned __int64>^ row = gcnew array<unsigned __int64>(TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i]);
				int pos = (int)TTYPEHEAPARRAYNELSPOS[ttypeindex][1, i];
				array<unsigned char>^ ui64 = gcnew array<unsigned char>(8);

				for (int j = 0; j < row->Length; j++)
				{
					ui64[7] = HEAPDATA[pos];
					ui64[6] = HEAPDATA[pos + 1];
					ui64[5] = HEAPDATA[pos + 2];
					ui64[4] = HEAPDATA[pos + 3];
					ui64[3] = HEAPDATA[pos + 4];
					ui64[2] = HEAPDATA[pos + 5];
					ui64[1] = HEAPDATA[pos + 6];
					ui64[0] = HEAPDATA[pos + 7];
					row[j] = BitConverter::ToInt64(ui64, 0) + bzero;
					pos += 8;
				}
				arrya[i] = row;
			}
			return arrya;
		}

		case ::TypeCode::Int32:
		{
			array<array<__int32>^>^ arrya = gcnew array<array<__int32>^>(NAXIS2);
			//#pragma omp parallel for
			for (int i = 0; i < NAXIS2; i++)
			{
				array<__int32>^ row = gcnew array<__int32>(TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i]);
				int pos = (int)TTYPEHEAPARRAYNELSPOS[ttypeindex][1, i];
				array<unsigned char>^ int32 = gcnew array<unsigned char>(4);

				for (int j = 0; j < row->Length; j++)
				{					
					int32[3] = HEAPDATA[pos];
					int32[2] = HEAPDATA[pos + 1];
					int32[1] = HEAPDATA[pos + 2];
					int32[0] = HEAPDATA[pos + 3];
					row[j] = BitConverter::ToInt32(int32, 0);
					pos += 4;
				}
				arrya[i] = row;
			}
			return arrya;
		}

		case ::TypeCode::UInt32:
		{
			unsigned __int32 bzero = 2147483648;
			array<array<unsigned __int32>^>^ arrya = gcnew array<array<unsigned __int32>^>(NAXIS2);
			//#pragma omp parallel for
			for (int i = 0; i < NAXIS2; i++)
			{
				array<unsigned __int32>^ row = gcnew array<unsigned __int32>(TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i]);
				int pos = (int)TTYPEHEAPARRAYNELSPOS[ttypeindex][1, i];
				array<unsigned char>^ uint32 = gcnew array<unsigned char>(4);

				for (int j = 0; j < row->Length; j++)
				{
					uint32[3] = HEAPDATA[pos];
					uint32[2] = HEAPDATA[pos + 1];
					uint32[1] = HEAPDATA[pos + 2];
					uint32[0] = HEAPDATA[pos + 3];
					row[j] = BitConverter::ToInt32(uint32, 0) + bzero;
					pos += 4;
				}
				arrya[i] = row;
			}
			return arrya;
		}

		case ::TypeCode::Int16:
		{
			array<array<__int16>^>^ arrya = gcnew array<array<__int16>^>(NAXIS2);
			//#pragma omp parallel for
			for (int i = 0; i < NAXIS2; i++)
			{
				array<__int16>^ row = gcnew array<__int16>(TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i]);
				int pos = (int)TTYPEHEAPARRAYNELSPOS[ttypeindex][1, i];
				array<unsigned char>^ int16 = gcnew array<unsigned char>(2);

				for (int j = 0; j < row->Length; j++)
				{					
					int16[1] = HEAPDATA[pos];
					int16[0] = HEAPDATA[pos + 1];
					row[j] = BitConverter::ToInt16(int16, 0);
					pos += 2;
				}
				arrya[i] = row;
			}
			return arrya;
		}

		case ::TypeCode::UInt16:
		{
			unsigned __int16 bzero = 32768;
			array<array<unsigned __int16>^>^ arrya = gcnew array<array<unsigned __int16>^>(NAXIS2);
			//#pragma omp parallel for
			for (int i = 0; i < NAXIS2; i++)
			{
				array<unsigned __int16>^ row = gcnew array<unsigned __int16>(TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i]);
				int pos = (int)TTYPEHEAPARRAYNELSPOS[ttypeindex][1, i];
				array<unsigned char>^ int16 = gcnew array<unsigned char>(2);

				for (int j = 0; j < row->Length; j++)
				{
					int16[1] = HEAPDATA[pos];
					int16[0] = HEAPDATA[pos + 1];
					row[j] = BitConverter::ToInt16(int16, 0) + bzero;
					pos += 2;
				}
				arrya[i] = row;
			}
			return arrya;
		}

		case ::TypeCode::SByte:
		{
			array<array<__int8>^>^ arrya = gcnew array<array<__int8>^>(NAXIS2);
			//#pragma omp parallel for
			for (int i = 0; i < NAXIS2; i++)
			{
				array<__int8>^ row = gcnew array<__int8>(TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i]);
				int pos = (int)TTYPEHEAPARRAYNELSPOS[ttypeindex][1, i];

				for (int j = 0; j < row->Length; j++)
					row[j] = (__int8)HEAPDATA[pos + j];
				arrya[i] = row;
			}
			return arrya;
		}

		case ::TypeCode::Byte:
		{
			array<array<unsigned __int8>^>^ arrya = gcnew array<array<unsigned __int8>^>(NAXIS2);
			//#pragma omp parallel for
			for (int i = 0; i < NAXIS2; i++)
			{
				array<unsigned __int8>^ row = gcnew array<unsigned __int8>(TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i]);
				int pos = (int)TTYPEHEAPARRAYNELSPOS[ttypeindex][1, i];
				
				for (int j = 0; j < row->Length; j++)
					row[j] = (unsigned __int8)HEAPDATA[pos + j];
				arrya[i] = row;
			}
			return arrya;
		}

		case ::TypeCode::Boolean:
		{
			array<array<bool>^>^ arrya = gcnew array<array<bool>^>(NAXIS2);
			//#pragma omp parallel for
			for (int i = 0; i < NAXIS2; i++)
			{
				array<bool>^ row = gcnew array<bool>(TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i]);
				int pos = (int)TTYPEHEAPARRAYNELSPOS[ttypeindex][1, i];

				for (int j = 0; j < row->Length; j++)
					row[j] = (bool)HEAPDATA[pos + j];// Convert::ToBoolean(HEAPDATA[pos + j]);
				arrya[i] = row;
			}
			return arrya;
		}

		case ::TypeCode::Char:
		{
			array<String^>^ arrya = gcnew array<String^>(NAXIS2);
			//#pragma omp parallel for
			for (int i = 0; i < NAXIS2; i++)
				arrya[i] = System::Text::Encoding::ASCII->GetString(HEAPDATA, TTYPEHEAPARRAYNELSPOS[ttypeindex][1, i], TTYPEHEAPARRAYNELSPOS[ttypeindex][0, i]);

			return arrya;
		}

		default:
		{
			throw gcnew Exception("Unrecognized TypeCode: '" + objectTypeCode.ToString() + "'");
			return nullptr;
		}
	}

	return nullptr;
}

void JPFITS::FITSBinTable::REMOVEHEAPTTYPE(int ttypeindex)
{
	int startpos = TTYPEHEAPARRAYNELSPOS[ttypeindex][1, 0];
	int endpos = TTYPEHEAPARRAYNELSPOS[ttypeindex][1, NAXIS2 - 1];
	if (TTYPEISCOMPLEX[ttypeindex])
		endpos += TTYPEHEAPARRAYNELSPOS[ttypeindex][0, NAXIS2 - 1] * TYPECODETONBYTES(HEAPTCODES[ttypeindex]) * 2;
	else
		endpos += TTYPEHEAPARRAYNELSPOS[ttypeindex][0, NAXIS2 - 1] * TYPECODETONBYTES(HEAPTCODES[ttypeindex]);

	array<unsigned char>^ prepend = gcnew array<unsigned char>(startpos);
	for (int i = 0; i < prepend->Length; i++)
		prepend[i] = HEAPDATA[i];

	array<unsigned char>^ append = gcnew array<unsigned char>(HEAPDATA->Length - endpos);
	for (int i = 0; i < append->Length; i++)
		append[i] = HEAPDATA[endpos + i];

	HEAPDATA = gcnew array<unsigned char>(prepend->Length + append->Length);
	for (int i = 0; i < prepend->Length; i++)
		HEAPDATA[i] = prepend[i];
	for (int i = prepend->Length; i < append->Length + prepend->Length; i++)
		HEAPDATA[i] = prepend[i - prepend->Length];
}

String^ JPFITS::FITSBinTable::GetTTypeEntryRow(String^ ttypeEntry, int rowindex)
{
	int ttypeindex = -1;
	for (int i = 0; i < TTYPES->Length; i++)
		if (TTYPES[i] == ttypeEntry)
		{
			ttypeindex = i;
			break;
		}

	if (ttypeindex == -1)
	{
		throw gcnew Exception("Extension Entry TTYPE Label wasn't found: '" + ttypeEntry + "'");
		return nullptr;
	}

	if (!TTYPEISHEAPARRAYDESC[ttypeindex])
	{
		int objectArrayRank;
		if (TREPEATS[ttypeindex] == 1)
			objectArrayRank = 1;
		else
			objectArrayRank = 2;

		int byteoffset = 0;
		for (int i = 0; i < ttypeindex; i++)
			byteoffset += TBYTES[i];
		int currentbyte = byteoffset + rowindex * NAXIS1;
		String^ str = "";

		switch (TCODES[ttypeindex])
		{
			case ::TypeCode::Double:
			{
				array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
				if (objectArrayRank == 1)
				{
					dbl[7] = BINTABLE[currentbyte];
					dbl[6] = BINTABLE[currentbyte + 1];
					dbl[5] = BINTABLE[currentbyte + 2];
					dbl[4] = BINTABLE[currentbyte + 3];
					dbl[3] = BINTABLE[currentbyte + 4];
					dbl[2] = BINTABLE[currentbyte + 5];
					dbl[1] = BINTABLE[currentbyte + 6];
					dbl[0] = BINTABLE[currentbyte + 7];
					return BitConverter::ToDouble(dbl, 0).ToString();
				}
				else
				{
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + rowindex * NAXIS1 + j * 8;
						dbl[7] = BINTABLE[currentbyte];
						dbl[6] = BINTABLE[currentbyte + 1];
						dbl[5] = BINTABLE[currentbyte + 2];
						dbl[4] = BINTABLE[currentbyte + 3];
						dbl[3] = BINTABLE[currentbyte + 4];
						dbl[2] = BINTABLE[currentbyte + 5];
						dbl[1] = BINTABLE[currentbyte + 6];
						dbl[0] = BINTABLE[currentbyte + 7];
						str += BitConverter::ToDouble(dbl, 0).ToString() + "; ";
					}
					return str;
				}
				break;
			}

			case ::TypeCode::Single:
			{
				array<unsigned char>^ sng = gcnew array<unsigned char>(4);
				if (objectArrayRank == 1)
				{
					sng[3] = BINTABLE[currentbyte];
					sng[2] = BINTABLE[currentbyte + 1];
					sng[1] = BINTABLE[currentbyte + 2];
					sng[0] = BINTABLE[currentbyte + 3];
					return BitConverter::ToSingle(sng, 0).ToString();
				}
				else
				{
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + rowindex * NAXIS1 + j * 4;
						sng[3] = BINTABLE[currentbyte];
						sng[2] = BINTABLE[currentbyte + 1];
						sng[1] = BINTABLE[currentbyte + 2];
						sng[0] = BINTABLE[currentbyte + 3];
						str += BitConverter::ToSingle(sng, 0).ToString() + "; ";
					}
					return str;
				}
				break;
			}

			case (::TypeCode::Int64):
			{
				array<unsigned char>^ i64 = gcnew array<unsigned char>(8);
				if (objectArrayRank == 1)
				{
					i64[7] = BINTABLE[currentbyte];
					i64[6] = BINTABLE[currentbyte + 1];
					i64[5] = BINTABLE[currentbyte + 2];
					i64[4] = BINTABLE[currentbyte + 3];
					i64[3] = BINTABLE[currentbyte + 4];
					i64[2] = BINTABLE[currentbyte + 5];
					i64[1] = BINTABLE[currentbyte + 6];
					i64[0] = BINTABLE[currentbyte + 7];
					return BitConverter::ToInt64(i64, 0).ToString();
				}
				else
				{
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + rowindex * NAXIS1 + j * 8;
						i64[7] = BINTABLE[currentbyte];
						i64[6] = BINTABLE[currentbyte + 1];
						i64[5] = BINTABLE[currentbyte + 2];
						i64[4] = BINTABLE[currentbyte + 3];
						i64[3] = BINTABLE[currentbyte + 4];
						i64[2] = BINTABLE[currentbyte + 5];
						i64[1] = BINTABLE[currentbyte + 6];
						i64[0] = BINTABLE[currentbyte + 7];
						str += BitConverter::ToInt64(i64, 0).ToString() + "; ";
					}
					return str;
				}
				break;
			}

			case (::TypeCode::UInt64):
			{
				unsigned __int64 bzero = 9223372036854775808;
				array<unsigned char>^ ui64 = gcnew array<unsigned char>(8);
				if (objectArrayRank == 1)
				{
					ui64[7] = BINTABLE[currentbyte];
					ui64[6] = BINTABLE[currentbyte + 1];
					ui64[5] = BINTABLE[currentbyte + 2];
					ui64[4] = BINTABLE[currentbyte + 3];
					ui64[3] = BINTABLE[currentbyte + 4];
					ui64[2] = BINTABLE[currentbyte + 5];
					ui64[1] = BINTABLE[currentbyte + 6];
					ui64[0] = BINTABLE[currentbyte + 7];
					return (BitConverter::ToInt64(ui64, 0) + bzero).ToString();
				}
				else
				{
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + rowindex * NAXIS1 + j * 8;
						ui64[7] = BINTABLE[currentbyte];
						ui64[6] = BINTABLE[currentbyte + 1];
						ui64[5] = BINTABLE[currentbyte + 2];
						ui64[4] = BINTABLE[currentbyte + 3];
						ui64[3] = BINTABLE[currentbyte + 4];
						ui64[2] = BINTABLE[currentbyte + 5];
						ui64[1] = BINTABLE[currentbyte + 6];
						ui64[0] = BINTABLE[currentbyte + 7];
						str += (BitConverter::ToInt64(ui64, 0) + bzero).ToString() + "; ";
					}
					return str;
				}
				break;
			}

			case ::TypeCode::UInt32:
			{
				unsigned __int32 bzero = 2147483648;
				array<unsigned char>^ uint32 = gcnew array<unsigned char>(4);
				if (objectArrayRank == 1)
				{
					uint32[3] = BINTABLE[currentbyte];
					uint32[2] = BINTABLE[currentbyte + 1];
					uint32[1] = BINTABLE[currentbyte + 2];
					uint32[0] = BINTABLE[currentbyte + 3];
					return (BitConverter::ToInt32(uint32, 0) + bzero).ToString();
				}
				else
				{
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + rowindex * NAXIS1 + j * 4;
						uint32[3] = BINTABLE[currentbyte];
						uint32[2] = BINTABLE[currentbyte + 1];
						uint32[1] = BINTABLE[currentbyte + 2];
						uint32[0] = BINTABLE[currentbyte + 3];
						str += (BitConverter::ToInt32(uint32, 0) + bzero).ToString() + "; ";
					}
					return str;
				}
				break;
			}

			case ::TypeCode::Int32:
			{
				array<unsigned char>^ int32 = gcnew array<unsigned char>(4);
				if (objectArrayRank == 1)
				{
					int32[3] = BINTABLE[currentbyte];
					int32[2] = BINTABLE[currentbyte + 1];
					int32[1] = BINTABLE[currentbyte + 2];
					int32[0] = BINTABLE[currentbyte + 3];
					return BitConverter::ToInt32(int32, 0).ToString();
				}
				else
				{

					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + rowindex * NAXIS1 + j * 4;
						int32[3] = BINTABLE[currentbyte];
						int32[2] = BINTABLE[currentbyte + 1];
						int32[1] = BINTABLE[currentbyte + 2];
						int32[0] = BINTABLE[currentbyte + 3];
						str += BitConverter::ToInt32(int32, 0).ToString() + "; ";
					}
					return str;
				}
				break;
			}

			case ::TypeCode::UInt16:
			{
				unsigned __int16 bzero = 32768;
				array<unsigned char>^ uint16 = gcnew array<unsigned char>(2);
				if (objectArrayRank == 1)
				{
					uint16[1] = BINTABLE[currentbyte];
					uint16[0] = BINTABLE[currentbyte + 1];
					return (BitConverter::ToInt16(uint16, 0) + bzero).ToString();
				}
				else
				{
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + rowindex * NAXIS1 + j * 2;
						uint16[1] = BINTABLE[currentbyte];
						uint16[0] = BINTABLE[currentbyte + 1];
						str += (BitConverter::ToInt16(uint16, 0) + bzero).ToString() + "; ";
					}
					return str;
				}
				break;
			}

			case ::TypeCode::Int16:
			{
				array<unsigned char>^ int16 = gcnew array<unsigned char>(2);
				if (objectArrayRank == 1)
				{

					int16[1] = BINTABLE[currentbyte];
					int16[0] = BINTABLE[currentbyte + 1];
					return BitConverter::ToInt16(int16, 0).ToString();
				}
				else
				{
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + rowindex * NAXIS1 + j * 2;
						int16[1] = BINTABLE[currentbyte];
						int16[0] = BINTABLE[currentbyte + 1];
						str += BitConverter::ToInt16(int16, 0).ToString() + "; ";
					}
					return str;
				}
				break;
			}

			case ::TypeCode::Byte:
			{
				if (objectArrayRank == 1)
				{
					unsigned __int8 ret = (unsigned __int8)BINTABLE[currentbyte];
					return ret.ToString();
				}
				else
				{
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + rowindex * NAXIS1 + j;
						unsigned __int8 ret = (unsigned __int8)BINTABLE[currentbyte];
						str += ret.ToString() + "; ";
					}
					return str;
				}
				break;
			}

			case ::TypeCode::SByte:
			{
				if (objectArrayRank == 1)
					return ((__int8)(BINTABLE[currentbyte])).ToString();
				else
				{
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + rowindex * NAXIS1 + j;
						str += ((__int8)(BINTABLE[currentbyte])).ToString() + "; ";
					}
					return str;
				}
				break;
			}

			case ::TypeCode::Boolean:
			{
				if (objectArrayRank == 1)
					return Convert::ToBoolean(BINTABLE[currentbyte]).ToString();
				else
				{
					for (int j = 0; j < TREPEATS[ttypeindex]; j++)
					{
						currentbyte = byteoffset + rowindex * NAXIS1 + j;
						str += Convert::ToBoolean(BINTABLE[currentbyte]).ToString() + "; ";
					}
					return str;
				}
				break;
			}

			case TypeCode::Char:
			{
				return System::Text::Encoding::ASCII->GetString(BINTABLE, currentbyte, TREPEATS[ttypeindex]);
			}

			default:
			{
				throw gcnew Exception("Unrecognized TypeCode: '" + TCODES[ttypeindex].ToString() + "'");
				return nullptr;
			}
		}
	}
	else
	{
		int currentbyte = TTYPEHEAPARRAYNELSPOS[ttypeindex][1, rowindex];
		String^ str = "";

		switch (HEAPTCODES[ttypeindex])
		{
			case ::TypeCode::Double:
			{
				array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
				for (int j = 0; j < TTYPEHEAPARRAYNELSPOS[ttypeindex][0, rowindex]; j++)
				{					
					dbl[7] = HEAPDATA[currentbyte];
					dbl[6] = HEAPDATA[currentbyte + 1];
					dbl[5] = HEAPDATA[currentbyte + 2];
					dbl[4] = HEAPDATA[currentbyte + 3];
					dbl[3] = HEAPDATA[currentbyte + 4];
					dbl[2] = HEAPDATA[currentbyte + 5];
					dbl[1] = HEAPDATA[currentbyte + 6];
					dbl[0] = HEAPDATA[currentbyte + 7];
					str += BitConverter::ToDouble(dbl, 0).ToString() + "; ";
					currentbyte += 8;
				}
				return str;
			}

			case ::TypeCode::Single:
			{
				array<unsigned char>^ sng = gcnew array<unsigned char>(4);
				for (int j = 0; j < TTYPEHEAPARRAYNELSPOS[ttypeindex][0, rowindex]; j++)
				{					
					sng[3] = HEAPDATA[currentbyte];
					sng[2] = HEAPDATA[currentbyte + 1];
					sng[1] = HEAPDATA[currentbyte + 2];
					sng[0] = HEAPDATA[currentbyte + 3];
					str += BitConverter::ToSingle(sng, 0).ToString() + "; ";
					currentbyte += 4;
				}
				return str;
			}

			case (::TypeCode::Int64):
			{
				array<unsigned char>^ i64 = gcnew array<unsigned char>(8);
				for (int j = 0; j < TTYPEHEAPARRAYNELSPOS[ttypeindex][0, rowindex]; j++)
				{					
					i64[7] = HEAPDATA[currentbyte];
					i64[6] = HEAPDATA[currentbyte + 1];
					i64[5] = HEAPDATA[currentbyte + 2];
					i64[4] = HEAPDATA[currentbyte + 3];
					i64[3] = HEAPDATA[currentbyte + 4];
					i64[2] = HEAPDATA[currentbyte + 5];
					i64[1] = HEAPDATA[currentbyte + 6];
					i64[0] = HEAPDATA[currentbyte + 7];
					str += BitConverter::ToInt64(i64, 0).ToString() + "; ";
					currentbyte += 8;
				}				
				return str;
			}

			case (::TypeCode::UInt64):
			{
				unsigned __int64 bzero = 9223372036854775808;
				array<unsigned char>^ ui64 = gcnew array<unsigned char>(8);
				for (int j = 0; j < TTYPEHEAPARRAYNELSPOS[ttypeindex][0, rowindex]; j++)
				{
					ui64[7] = HEAPDATA[currentbyte];
					ui64[6] = HEAPDATA[currentbyte + 1];
					ui64[5] = HEAPDATA[currentbyte + 2];
					ui64[4] = HEAPDATA[currentbyte + 3];
					ui64[3] = HEAPDATA[currentbyte + 4];
					ui64[2] = HEAPDATA[currentbyte + 5];
					ui64[1] = HEAPDATA[currentbyte + 6];
					ui64[0] = HEAPDATA[currentbyte + 7];
					str += (BitConverter::ToInt64(ui64, 0) + bzero).ToString() + "; ";
					currentbyte += 8;
				}
				return str;
			}

			case ::TypeCode::Int32:
			{
				array<unsigned char>^ int32 = gcnew array<unsigned char>(4);
				for (int j = 0; j < TTYPEHEAPARRAYNELSPOS[ttypeindex][0, rowindex]; j++)
				{					
					int32[3] = HEAPDATA[currentbyte];
					int32[2] = HEAPDATA[currentbyte + 1];
					int32[1] = HEAPDATA[currentbyte + 2];
					int32[0] = HEAPDATA[currentbyte + 3];
					str += BitConverter::ToInt32(int32, 0).ToString() + "; ";
					currentbyte += 4;
				}
				return str;
			}

			case ::TypeCode::UInt32:
			{
				unsigned __int32 bzero = 2147483648;
				array<unsigned char>^ uint32 = gcnew array<unsigned char>(4);
				for (int j = 0; j < TTYPEHEAPARRAYNELSPOS[ttypeindex][0, rowindex]; j++)
				{					
					uint32[3] = HEAPDATA[currentbyte];
					uint32[2] = HEAPDATA[currentbyte + 1];
					uint32[1] = HEAPDATA[currentbyte + 2];
					uint32[0] = HEAPDATA[currentbyte + 3];
					str += (BitConverter::ToInt32(uint32, 0) + bzero).ToString() + "; ";
					currentbyte += 4;
				}
				return str;
			}

			case ::TypeCode::Int16:
			{
				array<unsigned char>^ int16 = gcnew array<unsigned char>(2);
				for (int j = 0; j < TTYPEHEAPARRAYNELSPOS[ttypeindex][0, rowindex]; j++)
				{					
					int16[1] = HEAPDATA[currentbyte];
					int16[0] = HEAPDATA[currentbyte + 1];
					str += BitConverter::ToInt16(int16, 0).ToString() + "; ";
					currentbyte += 2;
				}
				return str;
			}

			case ::TypeCode::UInt16:
			{
				unsigned __int16 bzero = 32768;
				array<unsigned char>^ uint16 = gcnew array<unsigned char>(2);
				for (int j = 0; j < TTYPEHEAPARRAYNELSPOS[ttypeindex][0, rowindex]; j++)
				{
					uint16[1] = HEAPDATA[currentbyte];
					uint16[0] = HEAPDATA[currentbyte + 1];
					str += (BitConverter::ToInt16(uint16, 0) + bzero).ToString() + "; ";
					currentbyte += 2;
				}
				return str;
			}

			case ::TypeCode::SByte:
			{
				for (int j = 0; j < TTYPEHEAPARRAYNELSPOS[ttypeindex][0, rowindex]; j++)
					str += ((__int8)(HEAPDATA[currentbyte + j])).ToString() + "; ";
				return str;
			}

			case ::TypeCode::Byte:
			{
				for (int j = 0; j < TTYPEHEAPARRAYNELSPOS[ttypeindex][0, rowindex]; j++)
					str += ((unsigned __int8)HEAPDATA[currentbyte + j]).ToString() + "; ";
				return str;
			}
			
			case ::TypeCode::Boolean:
			{
				for (int j = 0; j < TTYPEHEAPARRAYNELSPOS[ttypeindex][0, rowindex]; j++)
					str += Convert::ToBoolean(HEAPDATA[currentbyte + j]).ToString() + "; ";
				return str;
			}

			case TypeCode::Char:
			{
				return System::Text::Encoding::ASCII->GetString(HEAPDATA, TTYPEHEAPARRAYNELSPOS[ttypeindex][1, rowindex], TTYPEHEAPARRAYNELSPOS[ttypeindex][0, rowindex]);
			}

			default:
			{
				throw gcnew Exception("Unrecognized TypeCode: '" + HEAPTCODES[ttypeindex].ToString() + "'");
				return nullptr;
			}
		}
	}
}

void JPFITS::FITSBinTable::RemoveTTYPEEntry(String^ ttypeEntry)
{
	int ttypeindex = -1;
	for (int i = 0; i < TTYPES->Length; i++)
		if (TTYPES[i] == ttypeEntry)
		{
			ttypeindex = i;
			break;
		}

	if (ttypeindex == -1)
	{
		throw gcnew Exception("Extension Entry TTYPE wasn't found: '" + ttypeEntry + "'");
		return;
	}

	array<Object^>^ newEntryDataObjs = gcnew array<Object^>(TFIELDS - 1);
	array<String^>^ newTTYPES = gcnew array<String^>(TFIELDS - 1);
	array<String^>^ newTFORMS = gcnew array<String^>(TFIELDS - 1);
	array<String^>^ newTUNITS = gcnew array<String^>(TFIELDS - 1);
	array<int>^ newTBYTES = gcnew array<int>(TFIELDS - 1);
	array<int>^ newTREPEATS = gcnew array<int>(TFIELDS - 1);
	array<TypeCode>^ newTCODES = gcnew array<TypeCode>(TFIELDS - 1);
	array<array<int>^>^ newTDIMS = gcnew array<array<int>^>(TFIELDS - 1);
	array<bool>^ newTTYPEISCOMPLEX = gcnew array<bool>(TFIELDS - 1);
	array<bool>^ newTTYPEISHEAPARRAYDESC = gcnew array<bool>(TFIELDS - 1);
	array<TypeCode>^ newHEAPTCODES = gcnew array<TypeCode>(TFIELDS - 1);
	array<array<int, 2>^>^ newTTYPEHEAPARRAYNELSPOS = gcnew array<array<int, 2>^>(TFIELDS - 1);

	int c = 0;
	for (int i = 0; i < TFIELDS; i++)
		if (i == ttypeindex)
		{
			if (TTYPEISHEAPARRAYDESC[ttypeindex])
				REMOVEHEAPTTYPE(ttypeindex);
		}
		else
		{
			TypeCode code;
			array<int>^ dimnelements;
			newEntryDataObjs[c] = this->GetTTYPEEntry(TTYPES[i], code, dimnelements);
			newTTYPES[c] = TTYPES[i];
			newTFORMS[c] = TFORMS[i];
			newTUNITS[c] = TUNITS[i];
			newTBYTES[c] = TBYTES[i];
			newTREPEATS[c] = TREPEATS[i];
			newTCODES[c] = TCODES[i];
			newTDIMS[c] = TDIMS[i];
			newTTYPEISCOMPLEX[c] = TTYPEISCOMPLEX[i];
			newTTYPEISHEAPARRAYDESC[c] = TTYPEISHEAPARRAYDESC[i];
			newHEAPTCODES[c] = HEAPTCODES[i];
			newTTYPEHEAPARRAYNELSPOS[c] = TTYPEHEAPARRAYNELSPOS[i];
			c++;
		}

	TFIELDS--;
	TTYPES = newTTYPES;
	TFORMS = newTFORMS;
	TUNITS = newTUNITS;
	TBYTES = newTBYTES;
	TREPEATS = newTREPEATS;
	TCODES = newTCODES;
	TDIMS = newTDIMS;
	TTYPEISCOMPLEX = newTTYPEISCOMPLEX;
	TTYPEISHEAPARRAYDESC = newTTYPEISHEAPARRAYDESC;
	HEAPTCODES = newHEAPTCODES;
	TTYPEHEAPARRAYNELSPOS = newTTYPEHEAPARRAYNELSPOS;

	MAKEBINTABLEBYTEARRAY(newEntryDataObjs);
}

void JPFITS::FITSBinTable::SetTTYPEEntries(array<String^>^ ttypeEntries, array<String^>^ entryUnits, array<Object^>^ entryArrays)
{
	for (int i = 0; i < entryArrays->Length; i++)
		if (((Array^)entryArrays[i])->Rank > 2)
		{
			throw gcnew Exception("Error: Do not use this function to add an n &gt; 2 dimensional array. Use AddTTYPEEntry.");
			return;
		}
		else if (((Array^)entryArrays[i])->Rank == 1 && Type::GetTypeCode(entryArrays[i]->GetType()->GetElementType()) == TypeCode::String)
		{
			array<String^>^ strarr = (array<String^>^)entryArrays[i];
			int nels = strarr[0]->Length;
			for (int j = 1; j < strarr->Length; j++)
				if (strarr[j]->Length != nels)
				{
					throw gcnew Exception("Error: String array entries '" + ttypeEntries[i] + "' are not all the same namber of characters (repeats) long.");
					return;
				}
		}

	bool equalnaxis2 = true;
	bool stringarrayrank2 = false;
	int naxis2;
	if (((Array^)entryArrays[0])->Rank == 1)
		naxis2 = ((Array^)entryArrays[0])->Length;
	else
		naxis2 = ((Array^)entryArrays[0])->GetLength(1);
	for (int i = 1; i < entryArrays->Length; i++)
		if (((Array^)entryArrays[i])->Rank == 1)
		{
			if (((Array^)entryArrays[i])->Length != naxis2)
			{
				equalnaxis2 = false;
				break;
			}
		}
		else
		{
			if (((Array^)entryArrays[i])->GetLength(1) != naxis2)
			{
				equalnaxis2 = false;
				break;
			}

			if (Type::GetTypeCode(entryArrays[i]->GetType()->GetElementType()) == TypeCode::String)
			{
				stringarrayrank2 = true;
				break;
			}
		}
	if (!equalnaxis2)
	{
		throw gcnew Exception("Error: all entry column heights, NAXIS2s, are not equal. Use an overload to add a variable length array to the heap.");
		return;
	}
	if (stringarrayrank2)
	{
		throw gcnew Exception("Error: Cannot pass a 2d String array. Only a 1D array of Strings is allowed.");
		return;
	}

	TFIELDS = entryArrays->Length;
	TTYPES = ttypeEntries;
	if (entryUnits != nullptr)
		TUNITS = entryUnits;
	else
		TUNITS = gcnew array<String^>(entryArrays->Length);
	TCODES = gcnew array<TypeCode>(entryArrays->Length);
	TREPEATS = gcnew array<int>(entryArrays->Length);
	TFORMS = gcnew array<String^>(entryArrays->Length);
	TBYTES = gcnew array<int>(entryArrays->Length);
	TTYPEISCOMPLEX = gcnew array<bool>(entryArrays->Length);
	TTYPEISHEAPARRAYDESC = gcnew array<bool>(entryArrays->Length);
	TDIMS = gcnew array<array<int>^>(entryArrays->Length);
	HEAPTCODES = gcnew array<TypeCode>(entryArrays->Length);
	TTYPEHEAPARRAYNELSPOS = gcnew array<array<int, 2>^>(entryArrays->Length);

	for (int i = 0; i < entryArrays->Length; i++)
	{
		TCODES[i] = Type::GetTypeCode((((Array^)entryArrays[i])->GetType())->GetElementType());
		if (TCODES[i] != TypeCode::String)
			if (((Array^)entryArrays[i])->Rank == 1)
				TREPEATS[i] = 1;
			else
				TREPEATS[i] = ((Array^)entryArrays[i])->GetLength(0);
		else
		{
			TCODES[i] = TypeCode::Char;
			TREPEATS[i] = ((array<String^>^)entryArrays[i])[0]->Length;
		}
		TFORMS[i] = TREPEATS[i].ToString() + TYPECODETFORM(TCODES[i]);
		TBYTES[i] = TYPECODETONBYTES(TCODES[i]) * TREPEATS[i];
	}	

	//new table, so these either need set for the first time, or updated
	BITPIX = 8;
	NAXIS = 2;
	NAXIS1 = 0;
	for (int i = 0; i < entryArrays->Length; i++)
		NAXIS1 += TBYTES[i];
	if (((Array^)entryArrays[0])->Rank == 1)
		NAXIS2 = ((Array^)entryArrays[0])->Length;
	else
		NAXIS2 = ((Array^)entryArrays[0])->GetLength(1);

	MAKEBINTABLEBYTEARRAY(entryArrays);
	HEAPDATA = nullptr;
}

void JPFITS::FITSBinTable::AddTTYPEEntry(String^ ttypeEntry, bool replaceIfExists, String^ entryUnits, Object^ entryArray)
{
	if (((Array^)entryArray)->Rank > 2)
	{
		throw gcnew Exception("Error: Do not use this function to add an n &gt; 2 dimensional array. Use an overload.");
		return;
	}

	AddTTYPEEntry(ttypeEntry, replaceIfExists, entryUnits, entryArray, nullptr, false, false);
}

void JPFITS::FITSBinTable::AddTTYPEEntry(String^ ttypeEntry, bool replaceIfExists, String^ entryUnits, Object^ entryArray, array<int>^ dimNElements, bool isComplex, bool addAsHeapVarRepeatArray)
{
	//int heapentryarrayLengthNRowsNAXIS2 = ((Array^)entryArray)->Length;
	//int heapentryarrayindexRowLengthRepeatsNels = ((Array^)(((Array^)entryArray)->GetValue(i)))->Length;//this is the variable repeat count
	//int heapentryarrayindexRowRank = ((Array^)(((Array^)entryArray)->GetValue(i)))->Rank;//all must be 1
	//TypeCode heapentryarrayTypeCode = Type::GetTypeCode((((Array^)entryArray)->GetValue(0))->GetType()->GetElementType());
	//bool isarray = (((Array^)entryArray)->GetValue(0))->GetType()->IsArray;//???????????????????
	//TypeCode heapentryarrayTypeCodeSTRING = Type::GetTypeCode(entryArray->GetType()->GetElementType());
	//return;

	int ttypeindex = -1;
	if (TTYPES != nullptr)
		for (int i = 0; i < TTYPES->Length; i++)
			if (TTYPES[i] == ttypeEntry)
			{
				ttypeindex = i;
				break;
			}

	if (ttypeindex != -1 && !replaceIfExists)
	{
		throw gcnew Exception("Extension Entry TTYPE '" + ttypeEntry + "' already exists, but was told to not overwrite it.");
		return;
	}

	bool isheapvarrepeatString = false;
	if (addAsHeapVarRepeatArray)
		if (Type::GetTypeCode(entryArray->GetType()->GetElementType()) == TypeCode::String)
			isheapvarrepeatString = true;

	if (isComplex && isheapvarrepeatString)
	{
		throw gcnew Exception("Adding a char String TTYPE but told it is numerical complex, which doesn't make sense: '" + ttypeEntry + ".");
		return;
	}

	if (isComplex && !addAsHeapVarRepeatArray)
	{
		if (((Array^)entryArray)->Rank == 1)
		{
			throw gcnew Exception("Extension Entry TTYPE '" + ttypeEntry + "' is supposed to be complex, but has a rank of 1. A complex array must have a rank of at least 2 for spatial and temporal pairings.");
			return;
		}
		if (Type::GetTypeCode((((Array^)entryArray)->GetType())->GetElementType()) != TypeCode::Double && Type::GetTypeCode((((Array^)entryArray)->GetType())->GetElementType()) != TypeCode::Single)
		{
			throw gcnew Exception("Extension Entry TTYPE '" + ttypeEntry + "' may only be single or double precision floating point if complex, but was " + Type::GetTypeCode((((Array^)entryArray)->GetType())->GetElementType()).ToString());
			return;
		}
		if (!JPMath::IsEven(((Array^)entryArray)->GetLength(0)))
		{
			throw gcnew Exception("Extension Entry TTYPE '" + ttypeEntry + "' is supposed to be complex, but is not an even pairing of spatial and temporal columns.");
			return;
		}			
	}
	
	if (isComplex && addAsHeapVarRepeatArray)
	{
		for (int i = 0; i < ((Array^)entryArray)->Length; i++)
			if (!JPMath::IsEven(((Array^)(((Array^)entryArray)->GetValue(i)))->Length))
			{
				throw gcnew Exception("Extension Entry TTYPE '" + ttypeEntry + "' is supposed to be complex, but is not an even pairing of spatial and temporal columns.");
				return;
			}

		if (Type::GetTypeCode((((Array^)entryArray)->GetValue(0))->GetType()->GetElementType()) != TypeCode::Double && Type::GetTypeCode((((Array^)entryArray)->GetValue(0))->GetType()->GetElementType()) != TypeCode::Single)
		{
			throw gcnew Exception("Extension Entry TTYPE '" + ttypeEntry + "' may only be single or double precision floating point if complex, but was " + Type::GetTypeCode((((Array^)entryArray)->GetValue(0))->GetType()->GetElementType()).ToString());
			return;
		}
	}

	if (addAsHeapVarRepeatArray && !isheapvarrepeatString)
		for (int i = 0; i < ((Array^)entryArray)->Length; i++)
			if (((Array^)(((Array^)entryArray)->GetValue(i)))->Rank != 1)
			{
				throw gcnew Exception("Extension Entry TTYPE '" + ttypeEntry + "' must be an array of rank = 1 arrays. Index '" + i.ToString() + "' is rank = " + ((Array^)(((Array^)entryArray)->GetValue(i)))->Rank);
				return;
			}

	if (!addAsHeapVarRepeatArray && Type::GetTypeCode((((Array^)entryArray)->GetType())->GetElementType()) == TypeCode::String)
		if (((Array^)entryArray)->Rank == 2)
		{
			throw gcnew Exception("Error: Cannot pass a 2d String array '" + ttypeEntry + "' . Only a 1D array of Strings is allowed.");
			return;
		}
		else
		{
			array<String^>^ strarr = (array<String^>^)entryArray;
			int nels = strarr[0]->Length;
			for (int j = 1; j < strarr->Length; j++)
				if (strarr[j]->Length != nels)
				{
					throw gcnew Exception("Error: String array entries '" + ttypeEntry + "' are not all the same namber of characters (repeats) long.");
					return;
				}
		}

	if (ttypeindex != -1)//then remove it
		this->RemoveTTYPEEntry(ttypeEntry);
	else
		ttypeindex = TFIELDS;//then put the entry at the last column of the table...NB this is a zero-based index...TFIELDS will increment by one below

	//either it was an add to a blank table, or a replacement, or an additional, so these either need set for the first time, or updated
	if (TFIELDS == 0)
		if (((Array^)entryArray)->Rank == 1)//true for heapentry too as array of arrays
			NAXIS2 = ((Array^)entryArray)->Length;
		else
			NAXIS2 = ((Array^)entryArray)->GetLength(1);
	else
		if (((Array^)entryArray)->Rank == 1 && ((Array^)entryArray)->Length != NAXIS2 || ((Array^)entryArray)->Rank > 1 && ((Array^)entryArray)->GetLength(1) != NAXIS2)
		{
			int naxis2;
			if (((Array^)entryArray)->Rank == 1)//true for heapentry too as array of arrays
				naxis2 = ((Array^)entryArray)->Length;
			else
				naxis2 = ((Array^)entryArray)->GetLength(1);
			throw gcnew Exception("Error: Existing NAXIS2 = " + NAXIS2 + "; new entryArray '" + ttypeEntry + "'  NAXIS2 = " + naxis2 + ".");
			return;
		}

	TFIELDS++;
	array<Object^>^ newEntryDataObjs = gcnew array<Object^>(TFIELDS);
	array<String^>^ newTTYPES = gcnew array<String^>(TFIELDS);
	array<String^>^ newTFORMS = gcnew array<String^>(TFIELDS);
	array<String^>^ newTUNITS = gcnew array<String^>(TFIELDS);
	array<int>^ newTBYTES = gcnew array<int>(TFIELDS);
	array<int>^ newTREPEATS = gcnew array<int>(TFIELDS);
	array<TypeCode>^ newTCODES = gcnew array<TypeCode>(TFIELDS);
	array<array<int>^>^ newTDIMS = gcnew array<array<int>^>(TFIELDS);
	array<bool>^ newTTYPEISCOMPLEX = gcnew array<bool>(TFIELDS);
	array<bool>^ newTTYPEISHEAPARRAYDESC = gcnew array<bool>(TFIELDS);
	array<TypeCode>^ newHEAPTCODES = gcnew array<TypeCode>(TFIELDS);
	array<array<int, 2>^>^ newTTYPEHEAPARRAYNELSPOS = gcnew array<array<int, 2>^>(TFIELDS);

	int c = 0;
	for (int i = 0; i < TFIELDS; i++)
		if (i == ttypeindex)
		{
			int instances = 1;
			if (((Array^)entryArray)->Rank > 1)
				instances = ((Array^)entryArray)->GetLength(0);

			if (!addAsHeapVarRepeatArray)
			{
				newTCODES[i] = Type::GetTypeCode((((Array^)entryArray)->GetType())->GetElementType());
				if (newTCODES[i] == TypeCode::String)
				{
					newTCODES[i] = TypeCode::Char;
					instances = ((array<String^>^)entryArray)[0]->Length;
				}
				newTBYTES[i] = TYPECODETONBYTES(newTCODES[i]) * instances;
				if (isComplex)
					if (newTCODES[i] == TypeCode::Double)
						newTFORMS[i] = (instances / 2).ToString() + "M";
					else
						newTFORMS[i] = (instances / 2).ToString() + "C";
				else
					newTFORMS[i] = instances.ToString() + TYPECODETFORM(newTCODES[i]);
				newTREPEATS[i] = instances;
			}
			else
			{
				newTCODES[i] = TypeCode::Int32;
				newTBYTES[i] = TYPECODETONBYTES(newTCODES[i]) * 2;
				if (isheapvarrepeatString)
					newHEAPTCODES[i] = TypeCode::Char;
				else
					newHEAPTCODES[i] = Type::GetTypeCode((((Array^)entryArray)->GetValue(0))->GetType()->GetElementType());
				newTTYPEISHEAPARRAYDESC[i] = addAsHeapVarRepeatArray;
				newTTYPEHEAPARRAYNELSPOS[i] = nullptr;//this gets set in MAKEBINTABLEBYTEARRAY......
				if (isComplex)
					if (newHEAPTCODES[i] == TypeCode::Double)
						newTFORMS[i] = "PM";
					else
						newTFORMS[i] = "PC";
				else
					newTFORMS[i] = "P" + TYPECODETFORM(newHEAPTCODES[i]);
				newTREPEATS[i] = 1;
				newTFORMS[i] = newTREPEATS[i] + newTFORMS[i];
			}

			newEntryDataObjs[i] = entryArray;
			newTTYPES[i] = ttypeEntry;			
			newTUNITS[i] = entryUnits;			
			newTDIMS[i] = dimNElements;
			newTTYPEISCOMPLEX[i] = isComplex;
		}
		else
		{
			TypeCode code;
			array<int>^ dimnelements;
			newEntryDataObjs[i] = this->GetTTYPEEntry(TTYPES[c], code, dimnelements);
			newTTYPES[i] = TTYPES[c];
			newTFORMS[i] = TFORMS[c];
			newTUNITS[i] = TUNITS[c];
			newTBYTES[i] = TBYTES[c];
			newTREPEATS[i] = TREPEATS[c];
			newTCODES[i] = TCODES[c];
			newTDIMS[i] = TDIMS[c];
			newTTYPEISCOMPLEX[i] = TTYPEISCOMPLEX[c];
			newTTYPEISHEAPARRAYDESC[i] = TTYPEISHEAPARRAYDESC[c];
			newHEAPTCODES[i] = HEAPTCODES[c];
			newTTYPEHEAPARRAYNELSPOS[i] = TTYPEHEAPARRAYNELSPOS[c];
			c++;
		}

	TTYPES = newTTYPES;
	TFORMS = newTFORMS;
	TUNITS = newTUNITS;
	TBYTES = newTBYTES;
	TREPEATS = newTREPEATS;
	TCODES = newTCODES;
	TDIMS = newTDIMS;
	TTYPEISCOMPLEX = newTTYPEISCOMPLEX;
	TTYPEISHEAPARRAYDESC = newTTYPEISHEAPARRAYDESC;
	HEAPTCODES = newHEAPTCODES;
	TTYPEHEAPARRAYNELSPOS = newTTYPEHEAPARRAYNELSPOS;

	//either it was an add to a blank table, or a replacement, or an additional, so these either need set for the first time, or updated
	BITPIX = 8;
	NAXIS = 2;
	NAXIS1 = 0;
	for (int i = 0; i < TBYTES->Length; i++)
		NAXIS1 += TBYTES[i];

	MAKEBINTABLEBYTEARRAY(newEntryDataObjs);
}

void JPFITS::FITSBinTable::AddExtraHeaderKey(String^ keyName, String^ keyValue, String^ keyComment)
{
	if (EXTRAKEYS == nullptr)
	{
		EXTRAKEYS = gcnew array<String^>(1) { keyName };
		EXTRAKEYVALS = gcnew array<String^>(1) { keyValue };
		EXTRAKEYCOMS = gcnew array<String^>(1) { keyComment };
	}
	else
	{
		array<String^>^ newkeys = gcnew array<String^>(EXTRAKEYS->Length + 1);
		array<String^>^ newvals = gcnew array<String^>(EXTRAKEYS->Length + 1);
		array<String^>^ newcoms = gcnew array<String^>(EXTRAKEYS->Length + 1);

		for (int i = 0; i < EXTRAKEYS->Length; i++)
		{
			newkeys[i] = EXTRAKEYS[i];
			newvals[i] = EXTRAKEYVALS[i];
			newcoms[i] = EXTRAKEYCOMS[i];
		}
		newkeys[EXTRAKEYS->Length] = keyName;
		newvals[EXTRAKEYS->Length] = keyValue;
		newcoms[EXTRAKEYS->Length] = keyComment;

		EXTRAKEYS = newkeys;
		EXTRAKEYVALS = newvals;
		EXTRAKEYCOMS = newcoms;
	}
}

String^ JPFITS::FITSBinTable::GetExtraHeaderKeyValue(String^ keyName)
{
	for (int i = 0; i < EXTRAKEYS->Length; i++)
		if (keyName->Equals(EXTRAKEYS[i]))
			return EXTRAKEYVALS[i];
	return "";
}

void JPFITS::FITSBinTable::RemoveExtraHeaderKey(String^ keyName, String^ keyValue)
{
	if (EXTRAKEYS == nullptr)
		return;
	
	int keyindex = -1;

	for (int i = 0; i < EXTRAKEYS->Length; i++)
		if (EXTRAKEYS[i] == keyName && EXTRAKEYVALS[i] == keyValue)
		{
			keyindex = i;
			break;
		}

	if (keyindex == -1)
		return;

	array<String^>^ newkeys = gcnew array<String^>(EXTRAKEYS->Length - 1);
	array<String^>^ newvals = gcnew array<String^>(EXTRAKEYS->Length - 1);
	array<String^>^ newcoms = gcnew array<String^>(EXTRAKEYS->Length - 1);
	int c = 0;
	for (int i = 0; i < EXTRAKEYS->Length; i++)
		if (i == keyindex)
			continue;
		else
		{
			newkeys[c] = EXTRAKEYS[i];
			newvals[c] = EXTRAKEYVALS[i];
			newcoms[c] = EXTRAKEYCOMS[i];
			c++;
		}
	EXTRAKEYS = newkeys;
	EXTRAKEYVALS = newvals;
	EXTRAKEYCOMS = newcoms;
}

void JPFITS::FITSBinTable::RemoveAllExtraHeaderKeys()
{
	EXTRAKEYS = nullptr;
	EXTRAKEYVALS = nullptr;
	EXTRAKEYCOMS = nullptr;
}

void JPFITS::FITSBinTable::RemoveExtension(String^ FileName, String^ ExtensionName)
{
	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	bool hasext = false;
	if (!FITSFILEOPS::SCANPRIMARYUNIT(fs, true, nullptr, hasext) || !hasext)
	{
		fs->Close();
		if (!hasext)
			throw gcnew Exception("File indicates no extensions present.");
		else
			throw gcnew Exception("File not formatted as FITS file.");
		return;
	}

	__int64 extensionstartposition, extensionendposition, tableendposition, pcount, theap;
	bool exists = FITSFILEOPS::SEEKEXTENSION(fs, "BINTABLE", ExtensionName, nullptr, extensionstartposition, extensionendposition, tableendposition, pcount, theap);
	if (!exists)
	{
		fs->Close();
		throw gcnew Exception("Could not find BINTABLE with name '" + ExtensionName + "'");
		return;
	}

	//if here then we found the extension and can excise it given its start and end position
	array<unsigned char>^ arrstart = gcnew array<unsigned char>((int)extensionstartposition);
	fs->Position = 0;
	fs->Read(arrstart, 0, arrstart->Length);

	array<unsigned char>^ arrend = gcnew array<unsigned char>(int(fs->Length - extensionendposition));
	fs->Position = extensionendposition;
	fs->Read(arrend, 0, arrend->Length);
	fs->Close();

	fs = gcnew FileStream(FileName, IO::FileMode::Create);
	fs->Write(arrstart, 0, arrstart->Length);
	fs->Write(arrend, 0, arrend->Length);
	fs->Close();
}

bool JPFITS::FITSBinTable::ExtensionExists(String^ FileName, String^ ExtensionName)
{
	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	bool hasext = false;
	if (!FITSFILEOPS::SCANPRIMARYUNIT(fs, true, nullptr, hasext) || !hasext)
	{
		fs->Close();
		if (!hasext)
			throw gcnew Exception("File indicates no extensions present.");
		else
			throw gcnew Exception("File not formatted as FITS file.");
		return false;
	}

	__int64 extensionstartposition, extensionendposition, tableendposition, pcount, theap;
	bool exists = FITSFILEOPS::SEEKEXTENSION(fs, "BINTABLE", ExtensionName, nullptr, extensionstartposition, extensionendposition, tableendposition, pcount, theap);
	fs->Close();
	return exists;
}

void JPFITS::FITSBinTable::Write(String^ FileName, String^ ExtensionName, bool OverWriteExtensionIfExists)
{
	EXTENSIONNAME = ExtensionName;
	FILENAME = FileName;

	if (!File::Exists(FILENAME))//then write a new file, otherwise check the existing file for existing table, etc.
	{
		JPFITS::FITSImage^ ff = gcnew FITSImage(FILENAME, true);
		ff->WriteImage(TypeCode::Double, true);
	}

	FileStream^ fs = gcnew FileStream(FILENAME, IO::FileMode::Open);

	bool hasext = false;
	if (!FITSFILEOPS::SCANPRIMARYUNIT(fs, true, nullptr, hasext))
	{
		fs->Close();
		throw gcnew Exception("File not formatted as FITS file. Use a new file.");
		return;
	}
	if (!hasext)
	{
		fs->Position = 0;
		FITSFILEOPS::SCANPRIMARYUNIT(fs, false, nullptr, hasext);
		array<unsigned char>^ primarydataarr = gcnew array<unsigned char>((int)(fs->Length - fs->Position));
		fs->Read(primarydataarr, 0, primarydataarr->Length);
		fs->Close();

		FITSImage^ ff = gcnew FITSImage(FILENAME, nullptr, true, false, false, false);
		int n = ff->Header->GetKeyIndex("NAXIS", false);
		if (n == -1)
		{
			throw gcnew Exception("File not formatted as FITS file (NAXIS not present). Use a new file.");
			return;
		}
		n = Convert::ToInt32(ff->Header->GetKeyValue("NAXIS"));
		if (n > 0)
		{
			n = ff->Header->GetKeyIndex("NAXIS" + n.ToString(), false);
			if (ff->Header->GetKeyIndex("BZERO", false) > n)
				n = ff->Header->GetKeyIndex("BZERO", false);
			if (ff->Header->GetKeyIndex("BSCALE", false) > n)
				n = ff->Header->GetKeyIndex("BSCALE", false);
		}
		else
			n = ff->Header->GetKeyIndex("NAXIS", false);
		ff->Header->SetKey("EXTEND", "T", "FITS file may contain extensions", true, n + 1);
		array<String^>^ HEADER = ff->Header->GetFormattedHeaderBlock(true, false);           // FITSHEADER::MAKEFORMATTEDIMAGEHEADER(ff->Header->HeaderKeys, ff->Header->HeaderKeyValues, ff->Header->HeaderKeyComments, false);

		array<unsigned char>^ headarr = gcnew array<unsigned char>(HEADER->Length * 80);
		for (int i = 0; i < HEADER->Length; i++)
			for (int j = 0; j < 80; j++)
				headarr[i * 80 + j] = (unsigned char)HEADER[i][j];

		fs = gcnew FileStream(FILENAME, IO::FileMode::Create);
		fs->Write(headarr, 0, headarr->Length);
		fs->Write(primarydataarr, 0, primarydataarr->Length);
		fs->Close();

		fs = gcnew FileStream(FILENAME, IO::FileMode::Open);
		FITSFILEOPS::SCANPRIMARYUNIT(fs, true, nullptr, hasext);
	}

	__int64 extensionstartposition, extensionendposition, tableendposition, pcount, theap;
	bool extensionfound = FITSFILEOPS::SEEKEXTENSION(fs, "BINTABLE", EXTENSIONNAME, nullptr, extensionstartposition, extensionendposition, tableendposition, pcount, theap);
	if (extensionfound && !OverWriteExtensionIfExists)
	{
		fs->Close();
		throw gcnew Exception("ExtensionName '" + EXTENSIONNAME + "' already exists and was told to not overwrite it...");
		return;
	}

	array<unsigned char>^ arr_prepend;
	array<unsigned char>^ arr_append;
	if (extensionfound)
	{
		arr_prepend = gcnew array<unsigned char>((int)extensionstartposition);
		fs->Position = 0;
		fs->Read(arr_prepend, 0, arr_prepend->Length);

		if (extensionendposition != fs->Length)//then this was not the end of the file...get the appendage data
		{
			fs->Position = extensionendposition;
			arr_append = gcnew array<unsigned char>(int(fs->Length - extensionendposition));
			fs->Read(arr_append, 0, arr_append->Length);
		}
		fs->Position = extensionstartposition;
	}
	else
	{
		arr_prepend = gcnew array<unsigned char>((int)fs->Length);
		fs->Position = 0;
		fs->Read(arr_prepend, 0, arr_prepend->Length);
	}
	fs->Close();

	fs = gcnew FileStream(FILENAME, IO::FileMode::Create);
	if (arr_prepend != nullptr)
		fs->Write(arr_prepend, 0, arr_prepend->Length);

	//format the header for writing
	array<String^>^ header = FORMATBINARYTABLEEXTENSIONHEADER();
	array<unsigned char>^ headerdata = gcnew array<unsigned char>(header->Length * 80);

	for (int i = 0; i < header->Length; i++)
		for (int j = 0; j < 80; j++)
			headerdata[i * 80 + j] = (unsigned char)header[i][j];

	fs->Write(headerdata, 0, headerdata->Length);
	fs->Write(BINTABLE, 0, BINTABLE->Length);

	if (HEAPDATA!= nullptr)
		fs->Write(HEAPDATA, 0, HEAPDATA->Length);
	int Tbytes = BINTABLE->Length;
	if (HEAPDATA != nullptr)
		Tbytes += HEAPDATA->Length;

	int Nfillbytes = int(Math::Ceiling(double(Tbytes) / 2880.0)) * 2880 - Tbytes;
	for (int i = 0; i < Nfillbytes; i++)
		fs->WriteByte(0);
	if (arr_append != nullptr)
		fs->Write(arr_append, 0, arr_append->Length);
	fs->Close();
}

void JPFITS::FITSBinTable::MAKEBINTABLEBYTEARRAY(array<Object^>^ ExtensionEntryData)
{
	for (int i = 0; i < ExtensionEntryData->Length; i++)
		if (!ExtensionEntryData[i]->GetType()->IsArray)
		{
			throw gcnew Exception("Error: Object at index '" + i + "' is not an array. Stopping write.");
			return;
		}

	MAKEHEAPBYTEARRAY(ExtensionEntryData);//will do nothing if there's no heap data

	int TBytes = NAXIS1 * NAXIS2;
	BINTABLE = gcnew array<unsigned char>(TBytes);

	int nthread = omp_get_max_threads();
	bool parallel_top = false;
	if (NAXIS2 >= nthread)
		parallel_top = true;
	bool exception = false;
	TypeCode exceptiontypecode;

	//now write the table data into the array
	#pragma omp parallel for if(parallel_top)
	for (int i = 0; i < NAXIS2; i++)
	{
		if (exception)
			break;

		for (int j = 0; j < ExtensionEntryData->Length; j++)
		{
			int jtps = 0;
			for (int jj = 0; jj < j; jj++)
				jtps += TBYTES[jj];

			int cc = i * NAXIS1 + jtps;

			if (TTYPEISHEAPARRAYDESC[j])
			{
				for (int ii = 0; ii < 2; ii++)
				{
					int nelpos = TTYPEHEAPARRAYNELSPOS[j][ii, i];
					BINTABLE[cc + ii * 4] = ((nelpos >> 24) & 0xff);
					BINTABLE[cc + ii * 4 + 1] = ((nelpos >> 16) & 0xff);
					BINTABLE[cc + ii * 4 + 2] = ((nelpos >> 8) & 0xff);
					BINTABLE[cc + ii * 4 + 3] = (nelpos & 0xff);
				}
				continue;
			}

			switch (TCODES[j])
			{
				case TypeCode::Double:
				{
					array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
					if (TREPEATS[j] == 1)
					{
						dbl = BitConverter::GetBytes(((array<double>^)ExtensionEntryData[j])[i]);
						BINTABLE[cc] = dbl[7];
						BINTABLE[cc + 1] = dbl[6];
						BINTABLE[cc + 2] = dbl[5];
						BINTABLE[cc + 3] = dbl[4];
						BINTABLE[cc + 4] = dbl[3];
						BINTABLE[cc + 5] = dbl[2];
						BINTABLE[cc + 6] = dbl[1];
						BINTABLE[cc + 7] = dbl[0];
					}
					else
						for (int ii = 0; ii < TREPEATS[j]; ii++)
						{
							dbl = BitConverter::GetBytes(((array<double, 2>^)ExtensionEntryData[j])[ii, i]);
							BINTABLE[cc + ii * 8] = dbl[7];
							BINTABLE[cc + ii * 8 + 1] = dbl[6];
							BINTABLE[cc + ii * 8 + 2] = dbl[5];
							BINTABLE[cc + ii * 8 + 3] = dbl[4];
							BINTABLE[cc + ii * 8 + 4] = dbl[3];
							BINTABLE[cc + ii * 8 + 5] = dbl[2];
							BINTABLE[cc + ii * 8 + 6] = dbl[1];
							BINTABLE[cc + ii * 8 + 7] = dbl[0];
						}
					break;
				}

				case TypeCode::Single:
				{
					array<unsigned char>^ flt = gcnew array<unsigned char>(4);
					if (TREPEATS[j] == 1)
					{
						flt = BitConverter::GetBytes(((array<float>^)ExtensionEntryData[j])[i]);
						BINTABLE[cc] = flt[3];
						BINTABLE[cc + 1] = flt[2];
						BINTABLE[cc + 2] = flt[1];
						BINTABLE[cc + 3] = flt[0];
					}
					else
						for (int ii = 0; ii < TREPEATS[j]; ii++)
						{
							flt = BitConverter::GetBytes(((array<float, 2>^)ExtensionEntryData[j])[ii, i]);
							BINTABLE[cc + ii * 4] = flt[3];
							BINTABLE[cc + ii * 4 + 1] = flt[2];
							BINTABLE[cc + ii * 4 + 2] = flt[1];
							BINTABLE[cc + ii * 4 + 3] = flt[0];
						}
					break;
				}

				case TypeCode::Int64:
				{
					__int64 val;
					if (TREPEATS[j] == 1)
					{
						val = ((array<__int64>^)ExtensionEntryData[j])[i];
						BINTABLE[cc] = ((val >> 56) & 0xff);
						BINTABLE[cc+ 1] = ((val >> 48) & 0xff);
						BINTABLE[cc + 2] = ((val >> 40) & 0xff);
						BINTABLE[cc + 3] = ((val >> 32) & 0xff);
						BINTABLE[cc + 4] = ((val >> 24) & 0xff);
						BINTABLE[cc + 5] = ((val >> 16) & 0xff);
						BINTABLE[cc + 6] = ((val >> 8) & 0xff);
						BINTABLE[cc + 7] = (val & 0xff);
					}
					else
						for (int ii = 0; ii < TREPEATS[j]; ii++)
						{
							val = ((array<__int64, 2>^)ExtensionEntryData[j])[ii, i];
							BINTABLE[cc + ii * 8] = ((val >> 56) & 0xff);
							BINTABLE[cc + ii * 8 + 1] = ((val >> 48) & 0xff);
							BINTABLE[cc + ii * 8 + 2] = ((val >> 40) & 0xff);
							BINTABLE[cc + ii * 8 + 3] = ((val >> 32) & 0xff);
							BINTABLE[cc + ii * 8 + 4] = ((val >> 24) & 0xff);
							BINTABLE[cc + ii * 8 + 5] = ((val >> 16) & 0xff);
							BINTABLE[cc + ii * 8 + 6] = ((val >> 8) & 0xff);
							BINTABLE[cc + ii * 8 + 7] = (val & 0xff);
						}
					break;
				}

				case TypeCode::UInt64:
				{
					unsigned __int64 bzero = 9223372036854775808;
					unsigned __int64 val;
					if (TREPEATS[j] == 1)
					{
						val = (((array<__int64>^)ExtensionEntryData[j])[i] - bzero);
						BINTABLE[cc] = ((val >> 56) & 0xff);
						BINTABLE[cc + 1] = ((val >> 48) & 0xff);
						BINTABLE[cc + 2] = ((val >> 40) & 0xff);
						BINTABLE[cc + 3] = ((val >> 32) & 0xff);
						BINTABLE[cc + 4] = ((val >> 24) & 0xff);
						BINTABLE[cc + 5] = ((val >> 16) & 0xff);
						BINTABLE[cc + 6] = ((val >> 8) & 0xff);
						BINTABLE[cc + 7] = (val & 0xff);
					}
					else
						for (int ii = 0; ii < TREPEATS[j]; ii++)
						{
							val = (((array<__int64, 2>^)ExtensionEntryData[j])[ii, i] - bzero);
							BINTABLE[cc + ii * 8] = ((val >> 56) & 0xff);
							BINTABLE[cc + ii * 8 + 1] = ((val >> 48) & 0xff);
							BINTABLE[cc + ii * 8 + 2] = ((val >> 40) & 0xff);
							BINTABLE[cc + ii * 8 + 3] = ((val >> 32) & 0xff);
							BINTABLE[cc + ii * 8 + 4] = ((val >> 24) & 0xff);
							BINTABLE[cc + ii * 8 + 5] = ((val >> 16) & 0xff);
							BINTABLE[cc + ii * 8 + 6] = ((val >> 8) & 0xff);
							BINTABLE[cc + ii * 8 + 7] = (val & 0xff);
						}
					break;
				}

				case TypeCode::Int32:
				{
					__int32 val;
					if (TREPEATS[j] == 1)
					{
						val = ((array<__int32>^)ExtensionEntryData[j])[i];
						BINTABLE[cc] = ((val >> 24) & 0xff);
						BINTABLE[cc + 1] = ((val >> 16) & 0xff);
						BINTABLE[cc + 2] = ((val >> 8) & 0xff);
						BINTABLE[cc + 3] = (val & 0xff);
					}
					else
						for (int ii = 0; ii < TREPEATS[j]; ii++)
						{
							val = ((array<__int32, 2>^)ExtensionEntryData[j])[ii, i];
							BINTABLE[cc + ii * 4] = ((val >> 24) & 0xff);
							BINTABLE[cc + ii * 4 + 1] = ((val >> 16) & 0xff);
							BINTABLE[cc + ii * 4 + 2] = ((val >> 8) & 0xff);
							BINTABLE[cc + ii * 4 + 3] = (val & 0xff);
						}
					break;
				}

				case TypeCode::UInt32:
				{
					unsigned __int32 bzero = 2147483648;
					unsigned __int32 val;
					if (TREPEATS[j] == 1)
					{
						val = (((array<__int32>^)ExtensionEntryData[j])[i] - bzero);
						BINTABLE[cc] = ((val >> 24) & 0xff);
						BINTABLE[cc + 1] = ((val >> 16) & 0xff);
						BINTABLE[cc + 2] = ((val >> 8) & 0xff);
						BINTABLE[cc + 3] = (val & 0xff);
					}
					else
						for (int ii = 0; ii < TREPEATS[j]; ii++)
						{
							val = (((array<__int32, 2>^)ExtensionEntryData[j])[ii, i] - bzero);
							BINTABLE[cc + ii * 4] = ((val >> 24) & 0xff);
							BINTABLE[cc + ii * 4 + 1] = ((val >> 16) & 0xff);
							BINTABLE[cc + ii * 4 + 2] = ((val >> 8) & 0xff);
							BINTABLE[cc + ii * 4 + 3] = (val & 0xff);
						}
					break;
				}

				case TypeCode::Int16:
				{
					__int16 val;
					if (TREPEATS[j] == 1)
					{
						val = ((array<__int16>^)ExtensionEntryData[j])[i];
						BINTABLE[cc] = ((val >> 8) & 0xff);
						BINTABLE[cc + 1] = (val & 0xff);
					}
					else
						for (int ii = 0; ii < TREPEATS[j]; ii++)
						{
							val = ((array<__int16, 2>^)ExtensionEntryData[j])[ii, i];
							BINTABLE[cc + ii * 2] = ((val >> 8) & 0xff);
							BINTABLE[cc + ii * 2 + 1] = (val & 0xff);
						}
					break;
				}

				case TypeCode::UInt16:
				{
					unsigned __int16 bzero = 32768;
					unsigned __int16 val;
					if (TREPEATS[j] == 1)
					{
						val = (((array<__int16>^)ExtensionEntryData[j])[i] - bzero);
						BINTABLE[cc] = ((val >> 8) & 0xff);
						BINTABLE[cc + 1] = (val & 0xff);
					}
					else
						for (int ii = 0; ii < TREPEATS[j]; ii++)
						{
							val = (((array<__int16, 2>^)ExtensionEntryData[j])[ii, i] - bzero);
							BINTABLE[cc + ii * 2] = ((val >> 8) & 0xff);
							BINTABLE[cc + ii * 2 + 1] = (val & 0xff);
						}
					break;
				}

				case TypeCode::SByte:
				{
					if (TREPEATS[j] == 1)
						BINTABLE[cc] = (((array<__int8>^)ExtensionEntryData[j])[i]);
					else
						for (int ii = 0; ii < TREPEATS[j]; ii++)
							BINTABLE[cc + ii] = (((array<__int8, 2>^)ExtensionEntryData[j])[ii, i]);
					break;
				}

				case TypeCode::Byte:
				{
					if (TREPEATS[j] == 1)
						BINTABLE[cc] = (((array<unsigned __int8>^)ExtensionEntryData[j])[i]);
					else
						for (int ii = 0; ii < TREPEATS[j]; ii++)
							BINTABLE[cc + ii] = (((array<unsigned __int8, 2>^)ExtensionEntryData[j])[ii, i]);
					break;
				}

				case TypeCode::Boolean:
				{
					if (TREPEATS[j] == 1)
						BINTABLE[cc] = ((array<bool>^)ExtensionEntryData[j])[i];
					else
						for (int ii = 0; ii < TREPEATS[j]; ii++)
							BINTABLE[cc + ii] = ((array<bool, 2>^)ExtensionEntryData[j])[ii, i];
					break;
				}

				case TypeCode::Char:
				{
					for (int ii = 0; ii < TREPEATS[j]; ii++)
						BINTABLE[cc + ii] = (unsigned char)(((array<String^>^)ExtensionEntryData[j])[i][ii]);
					break;
				}

				default:
				{
					exception = true;
					exceptiontypecode = TCODES[j];
					break;
				}
			}
		}
	}

	if (exception)
	{
		throw gcnew Exception("Data type not recognized for writing as FITS table: '" + exceptiontypecode.ToString() + "'");
		return;
	}
}

void JPFITS::FITSBinTable::MAKETTYPEHEAPARRAYNELSPOS(array<Object^>^ ExtensionEntryData, __int64 &totalBytes)
{
	int pos = 0, nels;
	totalBytes = 0;
	for (int i = 0; i < TTYPEISHEAPARRAYDESC->Length; i++)//same as extension entry data length
		if (TTYPEISHEAPARRAYDESC[i])
		{
			TTYPEHEAPARRAYNELSPOS[i] = gcnew array<int, 2>(2, NAXIS2);
			if (HEAPTCODES[i] == TypeCode::Char)
				for (int j = 0; j < NAXIS2; j++)
				{
					nels = ((String^)(((Array^)ExtensionEntryData[i])->GetValue(j)))->Length;
					TTYPEHEAPARRAYNELSPOS[i][0, j] = nels;
					TTYPEHEAPARRAYNELSPOS[i][1, j] = pos;
					pos += (nels * TYPECODETONBYTES(HEAPTCODES[i]));
				}
			else
				for (int j = 0; j < NAXIS2; j++)
				{
					nels = ((Array^)(((Array^)ExtensionEntryData[i])->GetValue(j)))->Length;
					TTYPEHEAPARRAYNELSPOS[i][0, j] = nels;
					TTYPEHEAPARRAYNELSPOS[i][1, j] = pos;
					pos += (nels * TYPECODETONBYTES(HEAPTCODES[i]));
				}
		}
	totalBytes = pos;
}

void JPFITS::FITSBinTable::MAKEHEAPBYTEARRAY(array<Object^>^ ExtensionEntryData)
{
	__int64 totalbytes = 0;
	MAKETTYPEHEAPARRAYNELSPOS(ExtensionEntryData, totalbytes);
	HEAPDATA = gcnew array<unsigned char>((int)totalbytes);

	int nthread = omp_get_max_threads();
	bool parallel = false;
	if (NAXIS2 >= nthread)
		parallel = true;

	for (int i = 0; i < TTYPEISHEAPARRAYDESC->Length; i++)//same as extension entry data length
	{
		if (!TTYPEISHEAPARRAYDESC[i])
			continue;

		switch (HEAPTCODES[i])
		{
			case TypeCode::Double:
			{
				//#pragma omp parallel for if(parallel)
				for (int y = 0; y < NAXIS2; y++)
				{
					array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
					int pos = TTYPEHEAPARRAYNELSPOS[i][1, y];

					for (int x = 0; x < TTYPEHEAPARRAYNELSPOS[i][0, y]; x++)
					{
						dbl = BitConverter::GetBytes(((array<double>^)(((array<array<double>^>^)(ExtensionEntryData[i]))[y]))[x]);
						HEAPDATA[pos] = dbl[7];
						HEAPDATA[pos + 1] = dbl[6];
						HEAPDATA[pos + 2] = dbl[5];
						HEAPDATA[pos + 3] = dbl[4];
						HEAPDATA[pos + 4] = dbl[3];
						HEAPDATA[pos + 5] = dbl[2];
						HEAPDATA[pos + 6] = dbl[1];
						HEAPDATA[pos + 7] = dbl[0];
						pos += 8;
					}
				}
				break;
			}

			case TypeCode::Single:
			{
				//#pragma omp parallel for if(parallel)
				for (int y = 0; y < NAXIS2; y++)
				{
					array<unsigned char>^ sng = gcnew array<unsigned char>(4);
					int pos = TTYPEHEAPARRAYNELSPOS[i][1, y];

					for (int x = 0; x < TTYPEHEAPARRAYNELSPOS[i][0, y]; x++)
					{
						sng = BitConverter::GetBytes(((array<float>^)(((array<array<float>^>^)(ExtensionEntryData[i]))[y]))[x]);
						HEAPDATA[pos] = sng[3];
						HEAPDATA[pos + 1] = sng[2];
						HEAPDATA[pos + 2] = sng[1];
						HEAPDATA[pos + 3] = sng[0];
						pos += 4;
					}
				}
				break;
			}

			case TypeCode::Int64:
			{
				//#pragma omp parallel for if(parallel)
				for (int y = 0; y < NAXIS2; y++)
				{
					__int64	val;
					int pos = TTYPEHEAPARRAYNELSPOS[i][1, y];

					for (int x = 0; x < TTYPEHEAPARRAYNELSPOS[i][0, y]; x++)
					{
						val = ((array<__int64>^)(((array<array<__int64>^>^)(ExtensionEntryData[i]))[y]))[x];
						HEAPDATA[pos] = ((val >> 56) & 0xff);
						HEAPDATA[pos + 1] = ((val >> 48) & 0xff);
						HEAPDATA[pos + 2] = ((val >> 40) & 0xff);
						HEAPDATA[pos + 3] = ((val >> 32) & 0xff);
						HEAPDATA[pos + 4] = ((val >> 24) & 0xff);
						HEAPDATA[pos + 5] = ((val >> 16) & 0xff);
						HEAPDATA[pos + 6] = ((val >> 8) & 0xff);
						HEAPDATA[pos + 7] = (val & 0xff);
						pos += 8;
					}
				}
				break;
			}

			case TypeCode::UInt64:
			{
				unsigned __int64 bzero = 9223372036854775808;
				//#pragma omp parallel for if(parallel)
				for (int y = 0; y < NAXIS2; y++)
				{
					unsigned __int64 val;
					int pos = TTYPEHEAPARRAYNELSPOS[i][1, y];

					for (int x = 0; x < TTYPEHEAPARRAYNELSPOS[i][0, y]; x++)
					{
						val = (((array<__int64>^)(((array<array<__int64>^>^)(ExtensionEntryData[i]))[y]))[x] - bzero);
						HEAPDATA[pos] = ((val >> 56) & 0xff);
						HEAPDATA[pos + 1] = ((val >> 48) & 0xff);
						HEAPDATA[pos + 2] = ((val >> 40) & 0xff);
						HEAPDATA[pos + 3] = ((val >> 32) & 0xff);
						HEAPDATA[pos + 4] = ((val >> 24) & 0xff);
						HEAPDATA[pos + 5] = ((val >> 16) & 0xff);
						HEAPDATA[pos + 6] = ((val >> 8) & 0xff);
						HEAPDATA[pos + 7] = (val & 0xff);
						pos += 8;
					}
				}
				break;
			}

			case TypeCode::Int32:
			{
				//#pragma omp parallel for if(parallel)
				for (int y = 0; y < NAXIS2; y++)
				{
					__int32 val;
					int pos = TTYPEHEAPARRAYNELSPOS[i][1, y];	

					for (int x = 0; x < TTYPEHEAPARRAYNELSPOS[i][0, y]; x++)
					{						
						val = ((array<__int32>^)(((array<array<__int32>^>^)(ExtensionEntryData[i]))[y]))[x];						
						HEAPDATA[pos] = ((val >> 24) & 0xff);
						HEAPDATA[pos + 1] = ((val >> 16) & 0xff);
						HEAPDATA[pos + 2] = ((val >> 8) & 0xff);
						HEAPDATA[pos + 3] = (val & 0xff);
						pos += 4;
					}
				}
				break;
			}

			case TypeCode::UInt32:
			{
				unsigned __int32 bzero = 2147483648;
				//#pragma omp parallel for if(parallel)
				for (int y = 0; y < NAXIS2; y++)
				{
					unsigned __int32 val;
					int pos = TTYPEHEAPARRAYNELSPOS[i][1, y];

					for (int x = 0; x < TTYPEHEAPARRAYNELSPOS[i][0, y]; x++)
					{
						val = (((array<__int32>^)(((array<array<__int32>^>^)(ExtensionEntryData[i]))[y]))[x] - bzero);
						HEAPDATA[pos] = ((val >> 24) & 0xff);
						HEAPDATA[pos + 1] = ((val >> 16) & 0xff);
						HEAPDATA[pos + 2] = ((val >> 8) & 0xff);
						HEAPDATA[pos + 3] = (val & 0xff);
						pos += 4;
					}
				}
				break;
			}

			case TypeCode::Int16:
			{
				//#pragma omp parallel for if(parallel)
				for (int y = 0; y < NAXIS2; y++)
				{
					__int16 val;
					int pos = TTYPEHEAPARRAYNELSPOS[i][1, y];

					for (int x = 0; x < TTYPEHEAPARRAYNELSPOS[i][0, y]; x++)
					{
						val = ((array<__int16>^)(((array<array<__int16>^>^)(ExtensionEntryData[i]))[y]))[x];
						HEAPDATA[pos] = ((val >> 8) & 0xff);
						HEAPDATA[pos + 1] = (val & 0xff);
						pos += 2;
					}
				}
				break;
			}

			case TypeCode::UInt16:
			{
				unsigned __int16 bzero = 32768;
				//#pragma omp parallel for if(parallel)
				for (int y = 0; y < NAXIS2; y++)
				{
					unsigned __int16 val;
					int pos = TTYPEHEAPARRAYNELSPOS[i][1, y];

					for (int x = 0; x < TTYPEHEAPARRAYNELSPOS[i][0, y]; x++)
					{
						val = (((array<__int16>^)(((array<array<__int16>^>^)(ExtensionEntryData[i]))[y]))[x] - bzero);
						HEAPDATA[pos] = ((val >> 8) & 0xff);
						HEAPDATA[pos + 1] = (val & 0xff);
						pos += 2;
					}
				}
				break;
			}

			case TypeCode::SByte:
			{
				//#pragma omp parallel for if(parallel)
				for (int y = 0; y < NAXIS2; y++)
					for (int x = 0; x < TTYPEHEAPARRAYNELSPOS[i][0, y]; x++)
						HEAPDATA[TTYPEHEAPARRAYNELSPOS[i][1, y] + x] = ((array<__int8>^)(((array<array<__int8>^>^)(ExtensionEntryData[i]))[y]))[x];
				break;
			}

			case TypeCode::Byte:
			{
				//#pragma omp parallel for if(parallel)
				for (int y = 0; y < NAXIS2; y++)
					for (int x = 0; x < TTYPEHEAPARRAYNELSPOS[i][0, y]; x++)
						HEAPDATA[TTYPEHEAPARRAYNELSPOS[i][1, y] + x] = ((array<unsigned __int8>^)(((array<array<unsigned __int8>^>^)(ExtensionEntryData[i]))[y]))[x];
				break;
			}

			case TypeCode::Boolean:
			{
				//#pragma omp parallel for if(parallel)
				for (int y = 0; y < NAXIS2; y++)
					for (int x = 0; x < TTYPEHEAPARRAYNELSPOS[i][0, y]; x++)
						HEAPDATA[TTYPEHEAPARRAYNELSPOS[i][1, y] + x] = ((array<bool>^)(((array<array<bool>^>^)(ExtensionEntryData[i]))[y]))[x];
				break;
			}

			case TypeCode::Char:
			{
				//#pragma omp parallel for if(parallel)
				for (int y = 0; y < NAXIS2; y++)
					for (int x = 0; x < TTYPEHEAPARRAYNELSPOS[i][0, y]; x++)
						HEAPDATA[TTYPEHEAPARRAYNELSPOS[i][1, y] + x] = (unsigned char)((String^)(((array<String^>^)(ExtensionEntryData[i]))[y]))[x];
				break;
			}

			default:
			{
				throw gcnew Exception("Data type not recognized for writing as FITS table: '" + HEAPTCODES[i].ToString() + "'");
				return;
			}
		}
	}	
}

array<String^>^ JPFITS::FITSBinTable::FORMATBINARYTABLEEXTENSIONHEADER()
{
	ArrayList^ hkeyslist = gcnew ArrayList();
	ArrayList^ hvalslist = gcnew ArrayList();
	ArrayList^ hcomslist = gcnew ArrayList();

	hkeyslist->Add("XTENSION");
	hvalslist->Add("BINTABLE");
	hcomslist->Add("binary table extension");
	hkeyslist->Add("BITPIX");
	hvalslist->Add("8");
	hcomslist->Add("8-bit bytes");
	hkeyslist->Add("NAXIS");
	hvalslist->Add("2");
	hcomslist->Add("2-dimensional binary table");
	hkeyslist->Add("NAXIS1");
	hvalslist->Add(NAXIS1.ToString());
	hcomslist->Add("width of table in bytes");
	hkeyslist->Add("NAXIS2");
	hvalslist->Add(NAXIS2.ToString());
	hcomslist->Add("number of rows in table");
	hkeyslist->Add("PCOUNT");
	if (HEAPDATA == nullptr)
		hvalslist->Add("0");
	else
		hvalslist->Add(HEAPDATA->Length.ToString());//we do not write the heap with a gap, it comes right after the bintable, hence PCOUNT is simply the heap size, and the THEAP gap is zero
	hcomslist->Add("size of heap data area (bytes)");
	hkeyslist->Add("GCOUNT");
	hvalslist->Add("1");
	hcomslist->Add("one data group");
	hkeyslist->Add("TFIELDS");
	hvalslist->Add(TFIELDS.ToString());
	hcomslist->Add("number of fields in each row");
	if (EXTENSIONNAME != "")
	{
		hkeyslist->Add("EXTNAME");
		hvalslist->Add(EXTENSIONNAME);
		hcomslist->Add("name of this binary table extension");
	}

	//KEY formats
	for (int i = 0; i < TTYPES->Length; i++)
	{
		//TFORM
		hkeyslist->Add("TFORM" + (i + 1).ToString());
		if (!TTYPEISHEAPARRAYDESC[i])
			hvalslist->Add(TFORMS[i]);
		else
		{
			int max = 0;
			for (int j = 0; j < NAXIS2; j++)
				if (TTYPEHEAPARRAYNELSPOS[i][0, j] > max)
					max = TTYPEHEAPARRAYNELSPOS[i][0, j];
			hvalslist->Add(TFORMS[i] + "(" + max.ToString() + ")");
		}
		if (TTYPEISHEAPARRAYDESC[i])
			if (!TTYPEISCOMPLEX[i])
				hcomslist->Add((2 * TYPECODETONBYTES(TCODES[i])).ToString() + "-byte " + TYPECODESTRING(TCODES[i]) + " heap descriptor for " + TYPECODESTRING(HEAPTCODES[i]));
			else
				hcomslist->Add((2 * TYPECODETONBYTES(TCODES[i])).ToString() + "-byte " + TYPECODESTRING(TCODES[i]) + " heap descriptor for " + TYPECODESTRING(HEAPTCODES[i]) + " complex pair");
		else if (TTYPEISCOMPLEX[i])
			hcomslist->Add((2 * TYPECODETONBYTES(TCODES[i])).ToString() + "-byte " + TYPECODESTRING(TCODES[i]) + " complex pair"); 
		else
			hcomslist->Add(TYPECODETONBYTES(TCODES[i]).ToString() + "-byte " + TYPECODESTRING(TCODES[i]));

		//TTYPE
		hkeyslist->Add("TTYPE" + (i + 1).ToString());
		hvalslist->Add(TTYPES[i]);
		hcomslist->Add("label for field " + (i + 1).ToString());

		//TZERO and TSCAL
		if (!TTYPEISHEAPARRAYDESC[i] && TCODES[i] == TypeCode::SByte || HEAPTCODES[i] == TypeCode::SByte)
		{
			hkeyslist->Add("TZERO" + (i + 1).ToString());
			hvalslist->Add("-128");
			hcomslist->Add("offset for signed 8-bit integers");

			hkeyslist->Add("TSCAL" + (i + 1).ToString());
			hvalslist->Add("1");
			hcomslist->Add("data are not scaled");
		}
		else if (!TTYPEISHEAPARRAYDESC[i] && TCODES[i] == TypeCode::UInt16 || HEAPTCODES[i] == TypeCode::UInt16)
		{
			hkeyslist->Add("TZERO" + (i + 1).ToString());
			hvalslist->Add("32768");
			hcomslist->Add("offset for unsigned 16-bit integers");

			hkeyslist->Add("TSCAL" + (i + 1).ToString());
			hvalslist->Add("1");
			hcomslist->Add("data are not scaled");
		}
		else if (!TTYPEISHEAPARRAYDESC[i] && TCODES[i] == TypeCode::UInt32 || HEAPTCODES[i] == TypeCode::UInt32)
		{
			hkeyslist->Add("TZERO" + (i + 1).ToString());
			hvalslist->Add("2147483648");
			hcomslist->Add("offset for unsigned 32-bit integers");

			hkeyslist->Add("TSCAL" + (i + 1).ToString());
			hvalslist->Add("1");
			hcomslist->Add("data are not scaled");
		}
		else if (!TTYPEISHEAPARRAYDESC[i] && TCODES[i] == TypeCode::UInt64 || HEAPTCODES[i] == TypeCode::UInt64)
		{
			hkeyslist->Add("TZERO" + (i + 1).ToString());
			hvalslist->Add("9223372036854775808");
			hcomslist->Add("offset for unsigned 64-bit integers");

			hkeyslist->Add("TSCAL" + (i + 1).ToString());
			hvalslist->Add("1");
			hcomslist->Add("data are not scaled");
		}

		//TUNIT
		if (TUNITS != nullptr && TUNITS[i] != nullptr && TUNITS[i] != "")
		{
			hkeyslist->Add("TUNIT" + (i + 1).ToString());
			hvalslist->Add(TUNITS[i]);
			hcomslist->Add("physical unit of field");
		}

		//TDIM
		if (TDIMS[i] != nullptr)//then it is a multi D array, and the dims should exist for this entry
		{
			hkeyslist->Add("TDIM" + (i + 1).ToString());
			String^ dim = "(";
			for (int j = 0; j < TDIMS[i]->Length; j++)
				dim += (TDIMS[i][j].ToString() + ",");
			dim = dim->Remove(dim->Length - 1) + ")";
			hvalslist->Add(dim);
			hcomslist->Add("N-dim array dimensions");
		}
	}

	//EXTRAKEYS
	if (EXTRAKEYS != nullptr)
		for (int i = 0; i < EXTRAKEYS->Length; i++)
		{
			hkeyslist->Add(EXTRAKEYS[i]->ToUpper());
			hvalslist->Add(EXTRAKEYVALS[i]);
			hcomslist->Add(EXTRAKEYCOMS[i]);
		}

	hkeyslist->Add("END     ");
	hvalslist->Add("");
	hcomslist->Add("");

	int NKeys = hkeyslist->Count;
	int NCards = (NKeys - 1) / 36;
	int NBlankKeys = (NCards + 1) * 36 - NKeys;
	array<String^>^ headerkeys = gcnew array<String^>((NCards + 1) * 36);
	array<String^>^ headerkeyvals = gcnew array<String^>((NCards + 1) * 36);
	array<String^>^ headerkeycoms = gcnew array<String^>((NCards + 1) * 36);
	array<String^>^ header = gcnew array<String^>((NCards + 1) * 36);

	for (int i = 0; i < NKeys; i++)
	{
		headerkeys[i] = (String^)hkeyslist[i];
		headerkeyvals[i] = (String^)hvalslist[i];
		headerkeycoms[i] = (String^)hcomslist[i];
	}

	String^ key;
	String^ value;
	String^ comment;
	for (int i = 0; i < NKeys - 1; i++)
	{
		key = headerkeys[i];
		key = key->Trim();//in case some idiot put spaces in the front or back...if in the middle??
		int L = key->Length;
		if (key == "COMMENT")
			key = "COMMENT ";//8 long, comment follows; key now = "COMMENT "..., no "=" needed
		else if (L >= 8)
		{
			key = key->Substring(0, 8);
			key += "= ";//key formatting done
		}
		else if (L < 8)
		{
			for (int i = 0; i < 8 - L; i++)
				key += " ";//pad right
			key += "= ";//key formatting done
		}

		//do value formatting
		if (JPMath::IsNumeric(headerkeyvals[i]))//then we have a numeric key value
		{
			double val = ::Convert::ToDouble(headerkeyvals[i]);
			if (val == 9223372036854775808)
				value = "9223372036854775808";
			else
				value = val.ToString();
			L = value->Length;
			if (L > 20)
				value = value->Substring(0, 20);
			if (L < 20)
				for (int i = 0; i < 20 - L; i++)
					value = " " + value;//pad left
		}
		else//else it must be a string or comment.
		{
			value = headerkeyvals[i];
			L = value->Length;
			if (L >= 18)
			{
				value = value->Substring(0, 18);
				value = "'" + value + "'";
			}
			if (L < 18)
			{
				value = "'" + value + "'";
				for (int i = 0; i < 18 - L; i++)
					value += " ";//pad right
			}
			if (headerkeyvals[i]->Trim() == "T")
				value = "                   T";
			if (headerkeyvals[i]->Trim() == "F")
				value = "                   F";
		}
		//value formatting done

		//do comment formatting...always a string
		comment = headerkeycoms[i];
		L = comment->Length;
		if (L > 48)
			comment = comment->Substring(0, 48);
		if (L < 48)
			for (int i = 0; i < 48 - L; i++)
				comment += " ";//pad right
		comment = " /" + comment;//comment formatting done

		//check for COMMENT and reconfigure key line...it isn't the most efficient approach but I dont care.
		if (key == "COMMENT ")
		{
			comment = headerkeyvals[i] + headerkeycoms[i];
			value = "";
			L = comment->Length;
			if (L > 72)
				comment = comment->Substring(0, 72);
			if (L < 72)
				for (int i = 0; i < 72 - L; i++)
					comment += " ";//pad right
		}

		header[i] = key + value + comment;
	}

	header[NKeys - 1] = "END                                                                             ";

	for (int i = 0; i < NBlankKeys; i++)
		header[NKeys + i] = "                                                                                ";

	return header;
}

int JPFITS::FITSBinTable::TFORMTONBYTES(String^ tform, int& instances)
{
	int N = 1;
	if (tform->Length > 1)
		if (tform->Contains("Q") || tform->Contains("P"))////heap ttype
		{
			if (JPMath::IsNumeric(tform->Substring(0, 1)))
				N = Convert::ToInt32(tform->Substring(0, 1));//might be zero...one default
			if (tform->Contains("Q"))
				tform = "Q";
			else
				tform = "P";
		}
		else
			N = ::Convert::ToInt32(tform->Substring(0, tform->Length - 1));//bintable ttype
	instances = N;

	wchar_t f = Convert::ToChar(tform->Substring(tform->Length - 1));

	switch (f)
	{
		case 'M':
		case 'Q':
		{
			instances *= 2;
			return N *= 16;
		}

		case 'C':
		case 'P':
		{
			instances *= 2;
			return N *= 8;
		}
		
		case 'D':
		case 'K':
			return N *= 8;

		case 'J':
		case 'E':
			return N *= 4;

		case 'I':
			return N *= 2;

		case 'A':
		case 'B':
		case 'L':
			return N *= 1;

		case 'X':
			return int(Math::Ceiling(double(N) / 8));

		default:
			throw gcnew Exception("Unrecognized TFORM: '" + tform + "'");
	}
}

String^ JPFITS::FITSBinTable::TYPECODETFORM(TypeCode typecode)
{
	switch (typecode)
	{
		case TypeCode::Double:
			return "D";
		case TypeCode::UInt64:
		case TypeCode::Int64:
			return "K";
		case TypeCode::Single:
			return "E";
		case TypeCode::UInt32:
		case TypeCode::Int32:
			return "J";
		case TypeCode::UInt16:
		case TypeCode::Int16:
			return "I";
		case TypeCode::Byte:
		case TypeCode::SByte:
			return "B";
		case TypeCode::Boolean:
			return "L";
		case TypeCode::Char:
			return "A";
		default:
			throw gcnew Exception("Unrecognized typecode: '" + typecode.ToString() + "'");
	}
}

String^ JPFITS::FITSBinTable::TYPECODESTRING(TypeCode typecode)
{
	switch (typecode)
	{
		case TypeCode::Double:
			return "DOUBLE";

		case TypeCode::UInt64:
		case TypeCode::Int64:
		case TypeCode::UInt32:
		case TypeCode::Int32:
		case TypeCode::UInt16:
		case TypeCode::Int16:
			return "INTEGER";

		case TypeCode::Single:
			return "SINGLE";

		case TypeCode::Byte:
		case TypeCode::SByte:
			return "BYTE";

		case TypeCode::Boolean:
			return "LOGICAL";

		case TypeCode::Char:
			return "CHAR";

		default:
			throw gcnew Exception("Unrecognized typecode: '" + typecode.ToString() + "'");
	}
}

TypeCode JPFITS::FITSBinTable::TFORMTYPECODE(String^ tform)
{
	wchar_t c = Convert::ToChar(tform->Substring(tform->Length - 1));

	switch (c)
	{
		case 'L':
			return TypeCode::Boolean;

		case 'X':
		case 'B':
			return TypeCode::Byte;

		case 'I':
			return TypeCode::Int16;

		case 'J':
		case 'P':
			return TypeCode::Int32;

		case 'K':
		case 'Q':
			return TypeCode::Int64;

		case 'A':
			return TypeCode::Char;

		case 'E':
		case 'C':
			return TypeCode::Single;

		case 'D':
		case 'M':
			return TypeCode::Double;
	
		default:
			throw gcnew Exception("Unrecognized TFORM: '" + tform + "'");
	}
}

int JPFITS::FITSBinTable::TYPECODETONBYTES(TypeCode typecode)
{
	switch (typecode)
	{
		case TypeCode::Double:
		case TypeCode::UInt64:
		case TypeCode::Int64:
			return 8;

		case TypeCode::UInt32:
		case TypeCode::Int32:
		case TypeCode::Single:
			return 4;

		case TypeCode::UInt16:
		case TypeCode::Int16:
			return 2;

		case TypeCode::Byte:
		case TypeCode::SByte:
		case TypeCode::Boolean:
		case TypeCode::Char:
			return 1;

		default:
			throw gcnew Exception("Unrecognized typecode: '" + typecode.ToString() + "'");
	}
}

void JPFITS::FITSBinTable::EATRAWBINTABLEHEADER(ArrayList^ header)
{
	//reset
	BITPIX = 0, NAXIS = 0, NAXIS1 = 0, NAXIS2 = 0, TFIELDS = 0;

	ArrayList^ extras = gcnew ArrayList();//for possible extras

	HEADER = gcnew array<String^>(header->Count);
	String^ strheaderline;
	int ttypeindex = -1;

	for (int i = 0; i < header->Count; i++)
	{
		strheaderline = (String^)header[i];
		HEADER[i] = strheaderline;

		if (BITPIX == 0)
			if (strheaderline->Substring(0, 8)->Trim()->Equals("BITPIX"))
			{
				BITPIX = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				continue;
			}
		if (NAXIS == 0)
			if (strheaderline->Substring(0, 8)->Trim()->Equals("NAXIS"))
			{
				NAXIS = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				continue;
			}
		if (NAXIS1 == 0)
			if (strheaderline->Substring(0, 8)->Trim()->Equals("NAXIS1"))
			{
				NAXIS1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				continue;
			}
		if (NAXIS2 == 0)
			if (strheaderline->Substring(0, 8)->Trim()->Equals("NAXIS2"))
			{
				NAXIS2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				continue;
			}
		if (TFIELDS == 0)
			if (strheaderline->Substring(0, 8)->Trim()->Equals("TFIELDS"))
			{
				TFIELDS = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				TTYPES = gcnew array<String^>(TFIELDS);
				TFORMS = gcnew array<String^>(TFIELDS);
				TBYTES = gcnew array<int>(TFIELDS);
				TREPEATS = gcnew array<int>(TFIELDS);
				TCODES = gcnew array<::TypeCode>(TFIELDS);
				TUNITS = gcnew array<String^>(TFIELDS);
				TTYPEISCOMPLEX = gcnew array<bool>(TFIELDS);
				TTYPEISHEAPARRAYDESC = gcnew array<bool>(TFIELDS);
				HEAPTCODES = gcnew array<::TypeCode>(TFIELDS);
				TDIMS = gcnew array<array<int>^>(TFIELDS);
				TTYPEHEAPARRAYNELSPOS = gcnew array<array<int, 2>^>(TFIELDS);
				continue;
			}

		if (strheaderline->Substring(0, 8)->Trim()->Equals("TTYPE" + (ttypeindex + 2).ToString()) || strheaderline->Substring(0, 8)->Trim()->Equals("TFORM" + (ttypeindex + 2).ToString()))
			ttypeindex++;

		if (strheaderline->Substring(0, 8)->Trim()->Equals("TTYPE" + (ttypeindex + 1).ToString()))
		{
			int f = strheaderline->IndexOf("'");
			int l = strheaderline->LastIndexOf("'");
			TTYPES[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1)->Trim();
			continue;
		}

		if (strheaderline->Substring(0, 8)->Trim()->Equals("TFORM" + (ttypeindex + 1).ToString()))
		{
			int f = strheaderline->IndexOf("'");
			int l = strheaderline->LastIndexOf("'");
			TFORMS[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1)->Trim();
			int instances = 1;
			TBYTES[ttypeindex] = TFORMTONBYTES(TFORMS[ttypeindex], instances);
			TREPEATS[ttypeindex] = instances;	
			if (TFORMS[ttypeindex]->Contains("Q") || TFORMS[ttypeindex]->Contains("P"))//heap form
			{
				TTYPEHEAPARRAYNELSPOS[ttypeindex] = GETHEAPTTYPENELSPOS(ttypeindex);
				TTYPEISHEAPARRAYDESC[ttypeindex] = true;
				if (TFORMS[ttypeindex]->Contains("Q"))
				{
					TCODES[ttypeindex] = TFORMTYPECODE("Q");
					HEAPTCODES[ttypeindex] = TFORMTYPECODE(TFORMS[ttypeindex]->Substring(TFORMS[ttypeindex]->IndexOf("Q") + 1, 1));
				}
				else
				{
					TCODES[ttypeindex] = TFORMTYPECODE("P");
					HEAPTCODES[ttypeindex] = TFORMTYPECODE(TFORMS[ttypeindex]->Substring(TFORMS[ttypeindex]->IndexOf("P") + 1, 1));
				}
				if (HEAPTCODES[ttypeindex] == TypeCode::Double || HEAPTCODES[ttypeindex] == TypeCode::Single)
					TTYPEISCOMPLEX[ttypeindex] = (TFORMS[ttypeindex]->Contains("M") || TFORMS[ttypeindex]->Contains("C"));
			}
			else
			{
				TCODES[ttypeindex] = TFORMTYPECODE(TFORMS[ttypeindex]);
				if (TCODES[ttypeindex] == TypeCode::Double || TCODES[ttypeindex] == TypeCode::Single)
					TTYPEISCOMPLEX[ttypeindex] = (TFORMS[ttypeindex]->Contains("M") || TFORMS[ttypeindex]->Contains("C"));
			}
			continue;
		}

		if (strheaderline->Substring(0, 8)->Trim()->Equals("TUNIT" + (ttypeindex + 1).ToString()))
		{
			int f = strheaderline->IndexOf("'");
			int l = strheaderline->LastIndexOf("'");
			if (f != -1)
				TUNITS[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1)->Trim();
			else
				TUNITS[ttypeindex] = strheaderline->Substring(10, 20)->Trim();
			continue;
		}

		if (strheaderline->Substring(0, 8)->Trim()->Equals("TDIM" + (ttypeindex + 1).ToString()))
		{
			ArrayList^ dimslist = gcnew ArrayList();
			//TDIMn = '(5,4,3)'
			int f = strheaderline->IndexOf("'");
			int l = strheaderline->LastIndexOf("'");
			String^ dimline = strheaderline->Substring(f + 1, l - f - 1)->Trim();//(5,4,3)
			dimline = dimline->Substring(1);//5,4,3)
			dimline = dimline->Substring(0, dimline->Length - 1);//5,4,3
			int lastcommaindex = 0;
			int nextcommaindex = -1;
			while (dimline->IndexOf(",", lastcommaindex + 1) != -1)
			{
				nextcommaindex = dimline->IndexOf(",", lastcommaindex + 1);
				dimslist->Add(dimline->Substring(lastcommaindex, nextcommaindex - lastcommaindex));
				lastcommaindex = nextcommaindex + 1;				
			}
			dimslist->Add(dimline->Substring(lastcommaindex));
			TDIMS[ttypeindex] = gcnew array<int>(dimslist->Count);
			for (int i = 0; i < dimslist->Count; i++)
				TDIMS[ttypeindex][i] = Convert::ToInt32((String^)dimslist[i]);
			continue;
		}

		//need to determine if the TypeCode here is supposed to be for signed or unsigned IF the type is an integer (8, 16, 32, or 64 bit)
		//therefore find the TSCALE and TZERO for this entry...if they don't exist then it is signed, if they do exist
		//then it is whatever values they are, for either determined signed or unsigned
		//then set this current tcode[typeindex] to what it should be
		if (strheaderline->Substring(0, 8)->Trim()->Equals("TZERO" + (ttypeindex + 1).ToString()))
		{
			if (!TTYPEISHEAPARRAYDESC[ttypeindex])
				switch (TCODES[ttypeindex])
				{
					case TypeCode::Byte:
					{
						if (Convert::ToSByte(strheaderline->Substring(10, 20)->Trim()) == -128)//then it is a signed
							TCODES[ttypeindex] = TypeCode::SByte;
						break;
					}
					case TypeCode::Int16:
					{
						if (Convert::ToUInt16(strheaderline->Substring(10, 20)->Trim()) == 32768)//then it is an unsigned
							TCODES[ttypeindex] = TypeCode::UInt16;
						break;
					}
					case TypeCode::Int32:
					{
						if (Convert::ToUInt32(strheaderline->Substring(10, 20)->Trim()) == 2147483648)//then it is an unsigned
							TCODES[ttypeindex] = TypeCode::UInt32;
						break;
					}
					case TypeCode::Int64:
					{
						if (Convert::ToUInt64(strheaderline->Substring(10, 20)->Trim()) == 9223372036854775808)//then it is an unsigned
							TCODES[ttypeindex] = TypeCode::UInt64;
						break;
					}
					default:
					{
						throw gcnew Exception("Unrecognized TypeCode in EATRAWBINTABLEHEADER at TZERO analysis");
						break;
					}
				}
			else
				switch (HEAPTCODES[ttypeindex])
				{
					case TypeCode::Byte:
					{
						if (Convert::ToSByte(strheaderline->Substring(10, 20)->Trim()) == -128)//then it is a signed
							HEAPTCODES[ttypeindex] = TypeCode::SByte;
						break;
					}
					case TypeCode::Int16:
					{
						if (Convert::ToUInt16(strheaderline->Substring(10, 20)->Trim()) == 32768)//then it is an unsigned
							HEAPTCODES[ttypeindex] = TypeCode::UInt16;
						break;
					}
					case TypeCode::Int32:
					{
						if (Convert::ToUInt32(strheaderline->Substring(10, 20)->Trim()) == 2147483648)//then it is an unsigned
							HEAPTCODES[ttypeindex] = TypeCode::UInt32;
						break;
					}
					case TypeCode::Int64:
					{
						if (Convert::ToUInt64(strheaderline->Substring(10, 20)->Trim()) == 9223372036854775808)//then it is an unsigned
							HEAPTCODES[ttypeindex] = TypeCode::UInt64;
						break;
					}
					default:
					{
						throw gcnew Exception("Unrecognized TypeCode in EATRAWBINTABLEHEADER at TZERO analysis");
						break;
					}
				}
			continue;
		}

		if (strheaderline->Substring(0, 8)->Trim()->Equals("TSCAL" + (ttypeindex + 1).ToString()))//don't need to do anything with this
			continue;

		String^ key = strheaderline->Substring(0, 8)->Trim();
		if (key->Length > 0 && key->Substring(0, 1) == "T" && JPMath::IsNumeric(key->Substring(key->Length - 1)))//then likely it is some other T____n field which isn't explicitly coded above...
			continue;

		//should now only be where extra keys might remain...so add them etc
		extras->Add(header[i]);
	}

	if (extras->Count == 0)
		return;

	EXTRAKEYS = gcnew array<String^>(extras->Count);
	EXTRAKEYVALS = gcnew array<String^>(extras->Count);
	EXTRAKEYCOMS = gcnew array<String^>(extras->Count);

	for (int i = 0; i < extras->Count; i++)
	{
		String^ line = (String^)extras[i];
		EXTRAKEYS[i] = line->Substring(0, 8)->Trim();

		if (EXTRAKEYS[i] == "COMMENT")
		{
			EXTRAKEYVALS[i] = line->Substring(8, 18);
			EXTRAKEYCOMS[i] = line->Substring(26);
		}
		else
		{
			if (JPMath::IsNumeric(line->Substring(10, 20)))//this has to work if it is supposed to be a numeric value here
				EXTRAKEYVALS[i] = line->Substring(10, 20)->Trim();//get rid of leading and trailing white space
			else
			{
				String^ nock = "'";
				EXTRAKEYVALS[i] = line->Substring(10, 20)->Trim();
				EXTRAKEYVALS[i] = EXTRAKEYVALS[i]->Trim(nock->ToCharArray());
				EXTRAKEYVALS[i] = EXTRAKEYVALS[i]->Trim();
			}
			EXTRAKEYCOMS[i] = line->Substring(32)->Trim();
		}
	}
}

