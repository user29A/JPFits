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
		throw gcnew Exception("File not formatted as FITS file, or indicates no extensions present.");
		return;
	}

	ArrayList^ header = gcnew ArrayList();
	__int64 extensionstartposition, extensionendposition;

	if (!FITSFILEOPS::SEEKEXTENSION(fs, "BINTABLE", extensionName, header, extensionstartposition, extensionendposition))
	{
		fs->Close();
		throw gcnew Exception("Could not find BINTABLE with name '" + extensionName + "'");
		return;
	}

	FILENAME = fileName;
	EXTENSIONNAME = extensionName;
	
	EATRAWBINTABLEHEADER(header);

	BINTABLE = gcnew array<unsigned char>(int(extensionendposition - fs->Position));
	fs->Read(BINTABLE, 0, BINTABLE->Length);

	fs->Close();
}

array<String^>^ JPFITS::FITSBinTable::GetAllExtensionNames(String^ FileName)
{
	return FITSFILEOPS::GETALLEXTENSIONNAMES(FileName, "BINTABLE");
}

array<double>^ JPFITS::FITSBinTable::GetTTYPEEntry(String^ ExtensionEntryLabel)
{
	int w = 0, h = 0;

	return GetTTYPEEntry(ExtensionEntryLabel, w, h);
}

array<double>^ JPFITS::FITSBinTable::GetTTYPEEntry(String^ ttypeEntryLabel, int& width, int& height)
{
	TypeCode tcode;
	int rank;
	Object^ obj = GetTTYPEEntry(ttypeEntryLabel, tcode, rank);
	int naxis1 = 0, naxis2 = 0;
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
				for (int x = 0; x < width; x++)
					for (int y = 0; y < height; y++)
						result[x * height + y] = ((array<double, 2>^)(obj))[x, y];
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
				for (int x = 0; x < width; x++)
					for (int y = 0; y < height; y++)
						result[x * height + y] = (double)((array<__int64, 2>^)(obj))[x, y];
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
				for (int x = 0; x < width; x++)
					for (int y = 0; y < height; y++)
						result[x * height + y] = (double)((array<float, 2>^)(obj))[x, y];
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
				for (int x = 0; x < width; x++)
					for (int y = 0; y < height; y++)
						result[x * height + y] = (double)((array<unsigned __int32, 2>^)(obj))[x, y];
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
				for (int x = 0; x < width; x++)
					for (int y = 0; y < height; y++)
						result[x * height + y] = (double)((array<__int32, 2>^)(obj))[x, y];
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
				for (int x = 0; x < width; x++)
					for (int y = 0; y < height; y++)
						result[x * height + y] = (double)((array<unsigned __int16, 2>^)(obj))[x, y];
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
				for (int x = 0; x < width; x++)
					for (int y = 0; y < height; y++)
						result[x * height + y] = (double)((array<__int16, 2>^)(obj))[x, y];
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
				for (int x = 0; x < width; x++)
					for (int y = 0; y < height; y++)
						result[x * height + y] = (double)((array<unsigned __int8, 2>^)(obj))[x, y];
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
				for (int x = 0; x < width; x++)
					for (int y = 0; y < height; y++)
						result[x * height + y] = (double)((array<__int8, 2>^)(obj))[x, y];
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
				for (int x = 0; x < width; x++)
					for (int y = 0; y < height; y++)
						result[x * height + y] = (double)((array<bool, 2>^)(obj))[x, y];
			}
			break;
		}

		default:
			throw gcnew Exception("Unrecognized TypeCode: '" + tcode.ToString() + "'");
	}

	return result;
}

Object^ JPFITS::FITSBinTable::GetTTYPEEntry(String^ ttypeEntryLabel, TypeCode &objectTypeCode, int &objectArrayRank)
{
	int extensionentry_index = -1;
	for (int i = 0; i < TTYPES->Length; i++)
		if (TTYPES[i] == ttypeEntryLabel)
		{
			extensionentry_index = i;
			break;
		}

	if (extensionentry_index == -1)
	{
		throw gcnew Exception("Extension Entry TTYPE Label wasn't found: '" + ttypeEntryLabel + "'");
		return nullptr;
	}

	objectTypeCode = TCODES[extensionentry_index];
	if (TINSTANCES[extensionentry_index] == 1)
		objectArrayRank = 1;
	else
		objectArrayRank = 2;

	int byteoffset = 0;
	for (int i = 0; i < extensionentry_index; i++)
		byteoffset += TBYTES[i];
	int currentbyte;

	switch (TCODES[extensionentry_index])
	{
		case ::TypeCode::Double:
		{
			if (objectArrayRank == 1)
			{
				array<double>^ vector = gcnew array<double>(TINSTANCES[extensionentry_index] * NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
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
						vector[i * TINSTANCES[extensionentry_index] + j] = BitConverter::ToDouble(dbl, 0);
					}
				}
				return vector;
			}
			else
			{
				array<double, 2>^ arrya = gcnew array<double, 2>(TINSTANCES[extensionentry_index], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
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

		case (::TypeCode::Int64):
		{
			if (objectArrayRank == 1)
			{
				array<__int64>^ vector = gcnew array<__int64>(TINSTANCES[extensionentry_index] * NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ i64 = gcnew array<unsigned char>(8);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
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
						vector[i * TINSTANCES[extensionentry_index] + j] = BitConverter::ToInt64(i64, 0);
					}
				}
				return vector;
			}
			else
			{
				array<__int64, 2>^ arrya = gcnew array<__int64, 2>(TINSTANCES[extensionentry_index], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ i64 = gcnew array<unsigned char>(8);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
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

		case ::TypeCode::Single:
		{
			if (objectArrayRank == 1)
			{
				array<float>^ vector = gcnew array<float>(TINSTANCES[extensionentry_index] * NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ sng = gcnew array<unsigned char>(4);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 4;
						sng[3] = BINTABLE[currentbyte];
						sng[2] = BINTABLE[currentbyte + 1];
						sng[1] = BINTABLE[currentbyte + 2];
						sng[0] = BINTABLE[currentbyte + 3];
						vector[i * TINSTANCES[extensionentry_index] + j] = BitConverter::ToSingle(sng, 0);
					}
				}
				return vector;
			}
			else			
			{
				array<float, 2>^ arrya = gcnew array<float, 2>(TINSTANCES[extensionentry_index], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ sng = gcnew array<unsigned char>(4);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
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

		case ::TypeCode::UInt32:
		{
			if (objectArrayRank == 1)
			{
				//unsigned __int32 bzero = 2147483648;// + bzero somewhere or does bitcoverter do it all???
				array<unsigned __int32>^ vector = gcnew array<unsigned __int32>(TINSTANCES[extensionentry_index] * NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ uint32 = gcnew array<unsigned char>(4);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 4;
						uint32[3] = BINTABLE[currentbyte];
						uint32[2] = BINTABLE[currentbyte + 1];
						uint32[1] = BINTABLE[currentbyte + 2];
						uint32[0] = BINTABLE[currentbyte + 3];
						vector[i * TINSTANCES[extensionentry_index] + j] = BitConverter::ToUInt32(uint32, 0); // + bzero somewhere or does bitcoverter do it all???
					}
				}
				return vector;
			}
			else
			{
				array<unsigned __int32, 2>^ arrya = gcnew array<unsigned __int32, 2>(TINSTANCES[extensionentry_index], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ uint32 = gcnew array<unsigned char>(4);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 4;
						uint32[3] = BINTABLE[currentbyte];
						uint32[2] = BINTABLE[currentbyte + 1];
						uint32[1] = BINTABLE[currentbyte + 2];
						uint32[0] = BINTABLE[currentbyte + 3];
						arrya[j, i] = BitConverter::ToUInt32(uint32, 0); // + bzero somewhere or does bitcoverter do it all???
					}
				}
				return arrya;
			}
			break;
		}

		case ::TypeCode::Int32:
		{
			if (objectArrayRank == 1)
			{
				//unsigned __int32 bzero = 0;// + bzero somewhere or does bitcoverter do it all???
				array<__int32>^ vector = gcnew array<__int32>(TINSTANCES[extensionentry_index] * NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ int32 = gcnew array<unsigned char>(4);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 4;
						int32[3] = BINTABLE[currentbyte];
						int32[2] = BINTABLE[currentbyte + 1];
						int32[1] = BINTABLE[currentbyte + 2];
						int32[0] = BINTABLE[currentbyte + 3];
						vector[i * TINSTANCES[extensionentry_index] + j] = BitConverter::ToInt32(int32, 0); // + bzero somewhere or does bitcoverter do it all???
					}
				}
				return vector;
			}
			else			
			{
				array<__int32, 2>^ arrya = gcnew array<__int32, 2>(TINSTANCES[extensionentry_index], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ int32 = gcnew array<unsigned char>(4);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 4;
						int32[3] = BINTABLE[currentbyte];
						int32[2] = BINTABLE[currentbyte + 1];
						int32[1] = BINTABLE[currentbyte + 2];
						int32[0] = BINTABLE[currentbyte + 3];
						arrya[j, i] = BitConverter::ToInt32(int32, 0); // + bzero somewhere or does bitcoverter do it all???
					}
				}
				return arrya;
			}
			break;
		}

		case ::TypeCode::UInt16:
		{
			if (objectArrayRank == 1)
			{
				//unsigned __int32 bzero = 32768;// + bzero somewhere or does bitcoverter do it all???
				array<unsigned __int16>^ vector = gcnew array<unsigned __int16>(TINSTANCES[extensionentry_index] * NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ uint16 = gcnew array<unsigned char>(2);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 2;
						uint16[1] = BINTABLE[currentbyte];
						uint16[0] = BINTABLE[currentbyte + 1];
						vector[i * TINSTANCES[extensionentry_index] + j] = BitConverter::ToUInt16(uint16, 0); // + bzero somewhere or does bitcoverter do it all???
					}
				}
				return vector;
			}
			else			
			{
				array<unsigned __int16, 2>^ arrya = gcnew array<unsigned __int16, 2>(TINSTANCES[extensionentry_index], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ uint16 = gcnew array<unsigned char>(2);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 2;
						uint16[1] = BINTABLE[currentbyte];
						uint16[0] = BINTABLE[currentbyte + 1];
						arrya[j, i] = BitConverter::ToUInt16(uint16, 0); // + bzero somewhere or does bitcoverter do it all???
					}
				}
				return arrya;
			}
			break;
		}

		case ::TypeCode::Int16:
		{
			if (objectArrayRank == 1)
			{
				//unsigned __int32 bzero = 0;// + bzero somewhere or does bitcoverter do it all???
				array<__int16>^ vector = gcnew array<__int16>(TINSTANCES[extensionentry_index] * NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ int16 = gcnew array<unsigned char>(2);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 2;
						int16[1] = BINTABLE[currentbyte];
						int16[0] = BINTABLE[currentbyte + 1];
						vector[i * TINSTANCES[extensionentry_index] + j] = BitConverter::ToInt16(int16, 0); // + bzero somewhere or does bitcoverter do it all???
					}
				}
				return vector;
			}
			else			
			{
				array<__int16, 2>^ arrya = gcnew array<__int16, 2>(TINSTANCES[extensionentry_index], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
				{
					array<unsigned char>^ int16 = gcnew array<unsigned char>(2);
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j * 2;
						int16[1] = BINTABLE[currentbyte];
						int16[0] = BINTABLE[currentbyte + 1];
						arrya[j, i] = BitConverter::ToInt16(int16, 0); // + bzero somewhere or does bitcoverter do it all???
					}
				}
				return arrya;
			}
			break;
		}

		case ::TypeCode::Byte:
		{
			if (objectArrayRank == 1)
			{
				//unsigned __int8 bzero = 128;// + bzero somewhere or does bitcoverter do it all???
				array<unsigned __int8>^ vector = gcnew array<unsigned __int8>(TINSTANCES[extensionentry_index] * NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j;
						vector[i * TINSTANCES[extensionentry_index] + j] = (unsigned __int8)BINTABLE[currentbyte] + 128; // + bzero somewhere or does bitcoverter do it all???
					}
				return vector;
			}
			else			
			{
				array<unsigned __int8, 2>^ arrya = gcnew array<unsigned __int8, 2>(TINSTANCES[extensionentry_index], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j;
						arrya[j, i] = (unsigned __int8)BINTABLE[currentbyte] + 128; // + bzero somewhere or does bitcoverter do it all???
					}
				return arrya;
			}
			break;
		}

		case ::TypeCode::SByte:
		{
			if (objectArrayRank == 1)
			{
				//unsigned __int8 bzero = 0;// + bzero somewhere or does bitcoverter do it all???
				array<__int8>^ vector = gcnew array<__int8>(TINSTANCES[extensionentry_index] * NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j;
						vector[i * TINSTANCES[extensionentry_index] + j] = (__int8)BINTABLE[currentbyte]; // + bzero somewhere or does bitcoverter do it all???
					}
				return vector;
			}
			else			
			{
				array<__int8, 2>^ arrya = gcnew array<__int8, 2>(TINSTANCES[extensionentry_index], NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j;
						arrya[j, i] = (__int8)BINTABLE[currentbyte]; // + bzero somewhere or does bitcoverter do it all???
					}
				return arrya;
			}
			break;
		}

		case ::TypeCode::Boolean:
		{
			if (objectArrayRank == 1)
			{
				array<bool>^ vector = gcnew array<bool>(TINSTANCES[extensionentry_index] * NAXIS2);
				#pragma omp parallel for private(currentbyte)
				for (int i = 0; i < NAXIS2; i++)
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j;
						vector[i * TINSTANCES[extensionentry_index] + j] = Convert::ToBoolean(BINTABLE[currentbyte]);
					}
				return vector;
			}
			else			
			{
				array<bool, 2>^ arrya = gcnew array<bool, 2>(TINSTANCES[extensionentry_index], NAXIS2);
				#pragma omp parallel for
				for (int i = 0; i < NAXIS2; i++)
					for (int j = 0; j < TINSTANCES[extensionentry_index]; j++)
					{
						currentbyte = byteoffset + i * NAXIS1 + j;
						arrya[j, i] = Convert::ToBoolean(BINTABLE[currentbyte]);
					}
				return arrya;
			}
			break;
		}

		default:
		{
			throw gcnew Exception("Unrecognized TypeCode: '" + TCODES[extensionentry_index].ToString() + "'");
			return nullptr;
		}
	}
}

void JPFITS::FITSBinTable::RemoveTTYPEEntry(String^ ttypeEntryLabel)
{
	int extensionentry_index = -1;
	for (int i = 0; i < TTYPES->Length; i++)
		if (TTYPES[i] == ttypeEntryLabel)
		{
			extensionentry_index = i;
			break;
		}

	if (extensionentry_index == -1)
	{
		throw gcnew Exception("Extension Entry TTYPE Label wasn't found: '" + ttypeEntryLabel + "'");
		return;
	}

	array<Object^>^ newEntryDataObjs = gcnew array<Object^>(TFIELDS - 1);
	array<String^>^ newTTYPES = gcnew array<String^>(TFIELDS - 1);
	array<String^>^ newTFORMS = gcnew array<String^>(TFIELDS - 1);
	array<String^>^ newTUNITS = gcnew array<String^>(TFIELDS - 1);
	array<int>^ newTBYTES = gcnew array<int>(TFIELDS - 1);
	array<int>^ newTINSTANCES = gcnew array<int>(TFIELDS - 1);
	array<TypeCode>^ newTCODES = gcnew array<TypeCode>(TFIELDS - 1);

	int c = 0;
	for (int i = 0; i < TFIELDS; i++)
		if (i == extensionentry_index)
			continue;
		else
		{
			TypeCode code;
			int rank;
			newEntryDataObjs[c] = this->GetTTYPEEntry(TTYPES[i], code, rank);
			newTTYPES[c] = TTYPES[i];
			newTFORMS[c] = TFORMS[i];
			newTUNITS[c] = TUNITS[i];
			newTBYTES[c] = TBYTES[i];
			newTINSTANCES[c] = TINSTANCES[i];
			newTCODES[c] = TCODES[i];
			c++;
		}

	BINTABLE = MAKEBINTABLEBYTEARRAY(newEntryDataObjs);
	TFIELDS--;
	TTYPES = newTTYPES;
	TFORMS = newTFORMS;
	TUNITS = newTUNITS;
	TBYTES = newTBYTES;
	TINSTANCES = newTINSTANCES;
	TCODES = newTCODES;
}

void JPFITS::FITSBinTable::SetTTYPEEntries(array<String^>^ ttypeEntryLabels, array<String^>^ entryUnits, array<Object^>^ entryArrays)
{
	bool equalnaxis2 = true;
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
		}
	if (!equalnaxis2)
	{
		throw gcnew Exception("Error: all entry column heights, NAXIS2s, are not equal.");
		return;
	}

	BINTABLE = MAKEBINTABLEBYTEARRAY(entryArrays);

	TFIELDS = entryArrays->Length;
	TTYPES = ttypeEntryLabels;
	TUNITS = entryUnits;
	TCODES = gcnew array<TypeCode>(entryArrays->Length);
	TINSTANCES = gcnew array<int>(entryArrays->Length);
	TFORMS = gcnew array<String^>(entryArrays->Length);
	TBYTES = gcnew array<int>(entryArrays->Length);
	   
	for (int i = 0; i < entryArrays->Length; i++)
	{
		TCODES[i] = Type::GetTypeCode((((Array^)entryArrays[i])->GetType())->GetElementType());
		if (((Array^)entryArrays[i])->Rank == 1)
			TINSTANCES[i] = 1;
		else
			TINSTANCES[i] = ((Array^)entryArrays[i])->GetLength(0);
		TFORMS[i] = TYPECODETFORM(TCODES[i]);
		TBYTES[i] = TYPECODETONBYTES(TCODES[i]) * TINSTANCES[i];
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
}

void JPFITS::FITSBinTable::AddTTYPEEntry(String^ ttypeEntryLabel, bool replaceIfExists, String^ entryUnits, Object^ entryArray)
{
	int extensionentry_index = -1;
	if (TTYPES != nullptr)
		for (int i = 0; i < TTYPES->Length; i++)
			if (TTYPES[i] == ttypeEntryLabel)
			{
				extensionentry_index = i;
				break;
			}

	if (extensionentry_index != -1 && !replaceIfExists)
	{
		throw gcnew Exception("Extension Entry TTYPE Label '" + ttypeEntryLabel + "' already exists, but was told to not overwrite it.");
		return;
	}

	if (extensionentry_index != -1)//then remove it
		this->RemoveTTYPEEntry(ttypeEntryLabel);
	else
		extensionentry_index = TFIELDS;//then put the entry at the last column of the table...NB this is a zero-based index...TFIELDS will increment by one below

	this->TFIELDS++;
	array<Object^>^ newEntryDataObjs = gcnew array<Object^>(TFIELDS);
	array<String^>^ newTTYPES = gcnew array<String^>(TFIELDS);
	array<String^>^ newTFORMS = gcnew array<String^>(TFIELDS);
	array<String^>^ newTUNITS = gcnew array<String^>(TFIELDS);
	array<int>^ newTBYTES = gcnew array<int>(TFIELDS);
	array<int>^ newTINSTANCES = gcnew array<int>(TFIELDS);
	array<TypeCode>^ newTCODES = gcnew array<TypeCode>(TFIELDS);

	int c = 0;
	for (int i = 0; i < TFIELDS; i++)
		if (i == extensionentry_index)
		{
			int rank = ((Array^)entryArray)->Rank;
			int instances = 1;
			if (rank == 2)
				instances = ((Array^)entryArray)->GetLength(0);

			newEntryDataObjs[i] = entryArray;
			newTTYPES[i] = ttypeEntryLabel;
			newTCODES[i] = Type::GetTypeCode((((Array^)entryArray)->GetType())->GetElementType());
			newTFORMS[i] = TYPECODETFORM(newTCODES[i]);
			newTUNITS[i] = entryUnits;
			newTBYTES[i] = TYPECODETONBYTES(newTCODES[i]) * instances;
			newTINSTANCES[i] = instances;
		}
		else
		{
			TypeCode code;
			int rank;
			newEntryDataObjs[i] = this->GetTTYPEEntry(TTYPES[c], code, rank);
			newTTYPES[i] = TTYPES[c];
			newTFORMS[i] = TFORMS[c];
			newTUNITS[i] = TUNITS[c];
			newTBYTES[i] = TBYTES[c];
			newTINSTANCES[i] = TINSTANCES[c];
			newTCODES[i] = TCODES[c];
			c++;
		}

	BINTABLE = MAKEBINTABLEBYTEARRAY(newEntryDataObjs);
	TTYPES = newTTYPES;
	TFORMS = newTFORMS;
	TUNITS = newTUNITS;
	TBYTES = newTBYTES;
	TINSTANCES = newTINSTANCES;
	TCODES = newTCODES;

	//either it was an add to a blank table, or a replacement, or an additional, so these either need set for the first time, or updated
	BITPIX = 8;
	NAXIS = 2;
	NAXIS1 = 0;
	for (int i = 0; i < TBYTES->Length; i++)
		NAXIS1 += TBYTES[i];
	if (((Array^)entryArray)->Rank == 1)
		NAXIS2 = ((Array^)entryArray)->Length;
	else
		NAXIS2 = ((Array^)entryArray)->GetLength(1);
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
			newvals[i] = EXTRAKEYS[i];
			newcoms[i] = EXTRAKEYS[i];
		}
		newkeys[EXTRAKEYS->Length] = keyName;
		newvals[EXTRAKEYS->Length] = keyValue;
		newcoms[EXTRAKEYS->Length] = keyComment;

		EXTRAKEYS = newkeys;
		EXTRAKEYVALS = newvals;
		EXTRAKEYCOMS = newcoms;
	}
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
		throw gcnew Exception("File not formatted as FITS file, or indicates no extensions present.");
		fs->Close();
		return;
	}

	__int64 extensionstartposition, extensionendposition;
	bool exists = FITSFILEOPS::SEEKEXTENSION(fs, "BINTABLE", ExtensionName, nullptr, extensionstartposition, extensionendposition);
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
		throw gcnew Exception("File not formatted as FITS file, or indicates no extensions present.");
		fs->Close();
		return false;
	}

	__int64 start, end;
	bool exists = FITSFILEOPS::SEEKEXTENSION(fs, "BINTABLE", ExtensionName, nullptr, start, end);
	fs->Close();
	return exists;
}

void JPFITS::FITSBinTable::Write(String^ FileName, String^ ExtensionName, bool OverWriteExtensionIfExists)
{
	TypeCode code;
	int rank;
	array<Object^>^ writeobjarr = gcnew array<Object^>(this->TFIELDS);
	for (int i = 0; i < TFIELDS; i++)
		writeobjarr[i] = this->GetTTYPEEntry(TTYPES[i], code, rank);

	EXTENSIONNAME = ExtensionName;
	FILENAME = FileName;

	WRITEEXTENSION(FileName, ExtensionName, OverWriteExtensionIfExists, this->TTYPES, this->TUNITS, this->EXTRAKEYS, this->EXTRAKEYVALS, this->EXTRAKEYCOMS, writeobjarr);
}

void JPFITS::FITSBinTable::WRITEEXTENSION(String^ FileName, String^ ExtensionName, bool OverWriteExtensionIfExists, array<String^>^ ExtensionEntryLabels, array<String^>^ ExtensionEntryDataUnits, array<String^>^ ExtensionHeaderExtraKeys, array<String^>^ ExtensionHeaderExtraKeyValues, array<String^>^ ExtensionHeaderExtraKeyComments, array<Object^>^ ExtensionEntryData)
{
	if (!File::Exists(FileName))//then write a new file, otherwise check the existing file for existing table, etc.
	{
		JPFITS::FITSImage^ ff = gcnew FITSImage(FileName, true);
		ff->WriteImage(TypeCode::Double, true);
	}

	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);

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

		FITSImage^ ff = gcnew FITSImage(FileName, nullptr, true, false, false, false);
		int n = ff->GetKeyIndex("NAXIS");
		if (n == -1)
		{
			throw gcnew Exception("File not formatted as FITS file (NAXIS not present). Use a new file.");
			return;
		}
		n = Convert::ToInt32(ff->GetKeyValue("NAXIS"));
		if (n > 0)
		{
			n = ff->GetKeyIndex("NAXIS" + n.ToString());
			if (ff->GetKeyIndex("BZERO") > n)
				n = ff->GetKeyIndex("BZERO");
			if (ff->GetKeyIndex("BSCALE") > n)
				n = ff->GetKeyIndex("BSCALE");
		}
		else
			n = ff->GetKeyIndex("NAXIS");
		ff->SetKey("EXTEND", "T", "FITS file may contain extensions", true, n + 1);
		array<String^>^ HEADER = FITSFILEOPS::GETFORMATTEDIMAGEHEADER(ff->HeaderKeys, ff->HeaderKeyValues, ff->HeaderKeyComments, false);

		array<unsigned char>^ headarr = gcnew array<unsigned char>(HEADER->Length * 80);
		for (int i = 0; i < HEADER->Length; i++)
			for (int j = 0; j < 80; j++)
				headarr[i * 80 + j] = (unsigned char)HEADER[i][j];

		fs = gcnew FileStream(FileName, IO::FileMode::Create);
		fs->Write(headarr, 0, headarr->Length);
		fs->Write(primarydataarr, 0, primarydataarr->Length);
		fs->Close();

		fs = gcnew FileStream(FileName, IO::FileMode::Open);
		FITSFILEOPS::SCANPRIMARYUNIT(fs, true, nullptr, hasext);
	}

	__int64 startpos, endpos;
	bool extensionfound = FITSFILEOPS::SEEKEXTENSION(fs, "BINTABLE", ExtensionName, nullptr, startpos, endpos);
	if (extensionfound && !OverWriteExtensionIfExists)
	{
		fs->Close();
		throw gcnew Exception("ExtensionName '" + ExtensionName + "' already exists and was told to not overwrite it...");
		return;
	}
	array<unsigned char>^ arr_append;
	if (extensionfound && endpos != fs->Length)//then this was not the end of the file...get the appendage data
	{
		fs->Position = endpos;
		arr_append = gcnew array<unsigned char>(int(fs->Length - endpos));
		fs->Read(arr_append, 0, arr_append->Length);
	}
	if (extensionfound)
		fs->Position = startpos;

	//format the header for writing
	array<String^>^ header = FORMATBINARYTABLEEXTENSIONHEADER(ExtensionName, ExtensionEntryData, ExtensionEntryLabels, ExtensionEntryDataUnits, ExtensionHeaderExtraKeys, ExtensionHeaderExtraKeyValues, ExtensionHeaderExtraKeyComments);
	array<unsigned char>^ headerdata = gcnew array<unsigned char>(header->Length * 80);

	for (int i = 0; i < header->Length; i++)
		for (int j = 0; j < 80; j++)
			headerdata[i * 80 + j] = (unsigned char)header[i][j];

	array<unsigned char>^ data;
	try
	{
		data = MAKEBINTABLEBYTEARRAY(ExtensionEntryData);
	}
	catch (Exception^ e)
	{
		fs->Close();
		throw gcnew Exception(e->Message);
		return;
	}

	fs->Write(headerdata, 0, headerdata->Length);
	fs->Write(data, 0, data->Length);
	if (arr_append != nullptr)
		fs->Write(arr_append, 0, arr_append->Length);
	fs->Close();
}

array<unsigned char>^ JPFITS::FITSBinTable::MAKEBINTABLEBYTEARRAY(array<Object^>^ ExtensionEntryData)
{
	array<TypeCode>^ ExtensionEntryDataTypes = gcnew array<TypeCode>(ExtensionEntryData->Length);
	for (int i = 0; i < ExtensionEntryData->Length; i++)
	{
		if (!ExtensionEntryData[i]->GetType()->IsArray)
		{
			throw gcnew Exception("Error: Object at index '" + i + "' is not an array. Stopping write.");
			return nullptr;
		}
		ExtensionEntryDataTypes[i] = Type::GetTypeCode((((Array^)ExtensionEntryData[i])->GetType())->GetElementType());
	}

	array<int>^ ExtensionEntryDataTypeInstances = gcnew array<int>(ExtensionEntryData->Length);
	array<int>^ datanaxes2 = gcnew array<int>(ExtensionEntryData->Length);
	for (int i = 0; i < ExtensionEntryData->Length; i++)
	{
		datanaxes2[i] = ((Array^)ExtensionEntryData[i])->Length;
		int rank = ((Array^)ExtensionEntryData[i])->Rank;
		if (rank == 1)
			ExtensionEntryDataTypeInstances[i] = 1;
		else
		{
			ExtensionEntryDataTypeInstances[i] = ((Array^)ExtensionEntryData[i])->GetLength(0);
			datanaxes2[i] = ((Array^)ExtensionEntryData[i])->GetLength(1);
		}
	}

	//need to get the table width in bytes and height number of rows
	int naxis1 = 0;
	for (int i = 0; i < ExtensionEntryDataTypes->Length; i++)
		naxis1 += TYPECODETONBYTES(ExtensionEntryDataTypes[i]) * ExtensionEntryDataTypeInstances[i];
	int naxis2 = datanaxes2[0];
	for (int i = 1; i < datanaxes2->Length; i++)
		if (naxis2 != datanaxes2[i])
		{
			throw gcnew Exception("ExtensionEntryData do not all have the same number of rows, NAXES2. Cannot continue;");
			return nullptr;
		}

	int NBytesTot = int(Math::Ceiling(double(naxis1 * naxis2) / 2880.0)) * 2880;
	array<unsigned char>^ data = gcnew array<unsigned char>(NBytesTot);

	int nthread = omp_get_max_threads();
	bool parallel_top = false;
	if (naxis2 >= nthread)
		parallel_top = true;
	bool exception = false;
	TypeCode exceptiontypecode;

	//now write the table data into the array
	#pragma omp parallel for if(parallel_top)
	for (int i = 0; i < naxis2; i++)
	{
		if (exception)
			break;

		for (int j = 0; j < ExtensionEntryData->Length; j++)
		{
			int jtps = 0;
			for (int jj = 0; jj < j; jj++)
				jtps += TYPECODETONBYTES(ExtensionEntryDataTypes[jj]) * ExtensionEntryDataTypeInstances[jj];

			int cc = i * naxis1 + jtps;

			switch (ExtensionEntryDataTypes[j])
			{
				case TypeCode::Double:
				{
					array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
					for (int ii = 0; ii < ExtensionEntryDataTypeInstances[j]; ii++)
					{
						if (ExtensionEntryDataTypeInstances[j] > 1)
							dbl = BitConverter::GetBytes(((array<double, 2>^)ExtensionEntryData[j])[ii, i]);
						else
							dbl = BitConverter::GetBytes(((array<double>^)ExtensionEntryData[j])[i]);

						data[cc + ii * 8] = dbl[7];
						data[cc + ii * 8 + 1] = dbl[6];
						data[cc + ii * 8 + 2] = dbl[5];
						data[cc + ii * 8 + 3] = dbl[4];
						data[cc + ii * 8 + 4] = dbl[3];
						data[cc + ii * 8 + 5] = dbl[2];
						data[cc + ii * 8 + 6] = dbl[1];
						data[cc + ii * 8 + 7] = dbl[0];
					}
					break;
				}

				case TypeCode::Single:
				{
					array<unsigned char>^ flt = gcnew array<unsigned char>(4);
					for (int ii = 0; ii < ExtensionEntryDataTypeInstances[j]; ii++)
					{
						if (ExtensionEntryDataTypeInstances[j] > 1)
							flt = BitConverter::GetBytes(((array<float, 2>^)ExtensionEntryData[j])[ii, i]);
						else
							flt = BitConverter::GetBytes(((array<float>^)ExtensionEntryData[j])[i]);

						data[cc + ii * 4] = flt[3];
						data[cc + ii * 4 + 1] = flt[2];
						data[cc + ii * 4 + 2] = flt[1];
						data[cc + ii * 4 + 3] = flt[0];
					}
					break;
				}

				case TypeCode::Int64:
				{
					__int64 val;
					for (int ii = 0; ii < ExtensionEntryDataTypeInstances[j]; ii++)
					{
						if (ExtensionEntryDataTypeInstances[j] > 1)
							val = ((array<__int64, 2>^)ExtensionEntryData[j])[ii, i];
						else
							val = ((array<__int64>^)ExtensionEntryData[j])[i];

						data[cc + ii * 8] = ((val >> 56) & 0xff);
						data[cc + ii * 8 + 1] = ((val >> 48) & 0xff);
						data[cc + ii * 8 + 2] = ((val >> 40) & 0xff);
						data[cc + ii * 8 + 3] = ((val >> 32) & 0xff);
						data[cc + ii * 8 + 4] = ((val >> 24) & 0xff);
						data[cc + ii * 8 + 5] = ((val >> 16) & 0xff);
						data[cc + ii * 8 + 6] = ((val >> 8) & 0xff);
						data[cc + ii * 8 + 7] = (val & 0xff);
					}
					break;
				}

				case TypeCode::Int32:
				{
					__int32 val;
					for (int ii = 0; ii < ExtensionEntryDataTypeInstances[j]; ii++)
					{
						if (ExtensionEntryDataTypeInstances[j] > 1)
							val = ((array<__int32, 2>^)ExtensionEntryData[j])[ii, i];
						else
							val = ((array<__int32>^)ExtensionEntryData[j])[i];

						data[cc + ii * 4] = ((val >> 24) & 0xff);
						data[cc + ii * 4 + 1] = ((val >> 16) & 0xff);
						data[cc + ii * 4 + 2] = ((val >> 8) & 0xff);
						data[cc + ii * 4 + 3] = (val & 0xff);
					}
					break;
				}

				case TypeCode::UInt32:
				{
					__int32 val;
					for (int ii = 0; ii < ExtensionEntryDataTypeInstances[j]; ii++)
					{
						if (ExtensionEntryDataTypeInstances[j] > 1)
							val = ((array<__int32, 2>^)ExtensionEntryData[j])[ii, i] - 2147483648;
						else
							val = ((array<__int32>^)ExtensionEntryData[j])[i] - 2147483648;

						data[cc + ii * 4] = ((val >> 24) & 0xff);
						data[cc + ii * 4 + 1] = ((val >> 16) & 0xff);
						data[cc + ii * 4 + 2] = ((val >> 8) & 0xff);
						data[cc + ii * 4 + 3] = (val & 0xff);
					}
					break;
				}

				case TypeCode::Int16:
				{
					__int16 val;
					for (int ii = 0; ii < ExtensionEntryDataTypeInstances[j]; ii++)
					{
						if (ExtensionEntryDataTypeInstances[j] > 1)
							val = ((array<__int16, 2>^)ExtensionEntryData[j])[ii, i];
						else
							val = ((array<__int16>^)ExtensionEntryData[j])[i];

						data[cc + ii * 2] = ((val >> 8) & 0xff);
						data[cc + ii * 2 + 1] = (val & 0xff);
					}
					break;
				}

				case TypeCode::UInt16:
				{
					__int16 val;
					for (int ii = 0; ii < ExtensionEntryDataTypeInstances[j]; ii++)
					{
						if (ExtensionEntryDataTypeInstances[j] > 1)
							val = ((array<__int16, 2>^)ExtensionEntryData[j])[ii, i] - 32768;
						else
							val = ((array<__int16>^)ExtensionEntryData[j])[i] - 32768;

						data[cc + ii * 2] = ((val >> 8) & 0xff);
						data[cc + ii * 2 + 1] = (val & 0xff);
					}
					break;
				}

				case TypeCode::SByte:
				{
					char val;
					for (int ii = 0; ii < ExtensionEntryDataTypeInstances[j]; ii++)
					{
						if (ExtensionEntryDataTypeInstances[j] > 1)
							val = ((array<char, 2>^)ExtensionEntryData[j])[ii, i];
						else
							val = ((array<char>^)ExtensionEntryData[j])[i];

						data[cc + ii] = val;
					}
					break;
				}

				case TypeCode::Byte:
				{
					char val;
					for (int ii = 0; ii < ExtensionEntryDataTypeInstances[j]; ii++)
					{
						if (ExtensionEntryDataTypeInstances[j] > 1)
							val = ((array<char, 2>^)ExtensionEntryData[j])[ii, i] - 128;
						else
							val = ((array<char>^)ExtensionEntryData[j])[i] - 128;

						data[cc + ii] = val;
					}
					break;
				}

				case TypeCode::Boolean:
				{
					bool val;
					for (int ii = 0; ii < ExtensionEntryDataTypeInstances[j]; ii++)
					{
						if (ExtensionEntryDataTypeInstances[j] > 1)
							val = ((array<bool, 2>^)ExtensionEntryData[j])[ii, i];
						else
							val = ((array<bool>^)ExtensionEntryData[j])[i];

						data[cc + ii] = val;
					}
					break;
				}

				default:
				{
					exception = true;
					exceptiontypecode = ExtensionEntryDataTypes[j];
					break;
				}
			}
		}
	}

	if (exception)
	{
		throw gcnew Exception("Data type not recognized for writing as FITS table: '" + exceptiontypecode.ToString() + "'");
		return nullptr;
	}

	return data;
}

array<String^>^ JPFITS::FITSBinTable::FORMATBINARYTABLEEXTENSIONHEADER(String^ ExtensionName, array<Object^>^ ExtensionEntryData, array<String^>^ ExtensionEntryLabels, array<String^>^ ExtensionEntryDataUnits, array<String^>^ ExtensionHeaderExtraKeys, array<String^>^ ExtensionHeaderExtraKeyValues, array<String^>^ ExtensionHeaderExtraKeyComments)
{
	ArrayList^ hkeyslist = gcnew ArrayList();
	ArrayList^ hvalslist = gcnew ArrayList();
	ArrayList^ hcomslist = gcnew ArrayList();

	array<TypeCode>^ ExtensionEntryDataTypes = gcnew array<TypeCode>(ExtensionEntryData->Length);
	for (int i = 0; i < ExtensionEntryData->Length; i++)
		ExtensionEntryDataTypes[i] = Type::GetTypeCode((((Array^)ExtensionEntryData[i])->GetType())->GetElementType());

	array<int>^ ExtensionEntryDataTypeInstances = gcnew array<int>(ExtensionEntryData->Length);
	array<int>^ datanaxes2 = gcnew array<int>(ExtensionEntryData->Length);
	for (int i = 0; i < ExtensionEntryData->Length; i++)
	{
		datanaxes2[i] = ((Array^)ExtensionEntryData[i])->Length;
		int rank = ((Array^)ExtensionEntryData[i])->Rank;
		if (rank == 1)
			ExtensionEntryDataTypeInstances[i] = 1;
		else
		{
			ExtensionEntryDataTypeInstances[i] = ((Array^)ExtensionEntryData[i])->GetLength(0);
			datanaxes2[i] /= ExtensionEntryDataTypeInstances[i];
		}
	}

	//need to get the table width in bytes and height number of rows
	int naxis1 = 0;
	for (int i = 0; i < ExtensionEntryDataTypes->Length; i++)
		naxis1 += TYPECODETONBYTES(ExtensionEntryDataTypes[i]) * ExtensionEntryDataTypeInstances[i];
	int naxis2 = datanaxes2[0];
	for (int i = 1; i < datanaxes2->Length; i++)
		if (naxis2 != datanaxes2[i])
		{
			throw gcnew Exception("ExtensionEntryData do not all have the same number of rows, NAXES2. Cannot continue;");
			return nullptr;
		}

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
	hvalslist->Add(naxis1.ToString());
	hcomslist->Add("width of table in bytes");
	hkeyslist->Add("NAXIS2");
	hvalslist->Add(naxis2.ToString());
	hcomslist->Add("number of rows in table");
	hkeyslist->Add("PCOUNT");
	hvalslist->Add("0");
	hcomslist->Add("size of special data area");
	hkeyslist->Add("GCOUNT");
	hvalslist->Add("1");
	hcomslist->Add("one data group");
	hkeyslist->Add("TFIELDS");
	hvalslist->Add(ExtensionEntryLabels->Length.ToString());
	hcomslist->Add("number of fields in each row");
	if (ExtensionName != "")
	{
		hkeyslist->Add("EXTNAME");
		hvalslist->Add(ExtensionName);
		hcomslist->Add("name of this binary table extension");
	}

	//KEY formats
	for (int i = 0; i < ExtensionEntryLabels->Length; i++)
	{
		//TTYPE
		hkeyslist->Add("TTYPE" + (i + 1).ToString());
		hvalslist->Add(ExtensionEntryLabels[i]);
		hcomslist->Add("label for field " + (i + 1).ToString());

		//TFORM
		hkeyslist->Add("TFORM" + (i + 1).ToString());
		hvalslist->Add(ExtensionEntryDataTypeInstances[i].ToString() + TYPECODETFORM(ExtensionEntryDataTypes[i]));
		hcomslist->Add("data format of field: " + TYPECODETONBYTES(ExtensionEntryDataTypes[i]).ToString() + "-byte " + TYPECODESTRING(ExtensionEntryDataTypes[i]));

		//TZERO and TSCAL
		if (ExtensionEntryDataTypes[i] == TypeCode::Byte)
		{
			hkeyslist->Add("TZERO" + (i + 1).ToString());
			hvalslist->Add("128");
			hcomslist->Add("offset for unsigned 16-bit integers");

			hkeyslist->Add("TSCAL" + (i + 1).ToString());
			hvalslist->Add("1");
			hcomslist->Add("data are not scaled");
		}
		if (ExtensionEntryDataTypes[i] == TypeCode::UInt16)
		{
			hkeyslist->Add("TZERO" + (i + 1).ToString());
			hvalslist->Add("32768");
			hcomslist->Add("offset for unsigned 16-bit integers");

			hkeyslist->Add("TSCAL" + (i + 1).ToString());
			hvalslist->Add("1");
			hcomslist->Add("data are not scaled");
		}
		if (ExtensionEntryDataTypes[i] == TypeCode::UInt32)
		{
			hkeyslist->Add("TZERO" + (i + 1).ToString());
			hvalslist->Add("2147483648");
			hcomslist->Add("offset for unsigned 32-bit integers");

			hkeyslist->Add("TSCAL" + (i + 1).ToString());
			hvalslist->Add("1");
			hcomslist->Add("data are not scaled");
		}

		//TUNIT
		hkeyslist->Add("TUNIT" + (i + 1).ToString());
		hvalslist->Add(ExtensionEntryDataUnits[i]);
		hcomslist->Add("physical unit of field");
	}

	//EXTRAKEYS
	if (ExtensionHeaderExtraKeys != nullptr)
		if (ExtensionHeaderExtraKeys->Length > 0)
			for (int i = 0; i < ExtensionHeaderExtraKeys->Length; i++)
			{
				hkeyslist->Add(ExtensionHeaderExtraKeys[i]->ToUpper());
				hvalslist->Add(ExtensionHeaderExtraKeyValues[i]);
				hcomslist->Add(ExtensionHeaderExtraKeyComments[i]);
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
		N = ::Convert::ToInt32(tform->Substring(0, tform->Length - 1));
	instances = N;

	wchar_t f = Convert::ToChar(tform->Substring(tform->Length - 1));

	switch (f)
	{
	case 'M':
		return N *= 16;

	case 'D':
	case 'K':
	case 'C':
	case 'P':
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

	default:
		throw gcnew Exception("Unrecognized typecode: '" + typecode.ToString() + "'");
	}
}

TypeCode JPFITS::FITSBinTable::TFORMTYPECODE(String^ tform)
{
	wchar_t c = Convert::ToChar(tform->Substring(tform->Length - 1));

	switch (c)
	{
	case 'D':
		return TypeCode::Double;
	case 'K':
		return TypeCode::Int64;
	case 'E':
		return TypeCode::Single;
	case 'J':
		return TypeCode::Int32;
	case 'I':
		return TypeCode::Int16;
	case 'B':
		return TypeCode::SByte;
	case 'L':
		return TypeCode::Boolean;
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
		return 1;

	default:
		throw gcnew Exception("Unrecognized typecode: '" + typecode.ToString() + "'");
	}
}

void JPFITS::FITSBinTable::EATRAWBINTABLEHEADER(ArrayList^ header)
{
	ArrayList^ extras = gcnew ArrayList();//for possible extras

	HEADER = gcnew array<String^>(header->Count);
	String^ strheaderline;
	int ttypeindex = -1;

	for (int i = 0; i < header->Count; i++)
	{
		strheaderline = (String^)header[i];
		HEADER[i] = strheaderline;

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
		if (BITPIX == 0)
			if (strheaderline->Substring(0, 8)->Trim()->Equals("BITPIX"))
			{
				BITPIX = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				continue;
			}

		if (TFIELDS == 0)
			if (strheaderline->Substring(0, 8)->Trim()->Equals("TFIELDS"))
			{
				TFIELDS = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				TTYPES = gcnew array<String^>(TFIELDS);
				TFORMS = gcnew array<String^>(TFIELDS);
				TBYTES = gcnew array<int>(TFIELDS);
				TINSTANCES = gcnew array<int>(TFIELDS);
				TCODES = gcnew array<::TypeCode>(TFIELDS);
				TUNITS = gcnew array<String^>(TFIELDS);
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
			TBYTES[ttypeindex] = TFORMTONBYTES(TFORMS[ttypeindex], instances);//need to convert the tform to Nbytes
			TINSTANCES[ttypeindex] = instances;
			TCODES[ttypeindex] = TFORMTYPECODE(TFORMS[ttypeindex]);
			continue;
		}

		if (strheaderline->Substring(0, 8)->Trim()->Equals("TUNIT" + (ttypeindex + 1).ToString()))
		{
			int f = strheaderline->IndexOf("'");
			int l = strheaderline->LastIndexOf("'");
			TUNITS[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1)->Trim();
			continue;
		}

		//need to determine if the TypeCode here is supposed to be for signed or unsigned IF the type is an integer (8, 16 or 32 bit)
		//therefore find the TSCALE and TZERO for this entry...if they don't exist then it is signed, if they do exist
		//then it is whatever values they are, for either determined signed or unsigned
		//then set this current tcode[typeindex] to what it should be
		if (strheaderline->Substring(0, 8)->Trim()->Equals("TZERO" + (ttypeindex + 1).ToString()))
		{
			if (TCODES[ttypeindex] == TypeCode::SByte)//then get the value
				if (Convert::ToByte(strheaderline->Substring(10, 20)->Trim()) == 128)//then it is an unsigned
					TCODES[ttypeindex] = TypeCode::Byte;
			if (TCODES[ttypeindex] == TypeCode::Int16)//then get the value
				if (Convert::ToUInt16(strheaderline->Substring(10, 20)->Trim()) == 32768)//then it is an unsigned
					TCODES[ttypeindex] = TypeCode::UInt16;
			if (TCODES[ttypeindex] == TypeCode::Int32)//then get the value
				if (Convert::ToUInt32(strheaderline->Substring(10, 20)->Trim()) == 2147483648)//then it is an unsigned
					TCODES[ttypeindex] = TypeCode::UInt32;
			continue;
		}

		if (strheaderline->Substring(0, 8)->Trim()->Equals("TSCAL" + (ttypeindex + 1).ToString()))//don't need to do anything with this
			continue;

		String^ key = strheaderline->Substring(0, 8)->Trim();
		if (key->Substring(0, 1) == "T" && JPMath::IsNumeric(key->Substring(key->Length - 1)))//then likely it is some other T____N field which I haven't explicitly coded above...will need to check all essential keywords listing to do this properly
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
		String^ line = (String^)header[i];
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

