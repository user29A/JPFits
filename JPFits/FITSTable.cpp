///*Copyright 2021 Joseph Edwin Postma
//
//Contact email: joepostma@live.ca
//
//This file is part of JPFITS.
//
//JPFITS is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//JPFITS is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//GNU General Public License for more details.
//
//See http://www.gnu.org/licenses/. */
//
#include "stdafx.h"
//#include "JPFITS.h"
//
//using namespace JPFITS;
//
//JPFITS::FITSTable::FITSTable(String^ fileName, String^ extensionName)
//{
//	FileStream^ fs = gcnew FileStream(fileName, IO::FileMode::Open);
//	bool hasext = false;
//	if (!FITSFILEOPS::SCANPRIMARYUNIT(fs, true, nullptr, hasext) || !hasext)
//	{
//		fs->Close();
//		throw gcnew Exception("File not formatted as FITS file, or indicates no extensions present.");
//		return;
//	}
//
//	ArrayList^ header = gcnew ArrayList();
//	__int64 extensionstartposition, extensionendposition;
//
//	if (!FITSFILEOPS::SEEKEXTENSION(fs, "TABLE", extensionName, header, extensionstartposition, extensionendposition))
//	{
//		fs->Close();
//		throw gcnew Exception("Could not find TABLE with name '" + extensionName + "'");
//		return;
//	}
//
//	/*FILENAME = fileName;
//	EXTENSIONNAME = extensionName;
//
//	EXTENSIONPOSITIONSTART = extensionstartposition;
//	EXTENSIONPOSITIONDATA = fs->Position;
//	EXTENSIONPOSITIONEND = extensionendposition;
//
//	EATRAWBINTABLEHEADER(header);
//
//	BINTABLE = gcnew array<unsigned char>(int(EXTENSIONPOSITIONEND - EXTENSIONPOSITIONDATA));
//	fs->Read(BINTABLE, 0, BINTABLE->Length);*/
//
//	fs->Close();
//}
//
//void JPFITS::FITSTable::EATRAWTABLEHEADER(ArrayList^ header)
//{
//	//HEADER = gcnew array<String^>(header->Count);
//	//String^ strheaderline;
//	//int ttypeindex = -1;
//
//	//for (int i = 0; i < header->Count; i++)
//	//{
//	//	strheaderline = (String^)header[i];
//	//	HEADER[i] = strheaderline;
//	//	//do keys, values, and comments???????
//
//	//	if (NAXIS == 0)
//	//		if (strheaderline->Substring(0, 8) == "NAXIS   ")
//	//		{
//	//			NAXIS = ::Convert::ToInt32(strheaderline->Substring(10, 20));
//	//			continue;
//	//		}
//	//	if (NAXIS1 == 0)
//	//		if (strheaderline->Substring(0, 8) == "NAXIS1  ")
//	//		{
//	//			NAXIS1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
//	//			continue;
//	//		}
//	//	if (NAXIS2 == 0)
//	//		if (strheaderline->Substring(0, 8) == "NAXIS2  ")
//	//		{
//	//			NAXIS2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
//	//			continue;
//	//		}
//	//	if (BITPIX == 0)
//	//		if (strheaderline->Substring(0, 8) == "BITPIX  ")
//	//		{
//	//			BITPIX = ::Convert::ToInt32(strheaderline->Substring(10, 20));
//	//			continue;
//	//		}
//
//	//	if (TFIELDS == 0)
//	//		if (strheaderline->Substring(0, 8) == "TFIELDS ")
//	//		{
//	//			TFIELDS = ::Convert::ToInt32(strheaderline->Substring(10, 20));
//	//			TTYPES = gcnew array<String^>(TFIELDS);
//	//			TFORMS = gcnew array<String^>(TFIELDS);
//	//			TBYTES = gcnew array<int>(TFIELDS);
//	//			TREPEATS = gcnew array<int>(TFIELDS);
//	//			TCODES = gcnew array<::TypeCode>(TFIELDS);
//	//			TUNITS = gcnew array<String^>(TFIELDS);
//	//			continue;
//	//		}
//
//	//	if (strheaderline->Substring(0, 8)->Contains("TTYPE" + (ttypeindex + 2).ToString()) || strheaderline->Substring(0, 8)->Contains("TFORM" + (ttypeindex + 2).ToString()))
//	//		ttypeindex++;
//
//	//	if (strheaderline->Substring(0, 8)->Contains("TTYPE"))
//	//	{
//	//		int f = strheaderline->IndexOf("'");
//	//		int l = strheaderline->LastIndexOf("'");
//	//		TTYPES[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1)->Trim();
//	//		continue;
//	//	}
//
//	//	if (strheaderline->Substring(0, 8)->Contains("TFORM"))
//	//	{
//	//		int f = strheaderline->IndexOf("'");
//	//		int l = strheaderline->LastIndexOf("'");
//	//		TFORMS[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1)->Trim();
//	//		int instances = 1;
//	//		TBYTES[ttypeindex] = TFORMTONBYTES(TFORMS[ttypeindex], instances);//need to convert the tform to Nbytes
//	//		TREPEATS[ttypeindex] = instances;
//	//		TCODES[ttypeindex] = TFORMTYPECODE(TFORMS[ttypeindex]);
//	//		continue;
//	//	}
//
//	//	if (strheaderline->Substring(0, 8)->Contains("TUNIT"))
//	//	{
//	//		int f = strheaderline->IndexOf("'");
//	//		int l = strheaderline->LastIndexOf("'");
//	//		TUNITS[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1)->Trim();
//	//		continue;
//	//	}
//
//	//	//need to determine if the TypeCode here is supposed to be for signed or unsigned IF the type is an integer (8, 16 or 32 bit)
//	//	//therefore find the TSCALE and TZERO for this entry...if they don't exist then it is signed, if they do exist
//	//	//then it is whatever values they are, for either determined signed or unsigned
//	//	//then set this current tcode[typeindex] to what it should be
//	//	if (strheaderline->Substring(0, 8)->Contains("TZERO"))
//	//	{
//	//		if (TCODES[ttypeindex] == TypeCode::SByte)//then get the value
//	//			if (Convert::ToByte(strheaderline->Substring(10, 20)->Trim()) == 128)//then it is an unsigned
//	//				TCODES[ttypeindex] = TypeCode::Byte;
//	//		if (TCODES[ttypeindex] == TypeCode::Int16)//then get the value
//	//			if (Convert::ToUInt16(strheaderline->Substring(10, 20)->Trim()) == 32768)//then it is an unsigned
//	//				TCODES[ttypeindex] = TypeCode::UInt16;
//	//		if (TCODES[ttypeindex] == TypeCode::Int32)//then get the value
//	//			if (Convert::ToUInt32(strheaderline->Substring(10, 20)->Trim()) == 2147483648)//then it is an unsigned
//	//				TCODES[ttypeindex] = TypeCode::UInt32;
//	//	}
//	//}
//}
//
