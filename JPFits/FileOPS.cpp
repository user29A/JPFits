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

bool JPFITS::FITSFILEOPS::SCANPRIMARYUNIT(FileStream^ fs, bool scanpastprimarydata, ArrayList^ header_return, bool &has_extensions)
{
	array<unsigned char>^ c = gcnew array<unsigned char>(2880);
	int naxis = -1, bitpix = -1, Nnaxisn = 1;
	array<int>^ naxisn;

	//read primary header
	bool endheader = false, FITSformat = false, extend = false;

	while (!endheader)
	{
		//read 2880 block
		fs->Read(c, 0, 2880);//stream will be placed at end of header block, so that it can begin reading data values or 2ndry header
		
		for (int i = 0; i < 36; i++)
		{
			String^ line = System::Text::Encoding::ASCII->GetString(c, i * 80, 80);

			if (header_return != nullptr)
				header_return->Add(line);

			if (!FITSformat && i == 0)
				if (line->Substring(0, 8)->Trim() == "SIMPLE")
					if (line->Substring(10, 20)->Trim() == "T")
					{
						FITSformat = true;//if it doesn't exist, then it needs to be added
						continue;
					}
					else
					{
						endheader = true;
						break;
					}
				else
				{
					endheader = true;
					break;
				}
			
			if (bitpix == -1)
				if (line->Substring(0, 8)->Trim() == "BITPIX")
				{
					bitpix = ::Convert::ToInt32(line->Substring(10, 20));
					continue;
				}

			if (naxis == -1)
				if (line->Substring(0, 8)->Trim() == "NAXIS")
				{
					naxis = ::Convert::ToInt32(line->Substring(10, 20));
					naxisn = gcnew array<int>(naxis);
					continue;
				}

			if (Nnaxisn <= naxis)
				if (line->Substring(0, 8)->Trim() == ("NAXIS" + Nnaxisn.ToString()))
				{
					naxisn[Nnaxisn - 1] = ::Convert::ToInt32(line->Substring(10, 20));
					Nnaxisn++;
					continue;
				}

			if (!extend)
				if (line->Substring(0, 8)->Trim() == "EXTEND")
					if (line->Substring(10, 20)->Trim() == "T")
					{
						extend = true;
						has_extensions = extend;
						continue;
					}

			if (line->Substring(0, 8)->Trim() == "END")
			{
				endheader = true;
				break;
			}			
		}

		if (fs->Position >= fs->Length && !endheader)
		{
			FITSformat = false;
			break;
		}
	}
	//now at end of primary header block

	if (!FITSformat)
		return false;

	//if primary header has image data, may skip past it
	if (scanpastprimarydata)
		if (naxis != 0)
		{
			__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
			for (int i = 0; i < naxisn->Length; i++)
				NBytes *= naxisn[i];
			fs->Seek(__int64(Math::Ceiling(double(NBytes) / 2880) * 2880), ::SeekOrigin::Current);
		}

	return true;
	//should now be at the 1st extension header
}

array<String^>^ JPFITS::FITSFILEOPS::GETALLEXTENSIONNAMES(String^ FileName, String^ extension_type)
{
	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	bool hasext = false;
	if (!FITSFILEOPS::SCANPRIMARYUNIT(fs, true, nullptr, hasext) || !hasext)
	{
		fs->Close();
		throw gcnew Exception("File not formatted as FITS file, or indicates no extensions present.");
		return nullptr;
	}

	array<unsigned char>^ charheaderblock = gcnew array<unsigned char>(2880);
	int naxis = 0, naxis1 = 0, naxis2 = 0, bitpix = 8;
	bool endheader = false, extensionnamefound = false, extensiontypefound = false, endfile = false, extnamekeyexists = false;
	String^ strheaderline;
	ArrayList^ namelist = gcnew ArrayList();
	String^ extname = "";

	if (fs->Position >= fs->Length)
		endfile = true;

	while (!endfile)
	{
		//reset
		extname = "";
		endheader = false;
		extnamekeyexists = false;
		extensiontypefound = false;
		naxis = 0, naxis1 = 0, naxis2 = 0;
		while (!endheader)
		{
			fs->Read(charheaderblock, 0, 2880);

			for (int i = 0; i < 36; i++)
			{
				strheaderline = System::Text::Encoding::ASCII->GetString(charheaderblock, i * 80, 80);

				if (!extensiontypefound)
					if (strheaderline->Substring(0, 8) == "XTENSION")
					{
						int f = strheaderline->IndexOf("'");
						int l = strheaderline->LastIndexOf("'");
						if (strheaderline->Substring(f + 1, l - f - 1) == extension_type)
							extensiontypefound = true;
					}

				if (!extnamekeyexists)
					if (strheaderline->Substring(0, 8) == "EXTNAME ")
					{
						extnamekeyexists = true;
						int f = strheaderline->IndexOf("'");
						int l = strheaderline->LastIndexOf("'");
						extname = strheaderline->Substring(f + 1, l - f - 1)->Trim();
					}

				if (naxis == 0)
					if (strheaderline->Substring(0, 8) == "NAXIS   ")
						naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				if (naxis1 == 0)
					if (strheaderline->Substring(0, 8) == "NAXIS1  ")
						naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				if (naxis2 == 0)
					if (strheaderline->Substring(0, 8) == "NAXIS2  ")
						naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				if (bitpix == 0)
					if (strheaderline->Substring(0, 8) == "BITPIX  ")
						bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));

				if (strheaderline->Substring(0, 8) == "END     ")//check if we're at the end of the header keys
				{
					if (extensiontypefound)
						namelist->Add(extname);
					endheader = true;
					break;
				}
			}
		}

		if (naxis != 0)
		{
			__int64 NBytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
			fs->Seek(__int64(Math::Ceiling(double(NBytes) / 2880) * 2880), ::SeekOrigin::Current);
		}

		if (fs->Position >= fs->Length)
			endfile = true;
	}

	fs->Close();

	array<String^>^ list = gcnew array<String^>(namelist->Count);
	for (int i = 0; i < namelist->Count; i++)
		list[i] = (String^)namelist[i];

	return list;
}

bool JPFITS::FITSFILEOPS::SEEKEXTENSION(FileStream^ fs, String^ extension_type, String^ extension_name, ArrayList^ header_return, __int64 &startposition, __int64 &endposition)
{
	if (fs->Position == 0)
	{
		throw gcnew Exception("SEEKEXTENSION only after SCANPRIMARYUNIT with scanpastprimarydata");
		return false;
	}

	array<unsigned char>^ charheaderblock = gcnew array<unsigned char>(2880);
	int naxis = 0, naxis1 = 0, naxis2 = 0, bitpix = 0;
	bool endheader = false, extensionnamefound = false, extensiontypefound = false, endfile = false, extnamekeyexists = false, extensionfound = false;
	String^ strheaderline;

	if (fs->Position >= fs->Length)
		endfile = true;

	while (!extensionfound && !endfile)
	{
		//reset
		startposition = fs->Position;
		endheader = false;
		extnamekeyexists = false;
		extensiontypefound = false;
		naxis = 0, naxis1 = 0, naxis2 = 0, bitpix = 0;
		if (header_return != nullptr)
			header_return->Clear();
		while (!endheader)
		{
			fs->Read(charheaderblock, 0, 2880);

			for (int i = 0; i < 36; i++)
			{
				strheaderline = System::Text::Encoding::ASCII->GetString(charheaderblock, i * 80, 80);

				if (header_return != nullptr)
					header_return->Add(strheaderline);

				if (!extensiontypefound)
					if (strheaderline->Substring(0, 8) == "XTENSION")
					{
						int f = strheaderline->IndexOf("'");
						int l = strheaderline->LastIndexOf("'");
						if (strheaderline->Substring(f + 1, l - f - 1)->Trim() == extension_type)
							extensiontypefound = true;
					}

				if (!extnamekeyexists)
					if (strheaderline->Substring(0, 8) == "EXTNAME ")
					{
						extnamekeyexists = true;
						int f = strheaderline->IndexOf("'");
						int l = strheaderline->LastIndexOf("'");
						if (extension_name == strheaderline->Substring(f + 1, l - f - 1)->Trim())
							extensionnamefound = true;
					}

				if (naxis == 0)
					if (strheaderline->Substring(0, 8) == "NAXIS   ")
						naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				if (naxis1 == 0)
					if (strheaderline->Substring(0, 8) == "NAXIS1  ")
						naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				if (naxis2 == 0)
					if (strheaderline->Substring(0, 8) == "NAXIS2  ")
						naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				if (bitpix == 0)
					if (strheaderline->Substring(0, 8) == "BITPIX  ")
						bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));

				if (strheaderline->Substring(0, 8) == "END     ")//check if we're at the end of the header keys
				{
					endheader = true;
					break;
				}
			}
		}

		if (fs->Position >= fs->Length)
			endfile = true;

		if (extensiontypefound)
			if ((extnamekeyexists && extensionnamefound) || (!extnamekeyexists && extension_name == ""))
				extensionfound = true;

		endposition = fs->Position;
		if (naxis != 0)
		{
			__int64 NBytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
			if (!extensionfound)
				fs->Seek(__int64(Math::Ceiling(double(NBytes) / 2880) * 2880), ::SeekOrigin::Current);
			else
				endposition = fs->Position + __int64(Math::Ceiling(double(NBytes) / 2880) * 2880);
		}	
	}
	
	return extensionfound;
}

array<String^>^ JPFITS::FITSFILEOPS::MAKEPRIMARYDEFFORMATTEDHEADER(bool containsExtensions)
{
	array<String^>^ keys = gcnew array<String^>(5);
	array<String^>^ vals = gcnew array<String^>(5);
	array<String^>^ coms = gcnew array<String^>(5);

	keys[0] = "SIMPLE";
	keys[1] = "BITPIX";
	keys[2] = "NAXIS";
	keys[3] = "EXTEND";
	keys[4] = "END";
	
	vals[0] = "T";
	vals[1] = "-64";
	vals[2] = "0";
	if (containsExtensions)
		vals[3] = "T";
	else
		vals[3] = "F";
	vals[4] = "";
	
	coms[0] = "Valid Fits File";
	coms[1] = "Bits per Pixel";
	coms[2] = "Number of Axes";
	if (containsExtensions)
		coms[3] = "File may contain extensions";
	else
		coms[3] = "File does not contain extensions";
	coms[4] = "";

	return GETFORMATTEDIMAGEHEADER(keys, vals, coms, false);
}

array<String^>^ JPFITS::FITSFILEOPS::GETFORMATTEDIMAGEHEADER(array<String^>^ imageHeaderKeys, array<String^>^ imageHeaderKeyValues, array<String^>^ imageHeaderKeyComments, bool isExtension)
{
	if (isExtension)
	{
		imageHeaderKeys[0] = "XTENSION";
		imageHeaderKeyValues[0] = "IMAGE";
		imageHeaderKeyComments[0] = "Image extension";

		bool pcountkey = false, gcountkey = false;
		for (int i = 0; i < imageHeaderKeys->Length; i++)
			if (!pcountkey || !gcountkey)
			{
				if (imageHeaderKeys[i]->Trim() == "PCOUNT")
					pcountkey = true;
				if (imageHeaderKeys[i]->Trim() == "GCOUNT")
					pcountkey = true;
			}
		if (!pcountkey && !gcountkey)//they would BOTH not be present if things are being done correctly...need to add them
		{
			int naxis = -1;
			for (int i = 0; i < imageHeaderKeys->Length; i++)
				if (imageHeaderKeys[i]->Trim() == "NAXIS")
				{
					naxis = Convert::ToInt32(imageHeaderKeyValues[i]->Trim());
					break;
				}
			int naxisNindex = -1;
			for (int i = 0; i < imageHeaderKeys->Length; i++)
				if (imageHeaderKeys[i]->Trim() == "NAXIS" + naxis.ToString())
				{
					naxisNindex = i;
					break;
				}
			array<String^>^ keys = gcnew array<String^>(imageHeaderKeys->Length + 2);
			array<String^>^ vals = gcnew array<String^>(imageHeaderKeys->Length + 2);
			array<String^>^ coms = gcnew array<String^>(imageHeaderKeys->Length + 2);

			for (int i = 0; i <= naxisNindex; i++)
			{
				keys[i] = imageHeaderKeys[i];
				vals[i] = imageHeaderKeyValues[i];
				coms[i] = imageHeaderKeyComments[i];
			}
			keys[naxisNindex + 1] = "PCOUNT";
			vals[naxisNindex + 1] = "0";
			coms[naxisNindex + 1] = "";
			keys[naxisNindex + 2] = "GCOUNT";
			vals[naxisNindex + 2] = "1";
			coms[naxisNindex + 2] = "";

			for (int i = naxisNindex + 3; i < keys->Length; i++)
			{
				keys[i] = imageHeaderKeys[i - 2];
				vals[i] = imageHeaderKeyValues[i - 2];
				coms[i] = imageHeaderKeyComments[i - 2];
			}

			imageHeaderKeys = keys;
			imageHeaderKeyValues = vals;
			imageHeaderKeyComments = coms;
		}
	}
	else
	{
		imageHeaderKeys[0] = "SIMPLE";
		imageHeaderKeyValues[0] = "T";
		imageHeaderKeyComments[0] = "File conforms to FITS standard.";
	}

	int NKeys = imageHeaderKeyComments->Length;
	int NCards = (NKeys - 1) / 36;
	int NBlankKeys = (NCards + 1) * 36 - NKeys;
	array<String^>^ HEADER = gcnew array<String^>((NCards + 1) * 36);

	String^ key;
	String^ value;
	String^ comment;

	for (int i = 0; i < NKeys - 1; i++)
	{
		key = imageHeaderKeys[i];
		key = key->Trim();//in case some idiot put spaces in the front or back...
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
		if (JPMath::IsNumeric(imageHeaderKeyValues[i]))//then we have a numeric key value
		{
			double val = ::Convert::ToDouble(imageHeaderKeyValues[i]);

			if (Math::Abs(val) <= 1e-5 || Math::Abs(val) >= 1e13)
				value = val.ToString("0.00###########e+00");
			else
				value = val.ToString("G");
			if (val == 0)//change the exp to integer 0
				value = "0";

			L = value->Length;
			if (L > 20)
				value = value->Substring(0, 20);
			if (L < 20)
				for (int i = 0; i < 20 - L; i++)
					value = " " + value;//pad left
		}
		else//else it must be a string or comment.
		{
			value = imageHeaderKeyValues[i];
			L = value->Length;
			if (L > 18)
			{
				value = value->Substring(0, 18);
			}
			if (L < 18)
				for (int i = 0; i < 18 - L; i++)
					value += " ";//pad right
			if (value->Trim() == "T")
				value = "                   T";
			else if (value->Trim() == "F")
				value = "                   F";
			else
				value = "'" + value + "'";
		}
		//value formatting done

		//do comment formatting...always a string
		comment = imageHeaderKeyComments[i];
		L = comment->Length;
		if (L > 48)
			comment = comment->Substring(0, 48);
		if (L < 48)
			for (int i = 0; i < 48 - L; i++)
				comment += " ";//pad right
		comment = " /" + comment;//comment formatting done

		//check for COMMENT and reconfigure key line...
		if (key == "COMMENT ")
		{
			comment = imageHeaderKeyValues[i] + imageHeaderKeyComments[i];
			value = "";
			L = comment->Length;
			if (L > 72)
				comment = comment->Substring(0, 72);
			if (L < 72)
				for (int i = 0; i < 72 - L; i++)
					comment += " ";//pad right
		}

		HEADER[i] = key + value + comment;
	}

	HEADER[NKeys - 1] = "END                                                                             ";

	for (int i = 0; i < NBlankKeys; i++)
		HEADER[NKeys + i] = "                                                                                ";

	return HEADER;
}

