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

	//if primary header has image data, may skip past it, otherwise stay at primary data start
	if (scanpastprimarydata)
		if (naxis != 0)
		{
			__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
			for (int i = 0; i < naxisn->Length; i++)
				NBytes *= naxisn[i];
			fs->Seek(__int64(Math::Ceiling(double(NBytes) / 2880) * 2880), ::SeekOrigin::Current);//should now be at the 1st extension header
		}

	return true;
}

array<String^>^ JPFITS::FITSFILEOPS::GETALLEXTENSIONNAMES(String^ FileName, String^ extension_type)
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
		return nullptr;
	}

	array<unsigned char>^ charheaderblock = gcnew array<unsigned char>(2880);
	int naxis = 0, naxis1 = 0, naxis2 = 0, bitpix = 0;
	__int64 pcount = -1;
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
		naxis = 0, naxis1 = 0, naxis2 = 0, pcount = -1, bitpix = 0;
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
						if (strheaderline->Substring(f + 1, l - f - 1)->Trim() == extension_type)
							extensiontypefound = true;
						continue;
					}
				if (!extnamekeyexists)
					if (strheaderline->Substring(0, 8) == "EXTNAME ")
					{
						extnamekeyexists = true;
						int f = strheaderline->IndexOf("'");
						int l = strheaderline->LastIndexOf("'");
						extname = strheaderline->Substring(f + 1, l - f - 1)->Trim();
						continue;
					}
				if (naxis == 0)
					if (strheaderline->Substring(0, 8)->Trim()->Equals("NAXIS"))
					{
						naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
						continue;
					}
				if (naxis1 == 0)
					if (strheaderline->Substring(0, 8)->Trim()->Equals("NAXIS1"))
					{
						naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
						continue;
					}
				if (naxis2 == 0)
					if (strheaderline->Substring(0, 8)->Trim()->Equals("NAXIS2"))
					{
						naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
						continue;
					}
				if (bitpix == 0)
					if (strheaderline->Substring(0, 8)->Trim()->Equals("BITPIX"))
					{
						bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
						continue;
					}
				if (pcount == -1)
					if (strheaderline->Substring(0, 8)->Trim()->Equals("PCOUNT"))
					{
						pcount = ::Convert::ToInt64(strheaderline->Substring(10, 20));
						continue;
					}
				if (strheaderline->Substring(0, 8) == "END     ")//check if we're at the end of the header keys
				{
					if (extensiontypefound)
						namelist->Add(extname);
					endheader = true;
					break;
				}
			}
		}

		__int64 TableBytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
		fs->Seek(__int64(Math::Ceiling(double(TableBytes + pcount) / 2880) * 2880), ::SeekOrigin::Current);

		if (fs->Position >= fs->Length)
			endfile = true;
	}

	fs->Close();

	array<String^>^ list = gcnew array<String^>(namelist->Count);
	for (int i = 0; i < namelist->Count; i++)
		list[i] = (String^)namelist[i];

	return list;
}

bool JPFITS::FITSFILEOPS::SEEKEXTENSION(FileStream^ fs, String^ extension_type, int extension_number, ArrayList^ header_return, __int64& extensionStartPosition, __int64& extensionEndPosition, __int64& tableEndPosition, __int64& pcount, __int64& theap)
{
	return SEEKEXTENSION(fs, extension_type, "_FINDEXTNUM_" + extension_number.ToString(), header_return, extensionStartPosition, extensionEndPosition, tableEndPosition, pcount, theap);
}

bool JPFITS::FITSFILEOPS::SEEKEXTENSION(FileStream^ fs, String^ extension_type, String^ extension_name, ArrayList^ header_return, __int64 &extensionStartPosition, __int64 &extensionEndPosition, __int64 &tableEndPosition, __int64 &pcount, __int64 &theap)
{
	if (fs->Position == 0)
	{
		throw gcnew Exception("SEEKEXTENSION only after SCANPRIMARYUNIT with scanpastprimarydata");
		return false;
	}

	array<unsigned char>^ charheaderblock = gcnew array<unsigned char>(2880);
	int naxis = 0, naxis1 = 0, naxis2 = 0, bitpix = 0, extnum = -1, seekedextensionnum = 0;
	bool endheader = false, extensionnamefound = false, extensiontypefound = false, endfile = false, extnamekeyexists = false, extensionfound = false;
	String^ strheaderline;

	bool findextnum = extension_name->Contains("_FINDEXTNUM_");
	if (findextnum)
		extnum = Convert::ToInt32(extension_name->Substring(extension_name->LastIndexOf("_") + 1));

	if (fs->Position >= fs->Length)
		endfile = true;

	while (!extensionfound && !endfile)
	{
		//reset
		extensionStartPosition = fs->Position;
		endheader = false;
		extnamekeyexists = false;
		extensiontypefound = false;
		naxis = 0, naxis1 = 0, naxis2 = 0, bitpix = 0, pcount = -1, theap = -1;		
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
						{
							extensiontypefound = true;
							seekedextensionnum++;
						}
						continue;
					}
				if (!extnamekeyexists)
					if (strheaderline->Substring(0, 8) == "EXTNAME ")
					{
						extnamekeyexists = true;
						int f = strheaderline->IndexOf("'");
						int l = strheaderline->LastIndexOf("'");
						if (extension_name == strheaderline->Substring(f + 1, l - f - 1)->Trim())
							extensionnamefound = true;
						continue;
					}
				if (bitpix == 0)
					if (strheaderline->Substring(0, 8)->Trim()->Equals("BITPIX"))
					{
						bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
						continue;
					}
				if (naxis == 0)
					if (strheaderline->Substring(0, 8)->Trim()->Equals("NAXIS"))
					{
						naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
						continue;
					}
				if (naxis1 == 0)
					if (strheaderline->Substring(0, 8)->Trim()->Equals("NAXIS1"))
					{
						naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
						continue;
					}
				if (naxis2 == 0)
					if (strheaderline->Substring(0, 8)->Trim()->Equals("NAXIS2"))
					{
						naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
						theap = naxis1 * naxis2;
						continue;
					}				
				if (pcount == -1)
					if (strheaderline->Substring(0, 8)->Trim()->Equals("PCOUNT"))
					{
						pcount = ::Convert::ToInt64(strheaderline->Substring(10, 20));
						continue;
					}
				if (theap == naxis1 * naxis2)
					if (strheaderline->Substring(0, 8)->Trim()->Equals("THEAP"))
					{
						theap = ::Convert::ToInt64(strheaderline->Substring(10, 20));
						continue;
					}
				if (strheaderline->Substring(0, 8)->Trim()->Equals("END"))
				{
					endheader = true;
					break;
				}
			}
		}		

		if (extensiontypefound)
			if ((extnamekeyexists && extensionnamefound && !findextnum) || (!extnamekeyexists && extension_name == "" && !findextnum) || (findextnum && extnum == seekedextensionnum))
				extensionfound = true;

		__int64 TableBytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
		if (!extensionfound)
		{
			fs->Position += __int64(Math::Ceiling(double(TableBytes + pcount) / 2880) * 2880);
			if (fs->Position >= fs->Length)
				endfile = true;
		}
		else
		{
			extensionEndPosition = fs->Position + __int64(Math::Ceiling(double(TableBytes + pcount) / 2880) * 2880);
			tableEndPosition = fs->Position + TableBytes;
		}
	}
	
	return extensionfound;
}

