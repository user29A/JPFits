#include "stdafx.h"
#include "JPFITS.h"

using namespace JPFITS;

array<String^>^ JPFITS::FITSBinTable::GetExtensionHeader(String^ FileName, String^ ExtensionName)
{
	array<String^>^ header;
	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	array<unsigned char>^ charheaderblock;
	int naxis, naxis1, naxis2, bitpix, Nnaxisn = 1;
	array<int>^ naxisn;

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (strheaderline->Substring(0, 8) == "NAXIS   ")
		{
			naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			naxisn = gcnew array<int>(naxis);
		}
		if (strheaderline->Substring(0, 8)->Trim() == "NAXIS" + Nnaxisn.ToString() && naxis >= Nnaxisn)
		{
			naxisn[Nnaxisn - 1] = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			Nnaxisn++;
		}
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	//now at end of primary header; if primary header has image data, must skip past it
	if (naxis != 0)
	{
		__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
		for (int i = 0; i < naxisn->Length; i++)
			NBytes *= naxisn[i];
		__int64 rem;
		::Math::DivRem(NBytes, 2880, rem);
		fs->Seek(NBytes + (2880 - rem), ::SeekOrigin::Current);
	}
	//should now be at the 1st header extension
	//read it until it is at its end, and check to see if it is the extension we want

	endheader = false;//reset
	bool extensionentryfound = false;
	bool extensionfound = false;

	bool endfile = false;
	if (fs->Position >= fs->Length)
		endfile = true;

	while (!extensionfound && !endfile)
	{
		headerblocks = 0;
		headerlines = 0;
		endheader = false;//reset
		while (!endheader)
		{
			//read 80 byte line
			headerlines++;
			fs->Read(charheaderline, 0, 80);
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);
			if (strheaderline->Substring(0, 8) == "END     ")//check if we're at the end of the header keys
				endheader = true;
		}
		fs->Position -= (80 * headerlines);//go back to start of header
		headerblocks = (headerlines - 1) / 36 + 1;
		//read headernblocks 2880 byte blocks
		charheaderblock = gcnew array<unsigned char>(2880 * headerblocks);
		fs->Read(charheaderblock, 0, 2880 * headerblocks);//stream will be placed at end of header block, so that it can begin reading data values or 2ndry header later

		header = gcnew array<String^>(headerlines);

		for (int i = 0; i < headerlines; i++)
		{
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderblock, i * 80, 80);
			header[i] = strheaderline;

			if (ExtensionName == "" && strheaderline->Substring(0, 8) == "XTENSION")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (strheaderline->Substring(f + 1, l - f - 1) == "BINTABLE")
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "EXTNAME ")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (ExtensionName == strheaderline->Substring(f + 1, l - f - 1)->Trim())
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "NAXIS   ")
				naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS1  ")
				naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS2  ")
				naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "BITPIX  ")
				bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		}

		if (extensionfound == false)
		{
			//now at end of extension header, skip its data
			if (naxis != 0)
			{
				__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
				__int64 rem;
				::Math::DivRem(Nbytes, 2880, rem);
				fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
			}

			extensionentryfound = false;//re-find the table entry we want

			if (fs->Position >= fs->Length)
				endfile = true;
		}
	}

	fs->Close();

	if (!extensionfound)
		throw gcnew Exception("ExtensionName wasn't found: '" + ExtensionName + "'");

	return header;
}

int JPFITS::FITSBinTable::GetExtensionNumberOfEntries(String^ FileName, String^ ExtensionName)
{
	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	int naxis, naxis1, naxis2, bitpix, Nnaxisn = 1;
	array<int>^ naxisn;

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (strheaderline->Substring(0, 8) == "NAXIS   ")
		{
			naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			naxisn = gcnew array<int>(naxis);
		}
		if (strheaderline->Substring(0, 8)->Trim() == "NAXIS" + Nnaxisn.ToString() && naxis >= Nnaxisn)
		{
			naxisn[Nnaxisn - 1] = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			Nnaxisn++;
		}
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	//now at end of primary header; if primary header has image data, must skip past it
	if (naxis != 0)
	{
		__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
		for (int i = 0; i < naxisn->Length; i++)
			NBytes *= naxisn[i];
		__int64 rem;
		::Math::DivRem(NBytes, 2880, rem);
		fs->Seek(NBytes + (2880 - rem), ::SeekOrigin::Current);
	}
	//should now be at the 1st header extension
	//read it until it is at its end, and check to see if it is the extension we want

	int tfields;//number of field columns
	endheader = false;//reset
	bool extensionentryfound = false;
	bool extensionfound = false;

	bool endfile = false;
	if (fs->Position >= fs->Length)
		endfile = true;

	while (!extensionfound && !endfile)
	{
		headerblocks = 0;
		headerlines = 0;
		endheader = false;//reset
		while (!endheader)
		{
			headerlines++;
			fs->Read(charheaderline, 0, 80);
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

			if (ExtensionName == "" && strheaderline->Substring(0, 8) == "XTENSION")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (strheaderline->Substring(f + 1, l - f - 1) == "BINTABLE")
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "EXTNAME ")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (ExtensionName == strheaderline->Substring(f + 1, l - f - 1)->Trim())
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "NAXIS   ")
				naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS1  ")
				naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS2  ")
				naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "BITPIX  ")
				bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "TFIELDS ")
				tfields = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "END     ")
				endheader = true;
		}
		headerblocks = (headerlines - 1) / 36 + 1;
		fs->Position += (headerblocks * 36 - headerlines) * 80;

		if (extensionfound == false)
		{
			//now at end of extension header, skip its data
			if (naxis != 0)
			{
				__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
				__int64 rem;
				::Math::DivRem(Nbytes, 2880, rem);
				fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
			}

			extensionentryfound = false;//re-find the table entry we want

			if (fs->Position >= fs->Length)
				endfile = true;
		}
	}

	fs->Close();

	if (!extensionfound)
		throw gcnew Exception("ExtensionName wasn't found: '" + ExtensionName + "'");

	return tfields;
}

array<TypeCode>^ JPFITS::FITSBinTable::GetExtensionEntryDataTypes(String^ FileName, String^ ExtensionName)
{
	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	int naxis, naxis1, naxis2, bitpix, Nnaxisn = 1;
	array<int>^ naxisn;

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (strheaderline->Substring(0, 8) == "NAXIS   ")
		{
			naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			naxisn = gcnew array<int>(naxis);
		}
		if (strheaderline->Substring(0, 8)->Trim() == "NAXIS" + Nnaxisn.ToString() && naxis >= Nnaxisn)
		{
			naxisn[Nnaxisn - 1] = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			Nnaxisn++;
		}
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	//now at end of primary header; if primary header has image data, must skip past it
	if (naxis != 0)
	{
		__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
		for (int i = 0; i < naxisn->Length; i++)
			NBytes *= naxisn[i];
		__int64 rem;
		::Math::DivRem(NBytes, 2880, rem);
		fs->Seek(NBytes + (2880 - rem), ::SeekOrigin::Current);
	}
	//should now be at the 1st header extension
	//read it until it is at its end, and check to see if it is the extension we want

	int tfields;//number of field columns
	endheader = false;//reset
	bool extensionfound = false;
	array<String^>^ ttypes;
	array<String^>^ tforms;
	array<int>^ tbytes;
	array<int>^ tinstances;
	array<::TypeCode>^ tcodes;
	bool typesset = false;
	int ttypeindex = -1;

	bool endfile = false;
	if (fs->Position >= fs->Length)
		endfile = true;

	while (!extensionfound && !endfile)
	{
		headerblocks = 0;
		headerlines = 0;
		endheader = false;//reset
		while (!endheader)
		{
			headerlines++;
			fs->Read(charheaderline, 0, 80);
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

			if (ExtensionName == "" && strheaderline->Substring(0, 8) == "XTENSION")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (strheaderline->Substring(f + 1, l - f - 1) == "BINTABLE")
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "EXTNAME ")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (ExtensionName == strheaderline->Substring(f + 1, l - f - 1)->Trim())
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "NAXIS   ")
				naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS1  ")
				naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS2  ")
				naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "BITPIX  ")
				bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "END     ")
				endheader = true;

			if (strheaderline->Substring(0, 8) == "TFIELDS ")
			{
				tfields = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				ttypes = gcnew array<String^>(tfields);
				tforms = gcnew array<String^>(tfields);
				tbytes = gcnew array<int>(tfields);
				tinstances = gcnew array<int>(tfields);
				tcodes = gcnew array<::TypeCode>(tfields);
			}

			if (strheaderline->Substring(0, 8)->Contains("TTYPE" + (ttypeindex + 2).ToString()) || strheaderline->Substring(0, 8)->Contains("TFORM" + (ttypeindex + 2).ToString()))
				ttypeindex++;

			if (strheaderline->Substring(0, 8)->Contains("TTYPE"))
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				ttypes[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1);
			}

			if (strheaderline->Substring(0, 8)->Contains("TFORM"))
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				tforms[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1)->Trim();
				int instances = 1;
				tbytes[ttypeindex] = TFORMTONBYTES(tforms[ttypeindex], instances);//need to convert the tform to Nbytes
				tinstances[ttypeindex] = instances;
				tcodes[ttypeindex] = TFORMTYPECODE(tforms[ttypeindex]);
			}

			//need to determine if the TypeCode here is supposed to be for signed or unsigned IF the type is an integer (8, 16 or 32 bit)
			//therefore find the TSCALE and TZERO for this entry...if they don't exist then it is signed, if they do exist
			//then it is whatever values they are, for either determined signed or unsigned
			//then set this current tcode[typeindex] to what it should be
			if (strheaderline->Substring(0, 8)->Contains("TZERO") && tcodes[ttypeindex] == TypeCode::SByte)//then get the value
				if (Convert::ToByte(strheaderline->Substring(10, 20)->Trim()) == 128)//then it is an unsigned
					tcodes[ttypeindex] = TypeCode::Byte;
			if (strheaderline->Substring(0, 8)->Contains("TZERO") && tcodes[ttypeindex] == TypeCode::Int16)//then get the value
				if (Convert::ToUInt16(strheaderline->Substring(10, 20)->Trim()) == 32768)//then it is an unsigned
					tcodes[ttypeindex] = TypeCode::UInt16;
			if (strheaderline->Substring(0, 8)->Contains("TZERO") && tcodes[ttypeindex] == TypeCode::Int32)//then get the value
				if (Convert::ToUInt32(strheaderline->Substring(10, 20)->Trim()) == 2147483648)//then it is an unsigned
					tcodes[ttypeindex] = TypeCode::UInt32;
		}
		headerblocks = (headerlines - 1) / 36 + 1;
		fs->Position += (headerblocks * 36 - headerlines) * 80;

		if (extensionfound == false)
		{
			//now at end of extension header, skip its data
			if (naxis != 0)
			{
				__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
				__int64 rem;
				::Math::DivRem(Nbytes, 2880, rem);
				fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
			}

			endheader = false;//reset to scan next header
			typesset = false;//reset the type arrays
			ttypeindex = -1;

			if (fs->Position >= fs->Length)
				endfile = true;
		}
	}

	fs->Close();

	if (!extensionfound)
		throw gcnew Exception("ExtensionName wasn't found: '" + ExtensionName + "'");

	return tcodes;
}

array<int>^ JPFITS::FITSBinTable::GetExtensionEntryDataInstances(String^ FileName, String^ ExtensionName)
{
	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	int naxis, naxis1, naxis2, bitpix, Nnaxisn = 1;
	array<int>^ naxisn;

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (strheaderline->Substring(0, 8) == "NAXIS   ")
		{
			naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			naxisn = gcnew array<int>(naxis);
		}
		if (strheaderline->Substring(0, 8)->Trim() == "NAXIS" + Nnaxisn.ToString() && naxis >= Nnaxisn)
		{
			naxisn[Nnaxisn - 1] = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			Nnaxisn++;
		}
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	//at end of primary header; if primary header has image data, must skip past it
	if (naxis != 0)
	{
		__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
		for (int i = 0; i < naxisn->Length; i++)
			NBytes *= naxisn[i];
		__int64 rem;
		::Math::DivRem(NBytes, 2880, rem);
		fs->Seek(NBytes + (2880 - rem), ::SeekOrigin::Current);
	}
	//should now be at the 1st header extension
	//read it until it is at its end, and check to see if it is the extension we want

	int tfields;//number of field columns
	endheader = false;//reset
	bool extensionfound = false;
	array<String^>^ ttypes;
	array<String^>^ tforms;
	array<int>^ tbytes;
	array<int>^ tinstances;
	array<::TypeCode>^ tcodes;
	bool typesset = false;
	int ttypeindex = -1;

	bool endfile = false;
	if (fs->Position >= fs->Length)
		endfile = true;

	while (!extensionfound && !endfile)
	{
		headerblocks = 0;
		headerlines = 0;
		endheader = false;//reset
		while (!endheader)
		{
			headerlines++;
			fs->Read(charheaderline, 0, 80);
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

			if (ExtensionName == "" && strheaderline->Substring(0, 8) == "XTENSION")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (strheaderline->Substring(f + 1, l - f - 1) == "BINTABLE")
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "EXTNAME ")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (ExtensionName == strheaderline->Substring(f + 1, l - f - 1)->Trim())
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "NAXIS   ")
				naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS1  ")
				naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS2  ")
				naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "BITPIX  ")
				bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "END     ")
				endheader = true;

			if (strheaderline->Substring(0, 8) == "TFIELDS ")
			{
				tfields = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				ttypes = gcnew array<String^>(tfields);
				tforms = gcnew array<String^>(tfields);
				tbytes = gcnew array<int>(tfields);
				tinstances = gcnew array<int>(tfields);
				tcodes = gcnew array<::TypeCode>(tfields);
			}

			if (strheaderline->Substring(0, 8)->Contains("TTYPE" + (ttypeindex + 2).ToString()) || strheaderline->Substring(0, 8)->Contains("TFORM" + (ttypeindex + 2).ToString()))
				ttypeindex++;

			if (strheaderline->Substring(0, 8)->Contains("TTYPE"))
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				ttypes[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1);
			}

			if (strheaderline->Substring(0, 8)->Contains("TFORM"))
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				tforms[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1)->Trim();
				int instances = 1;
				tbytes[ttypeindex] = TFORMTONBYTES(tforms[ttypeindex], instances);//need to convert the tform to Nbytes
				tinstances[ttypeindex] = instances;
				tcodes[ttypeindex] = TFORMTYPECODE(tforms[ttypeindex]);
			}

			//need to determine if the TypeCode here is supposed to be for signed or unsigned IF the type is an integer (8, 16 or 32 bit)
			//therefore find the TSCALE and TZERO for this entry...if they don't exist then it is signed, if they do exist
			//then it is whatever values they are, for either determined signed or unsigned
			//then set this current tcode[typeindex] to what it should be
			if (strheaderline->Substring(0, 8)->Contains("TZERO") && tcodes[ttypeindex] == TypeCode::SByte)//then get the value
				if (Convert::ToByte(strheaderline->Substring(10, 20)->Trim()) == 128)//then it is an unsigned
					tcodes[ttypeindex] = TypeCode::Byte;
			if (strheaderline->Substring(0, 8)->Contains("TZERO") && tcodes[ttypeindex] == TypeCode::Int16)//then get the value
				if (Convert::ToUInt16(strheaderline->Substring(10, 20)->Trim()) == 32768)//then it is an unsigned
					tcodes[ttypeindex] = TypeCode::UInt16;
			if (strheaderline->Substring(0, 8)->Contains("TZERO") && tcodes[ttypeindex] == TypeCode::Int32)//then get the value
				if (Convert::ToUInt32(strheaderline->Substring(10, 20)->Trim()) == 2147483648)//then it is an unsigned
					tcodes[ttypeindex] = TypeCode::UInt32;
		}
		headerblocks = (headerlines - 1) / 36 + 1;
		fs->Position += (headerblocks * 36 - headerlines) * 80;

		if (extensionfound == false)
		{
			//now at end of extension header, skip its data
			if (naxis != 0)
			{
				__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
				__int64 rem;
				::Math::DivRem(Nbytes, 2880, rem);
				fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
			}

			endheader = false;//reset to scan next header
			typesset = false;//reset the type arrays
			ttypeindex = -1;

			if (fs->Position >= fs->Length)
				endfile = true;
		}
	}

	fs->Close();

	if (!extensionfound)
		throw gcnew Exception("ExtensionName wasn't found: '" + ExtensionName + "'");

	return tinstances;
}

array<String^>^ JPFITS::FITSBinTable::GetExtensionEntryLabels(String^ FileName, String^ ExtensionName)
{
	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	int naxis, naxis1, naxis2, bitpix, Nnaxisn = 1;
	array<int>^ naxisn;

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (strheaderline->Substring(0, 8) == "NAXIS   ")
		{
			naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			naxisn = gcnew array<int>(naxis);
		}
		if (strheaderline->Substring(0, 8)->Trim() == "NAXIS" + Nnaxisn.ToString() && naxis >= Nnaxisn)
		{
			naxisn[Nnaxisn - 1] = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			Nnaxisn++;
		}
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	//now at end of primary header; if primary header has image data, must skip past it
	if (naxis != 0)
	{
		__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
		for (int i = 0; i < naxisn->Length; i++)
			NBytes *= naxisn[i];
		__int64 rem;
		::Math::DivRem(NBytes, 2880, rem);
		fs->Seek(NBytes + (2880 - rem), ::SeekOrigin::Current);
	}
	//should now be at the 1st header extension
	//read it until it is at its end, and check to see if it is the extension we want

	int tfields;//number of field columns
	endheader = false;//reset
	bool extensionfound = false;
	int extensionentry_index = -1;
	array<String^>^ ttypes;
	bool typesset = false;
	int ttypeindex = -1;

	bool endfile = false;
	if (fs->Position >= fs->Length)
		endfile = true;

	while (!extensionfound && !endfile)
	{
		headerblocks = 0;
		headerlines = 0;
		endheader = false;//reset
		while (!endheader)
		{
			headerlines++;
			fs->Read(charheaderline, 0, 80);
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

			if (ExtensionName == "" && strheaderline->Substring(0, 8) == "XTENSION")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (strheaderline->Substring(f + 1, l - f - 1) == "BINTABLE")
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "EXTNAME ")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (ExtensionName == strheaderline->Substring(f + 1, l - f - 1)->Trim())
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "NAXIS   ")
				naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS1  ")
				naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS2  ")
				naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "BITPIX  ")
				bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "END     ")
				endheader = true;

			if (strheaderline->Substring(0, 8) == "TFIELDS ")
			{
				tfields = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				ttypes = gcnew array<String^>(tfields);
				typesset = true;
			}

			if (strheaderline->Substring(0, 8)->Contains("TTYPE" + (ttypeindex + 2).ToString()) || strheaderline->Substring(0, 8)->Contains("TFORM" + (ttypeindex + 2).ToString()))
				ttypeindex++;

			if (strheaderline->Substring(0, 8)->Contains("TTYPE"))
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				ttypes[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1)->Trim();
			}
		}
		headerblocks = (headerlines - 1) / 36 + 1;
		fs->Position += (headerblocks * 36 - headerlines) * 80;

		if (extensionfound == false)
		{
			//now at end of extension header, skip its data
			if (naxis != 0)
			{
				__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
				__int64 rem;
				::Math::DivRem(Nbytes, 2880, rem);
				fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
			}

			endheader = false;//reset to scan next header
			extensionentry_index = -1;//re-find its index
			typesset = false;//reset the type arrays
			ttypeindex = -1;

			if (fs->Position >= fs->Length)
				endfile = true;
		}
	}
	fs->Close();

	if (!extensionfound)
	{
		throw gcnew Exception("ExtensionName wasn't found: '" + ExtensionName + "'");
	}

	return ttypes;
}

array<String^>^ JPFITS::FITSBinTable::GetExtensionEntryUnits(String^ FileName, String^ ExtensionName)
{
	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	int naxis, naxis1, naxis2, bitpix, Nnaxisn = 1;
	array<int>^ naxisn;

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (strheaderline->Substring(0, 8) == "NAXIS   ")
		{
			naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			naxisn = gcnew array<int>(naxis);
		}
		if (strheaderline->Substring(0, 8)->Trim() == "NAXIS" + Nnaxisn.ToString() && naxis >= Nnaxisn)
		{
			naxisn[Nnaxisn - 1] = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			Nnaxisn++;
		}
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	//now at end of primary header; if primary header has image data, must skip past it
	if (naxis != 0)
	{
		__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
		for (int i = 0; i < naxisn->Length; i++)
			NBytes *= naxisn[i];
		__int64 rem;
		::Math::DivRem(NBytes, 2880, rem);
		fs->Seek(NBytes + (2880 - rem), ::SeekOrigin::Current);
	}
	//should now be at the 1st header extension
	//read it until it is at its end, and check to see if it is the extension we want

	int tfields;//number of field columns
	endheader = false;//reset
	bool extensionfound = false;
	int extensionentry_index = -1;
	array<String^>^ tunits;
	bool typesset = false;
	int ttypeindex = -1;

	bool endfile = false;
	if (fs->Position >= fs->Length)
		endfile = true;

	while (!extensionfound && !endfile)
	{
		headerblocks = 0;
		headerlines = 0;
		endheader = false;//reset
		while (!endheader)
		{
			headerlines++;
			fs->Read(charheaderline, 0, 80);
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

			if (ExtensionName == "" && strheaderline->Substring(0, 8) == "XTENSION")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (strheaderline->Substring(f + 1, l - f - 1) == "BINTABLE")
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "EXTNAME ")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (ExtensionName == strheaderline->Substring(f + 1, l - f - 1)->Trim())
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "NAXIS   ")
				naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS1  ")
				naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS2  ")
				naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "BITPIX  ")
				bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "END     ")
				endheader = true;

			if (strheaderline->Substring(0, 8) == "TFIELDS ")
			{
				tfields = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				tunits = gcnew array<String^>(tfields);
				typesset = true;
			}

			if (strheaderline->Substring(0, 8)->Contains("TTYPE" + (ttypeindex + 2).ToString()) || strheaderline->Substring(0, 8)->Contains("TFORM" + (ttypeindex + 2).ToString()))
				ttypeindex++;

			if (strheaderline->Substring(0, 8)->Contains("TUNIT"))
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				tunits[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1)->Trim();
			}
		}
		headerblocks = (headerlines - 1) / 36 + 1;
		fs->Position += (headerblocks * 36 - headerlines) * 80;

		if (extensionfound == false)
		{
			//now at end of extension header, skip its data
			if (naxis != 0)
			{
				__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
				__int64 rem;
				::Math::DivRem(Nbytes, 2880, rem);
				fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
			}

			endheader = false;//reset to scan next header
			extensionentry_index = -1;//re-find its index
			typesset = false;//reset the type arrays
			ttypeindex = -1;

			if (fs->Position >= fs->Length)
				endfile = true;
		}
	}
	fs->Close();

	if (!extensionfound)
	{
		throw gcnew Exception("ExtensionName wasn't found: '" + ExtensionName + "'");
	}

	return tunits;
}

array<String^>^ JPFITS::FITSBinTable::GetAllExtensionNames(String^ FileName)
{
	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	int naxis, naxis1, naxis2, bitpix, Nnaxisn = 1;
	array<int>^ naxisn;

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (strheaderline->Substring(0, 8) == "NAXIS   ")
		{
			naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			naxisn = gcnew array<int>(naxis);
		}
		if (strheaderline->Substring(0, 8)->Trim() == "NAXIS" + Nnaxisn.ToString() && naxis >= Nnaxisn)
		{
			naxisn[Nnaxisn - 1] = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			Nnaxisn++;
		}
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	//now at end of primary header; if primary header has image data, must skip past it
	if (naxis != 0)
	{
		__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
		for (int i = 0; i < naxisn->Length; i++)
			NBytes *= naxisn[i];
		__int64 rem;
		::Math::DivRem(NBytes, 2880, rem);
		fs->Seek(NBytes + (2880 - rem), ::SeekOrigin::Current);
	}
	//should now be at the 1st header extension
	//read it until it is at its end, and check to see if it is the extension we want

	endheader = false;//reset
	bool extensionfound = false;
	bool endfile = false;
	ArrayList^ namelist = gcnew ArrayList();
	String^ extname = "";
	if (fs->Position >= fs->Length)
		endfile = true;

	while (!endfile)
	{
		headerblocks = 0;
		headerlines = 0;
		endheader = false;//reset
		while (!endheader)
		{
			headerlines++;
			fs->Read(charheaderline, 0, 80);
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

			if (strheaderline->Substring(0, 8) == "XTENSION")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (strheaderline->Substring(f + 1, l - f - 1) == "BINTABLE")
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "EXTNAME ")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				extname = strheaderline->Substring(f + 1, l - f - 1)->Trim();
			}

			if (strheaderline->Substring(0, 8) == "NAXIS   ")
				naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS1  ")
				naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS2  ")
				naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "BITPIX  ")
				bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));

			if (strheaderline->Substring(0, 8) == "END     ")//check if we're at the end of the header keys
			{
				if (extensionfound)
					namelist->Add(extname);
				extname = "";
				extensionfound = false;
				endheader = true;
				break;
			}
		}
		headerblocks = (headerlines - 1) / 36 + 1;
		fs->Position += (headerblocks * 36 - headerlines) * 80;

		//now at end of extension header, skip its data
		if (naxis != 0)
		{
			__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
			__int64 rem;
			::Math::DivRem(Nbytes, 2880, rem);
			fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
		}

		endheader = false;//reset to scan next header

		if (fs->Position >= fs->Length)
			endfile = true;
	}
	fs->Close();

	array<String^>^ list = gcnew array<String^>(namelist->Count);
	for (int i = 0; i < namelist->Count; i++)
		list[i] = (String^)namelist[i];

	return list;
}

array<unsigned char>^ JPFITS::FITSBinTable::GetExtensionAsByteArray(String^ FileName, String^ ExtensionName, int& Width, int& Height)
{
	array<unsigned char>^ bytearray;

	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	int naxis, naxis1, naxis2, bitpix, Nnaxisn = 1;
	array<int>^ naxisn;

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (strheaderline->Substring(0, 8) == "NAXIS   ")
		{
			naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			naxisn = gcnew array<int>(naxis);
		}
		if (strheaderline->Substring(0, 8)->Trim() == "NAXIS" + Nnaxisn.ToString() && naxis >= Nnaxisn)
		{
			naxisn[Nnaxisn - 1] = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			Nnaxisn++;
		}
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	//now at end of primary header; if primary header has image data, must skip past it
	if (naxis != 0)
	{
		__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
		for (int i = 0; i < naxisn->Length; i++)
			NBytes *= naxisn[i];
		__int64 rem;
		::Math::DivRem(NBytes, 2880, rem);
		fs->Seek(NBytes + (2880 - rem), ::SeekOrigin::Current);
	}
	//should now be at the 1st header extension
	//read it until it is at its end, and check to see if it is the extension we want

	endheader = false;//reset
	bool extensionfound = false;

	bool endfile = false;
	if (fs->Position >= fs->Length)
		endfile = true;

	while (!extensionfound && !endfile)
	{
		headerblocks = 0;
		headerlines = 0;
		endheader = false;//reset
		while (!endheader)
		{
			headerlines++;
			fs->Read(charheaderline, 0, 80);
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

			if (ExtensionName == "" && strheaderline->Substring(0, 8) == "XTENSION")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (strheaderline->Substring(f + 1, l - f - 1) == "BINTABLE")
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "EXTNAME ")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (ExtensionName == strheaderline->Substring(f + 1, l - f - 1)->Trim())
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "NAXIS   ")
				naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
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

		if (extensionfound == false)
		{
			//now at end of extension header, skip its data
			if (naxis != 0)
			{
				__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
				__int64 rem;
				::Math::DivRem(Nbytes, 2880, rem);
				fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
			}

			endheader = false;//reset to scan next header

			if (fs->Position >= fs->Length)
				endfile = true;
		}
	}

	if (!extensionfound)
	{
		fs->Close();
		throw gcnew Exception("ExtensionName wasn't found: '" + ExtensionName + "'");
	}

	//now at end of found extension header, can read the table
	int NBytes = naxis1 * naxis2;
	Width = naxis1;
	Height = naxis2;
	bytearray = gcnew array<unsigned char>(NBytes);
	fs->Read(bytearray, 0, NBytes);
	fs->Close();

	return bytearray;
}

array<unsigned char, 2>^ JPFITS::FITSBinTable::GetExtensionEntryAsByteArray(String^ FileName, String^ ExtensionName, String^ ExtensionEntryLabel)
{
	array<unsigned char, 2>^ bytearray;

	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	int naxis, naxis1, naxis2, bitpix, Nnaxisn = 1;
	array<int>^ naxisn;

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (strheaderline->Substring(0, 8) == "NAXIS   ")
		{
			naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			naxisn = gcnew array<int>(naxis);
		}
		if (strheaderline->Substring(0, 8)->Trim() == "NAXIS" + Nnaxisn.ToString() && naxis >= Nnaxisn)
		{
			naxisn[Nnaxisn - 1] = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			Nnaxisn++;
		}
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	//now at end of primary header; if primary header has image data, must skip past it
	if (naxis != 0)
	{
		__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
		for (int i = 0; i < naxisn->Length; i++)
			NBytes *= naxisn[i];
		__int64 rem;
		::Math::DivRem(NBytes, 2880, rem);
		fs->Seek(NBytes + (2880 - rem), ::SeekOrigin::Current);
	}
	//should now be at the 1st header extension
	//read it until it is at its end, and check to see if it is the extension we want

	int tfields;//number of field columns
	endheader = false;//reset
	bool extensionentryfound = false;
	bool extensionfound = false;
	int extensionentry_index = -1;
	array<String^>^ ttypes;
	array<String^>^ tforms;
	array<int>^ tbytes;
	array<int>^ tinstances;
	array<::TypeCode>^ tcodes;
	bool typesset = false;
	int ttypeindex = -1;

	bool endfile = false;
	if (fs->Position >= fs->Length)
		endfile = true;

	while (!extensionfound && !endfile)
	{
		headerblocks = 0;
		headerlines = 0;
		endheader = false;//reset
		while (!endheader)
		{
			headerlines++;
			fs->Read(charheaderline, 0, 80);
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

			if (ExtensionName == "" && strheaderline->Substring(0, 8) == "XTENSION")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (strheaderline->Substring(f + 1, l - f - 1) == "BINTABLE")
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "EXTNAME ")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (ExtensionName == strheaderline->Substring(f + 1, l - f - 1)->Trim())
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "NAXIS   ")
				naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS1  ")
				naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS2  ")
				naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "BITPIX  ")
				bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "END     ")
				endheader = true;

			if (strheaderline->Substring(0, 8) == "TFIELDS ")
			{
				tfields = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				ttypes = gcnew array<String^>(tfields);
				tforms = gcnew array<String^>(tfields);
				tbytes = gcnew array<int>(tfields);
				tinstances = gcnew array<int>(tfields);
				tcodes = gcnew array<::TypeCode>(tfields);
			}

			if (strheaderline->Substring(0, 8)->Contains("TTYPE" + (ttypeindex + 2).ToString()) || strheaderline->Substring(0, 8)->Contains("TFORM" + (ttypeindex + 2).ToString()))
				ttypeindex++;

			if (strheaderline->Substring(0, 8)->Contains("TTYPE"))
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				ttypes[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1);
			}

			if (strheaderline->Substring(0, 8)->Contains("TFORM"))
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				tforms[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1)->Trim();
				int instances = 1;
				tbytes[ttypeindex] = TFORMTONBYTES(tforms[ttypeindex], instances);//need to convert the tform to Nbytes
				tinstances[ttypeindex] = instances;
				tcodes[ttypeindex] = TFORMTYPECODE(tforms[ttypeindex]);
			}

			if (!extensionentryfound)
			{
				if (strheaderline->Substring(0, 8)->Contains("TTYPE"))
				{
					extensionentry_index++;
					int f = strheaderline->IndexOf("'");
					int l = strheaderline->LastIndexOf("'");
					if (ExtensionEntryLabel == strheaderline->Substring(f + 1, l - f - 1)->Trim())
						extensionentryfound = true;
				}
			}
		}
		headerblocks = (headerlines - 1) / 36 + 1;
		fs->Position += (headerblocks * 36 - headerlines) * 80;

		if (extensionfound == false)
		{
			//now at end of extension header, skip its data
			if (naxis != 0)
			{
				__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
				__int64 rem;
				::Math::DivRem(Nbytes, 2880, rem);
				fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
			}

			endheader = false;//reset to scan next header
			extensionentryfound = false;//re-find the table entry we want
			extensionentry_index = -1;//re-find its index
			typesset = false;//reset the type arrays
			ttypeindex = -1;

			if (fs->Position >= fs->Length)
				endfile = true;
		}
	}

	if (!extensionfound)
	{
		fs->Close();
		throw gcnew Exception("ExtensionName wasn't found: '" + ExtensionName + "'");
	}

	if (!extensionentryfound)
	{
		fs->Close();
		throw gcnew Exception("ExtensionEntryLabel wasn't found: '" + ExtensionEntryLabel + "'");
	}

	TypeCode extensionentrycode = tcodes[extensionentry_index];
	if (extensionentrycode != TypeCode::SByte)
	{
		fs->Close();
		throw gcnew Exception("Requested entry '" + ExtensionEntryLabel + "'" + " is not of type BYTE. Use an overload instead for numeric value types.");
	}

	//now at end of found extension header, can read the table
	int NBytes = naxis1 * naxis2;
	int extensionentryinstances = tinstances[extensionentry_index];
	bytearray = gcnew array<unsigned char, 2>(extensionentryinstances, naxis2);
	array<unsigned char>^ arr = gcnew array<unsigned char>(NBytes);
	fs->Read(arr, 0, NBytes);
	fs->Close();

	int byteoffset = 0;
	for (int i = 0; i < extensionentry_index; i++)
		byteoffset += tbytes[i];
	int currentbyte;

	#pragma omp parallel for private(currentbyte)
	for (int i = 0; i < naxis2; i++)
		for (int j = 0; j < extensionentryinstances; j++)
		{
			currentbyte = byteoffset + i * naxis1 + j;
			bytearray[j, i] = arr[currentbyte];
		}

	return bytearray;
}

array<double, 2>^ JPFITS::FITSBinTable::GetExtensionEntries(String^ FileName, String^ ExtensionName, array<String^>^ ExtensionEntryLabels)
{
	int w = 0, h = 0;
	array<double>^ instance = GetExtensionEntry(FileName, ExtensionName, ExtensionEntryLabels[0], w, h);
	array<double, 2>^ result = gcnew array<double, 2>(ExtensionEntryLabels->Length, instance->Length);
	int i = 0;
	{
		for (i = 0; i < instance->Length; i++)
			result[0, i] = instance[i];

		for (i = 1; i < ExtensionEntryLabels->Length; i++)
		{
			instance = GetExtensionEntry(FileName, ExtensionName, ExtensionEntryLabels[i], w, h);
			for (int j = 0; j < instance->Length; j++)
				result[i, j] = instance[j];
		}
		return result;
	}
}

array<double>^ JPFITS::FITSBinTable::GetExtensionEntry(String^ FileName, String^ ExtensionName, String^ ExtensionEntryLabel)
{
	int w = 0, h = 0;

	return GetExtensionEntry(FileName, ExtensionName, ExtensionEntryLabel, w, h);
}

array<double>^ JPFITS::FITSBinTable::GetExtensionEntry(String^ FileName, String^ ExtensionName, String^ ExtensionEntryLabel, int& Width, int& Height)
{
	array<double>^ result;

	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	int naxis, naxis1, naxis2, bitpix, Nnaxisn = 1;
	array<int>^ naxisn;

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

		if (strheaderline->Substring(0, 8) == "NAXIS   ")
		{
			naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			naxisn = gcnew array<int>(naxis);
		}
		if (strheaderline->Substring(0, 8)->Trim() == "NAXIS" + Nnaxisn.ToString() && naxis >= Nnaxisn)
		{
			naxisn[Nnaxisn - 1] = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			Nnaxisn++;
		}
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	//at end of primary header; if primary header has image data, must skip past it
	if (naxis != 0)
	{
		__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
		for (int i = 0; i < naxisn->Length; i++)
			NBytes *= naxisn[i];
		__int64 rem;
		::Math::DivRem(NBytes, 2880, rem);
		fs->Seek(NBytes + (2880 - rem), ::SeekOrigin::Current);
	}
	//should now be at the 1st header extension
	//read it until it is at its end, and check to see if it is the extension we want

	int tfields;//number of field columns
	bool extensionentryfound = false;
	bool extensionfound = false;
	int extensionentry_index = -1;
	array<String^>^ ttypes;
	array<String^>^ tforms;
	array<int>^ tbytes;
	array<int>^ tinstances;
	array<::TypeCode>^ tcodes;
	bool typesset = false;
	int ttypeindex = -1;

	bool endfile = false;
	if (fs->Position >= fs->Length)
		endfile = true;

	while (!extensionfound && !endfile)
	{
		headerblocks = 0;
		headerlines = 0;
		endheader = false;//reset
		while (!endheader)
		{
			headerlines++;
			fs->Read(charheaderline, 0, 80);
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

			if (ExtensionName == "" && strheaderline->Substring(0, 8) == "XTENSION")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (strheaderline->Substring(f + 1, l - f - 1) == "BINTABLE")
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "EXTNAME ")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (ExtensionName == strheaderline->Substring(f + 1, l - f - 1)->Trim())
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "NAXIS   ")
				naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS1  ")
				naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS2  ")
				naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "BITPIX  ")
				bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "END     ")
				endheader = true;

			if (strheaderline->Substring(0, 8) == "TFIELDS ")
			{
				tfields = ::Convert::ToInt32(strheaderline->Substring(10, 20));
				ttypes = gcnew array<String^>(tfields);
				tforms = gcnew array<String^>(tfields);
				tbytes = gcnew array<int>(tfields);
				tinstances = gcnew array<int>(tfields);
				tcodes = gcnew array<::TypeCode>(tfields);
			}

			if (strheaderline->Substring(0, 8)->Contains("TTYPE" + (ttypeindex + 2).ToString()) || strheaderline->Substring(0, 8)->Contains("TFORM" + (ttypeindex + 2).ToString()))
				ttypeindex++;

			if (strheaderline->Substring(0, 8)->Contains("TTYPE"))
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				ttypes[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1);

				if (!extensionentryfound)
				{
					extensionentry_index++;
					if (ExtensionEntryLabel == strheaderline->Substring(f + 1, l - f - 1)->Trim())
						extensionentryfound = true;
				}
			}

			if (strheaderline->Substring(0, 8)->Contains("TFORM"))
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				tforms[ttypeindex] = strheaderline->Substring(f + 1, l - f - 1)->Trim();
				int instances = 1;
				tbytes[ttypeindex] = TFORMTONBYTES(tforms[ttypeindex], instances);//need to convert the tform to Nbytes
				tinstances[ttypeindex] = instances;
				tcodes[ttypeindex] = TFORMTYPECODE(tforms[ttypeindex]);
			}

			//need to determine if the TypeCode here is supposed to be for signed or unsigned IF the type is an integer (8, 16 or 32 bit)
			//therefore find the TSCALE and TZERO for this entry...if they don't exist then it is signed, if they do exist
			//then it is whatever values they are, for either determined signed or unsigned
			//then set this current tcode[typeindex] to what it should be
			if (strheaderline->Substring(0, 8)->Contains("TZERO") && tcodes[ttypeindex] == TypeCode::SByte)//then get the value
				if (Convert::ToByte(strheaderline->Substring(10, 20)->Trim()) == 128)//then it is an unsigned
					tcodes[ttypeindex] = TypeCode::Byte;
			if (strheaderline->Substring(0, 8)->Contains("TZERO") && tcodes[ttypeindex] == TypeCode::Int16)//then get the value
				if (Convert::ToUInt16(strheaderline->Substring(10, 20)->Trim()) == 32768)//then it is an unsigned
					tcodes[ttypeindex] = TypeCode::UInt16;
			if (strheaderline->Substring(0, 8)->Contains("TZERO") && tcodes[ttypeindex] == TypeCode::Int32)//then get the value
				if (Convert::ToUInt32(strheaderline->Substring(10, 20)->Trim()) == 2147483648)//then it is an unsigned
					tcodes[ttypeindex] = TypeCode::UInt32;
		}
		headerblocks = (headerlines - 1) / 36 + 1;
		fs->Position += (headerblocks * 36 - headerlines) * 80;

		if (extensionfound == false)
		{
			//now at end of extension header, skip its data
			if (naxis != 0)
			{
				__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
				__int64 rem;
				::Math::DivRem(Nbytes, 2880, rem);
				fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
			}

			endheader = false;//reset to scan next header
			extensionentryfound = false;//re-find the table entry we want
			extensionentry_index = -1;//re-find its index
			typesset = false;//reset the type arrays
			ttypeindex = -1;

			if (fs->Position >= fs->Length)
				endfile = true;
		}
	}

	if (!extensionfound)
	{
		fs->Close();
		throw gcnew Exception("ExtensionName wasn't found: '" + ExtensionName + "'");
	}

	if (!extensionentryfound)
	{
		fs->Close();
		throw gcnew Exception("ExtensionEntryLabel wasn't found: '" + ExtensionEntryLabel + "'");
	}

	//now at end of found extension header, can read the table
	int NBytes = naxis1 * naxis2;
	int extensionentryinstances = tinstances[extensionentry_index];
	Width = extensionentryinstances;
	Height = naxis2;
	result = gcnew array<double>(naxis2 * extensionentryinstances);
	array<unsigned char>^ arr = gcnew array<unsigned char>(NBytes);
	fs->Read(arr, 0, NBytes);
	fs->Close();

	TypeCode extensionentrycode = tcodes[extensionentry_index];
	int byteoffset = 0;
	for (int i = 0; i < extensionentry_index; i++)
		byteoffset += tbytes[i];
	int currentbyte;

	switch (extensionentrycode)
	{
		case ::TypeCode::Double:
		{
			#pragma omp parallel for private(currentbyte)
			for (int i = 0; i < naxis2; i++)
			{
				array<unsigned char>^ dbl = gcnew array<unsigned char>(8);
				for (int j = 0; j < extensionentryinstances; j++)
				{
					currentbyte = byteoffset + i * naxis1 + j * 8;
					dbl[7] = arr[currentbyte];
					dbl[6] = arr[currentbyte + 1];
					dbl[5] = arr[currentbyte + 2];
					dbl[4] = arr[currentbyte + 3];
					dbl[3] = arr[currentbyte + 4];
					dbl[2] = arr[currentbyte + 5];
					dbl[1] = arr[currentbyte + 6];
					dbl[0] = arr[currentbyte + 7];
					result[i * extensionentryinstances + j] = BitConverter::ToDouble(dbl, 0);
				}
			}
			break;
		}

		case (::TypeCode::Int64):
		{
			#pragma omp parallel for private(currentbyte)
			for (int i = 0; i < naxis2; i++)
			{
				array<unsigned char>^ i64 = gcnew array<unsigned char>(8);
				for (int j = 0; j < extensionentryinstances; j++)
				{
					currentbyte = byteoffset + i * naxis1 + j * 8;
					i64[7] = arr[currentbyte];
					i64[6] = arr[currentbyte + 1];
					i64[5] = arr[currentbyte + 2];
					i64[4] = arr[currentbyte + 3];
					i64[3] = arr[currentbyte + 4];
					i64[2] = arr[currentbyte + 5];
					i64[1] = arr[currentbyte + 6];
					i64[0] = arr[currentbyte + 7];
					result[i * extensionentryinstances + j] = (double)BitConverter::ToInt64(i64, 0);
				}
			}
			break;
		}

		case ::TypeCode::Single:
		{
			#pragma omp parallel for private(currentbyte)
			for (int i = 0; i < naxis2; i++)
			{
				array<unsigned char>^ sng = gcnew array<unsigned char>(4);
				for (int j = 0; j < extensionentryinstances; j++)
				{
					currentbyte = byteoffset + i * naxis1 + j * 4;
					sng[3] = arr[currentbyte];
					sng[2] = arr[currentbyte + 1];
					sng[1] = arr[currentbyte + 2];
					sng[0] = arr[currentbyte + 3];
					result[i * extensionentryinstances + j] = (double)BitConverter::ToSingle(sng, 0);
				}
			}
			break;
		}

		case ::TypeCode::UInt32:
		case ::TypeCode::Int32:
		{
			double bzero = 2147483648;
			if (extensionentrycode == ::TypeCode::Int32)
				bzero = 0;
			__int32 val;

			#pragma omp parallel for private(val, currentbyte)
			for (int i = 0; i < naxis2; i++)
			{
				for (int j = 0; j < extensionentryinstances; j++)
				{
					currentbyte = byteoffset + i * naxis1 + j * 4;
					val = (arr[currentbyte] << 24) | (arr[currentbyte + 1] << 16) | (arr[currentbyte + 2] << 8) | arr[currentbyte + 3];
					result[i * extensionentryinstances + j] = double(val) + bzero;
				}
			}
			break;
		}

		case ::TypeCode::UInt16:
		case ::TypeCode::Int16:
		{
			double bzero = 32768;
			if (extensionentrycode == ::TypeCode::Int16)
				bzero = 0;
			__int16 val;

			#pragma omp parallel for private(val, currentbyte)
			for (int i = 0; i < naxis2; i++)
			{
				for (int j = 0; j < extensionentryinstances; j++)
				{
					currentbyte = byteoffset + i * naxis1 + j * 2;
					val = (arr[currentbyte] << 8) | arr[currentbyte + 1];
					result[i * extensionentryinstances + j] = double(val) + bzero;
				}
			}
			break;
		}

		case ::TypeCode::Byte:
		case ::TypeCode::SByte:
		{
			double bzero = 128;
			if (extensionentrycode == ::TypeCode::SByte)
				bzero = 0;

			#pragma omp parallel for private(currentbyte)
			for (int i = 0; i < naxis2; i++)
			{
				for (int j = 0; j < extensionentryinstances; j++)
				{
					currentbyte = byteoffset + i * naxis1 + j;
					result[i * extensionentryinstances + j] = double(arr[currentbyte]) + bzero;
				}
			}
			break;
		}

		case ::TypeCode::Boolean:
		{
			#pragma omp parallel for private(currentbyte)
			for (int i = 0; i < naxis2; i++)
			{
				for (int j = 0; j < extensionentryinstances; j++)
				{
					currentbyte = byteoffset + i * naxis1 + j;
					result[i * extensionentryinstances + j] = double((bool)arr[currentbyte]);
				}
			}
			break;
		}

		default:
			throw gcnew Exception("Unrecognized TypeCode: '" + extensionentrycode.ToString() + "'");
	}

	return result;
}

void JPFITS::FITSBinTable::RemoveExtension(String^ FileName, String^ ExtensionName)
{
	if (!File::Exists(FileName))
		throw gcnew Exception("FileName '" + FileName + "'not found!");

	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	array<unsigned char>^ charheaderblock;
	int naxis, naxis1, naxis2, bitpix, Nnaxisn = 1;
	array<int>^ naxisn;

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (strheaderline->Substring(0, 8) == "NAXIS   ")
		{
			naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			naxisn = gcnew array<int>(naxis);
		}
		if (strheaderline->Substring(0, 8)->Trim() == "NAXIS" + Nnaxisn.ToString() && naxis >= Nnaxisn)
		{
			naxisn[Nnaxisn - 1] = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			Nnaxisn++;
		}
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	//now at end of primary header; if primary header has image data, must skip past it
	if (naxis != 0)
	{
		__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
		for (int i = 0; i < naxisn->Length; i++)
			NBytes *= naxisn[i];
		__int64 rem;
		::Math::DivRem(NBytes, 2880, rem);
		fs->Seek(NBytes + (2880 - rem), ::SeekOrigin::Current);
	}
	//should now be at the 1st header extension
	//read it until it is at its end, and check to see if it is the extension we want

	endheader = false;//reset
	bool extensionfound = false;

	bool endfile = false;
	if (fs->Position >= fs->Length)
		endfile = true;

	__int64 extensionstartposition = fs->Position;
	__int64 extensionendposition;

	while (!extensionfound && !endfile)
	{
		extensionstartposition = fs->Position;

		headerblocks = 0;
		headerlines = 0;
		endheader = false;//reset
		while (!endheader)
		{
			//read 80 byte line
			headerlines++;
			fs->Read(charheaderline, 0, 80);
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);
			if (strheaderline->Substring(0, 8) == "END     ")//check if we're at the end of the header keys
				endheader = true;
		}
		fs->Position -= (80 * headerlines);//go back to start of header
		headerblocks = (headerlines - 1) / 36 + 1;
		//read headernblocks 2880 byte blocks
		charheaderblock = gcnew array<unsigned char>(2880 * headerblocks);
		fs->Read(charheaderblock, 0, 2880 * headerblocks);//stream will be placed at end of header block, so that it can begin reading data values or 2ndry header later

		for (int i = 0; i < headerlines; i++)
		{
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderblock, i * 80, 80);

			if (ExtensionName == "" && strheaderline->Substring(0, 8) == "XTENSION")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (strheaderline->Substring(f + 1, l - f - 1) == "BINTABLE")
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "EXTNAME ")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (ExtensionName == strheaderline->Substring(f + 1, l - f - 1)->Trim())
					extensionfound = true;
			}

			if (strheaderline->Substring(0, 8) == "NAXIS   ")
				naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS1  ")
				naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS2  ")
				naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "BITPIX  ")
				bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		}

		//now at end of extension header, skip its data
		if (naxis != 0)
		{
			__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
			__int64 rem;
			::Math::DivRem(Nbytes, 2880, rem);
			fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
		}

		extensionendposition = fs->Position;

		if (fs->Position >= fs->Length)
			endfile = true;
	}

	if (!extensionfound)
	{
		fs->Close();
		throw gcnew Exception("ExtensionName wasn't found: '" + ExtensionName + "'");
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
	if (!File::Exists(FileName))
		throw gcnew Exception("FileName '" + FileName + "'not found!");

	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);
	array<unsigned char>^ charheaderline = gcnew array<unsigned char>(80);
	array<unsigned char>^ charheaderblock;
	int naxis, naxis1, naxis2, bitpix, Nnaxisn = 1;
	array<int>^ naxisn;

	//read primary header
	bool endheader = false;
	int headerlines = 0, headerblocks = 0;
	String^ strheaderline;
	while (!endheader)
	{
		headerlines++;
		fs->Read(charheaderline, 0, 80);
		strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);

		if (strheaderline->Substring(0, 8) == "NAXIS   ")
		{
			naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			naxisn = gcnew array<int>(naxis);
		}
		if (strheaderline->Substring(0, 8)->Trim() == "NAXIS" + Nnaxisn.ToString() && naxis >= Nnaxisn)
		{
			naxisn[Nnaxisn - 1] = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			Nnaxisn++;
		}
		if (strheaderline->Substring(0, 8) == "BITPIX  ")
			bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		if (strheaderline->Substring(0, 8) == "END     ")
			endheader = true;
	}
	headerblocks = (headerlines - 1) / 36 + 1;
	fs->Position += (headerblocks * 36 - headerlines) * 80;

	//now at end of primary header; if primary header has image data, must skip past it
	if (naxis != 0)
	{
		__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
		for (int i = 0; i < naxisn->Length; i++)
			NBytes *= naxisn[i];
		__int64 rem;
		::Math::DivRem(NBytes, 2880, rem);
		fs->Seek(NBytes + (2880 - rem), ::SeekOrigin::Current);
	}
	//should now be at the 1st header extension
	//read it until it is at its end, and check to see if it is the extension we want

	endheader = false;//reset
	bool extensionfound = false;

	bool endfile = false;
	if (fs->Position >= fs->Length)
		endfile = true;

	while (!extensionfound && !endfile)
	{
		headerblocks = 0;
		headerlines = 0;
		endheader = false;//reset
		while (!endheader)
		{
			//read 80 byte line
			headerlines++;
			fs->Read(charheaderline, 0, 80);
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderline);
			if (strheaderline->Substring(0, 8) == "END     ")//check if we're at the end of the header keys
				endheader = true;
		}
		fs->Position -= (80 * headerlines);//go back to start of header
		headerblocks = (headerlines - 1) / 36 + 1;
		//read headernblocks 2880 byte blocks
		charheaderblock = gcnew array<unsigned char>(2880 * headerblocks);
		fs->Read(charheaderblock, 0, 2880 * headerblocks);//stream will be placed at end of header block, so that it can begin reading data values or 2ndry header later

		for (int i = 0; i < headerlines; i++)
		{
			strheaderline = System::Text::Encoding::ASCII->GetString(charheaderblock, i * 80, 80);

			if (ExtensionName == "" && strheaderline->Substring(0, 8) == "XTENSION")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (strheaderline->Substring(f + 1, l - f - 1) == "BINTABLE")
					extensionfound = true;
				return extensionfound;
			}

			if (strheaderline->Substring(0, 8) == "EXTNAME ")
			{
				int f = strheaderline->IndexOf("'");
				int l = strheaderline->LastIndexOf("'");
				if (ExtensionName == strheaderline->Substring(f + 1, l - f - 1)->Trim())
					extensionfound = true;
				return extensionfound;
			}

			if (strheaderline->Substring(0, 8) == "NAXIS   ")
				naxis = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS1  ")
				naxis1 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "NAXIS2  ")
				naxis2 = ::Convert::ToInt32(strheaderline->Substring(10, 20));
			if (strheaderline->Substring(0, 8) == "BITPIX  ")
				bitpix = ::Convert::ToInt32(strheaderline->Substring(10, 20));
		}

		//now at end of extension header, skip its data
		if (naxis != 0)
		{
			__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
			__int64 rem;
			::Math::DivRem(Nbytes, 2880, rem);
			fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
		}

		if (fs->Position >= fs->Length)
			endfile = true;
	}

	return extensionfound;
}

void JPFITS::FITSBinTable::WriteExtension(String^ FileName, String^ ExtensionName, bool OverWriteExtensionIfExists, array<String^>^ ExtensionEntryLabels, array<String^>^ ExtensionEntryDataUnits, array<String^>^ ExtensionHeaderExtraKeys, array<String^>^ ExtensionHeaderExtraKeyValues, array<String^>^ ExtensionHeaderExtraKeyComments, array<Object^>^ ExtensionEntryData)
{
	if (!File::Exists(FileName))//then write a new file, otherwise check the existing file for existing table, etc.
	{
		JPFITS::FITSImage^ ff = gcnew FITSImage(FileName);
		ff->WriteFile(TypeCode::Double, true);
	}

	array<TypeCode>^ ExtensionEntryDataTypes = gcnew array<TypeCode>(ExtensionEntryData->Length);
	for (int i = 0; i < ExtensionEntryData->Length; i++)
		ExtensionEntryDataTypes[i] = Type::GetTypeCode((((Array^)ExtensionEntryData[i])->GetType())->GetElementType());

	array<int>^ ExtensionEntryDataTypeInstances = gcnew array<int>(ExtensionEntryData->Length);
	for (int i = 0; i < ExtensionEntryData->Length; i++)
	{
		int rank = ((Array^)ExtensionEntryData[i])->Rank;
		if (rank == 1)
			ExtensionEntryDataTypeInstances[i] = 1;
		else
			ExtensionEntryDataTypeInstances[i] = ((Array^)ExtensionEntryData[i])->GetLength(0);
	}

	added_extension:;

	FileStream^ fs = gcnew FileStream(FileName, IO::FileMode::Open);

	array<unsigned char>^ c = gcnew array<unsigned char>(2880);
	int naxis, naxis1, naxis2, bitpix, Nnaxisn = 1;
	array<int>^ naxisn;

	//read primary header
	bool endheader = false, extendkeyexists = false;
	int ncards = 0;
	while (endheader == false)
	{
		//read 2880 block
		ncards++;
		fs->Read(c, 0, 2880);//stream will be placed at end of header block, so that it can begin reading data values or 2ndry header
		for (int i = 0; i < 36; i++)
		{
			String^ line = System::Text::Encoding::ASCII->GetString(c, i * 80, 80);

			if (line->Substring(0, 8) == "EXTEND  ")
				extendkeyexists = true;//if it doesn't exist, then it needs to be added

			if (line->Substring(0, 8) == "NAXIS   ")
			{
				naxis = ::Convert::ToInt32(line->Substring(10, 20));
				naxisn = gcnew array<int>(naxis);
			}
			if (line->Substring(0, 8)->Trim() == "NAXIS" + Nnaxisn.ToString() && naxis >= Nnaxisn)
			{
				naxisn[Nnaxisn - 1] = ::Convert::ToInt32(line->Substring(10, 20));
				Nnaxisn++;
			}
			if (line->Substring(0, 8) == "BITPIX  ")
				bitpix = ::Convert::ToInt32(line->Substring(10, 20));
			if (line->Substring(0, 8) == "END     ")//check if we're at the end of the header keys
			{
				endheader = true;
				break;
			}
		}
	}
	//now at end of primary header

	//if the primary header didn't contain the EXTEND keyword, it now needs it
	if (!extendkeyexists)
	{
		array<unsigned char>^ arr = gcnew array<unsigned char>((int)fs->Length - ncards * 2880);
		fs->Read(arr, 0, arr->Length);
		fs->Close();

		FITSImage^ ff = gcnew FITSImage(FileName, nullptr, true, false, false, false);
		ff->AddKey("EXTEND", "T", "FITS file may contain extensions", -1);
		ff->FORMATHEADER();

		array<unsigned char>^ newarr = gcnew array<unsigned char>(arr->Length + ff->HEADER->Length * 80);

		for (int i = 0; i < ff->HEADER->Length; i++)
			for (int j = 0; j < 80; j++)
				newarr[i * 80 + j] = (unsigned char)ff->HEADER[i][j];

		Array::Copy(arr, 0, newarr, ff->HEADER->Length * 80, arr->Length);

		fs = gcnew FileStream(FileName, IO::FileMode::Create);
		fs->Write(newarr, 0, newarr->Length);
		fs->Close();

		goto added_extension;
	}

	//if primary header has image data, must skip past it
	if (naxis != 0)
	{
		__int64 NBytes = __int64(::Math::Abs(bitpix)) / 8;
		for (int i = 0; i < naxisn->Length; i++)
			NBytes *= naxisn[i];
		__int64 rem;
		::Math::DivRem(NBytes, 2880, rem);
		fs->Seek(NBytes + (2880 - rem), ::SeekOrigin::Current);
	}
	//should now be at the 1st header extension
	//read it until it is at its end, and check to see if it is the extension we want

	endheader = false;//reset
	bool extensionfound = false;
	bool endfile = false;
	__int64 extensionstartposition = fs->Position;
	__int64 extensionendposition = extensionstartposition;
	if (fs->Position >= fs->Length)
		endfile = true;

	while (!extensionfound && !endfile)
	{
		extensionstartposition = fs->Position;
		while (endheader == false)
		{
			//read 2880 block
			fs->Read(c, 0, 2880);//stream will be placed at end of header block, so that it can begin reading data values or 2ndry header
			for (int i = 0; i < 36; i++)
			{
				String^ line = System::Text::Encoding::ASCII->GetString(c, i * 80, 80);

				if (line->Substring(0, 8) == "EXTNAME ")
				{
					int f = line->IndexOf("'");
					int l = line->LastIndexOf("'");
					if (ExtensionName == line->Substring(f + 1, l - f - 1)->Trim())
						extensionfound = true;
				}
				if (line->Substring(0, 8) == "NAXIS   ")
					naxis = ::Convert::ToInt32(line->Substring(10, 20));
				if (line->Substring(0, 8) == "NAXIS1  ")
					naxis1 = ::Convert::ToInt32(line->Substring(10, 20));
				if (line->Substring(0, 8) == "NAXIS2  ")
					naxis2 = ::Convert::ToInt32(line->Substring(10, 20));
				if (line->Substring(0, 8) == "BITPIX  ")
					bitpix = ::Convert::ToInt32(line->Substring(10, 20));
				if (line->Substring(0, 8) == "END     ")//check if we're at the end of the header keys
				{
					endheader = true;
					break;
				}
			}
		}

		if (extensionfound == false)
		{
			//now at end of extension header, skip its data
			if (naxis != 0)
			{
				__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
				__int64 rem;
				::Math::DivRem(Nbytes, 2880, rem);
				fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
			}

			endheader = false;//reset to scan next header

			extensionendposition = fs->Position;

			if (fs->Position >= fs->Length)
				endfile = true;
		}
	}

	if (!OverWriteExtensionIfExists && extensionfound)
	{
		fs->Close();
		throw gcnew Exception("ExtensionName '" + ExtensionName + "' already exists and was told to not overwrite it...");
	}

	//if extension is found (already exists) then will need to write the updated extension data over top of it, and then append on whatever remained in the original file (other extensions) if any

	array<unsigned char>^ arr_append;
	if (extensionfound && !endfile)
	{
		//skip to end of this extension data cards and see if we're at the end of the file; if not, then there is data to append
		__int64 Nbytes = __int64(naxis1)*__int64(naxis2)*__int64(::Math::Abs(bitpix)) / 8;
		__int64 rem;
		::Math::DivRem(Nbytes, 2880, rem);
		fs->Seek(Nbytes + (2880 - rem), ::SeekOrigin::Current);
		if (fs->Position != fs->Length)//then this was not the end of the file...get the appendage data
		{
			arr_append = gcnew array<unsigned char>(int(fs->Length - fs->Position));
			fs->Read(arr_append, 0, arr_append->Length);
		}
	}

	//can now write the extension at the extension start position...append remaining data if needed afterwards
	if (!endfile)
		fs->Position = extensionstartposition;

	//format the header for writing
	array<String^>^ header = FORMATBINARYTABLEEXTENSIONHEADER(ExtensionName, ExtensionEntryData, ExtensionEntryLabels, ExtensionEntryDataTypes, ExtensionEntryDataTypeInstances, ExtensionEntryDataUnits, ExtensionHeaderExtraKeys, ExtensionHeaderExtraKeyValues, ExtensionHeaderExtraKeyComments);

	//now format the data array
	//need to get the table width in bytes and height number of rows
	int rownbytes = 0, nrows = 0;
	for (int i = 0; i < ExtensionEntryDataTypes->Length; i++)
		rownbytes += TYPECODETONBYTES(ExtensionEntryDataTypes[i]) * ExtensionEntryDataTypeInstances[i];
	array<int>^ dataNrows = EXTENSIONENTRYDATANROWS(ExtensionEntryData, ExtensionEntryDataTypes, ExtensionEntryDataTypeInstances);
	nrows = JPMath::Max(dataNrows, false);
	int NbytesExtension = nrows * rownbytes;//number of bytes in table
	int NbytesHead = header->Length * 80;//number of bytes in extension header...should always be multiple of 2880.
	int NbytesHeadExtension = NbytesHead + NbytesExtension;//this is the number of bytes in the header + extension data...but need to write file so multiple of 2880 bytes
	int Ncards = int(Math::Ceiling(double(NbytesHeadExtension) / 2880.0));
	int NBytesTot = Ncards * 2880;
	array<unsigned char>^ data = gcnew array<unsigned char>(NBytesTot);

	for (int i = 0; i < header->Length; i++)
		for (int j = 0; j < 80; j++)
			data[i * 80 + j] = (unsigned char)header[i][j];
	//header is now written into the data array

	int nthread = omp_get_max_threads();
	bool parallel_top = false;
	if (nrows >= nthread)
		parallel_top = true;
	bool exception = false;
	TypeCode exceptiontypecode;

	//now write the table data into the array
	#pragma omp parallel for if(parallel_top)
	for (int i = 0; i < nrows; i++)
	{
		if (exception)
			break;

		for (int j = 0; j < ExtensionEntryData->Length; j++)
		{
			int jtps = 0;
			for (int jj = 0; jj < j; jj++)
				jtps += TYPECODETONBYTES(ExtensionEntryDataTypes[jj]) * ExtensionEntryDataTypeInstances[jj];

			int cc = NbytesHead + i * rownbytes + jtps;

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
		fs->Close();
		throw gcnew Exception("Data type not recognized for writing as FITS table: '" + exceptiontypecode.ToString() + "'");
		return;
	}

	fs->Write(data, 0, NBytesTot);
	if (arr_append != nullptr)
		fs->Write(arr_append, 0, arr_append->Length);
	fs->Close();
}

void JPFITS::FITSBinTable::WriteExtension(String^ FileName, String^ ExtensionName, bool OverWriteExtensionIfExists, array<String^>^ ExtensionEntryLabels, array<TypeCode>^ ExtensionEntryDataTypes, array<String^>^ ExtensionEntryDataUnits, array<String^>^ ExtensionHeaderExtraKeys, array<String^>^ ExtensionHeaderExtraKeyValues, array<String^>^ ExtensionHeaderExtraKeyComments, array<double, 2>^ ExtensionEntryData)
{
	array<Object^>^ dataobj = gcnew array<Object^>(ExtensionEntryData->GetLength(0));
	
	/*array<int>^ ExtensionEntryDataTypeInstances = gcnew array<int>(dataobj->Length);
	for (int i = 0; i < dataobj->Length; i++)
		ExtensionEntryDataTypeInstances[i] = 1;*/

	int nthread = omp_get_max_threads();
	bool parallel_top = false;
	if (ExtensionEntryData->GetLength(0) >= nthread)
		parallel_top = true;
	bool exception = false;
	TypeCode exceptiontypecode;
	bool parallel_in = false;
	if (!parallel_top && ExtensionEntryData->GetLength(1) >= nthread)
		parallel_in = true;

	#pragma omp parallel for if(parallel_top)
	for (int i = 0; i < ExtensionEntryData->GetLength(0); i++)
	{
		if (exception)
			break;

		switch (ExtensionEntryDataTypes[i])
		{
			case TypeCode::Double:
			{
				array<double>^ datacol = gcnew array<double>(ExtensionEntryData->GetLength(1));

				#pragma omp parallel for if(parallel_in)
				for (int j = 0; j < datacol->Length; j++)
					datacol[j] = ExtensionEntryData[i, j];

				dataobj[i] = (Object^)datacol;
				break;
			}

			case TypeCode::Single:
			{
				array<float>^ datacol = gcnew array<float>(ExtensionEntryData->GetLength(1));

				#pragma omp parallel for if(parallel_in)
				for (int j = 0; j < datacol->Length; j++)
					datacol[j] = (float)ExtensionEntryData[i, j];

				dataobj[i] = (Object^)datacol;
				break;
			}

			case TypeCode::Int64:
			{
				array<__int64>^ datacol = gcnew array<__int64>(ExtensionEntryData->GetLength(1));
				
				#pragma omp parallel for if(parallel_in)
				for (int j = 0; j < datacol->Length; j++)
					datacol[j] = (__int64)ExtensionEntryData[i, j];

				dataobj[i] = (Object^)datacol;
				break;
			}

			case TypeCode::Int32:
			{
				array<__int32>^ datacol = gcnew array<__int32>(ExtensionEntryData->GetLength(1));
				
				#pragma omp parallel for if(parallel_in)
				for (int j = 0; j < datacol->Length; j++)
					datacol[j] = (__int32)ExtensionEntryData[i, j];

				dataobj[i] = (Object^)datacol;
				break;
			}

			case TypeCode::UInt32:
			{
				array<unsigned __int32>^ datacol = gcnew array<unsigned __int32>(ExtensionEntryData->GetLength(1));
				
				#pragma omp parallel for if(parallel_in)
				for (int j = 0; j < datacol->Length; j++)
					datacol[j] = (unsigned __int32)ExtensionEntryData[i, j];

				dataobj[i] = (Object^)datacol;
				break;
			}

			case TypeCode::Int16:
			{
				array<__int16>^ datacol = gcnew array<__int16>(ExtensionEntryData->GetLength(1));
				
				#pragma omp parallel for if(parallel_in)
				for (int j = 0; j < datacol->Length; j++)
					datacol[j] = (__int16)ExtensionEntryData[i, j];

				dataobj[i] = (Object^)datacol;
				break;
			}

			case TypeCode::UInt16:
			{
				array<unsigned __int16>^ datacol = gcnew array<unsigned __int16>(ExtensionEntryData->GetLength(1));
				
				#pragma omp parallel for if(parallel_in)
				for (int j = 0; j < datacol->Length; j++)
					datacol[j] = (unsigned __int16)ExtensionEntryData[i, j];

				dataobj[i] = (Object^)datacol;
				break;
			}

			case TypeCode::SByte:
			{
				array<char>^ datacol = gcnew array<char>(ExtensionEntryData->GetLength(1));
				
				#pragma omp parallel for if(parallel_in)
				for (int j = 0; j < datacol->Length; j++)
					datacol[j] = (char)ExtensionEntryData[i, j];

				dataobj[i] = (Object^)datacol;
				break;
			}

			case TypeCode::Byte:
			{
				array<unsigned char>^ datacol = gcnew array<unsigned char>(ExtensionEntryData->GetLength(1));
				
				#pragma omp parallel for if(parallel_in)
				for (int j = 0; j < datacol->Length; j++)
					datacol[j] = (unsigned char)ExtensionEntryData[i, j];

				dataobj[i] = (Object^)datacol;
				break;
			}

			case TypeCode::Boolean:
			{
				array<bool>^ datacol = gcnew array<bool>(ExtensionEntryData->GetLength(1));
				
				#pragma omp parallel for if(parallel_in)
				for (int j = 0; j < datacol->Length; j++)
					datacol[j] = (bool)ExtensionEntryData[i, j];

				dataobj[i] = (Object^)datacol;
				break;
			}

			default:
			{
				exception = true;
				exceptiontypecode = ExtensionEntryDataTypes[i];
				break;
			}
		}
	}

	if (exception)
	{
		throw gcnew Exception("Data type not recognized for writing as FITS table: '" + exceptiontypecode.ToString() + "'");
		return;
	}

	WriteExtension(FileName, ExtensionName, OverWriteExtensionIfExists, ExtensionEntryLabels, /*ExtensionEntryDataTypes, ExtensionEntryDataTypeInstances,*/ ExtensionEntryDataUnits, ExtensionHeaderExtraKeys, ExtensionHeaderExtraKeyValues, ExtensionHeaderExtraKeyComments, dataobj);
}

void JPFITS::FITSBinTable::WriteExtension(String^ FileName, String^ ExtensionName, bool OverWriteExtensionIfExists, String^ ExtensionEntryLabel, TypeCode ExtensionEntryDataType, String^ ExtensionEntryDataUnit, array<String^>^ ExtensionHeaderExtraKeys, array<String^>^ ExtensionHeaderExtraKeyValues, array<String^>^ ExtensionHeaderExtraKeyComments, array<double>^ ExtensionEntryData)
{
	array<Object^>^ dataobj = gcnew array<Object^>(1);
	array<int>^ ExtensionEntryDataTypeInstances = gcnew array<int>(1) { 1 };

	int nthread = omp_get_max_threads();
	TypeCode exceptiontypecode;
	bool parallel_in = false;
	if (ExtensionEntryData->Length >= nthread)
		parallel_in = true;
	bool exception = false;

	switch (ExtensionEntryDataType)
	{
		case TypeCode::Double:
		{
			dataobj[0] = (Object^)ExtensionEntryData;
			break;
		}

		case TypeCode::Single:
		{
			array<float>^ datacol = gcnew array<float>(ExtensionEntryData->Length);

			#pragma omp parallel for if(parallel_in)
			for (int j = 0; j < datacol->Length; j++)
				datacol[j] = (float)ExtensionEntryData[j];

			dataobj[0] = (Object^)datacol;
			break;
		}

		case TypeCode::Int64:
		{
			array<__int64>^ datacol = gcnew array<__int64>(ExtensionEntryData->Length);

			#pragma omp parallel for if(parallel_in)
			for (int j = 0; j < datacol->Length; j++)
				datacol[j] = (__int64)ExtensionEntryData[j];

			dataobj[0] = (Object^)datacol;
			break;
		}

		case TypeCode::Int32:
		{
			array<__int32>^ datacol = gcnew array<__int32>(ExtensionEntryData->Length);

			#pragma omp parallel for if(parallel_in)
			for (int j = 0; j < datacol->Length; j++)
				datacol[j] = (__int32)ExtensionEntryData[j];

			dataobj[0] = (Object^)datacol;
			break;
		}

		case TypeCode::UInt32:
		{
			array<unsigned __int32>^ datacol = gcnew array<unsigned __int32>(ExtensionEntryData->Length);

			#pragma omp parallel for if(parallel_in)
			for (int j = 0; j < datacol->Length; j++)
				datacol[j] = (unsigned __int32)ExtensionEntryData[j];

			dataobj[0] = (Object^)datacol;
			break;
		}

		case TypeCode::Int16:
		{
			array<__int16>^ datacol = gcnew array<__int16>(ExtensionEntryData->Length);

			#pragma omp parallel for if(parallel_in)
			for (int j = 0; j < datacol->Length; j++)
				datacol[j] = (__int16)ExtensionEntryData[j];

			dataobj[0] = (Object^)datacol;
			break;
		}

		case TypeCode::UInt16:
		{
			array<unsigned __int16>^ datacol = gcnew array<unsigned __int16>(ExtensionEntryData->Length);

			#pragma omp parallel for if(parallel_in)
			for (int j = 0; j < datacol->Length; j++)
				datacol[j] = (unsigned __int16)ExtensionEntryData[j];

			dataobj[0] = (Object^)datacol;
			break;
		}

		case TypeCode::SByte:
		{
			array<char>^ datacol = gcnew array<char>(ExtensionEntryData->Length);

			#pragma omp parallel for if(parallel_in)
			for (int j = 0; j < datacol->Length; j++)
				datacol[j] = (char)ExtensionEntryData[j];

			dataobj[0] = (Object^)datacol;
			break;
		}

		case TypeCode::Byte:
		{
			array<unsigned char>^ datacol = gcnew array<unsigned char>(ExtensionEntryData->Length);

			#pragma omp parallel for if(parallel_in)
			for (int j = 0; j < datacol->Length; j++)
				datacol[j] = (unsigned char)ExtensionEntryData[j];

			dataobj[0] = (Object^)datacol;
			break;
		}

		case TypeCode::Boolean:
		{
			array<bool>^ datacol = gcnew array<bool>(ExtensionEntryData->Length);

			#pragma omp parallel for if(parallel_in)
			for (int j = 0; j < datacol->Length; j++)
				datacol[j] = (bool)ExtensionEntryData[j];

			dataobj[0] = (Object^)datacol;
			break;
		}

		default:
		{
			exception = true;
			exceptiontypecode = ExtensionEntryDataType;
			break;
		}
	}

	if (exception)
	{
		throw gcnew Exception("Data type not recognized for writing as FITS table: '" + exceptiontypecode.ToString() + "'");
		return;
	}

	array<String^>^ ExtensionEntryLabels = gcnew array<String^>(1) { ExtensionEntryLabel };
	array<TypeCode>^ ExtensionEntryDataTypes = gcnew array<TypeCode>(1) { ExtensionEntryDataType };
	array<String^>^ ExtensionEntryDataUnits = gcnew array<String^>(1) { ExtensionEntryDataUnit };

	WriteExtension(FileName, ExtensionName, OverWriteExtensionIfExists, ExtensionEntryLabels, /*ExtensionEntryDataTypes, ExtensionEntryDataTypeInstances,*/ ExtensionEntryDataUnits, ExtensionHeaderExtraKeys, ExtensionHeaderExtraKeyValues, ExtensionHeaderExtraKeyComments, dataobj);
}


array<String^>^ JPFITS::FITSBinTable::FORMATBINARYTABLEEXTENSIONHEADER(String^ ExtensionName, array<Object^>^ ExtensionEntryData, array<String^>^ ExtensionEntryLabels, array<TypeCode>^ ExtensionEntryDataTypes, array<int>^ ExtensionEntryDataTypeInstances, array<String^>^ ExtensionEntryDataUnits, array<String^>^ ExtensionHeaderExtraKeys, array<String^>^ ExtensionHeaderExtraKeyValues, array<String^>^ ExtensionHeaderExtraKeyComments)
{
	ArrayList^ hkeyslist = gcnew ArrayList();
	ArrayList^ hvalslist = gcnew ArrayList();
	ArrayList^ hcomslist = gcnew ArrayList();

	//need to get the table width in bytes and height number of rows
	int rownbytes = 0, nrows = 0;
	for (int i = 0; i < ExtensionEntryDataTypes->Length; i++)
		rownbytes += TYPECODETONBYTES(ExtensionEntryDataTypes[i]) * ExtensionEntryDataTypeInstances[i];
	array<int>^ dataNrows = EXTENSIONENTRYDATANROWS(ExtensionEntryData, ExtensionEntryDataTypes, ExtensionEntryDataTypeInstances);
	nrows = JPMath::Max(dataNrows, false);

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
	hvalslist->Add(rownbytes.ToString());
	hcomslist->Add("width of table in bytes");
	hkeyslist->Add("NAXIS2");
	hvalslist->Add(nrows.ToString());
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
		if (ExtensionEntryDataTypes[i] == TypeCode::UInt16)
		{
			hkeyslist->Add("TZERO" + (i + 1).ToString());
			hvalslist->Add("32768");
			hcomslist->Add("offset for unsigned 16-bit integers");

			hkeyslist->Add("TSCAL" + (i + 1).ToString());
			hvalslist->Add("1");
			hcomslist->Add("data are not scaled");
		}
		if (ExtensionEntryDataTypes[i] == TypeCode::Int16)
		{
			hkeyslist->Add("TZERO" + (i + 1).ToString());
			hvalslist->Add("0");
			hcomslist->Add("offset for signed 16-bit integers");

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
		if (ExtensionEntryDataTypes[i] == TypeCode::Int16)
		{
			hkeyslist->Add("TZERO" + (i + 1).ToString());
			hvalslist->Add("0");
			hcomslist->Add("offset for signed 32-bit integers");

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
			if (value->Trim() == "T")
				value = "                   T";
			if (value->Trim() == "F")
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

array<int>^ JPFITS::FITSBinTable::EXTENSIONENTRYDATANROWS(array<Object^>^ ExtensionEntryData, array<TypeCode>^ ExtensionEntryDataTypes, array<int>^ ExtensionEntryDataTypeInstances)
{
	array<int>^ lengths = gcnew array<int>(ExtensionEntryData->Length);

	for (int i = 0; i < ExtensionEntryData->Length; i++)
	{
		switch (ExtensionEntryDataTypes[i])
		{
			case TypeCode::Double:
			{
				if (ExtensionEntryDataTypeInstances[i] == 1)
					lengths[i] = ((array<double>^)ExtensionEntryData[i])->Length;
				else
					lengths[i] = ((array<double, 2>^)ExtensionEntryData[i])->Length / ExtensionEntryDataTypeInstances[i];
				break;
			}

			case TypeCode::Single:
			{
				if (ExtensionEntryDataTypeInstances[i] == 1)
					lengths[i] = ((array<float>^)ExtensionEntryData[i])->Length;
				else
					lengths[i] = ((array<float, 2>^)ExtensionEntryData[i])->Length / ExtensionEntryDataTypeInstances[i];
				break;
			}

			case TypeCode::Int64:
			case TypeCode::UInt64:
			{
				if (ExtensionEntryDataTypeInstances[i] == 1)
					lengths[i] = ((array<__int64>^)ExtensionEntryData[i])->Length;
				else
					lengths[i] = ((array<__int64, 2>^)ExtensionEntryData[i])->Length / ExtensionEntryDataTypeInstances[i];
				break;
			}

			case TypeCode::Int32:
			case TypeCode::UInt32:
			{
				if (ExtensionEntryDataTypeInstances[i] == 1)
					lengths[i] = ((array<__int32>^)ExtensionEntryData[i])->Length;
				else
					lengths[i] = ((array<__int32, 2>^)ExtensionEntryData[i])->Length / ExtensionEntryDataTypeInstances[i];
				break;
			}

			case TypeCode::Int16:
			case TypeCode::UInt16:
			{
				if (ExtensionEntryDataTypeInstances[i] == 1)
					lengths[i] = ((array<__int16>^)ExtensionEntryData[i])->Length;
				else
					lengths[i] = ((array<__int16, 2>^)ExtensionEntryData[i])->Length / ExtensionEntryDataTypeInstances[i];
				break;
			}

			case TypeCode::Byte:
			case TypeCode::SByte:
			{
				if (ExtensionEntryDataTypeInstances[i] == 1)
					lengths[i] = ((array<char>^)ExtensionEntryData[i])->Length;
				else
					lengths[i] = ((array<char, 2>^)ExtensionEntryData[i])->Length / ExtensionEntryDataTypeInstances[i];
				break;
			}

			case TypeCode::Boolean:
			{
				if (ExtensionEntryDataTypeInstances[i] == 1)
					lengths[i] = ((array<bool>^)ExtensionEntryData[i])->Length;
				else
					lengths[i] = ((array<bool, 2>^)ExtensionEntryData[i])->Length / ExtensionEntryDataTypeInstances[i];
				break;
			}

			default:
				throw gcnew Exception("Unrecognized typecode: '" + ExtensionEntryDataTypes[i].ToString() + "'");
		}
	}

	return lengths;
}

