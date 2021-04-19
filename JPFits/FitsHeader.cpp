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

JPFITS::FITSImageHeader::FITSImageHeader(bool mayContainExtensions, array<double, 2>^ image)
{
	MAKE_DEFAULT_HEADER(mayContainExtensions, image);
}

JPFITS::FITSImageHeader::FITSImageHeader(ArrayList^ headerlines, bool populate_nonessential)
{
	if (populate_nonessential)
	{
		HEADERKEYS = gcnew array<String^>(headerlines->Count);
		HEADERKEYVALS = gcnew array<String^>(headerlines->Count);
		HEADERKEYCOMS = gcnew array<String^>(headerlines->Count);
		HEADERLINEISCOMMENT = gcnew array<bool>(headerlines->Count);
	}
	else
	{
		int N = 10;
		if (headerlines->Count < 10)
			N = headerlines->Count;
		HEADERKEYS = gcnew array<String^>(N);
		HEADERKEYVALS = gcnew array<String^>(N);
		HEADERKEYCOMS = gcnew array<String^>(N);
		HEADERLINEISCOMMENT = gcnew array<bool>(N);
	}

	//#pragma omp parallel for
	for (int i = 0; i < HEADERKEYS->Length; i++)
	{
		String^ line = (String^)headerlines[i];

		try
		{
			HEADERKEYS[i] = line->Substring(0, 8)->Trim();
			if (HEADERKEYS[i] == "END")
			{
				HEADERKEYVALS[i] = "";
				HEADERKEYCOMS[i] = "";
			}
			else if (line->Substring(8, 1) != "=")
			{
				HEADERLINEISCOMMENT[i] = true;
				HEADERKEYS[i] = "";
				HEADERKEYVALS[i] = "";
				HEADERKEYCOMS[i] = line->Trim();
			}
			else
			{
				if (JPMath::IsNumeric(line->Substring(10, 20)))//this has to work if it is supposed to be a numeric value here
					HEADERKEYVALS[i] = line->Substring(10, 20)->Trim();//get rid of leading and trailing white space
				else if (line->Substring(10, 20)->Trim() == "T" || line->Substring(10, 20)->Trim() == "F")
					HEADERKEYVALS[i] = line->Substring(10, 20)->Trim();
				else
				{
					if (line->IndexOf("'") == -1)
					{
						HEADERLINEISCOMMENT[i] = true;
						HEADERKEYS[i] = "";
						HEADERKEYVALS[i] = "";
						HEADERKEYCOMS[i] = line->Trim();
					}
					else
					{
						HEADERKEYVALS[i] = line->Substring(line->IndexOf("'") + 1, line->LastIndexOf("'") - line->IndexOf("'") - 1);
						HEADERKEYVALS[i] = HEADERKEYVALS[i]->Trim();
					}
				}

				int indx = line->IndexOf("/");
				if (indx != -1)
					HEADERKEYCOMS[i] = line->Substring(indx + 1)->Trim();
				else
					HEADERKEYCOMS[i] = "";
			}
		}
		catch (...)
		{
			throw gcnew Exception("Header line: \r\r'" + line + "'\r\r is very nastily formatted. Continuing should allow the image to be loaded but the header will not likely contain whatever was intended.");
		}
	}
}

String^ JPFITS::FITSImageHeader::GetKeyValue(String^ key)
{
	String^ result = "";
	bool brek = false;

	//#pragma omp parallel for
	for (int i = 0; i < HEADERKEYS->Length; i++)
		if (brek)
			break;
		else if (key == HEADERKEYS[i])
		{
			result = HEADERKEYVALS[i];
			brek = true;
		}

	return result;
}

String^ JPFITS::FITSImageHeader::GetKeyComment(String^ key)
{
	String^ result = "";
	bool brek = false;

	//#pragma omp parallel for
	for (int i = 0; i < HEADERKEYS->Length; i++)
		if (brek)
			break;
		else if (key == HEADERKEYS[i])
		{
			result = HEADERKEYCOMS[i];
			brek = true;
		}

	return result;
}

String^ JPFITS::FITSImageHeader::GetKeyName(int index)
{
	if (index < HEADERKEYS->Length)
		return HEADERKEYS[index];
	else
		return "";
}

String^ JPFITS::FITSImageHeader::GetKeyValue(int index)
{
	if (index < HEADERKEYVALS->Length)
		return HEADERKEYVALS[index];
	else
		return "";
}

String^ JPFITS::FITSImageHeader::GetKeyComment(int index)
{
	if (index < HEADERKEYCOMS->Length)
		return HEADERKEYCOMS[index];
	else
		return "";
}

int JPFITS::FITSImageHeader::GetKeyIndex(String^ key, bool KeyIsFullLineFormatted)
{
	int result = -1;
	bool brek = false;

	if (!KeyIsFullLineFormatted)
	{
		//#pragma omp parallel for
		for (int i = 0; i < HEADERKEYS->Length; i++)
			if (brek)
				break;
			else if (key == HEADERKEYS[i])
			{
				result = i;
				brek = true;
			}
	}
	else
	{
		//#pragma omp parallel for
		for (int i = 0; i < HEADERKEYS->Length; i++)
			if (brek)
				break;
			else if (key == this->HeaderLine[i])
			{
				result = i;
				brek = true;
			}
	}

	return result;
}

int JPFITS::FITSImageHeader::GetKeyIndex(String^ key, String^ keyvalue)
{
	int result = -1;
	bool brek = false;

	//#pragma omp parallel for
	for (int i = 0; i < HEADERKEYS->Length; i++)
		if (brek)
			break;
		else if (key == HEADERKEYS[i])
			if (keyvalue == HEADERKEYVALS[i])
				result = i;

	return result;
}

int JPFITS::FITSImageHeader::GetKeyIndex(String^ key, String^ keyvalue, String^ keycomment)
{
	int result = -1;
	bool brek = false;

	//#pragma omp parallel for
	for (int i = 0; i < HEADERKEYS->Length; i++)
		if (brek)
			break;
		else if (key == HEADERKEYS[i])
			if (keyvalue == HEADERKEYVALS[i])
				if (keycomment == HEADERKEYCOMS[i])
					result = i;

	return result;
}

void JPFITS::FITSImageHeader::SetKey(String^ Key, String^ Value, bool AddIfNotFound, int AddAtIndex)
{
	UPDATEDISPLAYHEADER = true;
	Key = Key->ToUpper();
	for (int i = 0; i < HEADERKEYS->Length; i++)
		if (Key == HEADERKEYS[i])
		{
			HEADERKEYVALS[i] = Value;
			return;
		}
	if (AddIfNotFound)
		AddKey(Key, Value, "", AddAtIndex);
}

void JPFITS::FITSImageHeader::SetKey(String^ Key, String^ Value, String^ Comment, bool AddIfNotFound, int AddAtIndex)
{
	UPDATEDISPLAYHEADER = true;
	Key = Key->ToUpper();
	for (int i = 0; i < HEADERKEYS->Length; i++)
		if (Key == HEADERKEYS[i])
		{
			HEADERKEYVALS[i] = Value;
			HEADERKEYCOMS[i] = Comment;
			return;
		}
	if (AddIfNotFound)
		AddKey(Key, Value, Comment, AddAtIndex);
}

void JPFITS::FITSImageHeader::SetKey(int index, String^ Key, String^ Value, String^ Comment)
{
	Key = Key->ToUpper();
	if (index < 0 || index >= HEADERKEYS->Length)
		return;

	HEADERKEYS[index] = Key;
	HEADERKEYVALS[index] = Value;
	HEADERKEYCOMS[index] = Comment;
	HEADERLINEISCOMMENT[index] = false;
	UPDATEDISPLAYHEADER = true;
}

void JPFITS::FITSImageHeader::AddKey(String^ NewKey, String^ NewValue, String^ NewComment, int KeyIndex)
{
	int c = 0;
	int L = HEADERKEYS->Length;
	if (KeyIndex < 0 || KeyIndex >= L)
		KeyIndex = L - 1;//add to end of header (before END)

	array<String^>^ headerkeys = gcnew array<String^>(L + 1);
	array<String^>^ headerkeyvals = gcnew array<String^>(L + 1);
	array<String^>^ headerkeycoms = gcnew array<String^>(L + 1);
	array<bool>^ headerlineiscomment = gcnew array<bool>(L + 1);
		
	for (int i = 0; i < L + 1; i++)
	{
		if (i == KeyIndex)
		{
			headerkeys[i] = NewKey;
			headerkeyvals[i] = NewValue;
			headerkeycoms[i] = NewComment;
			headerlineiscomment[i] = false;
			continue;
		}
		headerkeys[i] = HEADERKEYS[c];
		headerkeyvals[i] = HEADERKEYVALS[c];
		headerkeycoms[i] = HEADERKEYCOMS[c];
		headerlineiscomment[i] = HEADERLINEISCOMMENT[c];
		c++;
	}

	HEADERKEYS = headerkeys;
	HEADERKEYVALS = headerkeyvals;
	HEADERKEYCOMS = headerkeycoms;
	HEADERLINEISCOMMENT = headerlineiscomment;
	UPDATEDISPLAYHEADER = true;
}

void JPFITS::FITSImageHeader::AddCommentKeyLine(String^ commentKeyLine, int keyIndex)
{
	commentKeyLine = commentKeyLine->PadRight(80);//pad since might just be empty intended as blank line
	int nels = commentKeyLine->Length;
	int Nnewlines = (int)Math::Ceiling(double(nels) / 80);
	array<String^>^ strnewlines = gcnew array<String^>(Nnewlines);

	for (int i = 0; i < Nnewlines; i++)
		if ((i + 1) * 80 > nels)
			strnewlines[i] = commentKeyLine->Substring(i * 80);
		else
			strnewlines[i] = commentKeyLine->Substring(i * 80, 80);

	if (keyIndex < 0 || keyIndex >= HEADERKEYS->Length)
		keyIndex = HEADERKEYS->Length - 1;//add to end of header (before END)

	array<String^>^ headerkeys = gcnew array<String^>(HEADERKEYS->Length + Nnewlines);
	array<String^>^ headerkeyvals = gcnew array<String^>(HEADERKEYS->Length + Nnewlines);
	array<String^>^ headerkeycoms = gcnew array<String^>(HEADERKEYS->Length + Nnewlines);
	array<bool>^ headerlineiscomment = gcnew array<bool>(HEADERKEYS->Length + Nnewlines);

	int c = 0;
	for (int i = 0; i < headerkeys->Length; i++)
	{
		if (i >= keyIndex && i < keyIndex + Nnewlines)
		{
			if (commentKeyLine->Substring(0, 7) == "COMMENT")
			{
				headerkeys[i] = "COMMENT";
				headerkeycoms[i] = commentKeyLine->Substring(7);
			}
			else
			{
				headerkeys[i] = "";
				headerkeycoms[i] = commentKeyLine;
			}
			headerkeyvals[i] = "";
			headerlineiscomment[i] = true;
		}
		else
		{
			headerkeys[i] = HEADERKEYS[c];
			headerkeyvals[i] = HEADERKEYVALS[c];
			headerkeycoms[i] = HEADERKEYCOMS[c];
			headerlineiscomment[i] = HEADERLINEISCOMMENT[c];
			c++;
		}
	}

	HEADERKEYS = headerkeys;
	HEADERKEYVALS = headerkeyvals;
	HEADERKEYCOMS = headerkeycoms;
	HEADERLINEISCOMMENT = headerlineiscomment;
	UPDATEDISPLAYHEADER = true;
}

void JPFITS::FITSImageHeader::RemoveKey(int KeyIndex)
{
	if (KeyIndex < 0 || KeyIndex > HEADERKEYS->Length - 1)
		return;

	int c = -1;
	array<String^>^ keys = gcnew array<String^>(HEADERKEYS->Length - 1);
	array<String^>^ vals = gcnew array<String^>(HEADERKEYS->Length - 1);
	array<String^>^ coms = gcnew array<String^>(HEADERKEYS->Length - 1);
	array<bool>^ bols = gcnew array<bool>(HEADERKEYS->Length - 1);

	for (int i = 0; i < HEADERKEYS->Length; i++)
	{
		if (i == KeyIndex)
			continue;
		c++;
		keys[c] = HEADERKEYS[i];
		vals[c] = HEADERKEYVALS[i];
		coms[c] = HEADERKEYCOMS[i];
		bols[c] = HEADERLINEISCOMMENT[i];
	}
	HEADERKEYS = keys;
	HEADERKEYVALS = vals;
	HEADERKEYCOMS = coms;
	HEADERLINEISCOMMENT = bols;
	UPDATEDISPLAYHEADER = true;
}

void JPFITS::FITSImageHeader::RemoveKey(String^ Key)
{
	RemoveKey(GetKeyIndex(Key, false));
}

void JPFITS::FITSImageHeader::RemoveKey(String^ Key, String^ Value)
{
	RemoveKey(GetKeyIndex(Key, Value));
}

void JPFITS::FITSImageHeader::RemoveAllKeys(array<double, 2>^ image)
{
	MAKE_DEFAULT_HEADER(true, image);
	UPDATEDISPLAYHEADER = true;
}

void JPFITS::FITSImageHeader::CopyHeaderFrom(JPFITS::FITSImageHeader^ sourceHeader)
{
	array<bool>^ oktocopy = gcnew array<bool>(sourceHeader->HeaderKeys->Length);
	int ntocopy = 0;

	for (int i = 0; i < sourceHeader->HeaderKeys->Length; i++)
		if (VALIDKEYEDIT(sourceHeader->GetKeyName(i)))
		{
			ntocopy++;
			oktocopy[i] = true;
		}

	array<String^>^ newheaderkeys = gcnew array<String^>(HEADERKEYS->Length + ntocopy);
	array<String^>^ newheadervals = gcnew array<String^>(HEADERKEYS->Length + ntocopy);
	array<String^>^ newheadercoms = gcnew array<String^>(HEADERKEYS->Length + ntocopy);
	array<bool>^ newheaderbols = gcnew array<bool>(HEADERKEYS->Length + ntocopy);

	for (int i = 0; i < this->HeaderKeys->Length - 1; i++)//leave off END for now
	{
		newheaderkeys[i] = HEADERKEYS[i];
		newheadervals[i] = HEADERKEYVALS[i];
		newheadercoms[i] = HEADERKEYCOMS[i];
		newheaderbols[i] = HEADERLINEISCOMMENT[i];
	}

	int c = 0;
	for (int i = 0; i < sourceHeader->HeaderKeys->Length; i++)
		if (oktocopy[i])
		{
			newheaderkeys[HeaderKeys->Length - 1 + c] = sourceHeader->HeaderKeys[i];
			newheadervals[HeaderKeys->Length - 1 + c] = sourceHeader->HeaderKeyValues[i];
			newheadercoms[HeaderKeys->Length - 1 + c] = sourceHeader->HeaderKeyComments[i];
			newheaderbols[HeaderKeys->Length - 1 + c] = sourceHeader->HeaderLineIsComment[i];
			c++;
		}

	//END
	newheaderkeys[newheaderkeys->Length - 1] = "END";
	newheadervals[newheaderkeys->Length - 1] = "";
	newheadercoms[newheaderkeys->Length - 1] = "";
	newheaderbols[newheaderkeys->Length - 1] = false;

	HEADERKEYS = newheaderkeys;
	HEADERKEYVALS = newheadervals;
	HEADERKEYCOMS = newheadercoms;
	HEADERLINEISCOMMENT = newheaderbols;
	UPDATEDISPLAYHEADER = true;
}

void JPFITS::FITSImageHeader::SetBITPIXNAXISBSCZ(System::TypeCode precision, array<double, 2>^ image)
{
	UPDATEDISPLAYHEADER = true;
	if (image == nullptr || image->Length == 0)
	{
		SetKey("BITPIX", "8", false, 0);
		SetKey("NAXIS", "0", false, 0);
		RemoveKey("NAXIS1");
		RemoveKey("NAXIS2");
		RemoveKey("BZERO");
		RemoveKey("BSCALE");
		return;
	}
	else
	{
		SetKey("NAXIS", "2", "Number of image axes", true, 2);
		SetKey("NAXIS1", image->GetLength(0).ToString(), true, 3);
		SetKey("NAXIS2", image->GetLength(1).ToString(), true, 4);
	}

	switch (precision)
	{
		case TypeCode::SByte:
		{
			SetKey("BITPIX", "8", false, 0);
			SetKey("BZERO", "0", "Data Offset; pixel = pixel*BSCALE+BZERO", true, 5);
			SetKey("BSCALE", "1", "Data Scaling; pixel = pixel*BSCALE+BZERO", true, 6);
			break;
		}

		case TypeCode::Byte:
		{
			SetKey("BITPIX", "8", false, 0);
			SetKey("BZERO", "128", "Data Offset; pixel = pixel*BSCALE+BZERO", true, 5);
			SetKey("BSCALE", "1", "Data Scaling; pixel = pixel*BSCALE+BZERO", true, 6);
			break;
		}

		case TypeCode::Int16:
		{
			SetKey("BITPIX", "16", false, 0);
			SetKey("BZERO", "0", "Data Offset; pixel = pixel*BSCALE+BZERO", true, 5);
			SetKey("BSCALE", "1", "Data Scaling; pixel = pixel*BSCALE+BZERO", true, 6);
			break;
		}

		case TypeCode::UInt16:
		{
			SetKey("BITPIX", "16", false, 0);
			SetKey("BZERO", "32768", "Data Offset; pixel = pixel*BSCALE+BZERO", true, 5);
			SetKey("BSCALE", "1", "Data Scaling; pixel = pixel*BSCALE+BZERO", true, 6);
			break;
		}

		case TypeCode::Int32:
		{
			SetKey("BITPIX", "32", false, 0);
			SetKey("BZERO", "0", "Data Offset; pixel = pixel*BSCALE+BZERO", true, 5);
			SetKey("BSCALE", "1", "Data Scaling; pixel = pixel*BSCALE+BZERO", true, 6);
			break;
		}

		case TypeCode::UInt32:
		{
			SetKey("BITPIX", "32", false, 0);
			SetKey("BZERO", "2147483648", "Data Offset; pixel = pixel*BSCALE+BZERO", true, 5);
			SetKey("BSCALE", "1", "Data Scaling; pixel = pixel*BSCALE+BZERO", true, 6);
			break;
		}

		case TypeCode::Int64:
		{
			SetKey("BITPIX", "64", false, 0);
			SetKey("BZERO", "0", "Data Offset; pixel = pixel*BSCALE+BZERO", true, 5);
			SetKey("BSCALE", "1", "Data Scaling; pixel = pixel*BSCALE+BZERO", true, 6);
			break;
		}

		case TypeCode::UInt64:
		{
			SetKey("BITPIX", "64", false, 0);
			SetKey("BZERO", "9223372036854775808", "Data Offset; pixel = pixel*BSCALE+BZERO", true, 5);
			SetKey("BSCALE", "1", "Data Scaling; pixel = pixel*BSCALE+BZERO", true, 6);
			break;
		}

		case TypeCode::Single:
		{
			SetKey("BITPIX", "-32", false, 0);
			SetKey("BZERO", "0", "Data Offset; pixel = pixel*BSCALE+BZERO", true, 5);
			SetKey("BSCALE", "1", "Data Scaling; pixel = pixel*BSCALE+BZERO", true, 6);
			break;
		}

		case TypeCode::Double:
		{
			SetKey("BITPIX", "-64", false, 0);
			SetKey("BZERO", "0", "Data Offset; pixel = pixel*BSCALE+BZERO", true, 5);
			SetKey("BSCALE", "1", "Data Scaling; pixel = pixel*BSCALE+BZERO", true, 6);
			break;
		}

		default:
		{
			throw gcnew Exception("TypeCode '" + precision.ToString() + "' not recognized at SetBITPIXNAXISBSCZ.");
			return;
		}
	}
}

bool JPFITS::FITSImageHeader::VALIDKEYEDIT(String^ checkkey)
{
	for (int i = 0; i < INVALIDEDITKEYS->Length; i++)
		if (checkkey == INVALIDEDITKEYS[i])
			return false;

	return true;
}

void JPFITS::FITSImageHeader::MAKE_DEFAULT_HEADER(bool mayContainExtensions, array<double, 2>^ image)
{
	int nlines = 0;
	if (mayContainExtensions && image == nullptr)
		nlines = 5;
	else if (mayContainExtensions && image != nullptr)
		nlines = 9;
	else if (!mayContainExtensions && image == nullptr)
		nlines = 4;
	else  //(!mayContainExtensions && image != nullptr)
		nlines = 8;

	HEADERKEYS = gcnew array<String^>(nlines);
	HEADERKEYVALS = gcnew array<String^>(nlines);
	HEADERKEYCOMS = gcnew array<String^>(nlines);
	HEADERLINEISCOMMENT = gcnew array<bool>(nlines);

	__int64 BZERO = 0, BSCALE = 1;
	int BITPIX = -64, NAXIS = 0, NAXIS1 = 0, NAXIS2 = 0;
	if (image != nullptr)
	{
		NAXIS = 2;
		NAXIS1 = image->GetLength(0);
		NAXIS2 = image->GetLength(1);
	}

	HEADERKEYS[0] = "SIMPLE";
	HEADERKEYS[1] = "BITPIX";
	HEADERKEYS[2] = "NAXIS";
	HEADERKEYVALS[0] = "T";
	HEADERKEYVALS[1] = "-64";
	HEADERKEYVALS[2] = NAXIS.ToString();
	HEADERKEYCOMS[0] = "Valid Fits File";
	HEADERKEYCOMS[1] = "Bits per Pixel";
	HEADERKEYCOMS[2] = "Number of Axes";
	nlines = 3;
	if (image != nullptr)
	{
		HEADERKEYS[nlines] = "NAXIS1";
		HEADERKEYS[nlines + 1] = "NAXIS2";
		HEADERKEYS[nlines + 2] = "BZERO";
		HEADERKEYS[nlines + 3] = "BSCALE";
		HEADERKEYVALS[nlines] = NAXIS1.ToString();
		HEADERKEYVALS[nlines + 1] = NAXIS2.ToString();
		HEADERKEYVALS[nlines + 2] = "0";
		HEADERKEYVALS[nlines + 3] = "1";
		HEADERKEYCOMS[nlines] = "Width (No. of Columns)";
		HEADERKEYCOMS[nlines + 1] = "Height (No. of Rows)";
		HEADERKEYCOMS[nlines + 2] = "Data Offset";
		HEADERKEYCOMS[nlines + 3] = "Data Scaling";
		nlines += 4;
	}
	if (mayContainExtensions)
	{
		HEADERKEYS[nlines] = "EXTEND";
		HEADERKEYVALS[nlines] = "T";
		HEADERKEYCOMS[nlines] = "File may contain extensions";
		nlines += 1;
	}
	HEADERKEYS[nlines] = "END";
	HEADERKEYVALS[nlines] = "";
	HEADERKEYCOMS[nlines] = "";
	UPDATEDISPLAYHEADER = true;
}

String^ JPFITS::FITSImageHeader::GET_FORMATTED_HEADERLINE(int lineIndex)
{
	if (HEADERLINEISCOMMENT[lineIndex])
		if (HEADERKEYVALS[lineIndex] != "")
		{
			throw gcnew Exception("Error: Line was indicated as comment but the value field is not empty. key MAY equal COMMENT and comment string MUST be supplied in comment, and value must be empty if it is a comment line.");
			return "";
		}
		else
			if ((HEADERKEYS[lineIndex]->Length + HEADERKEYCOMS[lineIndex]->Length) > 80)
			{
				throw gcnew Exception("Error: key and comment are more than 80 elements: '" + HEADERKEYS[lineIndex] + HEADERKEYCOMS[lineIndex] + "' are " + (HEADERKEYS[lineIndex]->Length + HEADERKEYCOMS[lineIndex]->Length) + " elements.");
				return "";
			}
			else
				return (HEADERKEYS[lineIndex] + HEADERKEYCOMS[lineIndex])->PadRight(80);

	if (HEADERKEYS[lineIndex]->Length > 8)
	{
		throw gcnew Exception("Error: key was supplied with more than 8 characters: '" + HEADERKEYS[lineIndex] + "' is " + HEADERKEYS[lineIndex]->Length + " elements.");
		return "";
	}

	if (HEADERKEYS[lineIndex]->Trim() == "END")
		return HEADERKEYS[lineIndex]->Trim()->PadRight(80);

	String^ key = HEADERKEYS[lineIndex]->PadRight(8);
	key += "= ";
	//key formatting done

	//do value formatting
	String^ value;
	if (JPMath::IsNumeric(HEADERKEYVALS[lineIndex]))//then we have a numeric key value
	{
		double val = ::Convert::ToDouble(HEADERKEYVALS[lineIndex]);

		if (val == 9223372036854775808)//this is the bzero for unsigned int64...and is so large it will get converted to exponential notation which we do not want for this
			value = "9223372036854775808";
		else if (Math::Abs(val) <= 1e-6 || Math::Abs(val) >= 1e13)
			value = val.ToString("0.00###########e+00");
		else
			value = val.ToString("G");
		if (val == 0)//change the exp to integer 0
			value = "0";

		if (value->Length < 20)
			value = value->PadLeft(20);
	}
	else//else it must be a string
	{
		if (HEADERKEYVALS[lineIndex]->Trim() == "T" || HEADERKEYVALS[lineIndex]->Trim() == "F")
			value = HEADERKEYVALS[lineIndex]->Trim()->PadLeft(20);
		else
			value = "'" + HEADERKEYVALS[lineIndex]->PadRight(18) + "'";
	}
	//value formatting done

	String^ comment = HEADERKEYCOMS[lineIndex];

	comment = " / " + comment;//comment formatting done
	String^ line = key + value + comment;
	if (line->Length > 80)
		line = line->Substring(0, 80);
	else
		line = line->PadRight(80);

	return line;
}

array<String^>^ JPFITS::FITSImageHeader::GetFormattedHeaderBlock(bool isExtension, bool keysOnly)
{
	if (isExtension && !ISEXTENSION)
	{
		UPDATEDISPLAYHEADER = true;
		ISEXTENSION = true;

		this->RemoveKey("EXTEND");

		HEADERKEYS[0] = "XTENSION";
		HEADERKEYVALS[0] = "IMAGE";
		HEADERKEYCOMS[0] = "Image extension";

		bool pcountkey = false, gcountkey = false;
		for (int i = 0; i < HEADERKEYS->Length; i++)
			if (!pcountkey || !gcountkey)
			{
				if (HEADERKEYS[i]->Trim() == "PCOUNT")
					pcountkey = true;
				if (HEADERKEYS[i]->Trim() == "GCOUNT")
					pcountkey = true;
			}
		if (!pcountkey && !gcountkey)//they would BOTH not be present if things are being done correctly...need to add them
		{
			int naxis = -1;
			for (int i = 0; i < HEADERKEYS->Length; i++)
				if (HEADERKEYS[i]->Trim() == "NAXIS")
				{
					naxis = Convert::ToInt32(HEADERKEYVALS[i]->Trim());
					break;
				}
			int naxisNindex = -1;
			for (int i = 0; i < HEADERKEYS->Length; i++)
				if (HEADERKEYS[i]->Trim() == "NAXIS" + naxis.ToString())
				{
					naxisNindex = i;
					break;
				}
			array<String^>^ keys = gcnew array<String^>(HEADERKEYS->Length + 2);
			array<String^>^ vals = gcnew array<String^>(HEADERKEYS->Length + 2);
			array<String^>^ coms = gcnew array<String^>(HEADERKEYS->Length + 2);
			array<bool>^ bols = gcnew array<bool>(HEADERKEYS->Length + 2);

			for (int i = 0; i <= naxisNindex; i++)
			{
				keys[i] = HEADERKEYS[i];
				vals[i] = HEADERKEYVALS[i];
				coms[i] = HEADERKEYCOMS[i];
				bols[i] = HEADERLINEISCOMMENT[i];
			}
			keys[naxisNindex + 1] = "PCOUNT";
			vals[naxisNindex + 1] = "0";
			coms[naxisNindex + 1] = "number of bytes in heap area";
			bols[naxisNindex + 1] = false;
			keys[naxisNindex + 2] = "GCOUNT";
			vals[naxisNindex + 2] = "1";
			coms[naxisNindex + 2] = "single data table";
			bols[naxisNindex + 2] = false;

			for (int i = naxisNindex + 3; i < keys->Length; i++)
			{
				keys[i] = HEADERKEYS[i - 2];
				vals[i] = HEADERKEYVALS[i - 2];
				coms[i] = HEADERKEYCOMS[i - 2];
				bols[i] = HEADERLINEISCOMMENT[i - 2];
			}

			HEADERKEYS = keys;
			HEADERKEYVALS = vals;
			HEADERKEYCOMS = coms;
			HEADERLINEISCOMMENT = bols;
		}
	}
	else if (!isExtension && ISEXTENSION)
	{
		UPDATEDISPLAYHEADER = true;
		ISEXTENSION = false;
		HEADERKEYS[0] = "SIMPLE";
		HEADERKEYVALS[0] = "T";
		HEADERKEYCOMS[0] = "File conforms to FITS standard.";
	}

	if (!UPDATEDISPLAYHEADER && !keysOnly)
		return FORMATTEDHEADER;

	UPDATEDISPLAYHEADER = false;

	if (keysOnly)
		FORMATTEDHEADER = gcnew array<String^>(HEADERKEYS->Length);
	else
	{
		int NKeys = HEADERKEYS->Length;
		int NCards = (NKeys - 1) / 36;
		FORMATTEDHEADER = gcnew array<String^>((NCards + 1) * 36);
	}

	//#pragma omp parallel for
	for (int i = 0; i < HEADERKEYS->Length; i++)
		FORMATTEDHEADER[i] = FITSImageHeader::GET_FORMATTED_HEADERLINE(i);

	if (keysOnly)
		return FORMATTEDHEADER;

	String^ empty = "";
	//#pragma omp parallel for
	for (int i = HEADERKEYS->Length; i < FORMATTEDHEADER->Length; i++)
		FORMATTEDHEADER[i] = empty->PadLeft(80);

	return FORMATTEDHEADER;
}

