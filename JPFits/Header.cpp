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


void JPFITS::FITSImage::FORMATHEADER()//for writing to a file
{
	int NKeys = HEADERKEYS->Length;
	int NCards = (NKeys-1)/36;
	int NBlankKeys = (NCards+1)*36-NKeys;
	HEADER = gcnew array<String^>((NCards+1)*36);

	String^ key;
	String^ value;
	String^ comment;

	for (int i = 0; i < NKeys - 1; i++)
	{
		key = HEADERKEYS[i];
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
		if (JPMath::IsNumeric(HEADERKEYVALS[i]))//then we have a numeric key value
		{
			double val = ::Convert::ToDouble(HEADERKEYVALS[i]);
			
			if (VALIDKEYEDIT(HEADERKEYS[i]))
			{
				//value = val.ToString();
				if (Math::Abs(val) <= 1e-5 || Math::Abs(val) >= 1e7)
					value = val.ToString("0.00###########e+00");
				else if (JPMath::IsInteger(val))
					value = val.ToString("G");
				else
					value = val.ToString("G");
				if (val == 0)//change the exp to integer 0
					value = "0";
			}
			else
				value = HEADERKEYVALS[i];

			L = value->Length;
			if (L > 20)
				value = value->Substring(0, 20);
			if (L < 20)
				for (int i = 0; i < 20 - L; i++)
					value = " " + value;//pad left
		}
		else//else it must be a string or comment.
		{
			value = HEADERKEYVALS[i];
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
		comment = HEADERKEYCOMS[i];
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
			comment = HEADERKEYVALS[i] + HEADERKEYCOMS[i];
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

	HEADER[NKeys-1] = "END                                                                             ";

	for (int i = 0; i < NBlankKeys; i++)
		HEADER[NKeys + i] = "                                                                                ";
}

void JPFITS::FITSImage::MAKEDEFHEADER()
{
	HEADERKEYS = gcnew array<String^>(8);
	HEADERKEYVALS = gcnew array<String^>(8);
	HEADERKEYCOMS = gcnew array<String^>(8);
	BZERO = 0;
	BSCALE = 1;
	BITPIX = -64;
	NAXIS1 = DIMAGE->GetLength(0);
	NAXIS2 = DIMAGE->GetLength(1);

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
		KeyIndex = L-1;//add to end of header (before END)
	array<String^>^ headerkeys = gcnew array<String^>(L+1);
	array<String^>^ headerkeyvals = gcnew array<String^>(L+1);
	array<String^>^ headerkeycoms = gcnew array<String^>(L+1);
	int c = 0;
	for (int i = 0; i < L+1; i++)
	{
		if (i == KeyIndex)
		{
			headerkeys[i] =	   NewKey;
			headerkeyvals[i] = NewValue;
			headerkeycoms[i] = NewComment;
			continue;
		}
		headerkeys[i] =	   HEADERKEYS[c];
		headerkeyvals[i] = HEADERKEYVALS[c];
		headerkeycoms[i] = HEADERKEYCOMS[c];
		c++;
	}
	HEADERKEYS = gcnew array<String^>(L+1);
	HEADERKEYVALS = gcnew array<String^>(L+1);
	HEADERKEYCOMS = gcnew array<String^>(L+1);
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
	if (KeyIndex > HEADERKEYS->Length-1 || KeyIndex < 0)
		return;

	int c = -1;
	array<String^>^ keys = gcnew array<String^>(HEADERKEYS->Length-1);
	array<String^>^ vals = gcnew array<String^>(HEADERKEYS->Length-1);
	array<String^>^ coms = gcnew array<String^>(HEADERKEYS->Length-1);
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
	RemoveKey(GetKeyIndex(Key,Value));
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

		this->AddKey(source->HeaderKeys[i],source->HeaderKeyValues[i],source->HeaderKeyComments[i],-1);
	}
}

