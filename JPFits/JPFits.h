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

#pragma once
using namespace System;
using namespace System::IO;
using namespace System::Collections::Specialized;
using namespace System::Collections;
using namespace System::Text;
using namespace System::ComponentModel;
using namespace System::Windows::Forms;
using namespace Microsoft::Win32;
using namespace System::Drawing;
using namespace System::Drawing::Drawing2D;
using namespace System::Drawing::Imaging;
#include "String.h"
#include "omp.h"
#include <fstream>
#include <algorithm>

/// <summary> JPFITS class to create, read, interact with, modify components of, and write FITS-format files.</summary>
namespace JPFITS
{
	static void SetReg(System::String^ ProgramName, System::String^ KeyName, System::Object^ KeyValue)
	{
		RegistryKey^ User = Registry::CurrentUser;
		RegistryKey^ SW = User->OpenSubKey("Software", true);
		RegistryKey^ AstroWerks = SW->CreateSubKey("AstroWerks");
		RegistryKey^ SUBKEY = AstroWerks->CreateSubKey(ProgramName);
		SUBKEY->SetValue(KeyName, KeyValue);
	}
	static Object^ GetReg(System::String^ ProgramName, System::String^ KeyName)
	{
		RegistryKey^ User = Registry::CurrentUser;
		RegistryKey^ SW = User->OpenSubKey("Software", true);
		RegistryKey^ AstroWerks = SW->CreateSubKey("AstroWerks");
		RegistryKey^ SUBKEY = AstroWerks->CreateSubKey(ProgramName);
		Object ^ result = SUBKEY->GetValue(KeyName);
		return result;
	}

	private ref class FITSFILEOPS
	{
		public:

		/// <summary>Scans the primary header unit and data of a FITS file. Returns false if the file is not a FITS file.</summary>
		/// <param name="fs">The FileStream of the FITS file.</param>
		/// <param name="scanpastprimarydata">True to set the FileStream fs position to the end of the data block, otherwise the fs position will be at the end of the primary header block, i.e. at the beginning of the primary data.</param>
		/// <param name="header_return">Returns the header of the extension as an ArrayList with each 80-character header line being a String^ member of this list. Pass nullptr if not required.</param>
		/// <param name="has_extensions">Returns whether or not the FITS file may contain extensions.</param>
		static bool SCANPRIMARYUNIT(FileStream^ fs, bool scanpastprimarydata, ArrayList^ header_return, bool &has_extensions);

		/// <summary>Find the FITS extension table of the given type and name. Returns false if the XTENSION type of the specified EXTNAME is not found. If extension_name is found the FileStream fs will be placed at the beginning of the extension's main data table block.</summary>
		/// <param name="fs">The FileStream of the FITS file.
		/// <para>If EXTNAME is found the FileStream fs will be placed at the beginning of the extension's main data table block.</para>
		/// <para>If EXTNAME is NOT found it will be at the end of the file.</para></param>
		/// <param name="extension_type">The XTENSION extension type, either: &quot;BINTABLE&quot;, &quot;TABLE&quot;, or &quot;IMAGE&quot;.</param>
		/// <param name="extension_name">The EXTNAME extension name. If the extension is known to have no EXTNAME keyword and name, then pass an empty String and the first nameless extension of the specified type will be seeked.</param>
		/// <param name="header_return">Returns the header of the extension as an ArrayList with each 80-character header line being a String^ member of this list. Pass nullptr if not required.</param>
		/// <param name="extensionStartPosition">Returns the start position within the FileStream of the extension...i.e. at the block boundary at the start of its header.</param>
		/// <param name="extensionEndPosition">Returns the end position within the FileStream of the extension, including after any heap, rounded up to a multiple of 2880 bytes at the last block boundary.</param>
		/// <param name="tableEndPosition">Returns the end position within the FileStream of the main data table, NOT rounded to a data block boundary.</param>
		/// <param name="pcount">Returns the number of bytes of any remaining fill plus supplemental heap data area after the main table endposition, IF any heap data exists. Does not represent fill bytes after the main table if no heap exists. Does not include fill bytes after the heap.</param>
		/// <param name="theap">Returns the position within the filestream of the beginning of the heap relative to the beginning of the main table. Nominally equal to NAXIS1 * NAXIS2 unless THEAP keyword specifies a larger value.</param>
		static bool SEEKEXTENSION(FileStream^ fs, String^ extension_type, String^ extension_name, ArrayList^ header_return, __int64 &extensionStartPosition, __int64 &extensionEndPosition, __int64 &tableEndPosition, __int64 &pcount, __int64 &theap);

		/// <summary>Find the FITS extension table of the given type and name. Returns false if the XTENSION type of the specified EXTNAME is not found. If extension_name is found the FileStream fs will be placed at the beginning of the extension's main data table block.</summary>
		/// <param name="fs">The FileStream of the FITS file.
		/// <para>If EXTNAME is found the FileStream fs will be placed at the beginning of the extension's main data table block.</para>
		/// <para>If EXTNAME is NOT found it will be at the end of the file.</para></param>
		/// <param name="extension_type">The XTENSION extension type, either: &quot;BINTABLE&quot;, &quot;TABLE&quot;, or &quot;IMAGE&quot;.</param>
		/// <param name="extension_number">The ONE-BASED extension number. This can be used when extensions aren't named with the EXTNAME keyword; alternatively if they are named this still returns the XTENSION extension type of the specified number.</param>
		/// <param name="header_return">Returns the header of the extension as an ArrayList with each 80-character header line being a String^ member of this list. Pass nullptr if not required.</param>
		/// <param name="extensionStartPosition">Returns the start position within the FileStream of the extension...i.e. at the block boundary at the start of its header.</param>
		/// <param name="extensionEndPosition">Returns the end position within the FileStream of the extension, including after any heap, rounded up to a multiple of 2880 bytes at the last block boundary.</param>
		/// <param name="tableEndPosition">Returns the end position within the FileStream of the main data table, NOT rounded to a data block boundary.</param>
		/// <param name="pcount">Returns the number of bytes of any remaining fill plus supplemental heap data area after the main table endposition, IF any heap data exists. Does not represent fill bytes after the main table if no heap exists. Does not include fill bytes after the heap.</param>
		/// <param name="theap">Returns the position within the filestream of the beginning of the heap relative to the beginning of the main table. Nominally equal to NAXIS1 * NAXIS2 unless THEAP keyword specifies a larger value.</param>
		static bool SEEKEXTENSION(FileStream^ fs, String^ extension_type, int extension_number, ArrayList^ header_return, __int64& extensionStartPosition, __int64& extensionEndPosition, __int64& tableEndPosition, __int64& pcount, __int64& theap);

		/// <summary>Gets all extension names of a specified extension type in the FITS file.</summary>
		/// <param name="FileName">The full file name to read from disk.</param>
		/// <param name="extension_type">The XTENSION extension type, either: &quot;BINTABLE&quot;, &quot;TABLE&quot;, or &quot;IMAGE&quot;.</param>
		static array<String^>^ GETALLEXTENSIONNAMES(String^ FileName, String^ extension_type);
	};

	/// <summary> FITSImageHeader class to create, read, interact with, modify components of, and write FITSImage headers. The FITSImage class contains an instance of this and is accessed there via FITSImage->Header.</summary>
	public ref class FITSImageHeader
	{
		public:

		#pragma region CONSTRUCTORS

		/// <summary>Constructor. Creates an instance of a FITSImageHeader, with options to indicate whether extensions are present, and sets essential keywords for a given image it will be the header for.
		/// <para>If the image is to be an extension, then use GetFormattedHeaderBlock to pull the header out with SIMPLE = T changed to XTENSION = IMAGE for writing.</para></summary>
		/// <param name="mayContainExtensions">If true, heyword EXTEND = T is added, otherwise it is left out.</param>
		/// <param name="image">If image is nullptr, then NAXIS = 0 and there are no NAXISn keywords or BSCALE or BZERO. Otherwise NAXIS, NAXISn, BSCALE and BZERO are set as per the image dimensions.
		/// <para>If the image will be saved at a different precision than double, use  SetBITPIXNAXISBSCZ(precision, image) at write time.</para></param>
		FITSImageHeader(bool mayContainExtensions, array<double, 2>^ image);

		/// <summary>Constructor. Creates an instance of a FITSImageHeader out of a list of header lines. Typically the headerlines would be returned from FITSFILEOPS::SCANPRIMARYUNIT.</summary>
		/// <param name="headerlines">A list of header lines to be extrated and formatted into keys, values, and comments, or as comment lines.</param>
		/// <param name="populate_nonessential">If false, non-essential key lines will be ignored. Saves a little bit of construction time if you don't need those, but they'll be lost if you re-write the file without them.</param>
		FITSImageHeader(ArrayList^ headerlines, bool populate_nonessential);
		#pragma endregion

		#pragma region HEADEROPS

		/// <summary>GetKeyName returns the key of the primary header line at index. Returns empty String if the index exceeds the number of header lines.</summary>
		/// <param name="index">The zero-based line number to get the key name from.</param>
		String^ GetKeyName(int index);

		/// <summary>GetKeyValue returns the value of the primary header line at index. Returns empty String if the index exceeds the number of header lines.</summary>
		/// <param name="index">The zero-based line number to get the key value from.</param>
		String^ GetKeyValue(int index);

		/// <summary>GetKeyComment returns the comment of the primary header line at index. Returns empty String if the index exceeds the number of header lines.</summary>
		/// <param name="index">The zero-based line number to get the key comment from.</param>
		String^ GetKeyComment(int index);

		/// <summary>GetKeyValue returns the value of the primary header key named Key. Returns empty String if the key is not found.</summary>
		/// <param name="Key">The header key to find the value of.</param>
		String^ GetKeyValue(String^ Key);

		/// <summary>GetKeyComment returns the comment of the primary header key named Key. Returns empty String if the key is not found.</summary>
		/// <param name="Key">The header key to find the comment of.</param>
		String^ GetKeyComment(String^ Key);		

		/// <summary>GetKeyIndex returns the zero-based index in the primary header of the key named Key. Returns -1 if the key is not found.</summary>
		/// <param name="Key">The header key to find the index of.</param>
		/// <param name="KeyIsFullLineFormatted">If true then the entire formatted 80-element long line is compared; helpful if multiple keys have the same name or are formatted as comment lines.</param>
		int	GetKeyIndex(String^ Key, bool KeyIsFullLineFormatted);

		/// <summary>GetKeyIndex returns the zero-based index in the primary header of the key with matching value. Returns -1 if the key and value combination is not found.</summary>
		/// <param name="Key">The header key to find the index of.</param>
		/// <param name="Value">The header key value to find the index of.</param>
		int	GetKeyIndex(String^ Key, String^ Value);

		/// <summary>GetKeyIndex returns the zero-based index in the primary header of the key with matching value and comment. Returns -1 if the key, value and comment combination is not found.</summary>
		/// <param name="Key">The header key to find the index of.</param>
		/// <param name="Value">The header key value to find the index of.</param>
		/// <param name="Comment">The header key comment to find the index of.</param>
		int	GetKeyIndex(String^ Key, String^ Value, String^ Comment);

		/// <summary>SetKey sets the value of the key. If the key already exists then the value will be replaced but the comment will remain the same.</summary>
		/// <param name="Key">The header key to access.</param>
		/// <param name="Value">The header key value to set.</param>
		/// <param name="AddIfNotFound">Optionally add the key to the header if it isn't found.</param>
		/// <param name="AddAtIndex">If the key wasn't found, add at this zero-based index. Use -1 to append to the end of the header (before END key).</param>
		void SetKey(String^ Key, String^ Value, bool AddIfNotFound, int AddAtIndex);

		/// <summary>SetKey sets the value and comment of the key. If the key already exists then the value and comment will be replaced.</summary>
		/// <param name="Key">The header key to access.</param>
		/// <param name="Value">The header key value to set.</param>
		/// <param name="Comment">The header key comment to set.</param>
		/// <param name="AddIfNotFound">Optionally add the key to the header if it isn't found.</param>
		/// <param name="AddAtIndex">If the key wasn't found, add at this zero-based index. Use -1 to append to the end of the header (before END key).</param>
		void SetKey(String^ Key, String^ Value, String^ Comment, bool AddIfNotFound, int AddAtIndex);

		/// <summary>SetKey sets the key, value and comment of the key at the given header index. This will overwrite whatever key exists at that index.</summary>
		/// <param name="index">The 0-based index of the header key to access. If the index does not occur within the header, then nothing happens.</param>
		/// <param name="Key">The header key to set.</param>
		/// <param name="Value">The header key value to set.</param>
		/// <param name="Comment">The header key comment to set.</param>
		void SetKey(int index, String^ Key, String^ Value, String^ Comment);

		/// <summary>AddKey adds a new key with value and comment to the primary header.</summary>
		/// <param name="NewKey">The header key to add.</param>
		/// <param name="NewValue">The header key value to add.</param>
		/// <param name="NewComment">The header key comment to add.</param>
		/// <param name="KeyIndex">Add at this zero-based index. Use -1 to append to the end of the header (before END key).</param>
		void AddKey(String^ NewKey, String^ NewValue, String^ NewComment, int KeyIndex);

		/// <summary>AddCommentKey adds a new key line formatted as a comment.
		/// <para>If the length of the commentKeyLine is more than 80 elements, the comment will be continued on subsequent lines until depleted.</para>
		/// <para>If the user wishes the line to begin with COMMENT, then write the input commentKeyLine beginning as such.</para>
		/// <para>If the user wishes the line to be blank, then pass commentKeyLine as an empty string or as only containing blanks (whitespace).</para></summary>
		/// <param name="commentKeyLine">The comment line.</param>
		/// <param name="keyIndex">Insert at this zero-based index. Use -1 to append to the end of the header (before END key). If keyIndex exceeds the header, the line is appended to the end of the header.</param>
		void AddCommentKeyLine(String^ commentKeyLine, int keyIndex);

		/// <summary>RemoveKey removes the key at the given index from the primary header.</summary>
		/// <param name="KeyIndex">The zero-based index of the key to remove. If the index is outside of the range of the header nothing happens.</param>
		void RemoveKey(int KeyIndex);

		/// <summary>RemoveKey removes the given key from the primary header. If there is more than one key with the given name, only the first occurence will be removed.</summary>
		/// <param name="Key">The name of the header key to remove.</param>
		void RemoveKey(String^ Key);

		/// <summary>RemoveKey removes the given key with matching value from the primary header. If there is more than one key with the given name and value, only the first occurence will be removed.</summary>
		/// <param name="Key">The name of the header key to remove.</param>
		/// <param name="Value">The corresponding header key value to remove.</param>
		void RemoveKey(String^ Key, String^ Value);

		/// <summary>RemoveAllKeys clears all keys from the primary header. Essential keywords will remain. image is supplied to re-create essential keywords, or pass nullptr to set NAXIS = 0.</summary>
		void RemoveAllKeys(array<double, 2>^ image);

		/// <summary>Copies a header from another FITSImageHeader into this one. Restricted keywords are neither copied nor overwritten.</summary>
		void CopyHeaderFrom(JPFITS::FITSImageHeader^ sourceHeader);

		/// <summary>This sets the BITPIX, NAXIS, NAXISn, BSCALE and BZERO keywords of the header given the TypeCode and the image. If the image is null then NAXIS = 0 and any NAXISn keywords are removed as well as BSCALE and BZERO.</summary>
		void SetBITPIXNAXISBSCZ(TypeCode precision, array<double, 2>^ image);

		/// <summary>Returns a formatted header block with the existing keys, and sets the first key to either SIMPLE = T or XTENSION = IMAGE. If a full 2880-multiple block is needed, set keysOnly to false.</summary>
		/// <param name="asExtension">If true then the first keyword is set to XTENSION = IMAGE, otherwise it is SIMPLE = T.</param>
		/// <param name="keysOnly">If true then only the existing keywords are returned formatted, otherwise if you need the entire 2880-multiple block pass false. True typically needed for display, false typically needed for writing.</param>
		array<String^>^ GetFormattedHeaderBlock(bool asExtension, bool keysOnly);

		/// <summary>ValidKeyEdit returns whether the given key is an essential key and shouldn't be user-modified.</summary>
		/// <param name="EditingKey">The name of the header key to remove.</param>
		static bool VALIDKEYEDIT(String^ EditingKey);

		#pragma endregion

		#pragma region HEADEPROPERTIES
		/// <summary>HeaderKeys accesses a string array of all keys of the primary header.
		/// <para>Individual String^ keys may be accessed by indexing HeaderKeys[i].</para></summary>
		property array<String^>^ HeaderKeys
		{
			array<String^>^ get() { return HEADERKEYS; }
		}

		/// <summary>HeaderKeyValues accesses a string array of all key values of the primary header.
		/// <para>Individual String^ values may be accessed by indexing HeaderKeyValues[i].</para></summary>
		property array<String^>^ HeaderKeyValues
		{
			array<String^>^ get() { return HEADERKEYVALS; }
		}

		/// <summary>HeaderKeyComments accesses a string array of all key comments of the primary header.
		/// <para>Individual String^ comments may be accessed by indexing HeaderKeyComments[i].</para></summary>
		property array<String^>^ HeaderKeyComments
		{
			array<String^>^ get() { return HEADERKEYCOMS; }
		}

		/// <summary>Returns a formatted header line from the existing key at the lineIndex</summary>
		property String^ HeaderLine[int]
		{
			String^ get(int i) { return GET_FORMATTED_HEADERLINE(i); }
		}

		/// <summary>Returns whether a header line at an index is formatted as a comment line.</summary>
		property bool HeaderLineIsComment[int]
		{
			bool get(int i) { return HEADERLINEISCOMMENT[i]; }
		}
		#pragma endregion

		#pragma region PRIVATEMEMBERS
		private:

		//Header things
		array<String^>^ HEADERKEYS;
		array<String^>^ HEADERKEYVALS;
		array<String^>^ HEADERKEYCOMS;
		array<String^>^ FORMATTEDHEADER;
		array<bool>^ HEADERLINEISCOMMENT;
		bool UPDATEDISPLAYHEADER = true;
		bool ISEXTENSION = false;

		static array<String^>^ INVALIDEDITKEYS = { "SIMPLE", "EXTEND", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "BZERO", "BSCALE", "END", "PCOUNT", "GCOUNT", "THEAP", "GROUPS", "XTENSION", "TFIELDS" };
		void MAKE_DEFAULT_HEADER(bool mayContainExtensions, array<double, 2>^ image);//make a default header
		String^ GET_FORMATTED_HEADERLINE(int lineIndex);//returns a formatted header line from the existing key at the lineIndex

		#pragma endregion
	};

	/// <summary>WorldCoordinateSolution class provides functionality for creation and interaction with World Coordinate Solutions for FITS files.</summary>
	public ref class WorldCoordinateSolution
	{
		public:
		~WorldCoordinateSolution();
		WorldCoordinateSolution();
		WorldCoordinateSolution(JPFITS::FITSImageHeader^ header);

		/// <summary>Gets or Sets the column-major CD matrix for this class instance.</summary>
		property array<double, 2>^ CD_Matrix
		{
			array<double, 2>^ get() { return CDMATRIX; }
			void set(array<double, 2>^ cdmatrix)
			{
				CDMATRIX = cdmatrix;
				CD1_1 = CDMATRIX[0, 0];
				CD1_2 = CDMATRIX[1, 0];
				CD2_1 = CDMATRIX[0, 1];
				CD2_2 = CDMATRIX[1, 1];
				SET_CDMATRIXINV();
			}
		}

		/// <summary>Gets the one-based row-major element from the CD matrix CDi_j[int i, int j], where i is the row index, and j is the column index.</summary>
		property double CDi_j[int, int]
		{
			double get(int i, int j)
			{
				return CDMATRIX[j - 1, i - 1];
			}
			void set(int i, int j, double val)
			{
				CDMATRIX[j - 1, i - 1] = val;
				CD1_1 = CDMATRIX[0, 0];
				CD1_2 = CDMATRIX[1, 0];
				CD2_1 = CDMATRIX[0, 1];
				CD2_2 = CDMATRIX[1, 1];
				SET_CDMATRIXINV();
			}
		}

		/// <summary>Gets the inverse of the CD matrix.</summary>
		property array<double, 2>^ CD_Matrix_Inverse
		{
			array<double, 2>^ get() { return CDMATRIXINV; }
		}

		/// <summary>Gets the array of coordinate values on one-based axis n (Coordinate_Values[n]) used for this World Coordinate Solution.</summary>
		property array<double>^ Coordinate_Values[int]
		{
			array<double>^ get(int Coordinate_Axis)
			{
				if (Coordinate_Axis == 1)
					return CVAL1;
				if (Coordinate_Axis == 2)
					return CVAL2;
				return nullptr;
			}
			void set(int Coordinate_Axis, array<double>^ cvals)
			{
				if (Coordinate_Axis == 1)
					CVAL1 = cvals;
				if (Coordinate_Axis == 2)
					CVAL2 = cvals;
			}
		}

		/// <summary>Gets the array of one-based coordinate pixels on one-based axis n (Coordinate_Pixels[n]) used for this World Coordinate Solution.</summary>
		property array<double>^ Coordinate_Pixels[int]
		{
			array<double>^ get(int Coordinate_Axis)
			{
				if (Coordinate_Axis == 1)
					return CPIX1;
				if (Coordinate_Axis == 2)
					return CPIX2;
				return nullptr;
			}
			void set(int Coordinate_Axis, array<double>^ cpixs)
			{
				if (Coordinate_Axis == 1)
					CPIX1 = cpixs;
				if (Coordinate_Axis == 2)
					CPIX2 = cpixs;
			}
		}

		/// <summary>Gets or Sets the Coordinate Reference Value for the one-based axis n: CRVALn[int n].</summary>
		property double CRVALn[int]
		{
			double get(int Coordinate_Axis)
			{
				return CRVALN[Coordinate_Axis - 1];
			}
			void set(int Coordinate_Axis, double val)
			{
				CRVALN[Coordinate_Axis - 1] = val;
			}
		}

		/// <summary>Gets or Sets the one-based Coordinate Reference Pixel for the one-based axis n: CRPIXn[int n].</summary>
		property double CRPIXn[int]
		{
			double get(int Coordinate_Axis)
			{
				return CRPIXN[Coordinate_Axis - 1];
			}
			void set(int Coordinate_Axis, double val)
			{
				CRPIXN[Coordinate_Axis - 1] = val;
			}
		}

		/// <summary>Gets the world coordinate solution plate scale (arcseconds per pixel) for one-based axis n: WCSSCALn[int n].</summary>
		property double CDELTn[int]
		{
			double get(int Coordinate_Axis)
			{
				return CDELTN[Coordinate_Axis - 1];
			}
		}

		/// <summary>Gets the world coordinate solution field rotation (degrees) for one-based axis n: WCSROTn[int n].</summary>
		property double CROTAn[int]
		{
			double get(int Coordinate_Axis)
			{
				return CROTAN[Coordinate_Axis - 1];
			}
		}

		/// <summary>Gets the world coordinate solution type for one-based axis n: CTYPEn[int n].</summary>
		property String^ CTYPEn[int]
		{
			String^ get(int Coordinate_Axis)
			{
				return CTYPEN[Coordinate_Axis - 1];
			}
		}

		property double WCSFitResidual_MeanPix
		{
			double get() { return CPIXRM; }
		}

		property double WCSFitResidual_StdvPix
		{
			double get() { return CPIXRS; }
		}

		property double WCSFitResidual_MeanSky
		{
			double get() { return CVALRM; }
		}

		property double WCSFitResidual_StdvSky
		{
			double get() { return CVALRS; }
		}

		/// <summary>Solves the projection parameters for a given list of pixel and coordinate values. Pass nullptr for FITS if writing WCS parameters to a primary header not required.</summary>
		/// <param name="WCS_Type">The world coordinate solution type. For example: TAN, for tangent-plane or Gnomic projection. Only TAN is currently supported.</param>
		/// <param name="X_pix">An array of the image x-axis pixel locations.</param>
		/// <param name="Y_pix">An array of the image y-axis pixel locations.</param>
		/// <param name="zero_based_pixels">A boolean to indicate if the X_Pix and Y_Pix are zero-based coordinates. They will be converted to one-based if true.</param>
		/// <param name="cval1">An array of coordinate values in degrees on coordinats axis 1.</param>
		/// <param name="cval2">An array of coordinate values in degrees on coordinats axis 2.</param>
		/// <param name="header">An FITSImageHeader^ instance to write the solution into. Pass nulltr if not required.</param>
		void Solve_WCS(String^ WCS_Type, array<double>^ X_pix, array<double>^ Y_pix, bool zero_based_pixels, array<double>^ cval1, array<double>^ cval2, JPFITS::FITSImageHeader^ header);

		/// <summary>Gets the image [x, y] pixel position for a given world coordinate in degrees at cval1 and cval2.</summary>
		void Get_Pixel(double cval1, double cval2, String^ WCS_Type, double &X_pix, double &Y_pix, bool return_zero_based_pixels);

		/// <summary>Gets arrays of image [x, y] pixel positions for a list of given world coordinates in degrees at cval1 and cval2.</summary>
		void Get_Pixels(array<double>^ cval1, array<double>^ cval2, String^ WCS_Type, array<double>^ &X_pix, array<double>^ &Y_pix, bool return_zero_based_pixels);

		/// <summary>Gets the cval1 and cval2 world coordinate in degrees for a given image [x, y] pixel position.</summary>
		void Get_Coordinate(double X_pix, double Y_pix, bool zero_based_pixels, String^ WCS_Type, double &cval1, double &cval2);

		/// <summary>Gets the cval1 and cval2 world coordinate in sexagesimal for a given image [x, y] pixel position.</summary>
		void Get_Coordinate(double X_pix, double Y_pix, bool zero_based_pixels, String^ WCS_Type, String^ &cval1_sxgsml, String^ &cval2_sxgsml);

		/// <summary>Gets the cval1 and cval2 world coordinate in degrees and sexagesimal for a given image [x, y] pixel position.</summary>
		void Get_Coordinate(double X_pix, double Y_pix, bool zero_based_pixels, String^ WCS_Type, double &cval1, double &cval2, String^ &cval1_sxgsml, String^ &cval2_sxgsml);

		/// <summary>Gets arrays of cval1 and cval2 world coordinates in degrees for a list of given image [x, y] pixel positions.</summary>
		void Get_Coordinates(array<double>^ X_pix, array<double>^ Y_pix, bool zero_based_pixels, String^ WCS_Type, array<double>^ &cval1, array<double>^ &cval2);

		/// <summary>Gets arrays of cval1 and cval2 world coordinates in sexagesimal for a list of given image [x, y] pixel positions.</summary>
		void Get_Coordinates(array<double>^ X_pix, array<double>^ Y_pix, bool zero_based_pixels, String^ WCS_Type, array<String^>^ &cval1_sxgsml, array<String^>^ &cval2_sxgsml);

		/// <summary>Gets arrays of cval1 and cval2 world coordinates in degrees and sexagesimal for a list of given image [x, y] pixel positions.</summary>
		void Get_Coordinates(array<double>^ X_pix, array<double>^ Y_pix, bool zero_based_pixels, String^ WCS_Type, array<double>^ &cval1, array<double>^ &cval2, array<String^>^ &cval1_sxgsml, array<String^>^ &cval2_sxgsml);

		void CopyFrom(JPFITS::WorldCoordinateSolution^ wcs_source);

		void CopyTo(JPFITS::FITSImageHeader^ header);

		void Clear();

		static void Clear(JPFITS::FITSImageHeader^ header);

		bool Exists() { return WCSEXISTS; }

		/// <summary>Checks to see if a complete WCS solution exists in the primary header of the given FITS object.</summary>
		static bool Exists(FITSImageHeader^ header, array<String^>^ WCS_CTYPEN);

		private:
		array<double, 2>^ CDMATRIX;
		array<double, 2>^ CDMATRIXINV;
		array<double>^ CVAL1;
		array<double>^ CVAL2;
		array<double>^ DVAL1;
		array<double>^ DVAL2;
		array<double>^ CPIX1;
		array<double>^ CPIX2;
		array<double>^ CRVALN;
		array<double>^ CRPIXN;
		array<double>^ CDELTN;
		array<double>^ CROTAN;
		array<String^>^ CTYPEN;
		double CD1_1, CD1_2, CD2_1, CD2_2, CPIX1RM, CPIX1RS, CVAL1RM, CVAL1RS, CPIX2RM, CPIX2RS, CVAL2RM, CVAL2RS, CPIXRM, CPIXRS, CVALRM, CVALRS, CCVALD1, CCVALD2;
		String^ CCVALS1;
		String^ CCVALS2;
		void SET_CDMATRIXINV();
		bool WCSEXISTS = false;

		void EATHEADERFORWCS(JPFITS::FITSImageHeader^ header);
	};

	/// <summary> FITSImage class to create, read, interact with, modify components of, and write FITS Primary image data and its Header.</summary>
	public ref class FITSImage
	{
		public:

		#pragma region CONSTRUCTORS
		/// <summary>Create a dummy FITSImage object with a simple primary header.</summary>
		/// <param name="FullFileName">File name. May be anything for this dummy object.</param>
		FITSImage(String^ FullFileName, bool mayContainExtensions);

		/// <summary>Create a FITSImage object from an array object containing existing data.
		/// <para>Image data is maintained at or converted to double precision.</para></summary>
		/// <param name="FullFileName">File name.</param>
		/// <param name="ImageData">The data array to use for the FITS image data. The precision and rank of the underlying array will be automatically determined. Vectors will be converted to an array with column-rank.</param>
		/// <param name="Do_Stats">Optionally perform the statistics to determine min, max, mean, median, and standard deviation of the image data - saves time if you don't need those.</param>
		/// <param name="do_parallel">Populate the FITSImage object ImageData and perform stats (if true) with parallelization.</param>
		FITSImage(String^ FullFileName, Object^ ImageData, bool Do_Stats, bool do_parallel);

		/// <summary>Create a FITSImage object with Primary image data loaded to RAM memory from disk.
		/// <para>Image data is loaded as double precision independent of storage precision on disk.</para></summary>
		/// <param name="FullFileName">File name.</param>
		/// <param name="Range">Range is ZERO based 1-D int array [xmin xmax ymin ymax].  Pass nullptr or Range[0] = -1 to default to full image size.</param>
		/// <param name="Populate_Header">Optionally populate the header - sometimes you just want the data, and can skip reading the non-essential header lines.</param>
		/// <param name="Populate_Data">Optionally populate the image data array - sometimes you just want the header and don't need the data.</param>
		/// <param name="Do_Stats">Optionally perform the statistics to determine min, max, mean, median, and standard deviation of the image data (if populated) - saves time if you don't need those.</param>
		/// <param name="do_parallel">Populate the FITSImage object ImageData and perform stats (if true) with parallelization.</param>
		FITSImage(String^ FullFileName, array<int, 1>^ Range, bool Populate_Header, bool Populate_Data, bool Do_Stats, bool do_parallel);

		/// <summary>Create a FITSImage object with extension image data loaded to RAM memory from disk.
		/// <para>Image data is loaded as double precision independent of storage precision on disk.</para></summary>
		/// <param name="FullFileName">File name.</param>
		/// <param name="extensionName">The EXTNAME extension name of the image. If an empty string is passed, the first nameless IMAGE extension will be read. Exception if no such extension exits.</param>
		/// <param name="Range">Range is ZERO based 1-D int array [xmin xmax ymin ymax].  Pass nullptr or Range[0] = -1 to default to full image size.</param>
		/// <param name="Populate_Header">Optionally populate the header - sometimes you just want the data, and can skip reading the non-essential header lines.</param>
		/// <param name="Populate_Data">Optionally populate the image data array - sometimes you just want the header and don't need the data.</param>
		/// <param name="Do_Stats">Optionally perform the statistics to determine min, max, mean, median, and standard deviation of the image data (if populated) - saves time if you don't need those.</param>
		/// <param name="do_parallel">Populate the FITSImage object ImageData and perform stats (if true) with parallelization.</param>
		FITSImage(String^ FullFileName, String^ extensionName, array<int, 1>^ Range, bool Populate_Header, bool Populate_Data, bool Do_Stats, bool do_parallel);

		/// <summary>Create a FITSImage object with extension image data loaded to RAM memory from disk.
		/// <para>Image data is loaded as double precision independent of storage precision on disk.</para></summary>
		/// <param name="FullFileName">File name.</param>
		/// <param name="extensionNumber">The ONE-BASED extension number of the image. Useful when extensions are not named with EXTNAME keyword. Will return the image at the extension number if they are named regardless.</param>
		/// <param name="Range">Range is ZERO based 1-D int array [xmin xmax ymin ymax].  Pass nullptr or Range[0] = -1 to default to full image size.</param>
		/// <param name="Populate_Header">Optionally populate the header - sometimes you just want the data, and can skip reading the non-essential header lines.</param>
		/// <param name="Populate_Data">Optionally populate the image data array - sometimes you just want the header and don't need the data.</param>
		/// <param name="Do_Stats">Optionally perform the statistics to determine min, max, mean, median, and standard deviation of the image data (if populated) - saves time if you don't need those.</param>
		/// <param name="do_parallel">Populate the FITSImage object ImageData and perform stats (if true) with parallelization.</param>
		FITSImage(String^ FullFileName, int extensionNumber, array<int, 1>^ Range, bool Populate_Header, bool Populate_Data, bool Do_Stats, bool do_parallel);

		/// <summary>Create a FITSImage object referencing raw UChar data on disk. Image data is loaded as double precision independent of storage precision on disk.</summary>
		/// <param name="FullFileName">File name for the FITS object.</param>
		/// <param name="DiskUCharBufferName">File name of the disk byte data.</param>
		/// <param name="Precision">Precision of the data stored in the disk char array.</param>
		/// <param name="NAxis1">Length of the 1st axis (x-axis)</param>
		/// <param name="NAxis2">Length of the 2nd axis (y-axis)</param>
		FITSImage(String^ FullFileName, String^ DiskUCharBufferName, TypeCode Precision, int NAxis1, int NAxis2);

		~FITSImage();
		#pragma endregion

		#pragma region WRITING
		/// <summary>Write a FITS image to disk as a primary header and primary image from the FITSImage object with its existing file name.
		/// <para>If the file name already exists on disk, the primary unit will be overwritten, and any existing extensions will be appended to conserve the data file.</para></summary>
		/// <param name="Precision">Byte precision at which to write the image data.</param>
		/// <param name="do_parallel">Populate the underlying byte arrays for writing with parallelization.</param>
		void WriteImage(TypeCode Precision, bool do_parallel);

		/// <summary>Write a FITS image to disk as a primary header and primary image from the FITSImage object with a given file name.
		/// <para>If the file name already exists on disk, the primary unit will be overwritten, and any existing extensions will be appended to conserve the data file.</para></summary>
		/// <param name="FullFileName">File name.</param>
		/// <param name="Precision">Byte precision at which to write the image data.</param>
		/// <param name="do_parallel">Populate the underlying byte arrays for writing with parallelization.</param>
		void WriteImage(String^ FullFileName, TypeCode Precision, bool do_parallel);

		/// <summary>Write a FITS image to disk as an extension from the FITSImage object with a given file name.</summary>
		/// <param name="FullFileName">File name. Pass the object&apos;s own FullFileName to write to its existing file name.
		/// <para>If the file doesn't yet exist on disk, then a new file will be created with an empty Primary Unit, and the image will be written as an extension.</para>
		/// <para>If the file does exist, then the extension will be written with the logic for the overwriteIfExists parameter, and</para>
		/// <para>the existing primary unit and any other extensions will be conserved to the file.</para></param>
		/// <param name="extensionName">The EXTNAME extension name of the IMAGE extension. 
		/// <para>If an empty string is passed, the first nameless IMAGE extension will be written to.</para>
		/// <para>If no such extension exists, the extension will be written as a new extension to the FITS file.</para></param>
		/// <param name="overwriteExtensionIfExists">If the image extension already exists it can be overwritten. If it exists and the option is given to not overwrite it, then an exception will be thrown.</param>
		/// <param name="Precision">Byte precision at which to write the image data.</param>
		/// <param name="do_parallel">Populate the underlying byte arrays for writing with parallelization.</param>
		void WriteImage(String^ FullFileName, String^ extensionName, bool overwriteExtensionIfExists, TypeCode Precision, bool do_parallel);

		void WriteFileFromDiskBuffer(bool DeleteOrigDiskBuffer);

		array<unsigned char>^ GetFormattedDataBlock(TypeCode precision, bool do_parallel);
		#pragma endregion

		#pragma region STATICFILEIO
		/// <summary>Convert a (possibly poorly formatted) whitespace-delimited text file to a double array.
		/// <para>If the text file is large (>2MB) the program may seem to hang...just let it run until control is returned.</para></summary>
		/// <param name="FullFileName">File name.</param>
		static array<double, 2>^ ConvertTxtToDblArray(String^ FullFileName, bool IsPoorlyFormatted);

		/// <summary>Return the primary image of the FITS file as a double 2-D array.</summary>
		/// <param name="file">The full file name to read from disk.</param>
		/// <param name="Range">Range is ZERO based 1-D int array [xmin xmax ymin ymax]. Pass nullptr or Range[0] = -1 to default to full image size.</param>
		static array<double, 2>^ ReadImageArrayOnly(String^ file, array<int, 1>^ Range, bool do_parallel);

		/// <summary>Return the primary image of the FITS file as a double 1-D array.</summary>
		/// <param name="file">The full file name to read from disk.</param>
		/// <param name="Range">Range is ZERO based 1-D int array [xmin xmax ymin ymax]. One of the axes ranges must be length equal to 1.
		/// <para> Pass nullptr or Range[0] = -1 to default to full image size, assuming the image data is a vector.</para></param>
		static array<double>^ ReadImageVectorOnly(String^ file, array<int, 1>^ Range, bool do_parallel);

		/// <summary>Reads an N-dimensional array and returns the results as a double array. User may reorginize the array based on the return variable axis lengths vector nAxisN.</summary>
		/// <param name="nAxisN">An initilized, but not instantiated, int vector to return the axis lengths for each axis.</param>
		static array<double>^ ReadPrimaryNDimensionalData(String^ fullFileName, array<int>^ &nAxisN);

		/// <summary>If a Primary data unit is saved as a layered image cube where each layer is unique, separate the layers into individual named extensions instead.</summary>
		/// <param name="sourceFullFileName">The file name of the FITS file with the layered primary data unit.</param>
		/// <param name="destFullFileName">The file name to write the extensions to. If it is the same name as the source, then the source will be completely overwritten, including any other existing extensions which that file may have had.</param>
		/// <param name="layerExtensionNames">The names for each layer extension. Must be equal in length to the number of layers to pull out of the primary data unit; all extenions must have a unique name.</param>
		static void ExtendizePrimaryImageLayerCube(String^ sourceFullFileName, String^ destFullFileName, array<String^>^ layerExtensionNames);

		/// <summary>Returns an array of all image table extension names in a FITS file. If there are no image table extensions, returns an empty array.</summary>
		/// <param name="FileName">The full file name to read from disk.</param>
		static array<String^>^ GetAllExtensionNames(String^ FileName)
		{
			return FITSFILEOPS::GETALLEXTENSIONNAMES(FileName, "IMAGE");
		}

		///// <summary>Convert a FITS image on disk to an image.</summary>
		//static void ConvertToImage(String^ Source_FullFileName, String^ Destination_FullFileName, String^ file_type, String^ contrast_scaling, bool invert_colormap, bool do_parallel);  
		#pragma endregion

		#pragma region IMAGEOPS
		/// <summary>StatsUpD updates the statistics for the primary image: maximum, minimum, mean, median, and standard deviation.</summary>
		void StatsUpD(bool do_parallel);

		/// <summary>Use SetImage to replace the existing double array for the FITSImage object with a new double array.</summary>
		/// <param name="ImageData">The 2-D double array to set for the FITSImage object.</param>
		/// <param name="Do_Stats">Optionally update the stats for the new array.</param>
		void SetImage(array<double, 2>^ ImageData, bool Do_Stats, bool do_parallel);

		/// <summary>Returns a double array of a subset of coordinates from the primary image.</summary>
		/// <param name="X_Center">The zero-based center-position of the primary x-axis of the subimage.</param>
		/// <param name="Y_Center">The zero-based center-position of the primary y-axis of the subimage.</param>
		/// <param name="X_HalfWidth">The +- half-width of the x-axis of the subimage.</param>
		/// <param name="Y_HalfWidth">The +- half-width of the y-axis of the subimage.</param>
		array<double, 2>^ GetSubImage(int X_Center, int Y_Center, int X_HalfWidth, int Y_HalfWidth);

		/// <summary>Returns a double array of a subset of coordinates from the primary image.</summary>
		/// <param name="X_Center">The zero-based center-position of the primary x-axis of the subimage.</param>
		/// <param name="Y_Center">The zero-based center-position of the primary y-axis of the subimage.</param>
		/// <param name="X_HalfWidth">The +- half-width of the x-axis of the subimage.</param>
		/// <param name="Y_HalfWidth">The +- half-width of the y-axis of the subimage.</param>
		/// <param name="xdata">The x-indices of the subimage.</param>
		/// <param name="ydata">The y-indices of the subimage.</param>
		array<double, 2>^ GetSubImage(int X_Center, int Y_Center, int X_HalfWidth, int Y_HalfWidth, array<int>^ &xdata, array<int>^& ydata);

		/// <summary>Returns a double array of a subset of coordinates from the primary image.</summary>
		/// <param name="Range">The zero-based start and end coordinates of the subimage in the primary image. Range is: [xmin xmax ymin ymax].</param>
		array<double, 2>^ GetSubImage(array<int, 1>^ Range);

		/// <summary>RotateCW rotates the primary image by 90 degrees.</summary>
		/// <param name="CW">True to rotate clockwise, false to rotate counter-clock-wise.</param>
		void RotateCW(bool CW);

		/// <summary>FlipVertical flips the image across the horizontal axis, i.e. up to down.</summary>
		void FlipVertical();

		/// <summary>FlipVertical flips the image across the vertical axis, i.e. left to right.</summary>
		void FlipHorizontal();
		#pragma endregion

		#pragma region OPERATORS
		static array<double, 2>^ operator +(FITSImage^ lhs_img, FITSImage^ rhs_img);
		static array<double, 2>^ operator +(FITSImage^ lhs_img, double scalar);
		static array<double, 2>^ operator +(double scalar, FITSImage^ rhs_img);
		static array<double, 2>^ operator -(FITSImage^ lhs_img, FITSImage^ rhs_img);
		static array<double, 2>^ operator -(FITSImage^ lhs_img, double scalar);
		static array<double, 2>^ operator -(double scalar, FITSImage^ rhs_img);
		static array<double, 2>^ operator /(FITSImage^ lhs_img, FITSImage^ rhs_img);
		static array<double, 2>^ operator /(FITSImage^ lhs_img, double scalar);
		static array<double, 2>^ operator /(double scalar, FITSImage^ rhs_img);
		static array<double, 2>^ operator *(FITSImage^ lhs_img, FITSImage^ rhs_img);
		static array<double, 2>^ operator *(FITSImage^ lhs_img, double scalar);
		static array<double, 2>^ operator *(double scalar, FITSImage^ rhs_img);
		static array<double, 2>^ operator ^(FITSImage^ lhs_img, FITSImage^ rhs_img);
		static array<double, 2>^ operator ^(FITSImage^ lhs_img, double scalar);
		static array<double, 2>^ operator ^(double scalar, FITSImage^ rhs_img);
		#pragma endregion

		#pragma region PROPERTIES

		/// <summary>Default indexer accesses the image element of the primary image of the FITSImage object.</summary>
		property double default[int]
		{
			double get(int x) { return DIMAGE[0, x]; }
			void set(int x, double val) { DIMAGE[0, x] = val; }
		}

		/// <summary>Default indexer accesses the image element of the primary image of the FITSImage object.</summary>
		property double default[int, int]
		{
			double get(int x, int y) { return DIMAGE[x, y]; }
			void set(int x, int y, double val) { DIMAGE[x, y] = val; }
		}

		/// <summary>Min returns the minimum of the FITS image data array.  Returns zero if there is no array loaded or if stats have not been performed.</summary>
		property double Min
		{
			double get() { return MIN; }
		}

		/// <summary>Max returns the maximum of the FITS image data array.  Returns zero if there is no array loaded or if stats have not been performed.</summary>
		property double Max
		{
			double get() { return MAX; }
		}

		/// <summary>Median returns the median of the FITS image data array.  Returns zero if there is no array loaded or if stats have not been performed.</summary>
		property double Median
		{
			double get() { return MEDIAN; }
		}

		/// <summary>Mean returns the average of the FITS image data array.  Returns zero if there is no array loaded or if stats have not been performed.</summary>
		property double Mean
		{
			double get() { return MEAN; }
		}

		/// <summary>Std returns the standard deviation of the FITS image data array.  Returns zero if there is no array loaded or if stats have not been performed.</summary>
		property double Std
		{
			double get() { return STD; }
		}

		/// <summary>Sum returns the sum of the FITS image data array.  Returns zero if there is no array loaded or if stats have not been performed.</summary>
		property double Sum
		{
			double get() { return SUM; }
		}

		/// <summary>Width returns the width of the FITS image data array.  Returns zero if there is no array loaded.</summary>
		property int Width
		{
			int get() { return NAXIS1; }
		}

		/// <summary>Height returns the height of the FITS image data array.  Returns zero if there is no array loaded.</summary>
		property int Height
		{
			int get() { return NAXIS2; }
		}

		/// <summary>Length returns the total number of elements of the FITS image data array.  Returns zero if there is no array loaded.</summary>
		property int Length
		{
			int get() { return NAXIS2 * NAXIS1; }
		}

		/// <summary>FileName accesses just the file name of the FITS object.</summary>
		property String^ FileName
		{
			String^ get() { return FILENAME; }
			void set(String^ newFileName)
			{
				FILENAME = newFileName;
				FULLFILENAME = FILEPATH + FILENAME;
			}
		}

		/// <summary>FilePath accesses just the file path of the FITS object.</summary>
		property String^ FilePath
		{
			String^ get() { return FILEPATH; }
			void set(String^ newFilePath)
			{
				FILEPATH = newFilePath + "\\";
				FULLFILENAME = FILEPATH + FILENAME;
			}
		}

		/// <summary>FullFileName accesses the full file path + name of the FITS object.</summary>
		property String^ FullFileName
		{
			String^ get() { return FULLFILENAME; }
			void set(String^ newFullFileName)
			{
				FULLFILENAME = newFullFileName;
				int index = FULLFILENAME->LastIndexOf("\\");
				FILENAME = FULLFILENAME->Substring(index + 1);
				FILEPATH = FULLFILENAME->Substring(0, index + 1);
			}
		}

		/// <summary>Image accesses the 2-D double array of the primary FITS object image.
		/// <para>Individual elements of the array can be accessed by indexing -&gt;Image[x,y].</para>
		/// <para>Property setter automatically performs image stats when Image is set.  Use -&gt;SetImage instead for option to not perform stats.</para></summary>
		property array<double, 2>^ Image
		{
			array<double, 2>^ get() { return DIMAGE; }
			void set(array<double, 2>^ img) { SetImage(img, true, true); }
		}

		/// <summary>Provides access to the image header.</summary>
		property JPFITS::FITSImageHeader^ Header
		{
			JPFITS::FITSImageHeader^ get() { return HEADER; }
		}

		/// <summary>Provides access to the image WCS.</summary>
		property JPFITS::WorldCoordinateSolution^ WCS
		{
			JPFITS::WorldCoordinateSolution^ get() { return WORLDCOORDINATESOLUTION; }
			void set(WorldCoordinateSolution^ WCS) { WORLDCOORDINATESOLUTION = WCS; }
		}

		#pragma endregion

		#pragma region PRIVATEMEMBERS
		private:
		//Image Conditions
		bool HEADER_POP;
		bool DATA_POP;
		bool STATS_POP;
		bool FROMDISKBUFFER;
		bool ISEXTENSION;
		String^ EXTNAME;
		bool EXTNAMEOVERWRITE;

		//Image
		array<double, 2>^ DIMAGE;//double precision image

		//Image Stats
		double MIN, MAX, MEAN, MEDIAN, STD, SUM;

		//Fits Info
		int NAXIS1 = -1, NAXIS2 = -1, BITPIX = -1, NAXIS = -1;
		double BZERO = -1, BSCALE = -1;
		array<int>^ NAXISN;

		//File Info
		String^ FILENAME;
		String^ FILEPATH;
		String^ FULLFILENAME;
		String^ DISKBUFFERFULLNAME;

		//File IO
		void READIMAGE(FileStream^ fs, array<int, 1>^ Range, bool do_parallel);
		void EATIMAGEHEADER();
		void WRITEIMAGE(TypeCode prec, bool do_parallel);

		//Header
		JPFITS::FITSImageHeader^ HEADER;

		//WCS
		JPFITS::WorldCoordinateSolution^ WORLDCOORDINATESOLUTION;

		#pragma endregion
	};

	/// <summary>FITSImageSet class is an ArrayList object to hold, manage, and perform operations on a set of FITSImage objects.</summary>
	public ref class FITSImageSet
	{
		public:

		/// <summary>Constructor for the FITSImageSet class.</summary>
		FITSImageSet();

		/// <summary>Loads FITS objects into the FITSImageSet. If the FITSImageSet already has members (not previously cleared), then the new memers are added (appended) to this FITSImageSet.</summary>
		/// <param name="files">The full path list of files to load into the set.</param>
		/// <param name="range">Range is ZERO based 1-D int array [xmin xmax ymin ymax].  Pass nullptr or Range[0] = -1 to default to full image size.</param>
		/// <param name="do_stats">Determine stats for each FITS object when loaded.</param>
		/// <param name="do_parallel">Load the FITS files in parallel.</param>
		/// <param name="show_waitbar">Optionally show a cancellable waitbar when loading. If cancelled, return value is false.</param>
		/// <param name="waitbar_message">Message to display on Waitbar progress if it is shown.</param>
		bool Load(array<String^>^ files, array<int>^ range, bool do_stats, bool do_parallel, bool show_waitbar, String^ waitbar_message);

		/// <summary>Write the FITSImage objects from the FITSImageSet to disk.</summary>
		/// <param name="precision">The precision at which to write the image data.</param>
		/// <param name="do_parallel">Write the images with parallelism. In the past with platter drives this would have been impossible, but fast solid state drives can handle it. If there's only a few images then don't bother, but useful when writing hundreds.</param>
		/// <param name="show_waitbar">Optionally show a cancellable waitbar when saving. If cancelled, return value is false.</param>
		/// <param name="waitbar_message">Message to display on Waitbar progress if it is shown.</param>
		bool Write(TypeCode precision, bool do_parallel, bool show_waitbar, String^ waitbar_message);

		/// <summary>Write the FITSImage objects from the FITSImageSet as extensions.</summary>
		/// <param name="fileName">The file name to write to.</param>
		/// <param name="firstAsPrimary">Option to write the first image in the set as the primary data block, otherwise all images to be written as extensions.</param>
		/// <param name="primaryHeader">IF the first image is NOT to be written as the primary data block, then a header may still be supplied for the primary block. Throws an exception if firstAsPrimary is true and primaryHeader is not nullptr.</param>
		/// <param name="extensionNames">The names of the extensions. If firstAsPrimary is true then the length of extensionNames should be ONE LESS than the number of images in the set. No elements may be empty strings; all elements must be unique. Pass nullptr for automatic naming as EXTENS_nnnnnn.</param>
		/// <param name="imagePrecisions">The precision at which to write the image data. If a single element array is passed then this precision is applied to all images.</param>
		bool WriteAsExtensions(String^ fileName, bool firstAsPrimary, JPFITS::FITSImageHeader^ primaryHeader, array<String^>^ extensionNames, array<TypeCode>^ imagePrecisions);

		/// <summary>Appends a FITSImage object to the ArrayList FITSImageSet object.</summary>
		void Add(FITSImage^ FITS)
		{
			FITSLIST->Add(FITS);
			CHECK_CODIMENSIONAL();
		}

		/// <summary>Inserts a FITSImage object to the ArrayList FITSImageSet object at a given index.
		/// <para>If index is larger than the FITSImageSet count, the FITS object will be appended to the end.</para></summary>
		void AddAt(int index, FITSImage^ FITS)
		{
			if (index >= FITSLIST->Count)
				FITSLIST->Add(FITS);
			else
				FITSLIST->Insert(index, FITS);
			CHECK_CODIMENSIONAL();
		}

		/// <summary>Removes the FITSImage object at index from the FITSImageSet.
		/// <para>If index is beyond the set size, nothing happens.</para></summary>
		void RemoveAt(int index)
		{
			if (index < FITSLIST->Count)
			{
				FITSLIST->RemoveAt(index);
				FITSLIST->TrimToSize();
				CHECK_CODIMENSIONAL();
			}
		}

		/// <summary>Removes the FITSImage objects starting at index from the FITSImageSet.
		/// <para>If index is beyond the set size, nothing happens.</para></summary>
		void RemoveFrom(int index)
		{
			if (index < FITSLIST->Count)
			{
				FITSLIST->RemoveRange(index, FITSLIST->Count - index);
				FITSLIST->TrimToSize();
				CHECK_CODIMENSIONAL();
			}
		}

		/// <summary>Removes the count range of FITSImage objects starting at index from the FITSImageSet.
		/// <para>If index is beyond the set size, nothing happens.</para>
		/// <para>If index plus count is beyond the set size, all elements from index are removed.</para></summary>
		void RemoveRange(int index, int count)
		{
			if (index < FITSLIST->Count && index + count <= FITSLIST->Count)
			{
				FITSLIST->RemoveRange(index, count);
				FITSLIST->TrimToSize();
				CHECK_CODIMENSIONAL();
				return;
			}

			if (index < FITSLIST->Count && index + count > FITSLIST->Count)
			{
				FITSLIST->RemoveRange(index, FITSLIST->Count - index);
				FITSLIST->TrimToSize();
				CHECK_CODIMENSIONAL();
			}
		}

		/// <summary>Gets the common directory of the FITSImage objects in the FITSImageSet based on their file paths.</summary>
		String^ GetCommonDirectory();

		/// <summary>Gets the common directory of a series of file names, based on their file paths.</summary>
		static String^ GetCommonDirectory(array<String^>^ filelist);

		/// <summary>Clears the ArrayList FITSImageSet object of all members.</summary>
		void Clear()
		{
			FITSLIST->Clear();
			System::GC::Collect();
		}

		/// <summary>Sort sorts the FITSImageSet list given the key. Returns -1 if there was an error with the sort.</summary>
		/// <param name="key">If key is &quot;filename&quot; then the FITSImageSet list is sorted according to the member file names.
		/// <para> For example if the file names are alphabetical or numeric then the FITSImageSet list will be sorted by increasing file name.</para>
		/// <para> Otherwise key is a primary header key and then their corresponding values will be used to sort the FITSImageSet list.</para></param>
		int Sort(String^ key);

		/// <summary>Create a FITSImage object with primary image that is the pixel-wise mean of the FITSImageSet primary images.</summary>
		/// <param name="ImageSet">The FITSImageSet object.</param>
		/// <param name="Do_Stats">Optionally perform the statistics to determine min, max, mean, median, and stdv of the FITSImage result - saves time if you don't need those.</param>
		/// <param name="Show_Waitbar">Optionally compute the function with a cancellable Waitbar. If cancelled, return value is nullptr.</param>
		static FITSImage^ Mean(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar);

		/// <summary>Create a FITSImage object with primary image that is the pixel-wise sigma-clipped mean of the FITSImageSet primary images.
		/// <para>The computation is iterative and may take a long time in some situations and so a cancellable WaitBar is mandatory.</para>
		/// <para>If the computation is cancelled the function will return with the most recent iteration of the sigma-clipped stack.</para></summary>
		/// <param name="ImageSet">The FITSImageSet object.</param>
		/// <param name="Do_Stats">Optionally perform the statistics to determine min, max, mean, median, and stdv of the FITSImage result - saves time if you don't need those.</param>
		/// <param name="sigma">The maximum standard deviation allowed for each pixel column; values beyond sigma are clipped and replaced with the median of the pixel column.</param>
		static FITSImage^ MeanClipped(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, double sigma);

		/// <summary>Create a FITSImage object with primary image that is the pixel-wise median of the FITSImageSet primary images.</summary>
		/// <param name="ImageSet">The FITSImageSet object.</param>
		/// <param name="Do_Stats">Optionally perform the statistics to determine min, max, mean, median, and stdv of the FITSImage result - saves time if you don't need those.</param>
		/// <param name="Show_Waitbar">Optionally compute the function with a cancellable Waitbar. If cancelled, return value is nullptr.</param>
		/// <param name="waitbar_message">Message to display on Waitbar progress if it is shown.</param>
		static FITSImage^ Median(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar, String^ waitbar_message);

		/// <summary>Create a FITSImage object with primary image that is the pixel-wise sum of the FITSImageSet primary images.</summary>
		/// <param name="ImageSet">The FITSImageSet object.</param>
		/// <param name="Do_Stats">Optionally perform the statistics to determine min, max, mean, median, and stdv of the FITSImage result - saves time if you don't need those.</param>
		/// <param name="Show_Waitbar">Optionally compute the function with a cancellable Waitbar. If cancelled, return value is nullptr.</param>
		static FITSImage^ Sum(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar);

		/// <summary>Create a FITSImage object with primary image that is the pixel-wise quadrature sum of the FITSImageSet primary images.</summary>
		/// <param name="ImageSet">The FITSImageSet object.</param>
		/// <param name="Do_Stats">Optionally perform the statistics to determine min, max, mean, median, and stdv of the FITSImage result - saves time if you don't need those.</param>
		/// <param name="Show_Waitbar">Optionally compute the function with a cancellable Waitbar. If cancelled, return value is nullptr.</param>
		static FITSImage^ Quadrature(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar);

		/// <summary>Create a FITSImage object with primary image that is the pixel-wise maximum of the FITSImageSet primary images.</summary>
		/// <param name="ImageSet">The FITSImageSet object.</param>
		/// <param name="Do_Stats">Optionally perform the statistics to determine min, max, mean, median, and stdv of the FITSImage result - saves time if you don't need those.</param>
		/// <param name="Show_Waitbar">Optionally compute the function with a cancellable Waitbar. If cancelled, return value is nullptr.</param>
		static FITSImage^ Max(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar);

		/// <summary>Create a FITSImage object with primary image that is the pixel-wise minimum of the FITSImageSet primary images.</summary>
		/// <param name="ImageSet">The FITSImageSet object.</param>
		/// <param name="Do_Stats">Optionally perform the statistics to determine min, max, mean, median, and stdv of the FITSImage result - saves time if you don't need those.</param>
		/// <param name="Show_Waitbar">Optionally compute the function with a cancellable Waitbar. If cancelled, return value is nullptr.</param>
		static FITSImage^ Min(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar);

		/// <summary>Create a FITSImage object with primary image that is the pixel-wise standard deviation of the FITSImageSet primary images.</summary>
		/// <param name="ImageSet">The FITSImageSet object.</param>
		/// <param name="Do_Stats">Optionally perform the statistics to determine min, max, mean, median, and stdv of the FITSImage result - saves time if you don't need those.</param>
		/// <param name="Show_Waitbar">Optionally compute the function with a cancellable Waitbar. If cancelled, return value is nullptr.</param>
		static FITSImage^ Stdv(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar);

		/// <summary>Auto-register non-rotational primary images from the FITSImageSet. Only works when there is no field rotation in the image set, only translational shifts, and the shifts are less than half of the field.</summary>
		/// <param name="ImageSet">The FITSImageSet object.</param>
		/// <param name="RefImgIndex">The index in the FitSet list of the reference image to register all the other images to.</param>
		/// <param name="Do_Stats">Optionally perform the statistics to determine min, max, mean, median, and stdv of the registered images - saves time if you don't need those.</param>
		static void Register(JPFITS::FITSImageSet^ ImageSet, int RefImgIndex, bool Do_Stats);

		/// <summary>Scans all primary FITS headers in the FITSImageSet for identical lines and copies such lines to the specified FITSImage destination primary header.
		/// <para>Usage is that perhaps you form the mean of the FITSImageSet as a new FITSImage, and this new FITSImage should contain all the primary header</para>
		/// <para> lines which are identical in the FITSImageSet.</para>
		/// <para>The existing primary header of the FITS_destination is cleared before the operation, except for essential keywords.</para></summary>
		static void GatherHeaders(JPFITS::FITSImageSet^ FITS_Set, FITSImage^ FITS_destination);

		/// <summary>Scans all primary FITS headers from the file names for identical lines and copies such lines to the specified FITSImage destination primary header.
		/// <para>Usage is that perhaps you form the mean of the FITSImageSet as a new FITSImage, and this new FITSImage should contain all the primary header</para>
		/// <para> lines which are identical in the file names.</para>
		/// <para>The existing primary header of the FITSImage is cleared before the operation, except for essential keywords.</para></summary>
		static void GatherHeaders(array<String^>^ filenames, JPFITS::FITSImage^ FITS_destination);

		/// <summary>FITSImageSet indexer accesses the FITSImage object in the FITSImageSet at a given index, i.e. FITSImage^ f = FITSImageSet[i].</summary>
		property FITSImage^ default[int]
		{
			FITSImage^ get(int i) { return ((FITSImage^)(FITSLIST[i])); }
			void set(int i, FITSImage^ FITS) { FITSLIST[i] = FITS; }
		}

		/// <summary>Returns the number of FITSImage objects currently held within the FITSImageSet.</summary>
		property int Count
		{
			int get() { return FITSLIST->Count; }
		}

		/// <summary>Returns whether all primary images in the current FITSImageSet have the same dimension.</summary>
		property bool CoDimensional
		{
			bool get() { CHECK_CODIMENSIONAL(); return CODIMENSIONAL; }
		}

		/// <summary>Returns a String array of the full file names (path + file name) of all FITSImage objects in the current FITSImageSet.</summary>
		property array<String^>^ FullFileNames
		{
			array<String^>^ get()
			{
				array<String^>^ names = gcnew array<String^>(FITSLIST->Count);
				for (int i = 0; i < FITSLIST->Count; i++)
					names[i] = ((FITSImage^)(FITSLIST[i]))->FullFileName;

				return names;
			}
		}

		/// <summary>Returns a String array of the file names (excluding file path) of all FITSImage objects in the current FITSImageSet.</summary>
		property array<String^>^ FileNames
		{
			array<String^>^ get()
			{
				array<String^>^ names = gcnew array<String^>(FITSLIST->Count);
				for (int i = 0; i < FITSLIST->Count; i++)
					names[i] = ((FITSImage^)(FITSLIST[i]))->FileName;

				return names;
			}
		}

		/// <summary>Returns a String array of the file paths (excluding file names) of all FITSImage objects in the current FITSImageSet.</summary>
		property array<String^>^ FilePaths
		{
			array<String^>^ get()
			{
				array<String^>^ names = gcnew array<String^>(FITSLIST->Count);
				for (int i = 0; i < FITSLIST->Count; i++)
					names[i] = ((FITSImage^)(FITSLIST[i]))->FilePath;

				return names;
			}
		}

		private:
		ArrayList^ FITSLIST;
		bool CODIMENSIONAL;
		void CHECK_CODIMENSIONAL();

		static JPWaitBar::WaitBar^ WAITBAR;
		static BackgroundWorker^ BGWRKR;
		static Object^ BGWRKR_RESULT;
		void BGWRKR_DoWork(System::Object^  sender, System::ComponentModel::DoWorkEventArgs^  e);
		void BGWRKR_ProgressChanged(System::Object^  sender, System::ComponentModel::ProgressChangedEventArgs^  e);
		void BGWRKR_RunWorkerCompleted(System::Object^  sender, System::ComponentModel::RunWorkerCompletedEventArgs^  e);

	};

	/// <summary> FITSBinTable class to create, read, interact with, modify components of, and write FITS BINTABLE binary table data extensions.</summary>
	public ref class FITSBinTable
	{
		public:

		/// <summary>Create an empty FITSBinTable object. TTYPE entries may be added later via SetTTYPEEntries or AddTTYPEEntry. An extension name can be added at writetime.</summary>
		FITSBinTable();

		/// <summary>Create a FITSBinTable object from an existing extension.</summary>
		/// <param name="fileName">The full path filename.</param>
		/// <param name="extensionName">The BINTABLE EXTNAME name of the extension. If an empty string is passed the first nameless extension will be found, if one exists.</param>
		FITSBinTable(String^ fileName, String^ extensionName);

		/// <summary>Check if a TTYPE entry exists within the bintable.</summary>
		/// <param name="ttypeEntry">The name of the binary table extension entry, i.e. the TTYPE value.</param>
		bool TTYPEEntryExists(String^ ttypeEntry);

		/// <summary>Return a binary table entry as a double 1-D array, assuming it is a single colunmn entry. If the entry has more than one column, use the overload function to get its dimensions.</summary>
		/// <param name="ttypeEntry">The name of the binary table extension entry, i.e. the TTYPE value.</param>
		array<double>^ GetTTYPEEntry(String^ ttypeEntry);

		/// <summary>Return a binary table entry as a double 1-D array.</summary>
		/// <param name="ttypeEntry">The name of the binary table extension entry, i.e. the TTYPE value.</param>
		/// <param name="dimNElements">A vector to return the number of elements along each dimension of the Object. 
		/// <para>Contains the TDIM key values for an n &gt; 2 dimensional array, otherwise contains the instances (repeats, i.e. columns) and NAXIS2. Its length gives the rank of the array Object. If rank = 1 then it contains only NAXIS2.</para></param>
		array<double>^ GetTTYPEEntry(String^ ttypeEntry, array<int>^ &dimNElements);

		/// <summary>Return a binary table entry as an Object. Its type and rank are given to the user. If you just need a double precision array to work on, use the overload for that.</summary>
		/// <param name="ttypeEntry">The name of the binary table extension entry, i.e. the TTYPE value.</param>
		/// <param name="objectTypeCode">The TypeCode precision of the underlying array in the object.</param>
		/// <param name="dimNElements">A vector to return the number of elements along each dimension of the Object. 
		/// <para>Contains the TDIM key values for an n &gt; 2 dimensional array, otherwise contains the instances (repeats, i.e. columns) and NAXIS2. Its length gives the rank of the array Object. If rank = 1 then it contains only NAXIS2.</para></param>
		Object^ GetTTYPEEntry(String^ ttypeEntry, TypeCode &objectTypeCode, array<int>^ &dimNElements);

		/// <summary>Use this to access individual elements of the table with a String return. Useful for looking at TTYPEs with multiple instances.</summary>
		/// <param name="ttypeEntry">The name of the binary table extension entry, i.e. the TTYPE value.</param>
		/// <param name="rowindex">The row index of the column.</param>
		String^ GetTTypeEntryRow(String^ ttypeEntry, int rowindex);

		/// <summary>Remove one of the entries from the binary table. Inefficient if the table has a very large number of entries with very large number of elements. Operates on heap-stored data if required.</summary>
		/// <param name="ttypeEntry">The name of the binary table extension entry, i.e. the TTYPE value.</param>
		void RemoveTTYPEEntry(String^ ttypeEntry);

		/// <summary>Add an entry to the binary table. Useful when dealing with a small table. Use SetTTYPEEntries for a large table to set all entries at once.</summary>
		/// <param name="ttypeEntry">The name of the binary table extension entry, i.e. the TTYPE value.</param>
		/// <param name="replaceIfExists">Replace the TTYPE entry if it already exists. If it already exists and the option is given to not replace, then an exception will be thrown.</param>
		/// <param name="entryUnits">The physical units of the values of the array. Pass empty string if not required.</param>
		/// <param name="entryArray">The vector or 2D array to enter into the table.</param>
		void AddTTYPEEntry(String^ ttypeEntry, bool replaceIfExists, String^ entryUnits, Object^ entryArray);

		/// <summary>Add an n &gt; 2 dimensional and/or complex entry to the binary table or heap area. If entries already exist then the user must have formatted the n &gt; 2 dimensional array to match the existing table height NAXIS2.
		/// <para>Otherwise it is recommended to create this table with ONLY the n &gt; 2 dimensional entry formatted simply as a vector, non-repeated instance. The height or NAXIS2 will then be the number of elements of the n &gt; 2 dimensional array.</para>
		/// <para>If dimensions need to be recorded then supply the dimNelements argument.</para>
		/// <para>If adding a complex number array to the binary table, the entryArray must be either single or double floating point.</para>
		/// <para>If complex the entryArray must be a factor of two columns repeats where the 1st and odd numbered columns are the spatial part, and the 2nd and even numbered columns are the temporal part.</para>
		/// <para>If it is a variable repeat heap array then the entry must be supplied as an array of arrays, or an array of Strings; if complex each subarray must contain an even pairing of values.</para></summary>
		/// <param name="ttypeEntry">The name of the binary table extension entry, i.e. the TTYPE value.</param>
		/// <param name="replaceIfExists">Replace the TTYPE entry if it already exists. If it already exists and the option is given to not replace, then an exception will be thrown.</param>
		/// <param name="entryUnits">The physical units of the values of the array. Pass empty string if not required.</param>
		/// <param name="entryArray">The array to enter into the table.</param>
		/// <param name="dimNElements">A vector giving the number of elements along each dimension of the array, to write as the TDIM key for the entry IF the entry is n &gt; 2 dimensional; pass nullptr if the entry is not n &gt; 2 dimensional.</param>
		/// <param name="isComplex">A boolean to set whether the array should be interpreted as complex value pairings.</param>
		/// <param name="addAsHeapVarRepeatArray">A boolean to set whether to save the array as a variable repeat array in the heap area. If true, the entryArray must be an array of arrays or an array of Strings.</param>
		void AddTTYPEEntry(String^ ttypeEntry, bool replaceIfExists, String^ entryUnits, Object^ entryArray, array<int>^ dimNElements, bool isComplex, bool addAsHeapVarRepeatArray);

		/// <summary>Set the bintable full of entries all at once. More efficient than adding a large number of entries once at a time. Useful to use with a brand new and empty FITSBinTable. NOTE: THIS CLEARS ANY EXISTING ENTRIES INCLUDING THE HEAP.
		/// <para>Do not use for n &gt; 2 dimensional and/or complex entries.</para></summary>
		/// <param name="ttypeEntries">The names of the binary table extension entries, i.e. the TTYPE values.</param>
		/// <param name="entryUnits">The physical units of the values of the arrays. Pass nullptr if not needed, or with null elements or empty elements where not required, etc.</param>
		/// <param name="entryArrays">An array of vectors or 2D arrays to enter into the table.</param>
		void SetTTYPEEntries(array<String^>^ ttypeEntries, array<String^>^ entryUnits, array<Object^>^ entryArrays);

		/// <summary>Add an extra key to the extension header. If it is to be a COMMENT, just fill the keyValue with eighteen characters, and the keyComment with 54 characters.</summary>
		/// <param name="keyName">The name of the key.</param>
		/// <param name="keyValue">The value of the key. Pass numeric types as a string.</param>
		/// <param name="keyComment">The comment of the key.</param>
		void AddExtraHeaderKey(String^ keyName, String^ keyValue, String^ keyComment);

		/// <summary>Get the value of an extra Key. If the key doesn't exist, an empty String is returned.</summary>
		/// <param name="keyName">The name of the key.</param>
		String^ GetExtraHeaderKeyValue(String^ keyName);

		/// <summary>Remove the extra header key with the given name and value.</summary>
		void RemoveExtraHeaderKey(String^ keyName, String^ keyValue);

		/// <summary>Clear all extra header keys.</summary>
		void RemoveAllExtraHeaderKeys();

		/// <summary>Returns an array of all binary table extension names in a FITS file. If there are no binary table extensions, returns an empty array.</summary>
		/// <param name="FileName">The full file name to read from disk.</param>
		static array<String^>^ GetAllExtensionNames(String^ FileName);

		/// <summary>Write the binary table into a new or existing FITS file. If the binary table already exists in an existing FITS file, it can optionally be replaced.</summary>
		/// <param name="FileName">The full file name to write the binary table into. The file can either be new or already exist.</param>
		/// <param name="ExtensionName">The EXTNAME name of the extension. Can be empty (unnamed) but this is poor practice.</param>
		/// <param name="OverWriteExtensionIfExists">If the binary table already exists it can be overwritten. If it exists and the option is given to not overwrite it, then an exception will be thrown.</param>
		void Write(String^ FileName, String^ ExtensionName, bool OverWriteExtensionIfExists);

		/// <summary>Remove a binary table extension from the given FITS file.</summary>
		/// <param name="FileName">The full-path file name.</param>
		/// <param name="ExtensionName">The name of the binary table extension. If the extension isn't found, an exception is thrown.</param>
		static void RemoveExtension(String^ FileName, String^ ExtensionName);

		/// <summary>Checks if the binary extension exists inside the given FITS file.</summary>
		/// <param name="FileName">The full-path file name.</param>
		/// <param name="ExtensionName">The name of the binary table extension.</param>
		static bool ExtensionExists(String^ FileName, String^ ExtensionName);

		#pragma region Class Properties
		public:
		/// <summary>NumberOfTableEntries reports the number of fields in the extension, i.e. the TFIELDS value.</summary>
		property int NumberOfTableEntriesTFIELDS
		{
			int get() { return TFIELDS; }
		}

		/// <summary>TableDataTypes reports the .NET typecodes for each entry in the table.</summary>
		property TypeCode TableDataTypes[int]
		{
			TypeCode get(int n) 
			{
				if (TTYPEISHEAPARRAYDESC[n])
					return HEAPTCODES[n];
				else
					return TCODES[n];
			}
		}

		/// <summary>Returns wheather the TTYPE entry at the given entry index is a variable repeat array.</summary>
		property bool TTYPEIsHeapVariableRepeatEntry[int]
		{
			bool get(int n) { return TTYPEISHEAPARRAYDESC[n]; }
		}

		/// <summary>Returns the number of elements (repeats) for a given heap entry at a given row.</summary>
		property int TTYPERowRepeatsHeapEntry[int, int]
		{
			int get(int ttypeIndex, int row) { return TTYPEHEAPARRAYNELSPOS[ttypeIndex][0, row]; }
		}

		/// <summary>TableDataTypes reports the number of columns or repeats in each table entry. Variable repeat entries only report 1...use </summary>
		property array<int>^ TableDataRepeats
		{
			array<int>^ get() { return TREPEATS; }
		}

		/// <summary>TableDataLabels reports the name of each table entry, i.e. the TTYPE values.</summary>
		property array<String^>^ TableDataLabelsTTYPE
		{
			array<String^>^ get() { return TTYPES; }
		}

		/// <summary>TableDataLabels reports the units of each table entry, i.e. the TUNITS values.</summary>
		property array<String^>^ GetExtensionEntryUnits
		{
			array<String^>^ get() { return TUNITS; }
		}

		/// <summary>Return the binary table header as an array of Strings for each line of the header.</summary>
		property array<String^>^ Header
		{
			array<String^>^ get() { return HEADER; }
		}

		property String^ ExtensionNameEXTNAME
		{
			String^ get() { return EXTENSIONNAME; }
		}

		/// <summary>Return the width, in bytes, of the table.</summary>
		property int Naxis1
		{
			int get() { return NAXIS1; }
		}

		/// <summary>Return the height, number of rows, of the table.</summary>
		property int Naxis2
		{
			int get() { return NAXIS2; }
		}

		/// <summary>Return the BINTABLE data block, excluding header, as a (unsigned) byte array.</summary>
		property array<unsigned char>^ BINTABLEByteArray
		{
			array<unsigned char>^ get() { return BINTABLE; }
		}
		#pragma endregion

		#pragma region PRIVATECLASSMEMBERS
		private:
		int BITPIX = 0, NAXIS = 0, NAXIS1 = 0, NAXIS2 = 0, TFIELDS = 0;
		array<String^>^ TTYPES;//names of each table entry
		array<String^>^ TFORMS;//FITS name for the table entry precisions
		array<bool>^ TTYPEISCOMPLEX;//for tracking complex singles and doubles
		array<bool>^ TTYPEISHEAPARRAYDESC;//for tracking array descriptor entries for heap area data
		array<array<int, 2>^>^ TTYPEHEAPARRAYNELSPOS;//for tracking array descriptor entries for heap area data
		array<TypeCode>^ HEAPTCODES;//.NET typcodes for each table entry
		array<array<int>^>^ TDIMS;//for tracking multidimensional (rank >= 3) arrays
		array<String^>^ TUNITS;//FITS name for the table entry units
		array<int>^ TBYTES;//number of total bytes for each table entry
		array<int>^ TREPEATS;//number of TFORM instances of the table entry
		array<TypeCode>^ TCODES;//.NET typcodes for each table entry
		array<String^>^ HEADER;
		String^ FILENAME;
		String^ EXTENSIONNAME;
		array<String^>^ EXTRAKEYS;
		array<String^>^ EXTRAKEYVALS;
		array<String^>^ EXTRAKEYCOMS;
		array<unsigned char>^ BINTABLE;
		array<unsigned char>^ HEAPDATA;

		void MAKEBINTABLEBYTEARRAY(array<Object^>^ ExtensionEntryData);
		void MAKEHEAPBYTEARRAY(array<Object^>^ ExtensionEntryData);
		void MAKETTYPEHEAPARRAYNELSPOS(array<Object^>^ ExtensionEntryData, __int64 &totalBytes);
		array<String^>^ FORMATBINARYTABLEEXTENSIONHEADER();
		int TFORMTONBYTES(String^ tform, int& instances);
		TypeCode TFORMTYPECODE(String^ tform);
		String^ TYPECODETFORM(TypeCode typecode);
		String^ TYPECODESTRING(TypeCode typecode);
		int TYPECODETONBYTES(TypeCode typecode);
		void EATRAWBINTABLEHEADER(ArrayList^ header);
		Object^ GETHEAPTTYPE(int ttypeindex, TypeCode &objectTypeCode, array<int>^ &dimNElements);
		array<int, 2>^ GETHEAPTTYPENELSPOS(int ttypeindex);
		void REMOVEHEAPTTYPE(int ttypeindex);
		#pragma endregion
	};

	/// <summary>JPMath class provides functionality for common mathematical operations.</summary>
	public ref class JPMath
	{
	public:

		ref class PointD
		{
		public:
			/// <summary>A double-precision point class.</summary>
			PointD(double x, double y, double value);

			/// <summary>The X-axis position.</summary>
			property double X
			{
				double get() { return POINTX; }
				void set(double x) { POINTX = x; }
			}

			/// <summary>The Y-axis position.</summary>
			property double Y
			{
				double get() { return POINTY; }
				void set(double y) { POINTY = y; }
			}

			/// <summary>The value of the point.</summary>
			property double Value
			{
				double get() { return POINTVAL; }
				void set(double val) { POINTVAL = val; }
			}

			/// <summary>Computes the distance from the current point to another point.</summary>
			double DistanceTo(JPMath::PointD^ other_point);

			static bool PointInPoly(JPMath::PointD^ P, array<JPMath::PointD^>^ V, int n);

			static void PolygonInteriorPointsRegion(array<bool, 2>^ region, array<JPMath::PointD^>^ polygon, int Xmin, int Ymin, int Xmax, int Ymax);

		private:
			double POINTX, POINTY, POINTVAL;
			static double ISLEFT(JPMath::PointD^ P0, JPMath::PointD^ P1, JPMath::PointD^ P2);

		};

		ref class Triangle
		{
		public:
			Triangle(JPMath::PointD^ point0, JPMath::PointD^ point1, JPMath::PointD^ point2);
			Triangle(array<JPMath::PointD^>^ points);

			property JPMath::PointD^ Vertex[int]
			{
				JPMath::PointD^ get(int i) { return POINTS[i]; }
			}

			property double VertexAngle[int]
			{
				double get(int i) { return VERTEXANGLES[i]; }
			}

			property array<JPMath::PointD^>^ Points
			{
				array<JPMath::PointD^>^ get() { return POINTS; }
				void set(array<JPMath::PointD^>^ points) { POINTS = points; }
			}

			property double SideLength[int]
			{
				double get(int i) { return SIDELENGTHS[i]; }
			}

			property JPMath::PointD^ FieldVector
			{
				JPMath::PointD^ get() { return FIELDVECTOR; }
			}

			property double FieldVectorRadAngle
			{
				double get() { return FIELDVECTORRADANGLE; }
			}

			property double VertexPointSum
			{
				double get() { return TOTALVERTEXPOINTVALUESUM; }
			}

		private:
			array<JPMath::PointD^>^ POINTS;
			array<double>^ VERTEXANGLES;
			array<double>^ SIDELENGTHS;
			void SORTTRIANGLE();
			void MAKEVERTEXANGLES();
			void MAKEFIELDVECTORS();
			JPMath::PointD^ FIELDVECTOR;
			double FIELDVECTORRADANGLE;
			double TOTALVERTEXPOINTVALUESUM;
		};

		/// <summary>Convert sexigesimal RA hh:mm:ss.s and DEC dd:am:as.as to degrees.</summary>
		/// <param name="ra_sex">Right ascension in sexigesimal format.</param>
		/// <param name="dec_sex">Declination in sexigesimal format.</param>
		/// <param name="separator">If known, specify the seperator. Ex. 00:00:00.0 - separator is :. If not known, pass empty string or nullptr and it will be automatically determined.</param>
		/// <param name="RA_deg">Right ascension in decimal format, return by reference.</param>
		/// <param name="DEC_deg">Declination in decimal format, return by reference.</param>
		static void RADecSexToDegree(String^ ra_sex, String^ dec_sex, String^ separator, double &RA_deg, double &DEC_deg);
		
		/// <summary>Convert RA and DEC degrees to sexigesimal RA hh:mm:ss.s and DEC dd:am:as.as.</summary>
		/// <param name="RA_deg">Right ascension in decimal format.</param>
		/// <param name="DEC_deg">Declination in decimal format.</param>
		/// <param name="separator">Specify the desired seperator. Ex. 00:00:00.0 - separator is :. Ex. 00-00-00.0 - separator is -.</param>
		/// <param name="ra_sex">Right ascension in sexigesimal format, return by reference.</param>
		/// <param name="dec_sex">Declination in sexigesimal format, return by reference.</param>
		static void RADecDegreeToSex(double RA_deg, double DEC_deg, String^ separator, String^ &ra_sex, String^ &dec_sex);

		/// <summary>Returns the angle between -PI to +PI radians following the CAST convention given the run (x) and rise (y) of the direction vector.</summary>
		/// <param name="run">The run (signed horizontal amplitude) of the vector.</param>
		/// <param name="rise">The rise (signed vertical amplitude) of the vector.</param>
		static double aTanAbsoluteAngle(double run, double rise);

		/// <summary>Returns the sum of all elements in the data array.</summary>
		/// <param name="vectorOrArray">A vecor or 2-D array.</param>
		static double Sum(Object^ vectorOrArray, bool do_parallel);

		/// <summary>Returns the sum of all elements in the data array.</summary>
		/// <param name="vectorOrArray">A vecor or 2-D array.</param>
		static double Mean(Object^ vectorOrArray, bool do_parallel);

		/*/// <summary>Returns the mean (average) of all elements in the data array.</summary>
		/// <param name="data">A 2-D double array.</param>
		static double Mean(array<double, 2>^ data, bool do_parallel);

		/// <summary>Returns the mean (average) of all elements in the data array.</summary>
		/// <param name="data">A 1-D double array.</param>
		static double Mean(array<double, 1>^ data, bool do_parallel);*/

		/*/// <summary>Returns the sum of all elements in the data array.</summary>
		/// <param name="data">A 2-D double array.</param>
		static double Sum(array<double, 2>^ data, bool do_parallel);

		/// <summary>Returns the sum of all elements in the data array.</summary>
		/// <param name="data">A 1-D double array.</param>
		static double Sum(array<double, 1>^ data, bool do_parallel);

		/// <summary>Returns the sum of all elements in the data array.</summary>
		/// <param name="data">A 2-D int array.</param>
		static int Sum(array<int, 2>^ data, bool do_parallel);

		/// <summary>Returns the sum of all elements in the data array.</summary>
		/// <param name="data">A 1-D int array.</param>
		static int Sum(array<int, 1>^ data, bool do_parallel);*/

		static array<double, 2>^ Pad(array<double, 2>^ data, array<int>^ padding, bool do_parallel);
		static array<double, 2>^ Crop(array<double, 2>^ data, array<int>^ cropping, bool do_parallel);
		static array<double, 2>^ Excise(array<double, 2>^ data, bool column, int X0, int halfWidth, bool do_parallel);

		/// <summary>Sum a 2-D array along one dimension, resulting in a 1-D vector array.</summary>
		/// <param name="data">A 2-D double array.</param>
		/// <param name="dim">The dimension along which to sum.  
		/// <para>0 (zero) sums along the horizontal axis, resulting in a vertical vector.</para>
		/// <para>1 (one) sums along the vertical axis, resulting in a horizontal vector.</para></param>
		static array<double, 1>^ Sum(array<double, 2>^ data, int dim, bool do_parallel);

		/// <summary>Sum a 2-D array along one dimension, resulting in a 1-D vector array.</summary>
		/// <param name="data">A 2-D int array.</param>
		/// <param name="dim">The dimension along which to sum.  
		/// <para>0 (zero) sums along the horizontal axis, resulting in a vertical vector.</para>
		/// <para>1 (one) sums along the vertical axis, resulting in a horizontal vector.</para></param>
		static array<int, 1>^ Sum(array<int, 2>^ data, int dim, bool do_parallel);		

		/// <summary>Average a 2-D array along one dimension, resulting in a 1-D vector array.</summary>
		/// <param name="data">A 2-D double array.</param>
		/// <param name="dim">The dimension along which to average.  
		/// <para>0 (zero) averages along the horizontal axis, resulting in a vertical vector.</para>
		/// <para>1 (one) averages along the vertical axis, resulting in a horizontal vector.</para></param>
		static array<double, 1>^ Mean(array<double, 2>^ data, int dim, bool do_parallel);

		/// <summary>Average a 2D array along one dimension, resulting in a 1D vector array.</summary>
		/// <param name="data">A 2D double array.</param>
		/// <param name="dim">The dimension along which to average.  
		/// <para>0 (zero) averages along the horizontal axis, resulting in a vertical vector.</para>
		/// <para>1 (one) averages along the vertical axis, resulting in a horizontal vector.</para></param>
		static array<double, 1>^ Stdv(array<double, 2>^ data, int dim, bool do_parallel);

		static double Mean_RobustClipped(array<double>^ data, double sigma);
		static double Mean_RobustClipped(array<double, 2>^ data, double sigma);

		/// <summary>Returns the median of all elements in the data array.</summary>
		/// <param name="data">A 2-D double array.</param>
		static double Median(array<double, 2>^ data);

		/// <summary>Returns the median of all elements in the data array.</summary>
		/// <param name="data">A 1-D double array.</param>
		static double Median(array<double>^ data);

		/// <summary>Returns the standard deviation of all elements in the data array.</summary>
		/// <param name="data">A 2-D double array.</param>
		static double Stdv(array<double, 2>^ data, bool do_parallel);

		/// <summary>Returns the standard deviation of all elements in the data array.</summary>
		/// <param name="data">A 2-D double array.</param>
		/// <param name="known_mean">If the mean of the data is already known, then save compute time by not having to calculate it first before the stdv is calculated.</param>
		static double Stdv(array<double, 2>^ data, double known_mean, bool do_parallel);

		/// <summary>Returns the standard deviation of all elements in the data array.</summary>
		/// <param name="data">A 1-D double array.</param>
		static double Stdv(array<double, 1>^ data, bool do_parallel);

		/// <summary>Returns the standard deviation of all elements in the data array.</summary>
		/// <param name="data">A 1-D double array.</param>
		/// <param name="known_mean">If the mean of the data is already known, then save compute time by not having to calculate it first before the stdv is calculated.</param>
		static double Stdv(array<double, 1>^ data, double known_mean, bool do_parallel);

		/// <summary>Returns the absolute values of all elements in the data array.</summary>
		/// <param name="data">A 2-D double array.</param>
		static array<double, 2>^ Abs(array<double, 2>^ data, bool do_parallel);

		/// <summary>Returns the absolute values of all elements in the data array.</summary>
		/// <param name="data">A 1-D double array.</param>
		static array<double, 1>^ Abs(array<double, 1>^ data, bool do_parallel);

		/// <summary>Returns the rounded values of all elements in the data array.</summary>
		/// <param name="data">A 2-D double array.</param>
		/// <param name="digits">The number of digits to which to round the data values.</param>
		static array<double, 2>^ Round(array<double, 2>^ data, int digits, bool do_parallel);

		/// <summary>Returns the floored (removed decimal parts, next lowest integer) values of all elements in the data array.</summary>
		/// <param name="data">A 2-D double array.</param>
		static array<double, 2>^ Floor(array<double, 2>^ data, bool do_parallel);

		/// <summary>Returns an array with all data values less than <i>clip_floor</i> replaced with <i>clip_floor</i>.</summary>
		/// <param name="data">A 2-D double array.</param>
		/// <param name="clip_floor">The value below which all data elements will be replaced with.</param>
		static array<double, 2>^ Floor(array<double, 2>^ data, double clip_floor, bool do_parallel);

		/// <summary>Returns the ceiling values (next higher integer) of all elements in the data array.</summary>
		/// <param name="data">A 2-D double array.</param>
		static array<double, 2>^ Ceil(array<double, 2>^ data, bool do_parallel);

		/// <summary>Applies an exponential power to a data array.</summary>
		/// <param name="data">A 2-D double array.</param>
		/// <param name="exponent">The exponent to apply to the array.</param>
		static array<double, 2>^ Power(array<double, 2>^ data, double exponent, bool do_parallel);

		/// <summary>Return the square root of an array.</summary>
		/// <param name="data">A 2-D double array.</param>
		static array<double, 2>^ Sqrt(array<double, 2>^ data, bool do_parallel);

		/// <summary>Return the base-10 logarithm of an array.</summary>
		/// <param name="data">A 2-D double array.</param>
		static array<double, 2>^ Log(array<double, 2>^ data, bool do_parallel);

		/// <summary>Return the custom-base logarithm of an array.</summary>
		/// <param name="data">A 2-D double array.</param>
		/// <param name="base">The base for the logarithm.</param>
		static array<double, 2>^ Log(array<double, 2>^ data, double base, bool do_parallel);

		/// <summary>Return the natural logarithm of an array.</summary>
		/// <param name="data">A 2-D double array.</param>
		static array<double, 2>^ Ln(array<double, 2>^ data, bool do_parallel);

		static array<double, 2>^ Exp(array<double, 2>^ data, bool do_parallel);
		static array<double, 2>^ Exp(array<double, 2>^ data, double base, bool do_parallel);

		/// <summary>Returns the maxima and their indices along a given dimension of the data array.
		/// <para> If dim = 0, the maxima are row-wise.</para>
		/// <para> If dim = 1, the maxima are column-wise.</para></summary>
		/// <param name="data">A 2-D double array.</param>
		/// <param name="dim">The dimension along which to reduce to maximums:  0 is x (rows), 1 is y (columns).</param>
		/// <param name="indices">An array passed to populate the indices at which the maxima appear along the dimension.</param>
		static array<double>^ Max(array<double, 2>^ data, int dim, array<int>^ &indices, bool do_parallel);

		/// <summary>Returns the global maximum and determines its [x, y] index in the 2-D array data.</summary>
		static double Max(array<double, 2>^ data, int& x, int& y, bool do_parallel);

		/// <summary>Returns the global maximum of the 2-D array data.</summary>
		static double Max(array<double, 2>^ data, bool do_parallel);

		/// <summary>Returns the global maximum and its index in the 1-D array data.</summary>
		static double Max(array<double>^ data, int& index, bool do_parallel);

		/// <summary>Returns the global maximum of the 1-D array data.</summary>
		static double Max(array<double>^ data, bool do_parallel);

		/// <summary>Returns the global maximum of the data array.</summary>
		static int Max(array<int>^ data, bool do_parallel);

		/// <summary>Returns the global maximum and its index within a subsection of the 1-D array data.</summary>
		/// <param name="data">A 1-D double array.</param>
		/// <param name="startIndex">The start index at which to begin checking for a maximum value.</param>
		/// <param name="endIndex">The end index within which to check for a maximum value.</param>
		/// <param name="maxIndex">The index in the array at which the maximum value occurs.</param>
		static double Max(array<double>^ data, int startIndex, int endIndex, int &maxIndex, bool do_parallel);

		/// <summary>Returns the global maximum of the data array.</summary>
		/// <param name="data">A 1-D int array.</param>
		/// <param name="index">The index at which the maximum occurs in the array.</param>
		static int Max(array<int>^ data, int &index, bool do_parallel);

		/// <summary>Returns the global maximum of the data array and its [x, y] indices.</summary>
		/// <param name="data">A 2-D int array.</param>
		static unsigned int Max(array<unsigned int, 2>^ data, int& x, int& y, bool do_parallel);

		/// <summary>Returns the minima and their indices along a given dimension of the data array.
		/// <para> If dim = 0, the minima are row-wise.</para>
		/// <para> If dim = 1, the minima are column-wise.</para></summary>
		/// <param name="data">A 2-D double array.</param>
		/// <param name="dim">The dimension along which to reduce to minimums:  0 is x (rows), 1 is y (columns).</param>
		/// <param name="indices">An array passed to populate the indices at which the minima appear along the dimension.</param>
		static array<double>^ Min(array<double, 2>^ data, int dim, array<int>^ &indices, bool do_parallel);

		/// <summary>Returns the global minimum and its [x, y] index in the 2-D array data.</summary>
		static double Min(array<double, 2>^ data, int& x, int& y, bool do_parallel);

		/// <summary>Returns the global minimum of the 2-D array data.</summary>
		static double Min(array<double, 2>^ data, bool do_parallel);

		/// <summary>Returns the global minimum and its index in the 1-D array data.</summary>
		static double Min(array<double>^ data, int& index, bool do_parallel);

		/// <summary>Returns the global minimum of the 1-D array data.</summary>
		static double Min(array<double>^ data, bool do_parallel);

		static void MinSetIndexes(array<double>^ data, array<int>^ indexes, bool do_parallel);

		/// <summary>Returns the global minimum and its indeces within a subsection of the 1-D array data.</summary>
		/// <param name="data">A 1-D double array.</param>
		/// <param name="startIndex">The start index at which to begin checking for a minimum value.</param>
		/// <param name="endIndex">The end index within which to check for a minimum value.</param>
		/// <param name="minIndex">The index in the array at which the minimum value occurs.</param>
		static double Min(array<double, 1>^ data, int startIndex, int endIndex, int &minIndex, bool do_parallel);

		/// <summary>Returns the global minimum and maximum of the 2-D array data.</summary>
		static void MinMax(array<double, 2>^ data, double &min, double &max, bool do_parallel);

		/// <summary>Returns the global minimum and maximum of the 1-D array data.</summary>
		static void MinMax(array<double>^ data, double &min, double &max, bool do_parallel);

		/// <summary>Returns a new array with all values at the given indeces replaced with the given value.</summary>
		/// <param name="data">A 2-D double array.</param>
		/// <param name="coords">An n x 2 array giving the row [n, 0] and column [n, 1] indices at which to replace the values.</param>
		/// <param name="val">The value with which to replace at the given indices.</param>
		static array<double, 2>^ Replace(array<double, 2>^ data, array<int, 2>^ coords, double val, bool do_parallel);

		/// <summary>Returns a new array with all values at the given indeces replaced with the given value.</summary>
		/// <param name="data">A 2-D double array.</param>
		/// <param name="coords">An array giving the indices at which to replace the values.</param>
		/// <param name="val">The value with which to replace at the given indices.</param>
		static array<double>^ Replace(array<double>^ data, array<int>^ coords, double val);

		/// <summary>Returns an array with the indeces at which the 2D data array satisfies the matching style for the given value.
		/// <para>The return array is an n x 2 array giving the row [n, 0] and column [n, 1] indices of the match.</para></summary>
		/// <param name="data">The data array to check for matches.</param>
		/// <param name="val">The value with which to check for a match in the data array.</param>
		/// <param name="style">The matching style can be &lt;, &lt;=, ==, &gt;=, &gt;, !=.</param>
		static array<int, 2>^ Find(array<double, 2>^ data, double val, String^ style, bool do_parallel);

		/// <summary>Returns an array with the indeces at which the 1D data array satisfies the matching style for the given value.</summary>
		/// <param name="data">The data array to check for matches.</param>
		/// <param name="val">The value with which to check for a match in the data array.</param>
		/// <param name="style">The matching style can be &lt;, &lt;=, ==, &gt;=, &gt;, !=.</param>
		static array<int, 1>^ Find(array<double, 1>^ data, double val, String^ style);

		/// <summary>Returns an array with the indeces at which the 1D data array satisfies the matching style for the given value.</summary>
		/// <param name="data">The data array to check for matches.</param>
		/// <param name="val">The value with which to check for a match in the data array.</param>
		/// <param name="style">The matching style can be &lt;, &lt;=, ==, &gt;=, &gt;, !=.</param>
		/// <param name="startindex">The starting index at which to begin checking for matches.</param>
		static array<int, 1>^ Find(array<double, 1>^ data, double val, String^ style, int startindex);

		/// <summary>Returns an array with the indeces at which the 1D data array satisfies the matching style for the given value.</summary>
		/// <param name="data">The data array to check for matches.</param>
		/// <param name="val">The value with which to check for a match in the data array.</param>
		/// <param name="style">The matching style can be &lt;, &lt;=, ==, &gt;=, &gt;, !=.</param>
		/// <param name="startindex">The starting index at which to begin checking for matches.</param>
		/// <param name="endindex">The ending index at which to stop checking for matches.</param>
		static array<int, 1>^ Find(array<double, 1>^ data, double val, String^ style, int startindex, int endindex);

		/// <summary>Returns either the first or last index in the data array that satisfies the match.</summary>
		/// <param name="data">The data array to check for matches.</param>
		/// <param name="val">The value with which to check for a match in the data array.</param>
		/// <param name="style">The matching style can be &lt;, &lt;=, ==, &gt;=, &gt;, !=.</param>
		/// <param name="return_first_true_last_false">Return first index of the match (true) or the last index (false).</param>
		static int Find(array<double, 1>^ data, double val, String^ style, bool return_first_true_last_false);

		/// <summary>Returns the cross correlation of two equal-length vectors and its lag shifts.</summary>
		/// <param name="reference">The reference data array against which to create the cross correlation.</param>
		/// <param name="relative">The comparison data array with which to create the cross correlation.</param>
		/// <param name="lags">An array passed (arrays pass by reference) to populate the cross correlation lags.</param>
		static array<double>^ XCorr(array<double>^ reference, array<double>^ relative, array<int>^ lags, bool do_parallel);

		/// <summary>Determines the cross correlation lags between two images using the image-reduction-to-vector method.</summary>
		/// <param name="reference">The reference data array against which to create the cross correlation.</param>
		/// <param name="COMPARISON">The comparison data array with which to create the cross correlation.</param>
		/// <param name="autoDeBias_refX">Option to automatically de-gradient the reference image along the x-dimension (horizontal degradient).</param>
		/// <param name="autoDeBias_refY">Option to automatically de-gradient the reference image along the y-dimension (vertical degradient).</param>
		/// <param name="autoDeBias_COMX">Option to automatically de-gradient the comparison image along the x-dimension (horizontal degradient).</param>
		/// <param name="autoDeBias_COMY">Option to automatically de-gradient the comparison image along the y-dimension (vertical degradient).</param>
		/// <param name="autoHanning_ref">Option to automatically Hanning-window the reference image.</param>
		/// <param name="autoHanning_COM">Option to automatically Hanning-window the comparison image.</param>
		/// <param name="xshift">The sub-integer x-shift of the comparison with respect to the reference, passed by reference.</param>
		/// <param name="yshift">The sub-integer y-shift of the comparison with respect to the reference, passed by reference.</param>
		/// <param name="do_parallel">Optionally perform all array operations in parallel. False when parallelizing upstream.</param>
		static void XCorrImageLagShifts(array<double, 2>^ reference, array<double, 2>^ COMPARISON, bool autoDeBias_refX, bool autoDeBias_refY, bool autoDeBias_COMX, bool autoDeBias_COMY, bool autoHanning_ref, bool autoHanning_COM, double& xshift, double& yshift, bool do_parallel);

		/// <summary>Determines the cross correlation lags between two images using the image-reduction-to-vector method where the reference image has already been reduced to X and Y vectors.</summary>
		/// <param name="referenceX">The reference horizontal vector against which to create the cross correlation.</param>
		/// <param name="referenceY">The reference vertical vector against which to create the cross correlation.</param>
		/// <param name="COMPARISON">The comparison data array with which to create the cross correlation.</param>
		/// <param name="autoDeBias_COMX">Option to automatically de-gradient the comparison image along the x-dimension (horizontal degradient).</param>
		/// <param name="autoDeBias_COMY">Option to automatically de-gradient the comparison image along the y-dimension (vertical degradient).</param>
		/// <param name="xshift">The sub-integer x-shift of the comparison with respect to the reference, passed by reference.</param>
		/// <param name="yshift">The sub-integer y-shift of the comparison with respect to the reference, passed by reference.</param>
		/// <param name="do_parallel">Optionally perform all array operations in parallel. False when parallelizing upstream.</param>
		static void XCorrImageLagShifts(array<double>^ referenceX, array<double>^ referenceY, array<double, 2>^ COMPARISON, bool autoDeBias_COMX, bool autoDeBias_COMY, bool autoHanning_COM, double& xshift, double& yshift, bool do_parallel);

		/// <summary>Determines the Curve of Growth photometry for a source centered in the ROI image.</summary>
		/// <param name="ROI">The region of interest image to determine the curve of growth for.</param>
		/// <param name="N_last_fit_pts">The number of tailing points to fit for linear slope - intercept of this line is source counts, slope is the background count per pixel.</param>
		/// <param name="N_points_COG">The number of points for each curve of growth point. Used as the abscissa against the return value.</param>
		/// <param name="background_signal_per_pix">The slope of the linear fit line to the tailing points, i.e., the counts per pixel background.</param>
		/// <param name="source_signal">The intercept of the linear fit line to the tailing points, i.e., the total central source counts.</param>
		static array<double>^ COG(array<double, 2>^ ROI, int N_last_fit_pts, array<double>^ &N_points_COG, double &background_signal_per_pix, double &source_signal);

		static double QuadFit3PtsCenterPos(array<double>^ x, array<double>^ y);
		static array<double>^ QuadFit3PtsParams(array<double>^ x, array<double>^ y);

		static array<double, 1>^ CosineBell(int length);
		static array<double, 1>^ Hanning(array<double, 1>^ data);
		static array<double, 2>^ Hanning(array<double, 2>^ image, bool do_parallel);

		static array<double, 2>^ Bin(array<double, 2>^ data, int Nx, int Ny, bool do_parallel);//remainders are dropped
		static array<unsigned int, 2>^ Bin(array<unsigned int, 2>^ data, int Nx, int Ny, bool do_parallel);//remainders are dropped
		static array<double>^ Bin(array<double>^ data, int Nx);//remainders are dropped

		static array<double>^ Smooth(array<double>^ data, int kernelsize, bool do_parallel);

		static array<double, 2>^ ShiftArrayInt(array<double, 2>^ data, int xshift, int yshift, bool do_parallel);

		/// <summary>Rotates an array about its center.</summary>
		/// <param name="data">The array to rotate.</param>
		/// <param name="radians">The angle to rotate the array, positive counter-clockwise.</param>
		/// <param name="x_center">The rotation center on the x-axis to rotate the array about. Pass Double.MaxValue for array center.</param>
		/// <param name="y_center">The rotation center on the y-axis to rotate the array about. Pass Double.MaxValue for array center.</param>
		/// <param name="style">&quot;nearest&quot; - nearest-neighbor pixel, or, &quot;bilinear&quot; - for 2x2 interpolation, or, &quot;lanc_n&quot; - for Lanczos interpolation of order n = 3, 4, 5.</param>
		static array<double, 2>^ RotateShiftArray(array<double, 2>^ data, double radians, double x_center, double y_center, String^ style, int xshift, int yshift, bool do_parallel);

		//normalized sinc function
		static double sinc_n(double x);

		//un-normalized sinc function
		static double sinc_un(double x);

		static double Lanczos(double x, int n);

		static double InterpolateBiLinear(array<double, 2>^ data, int width, int height, double x, double y);
		static double InterpolateLanczos(array<double, 2>^ data, int width, int height, double x, double y, int n);

		static unsigned __int64 Factorial(unsigned __int64 N);

		/// <summary>Calculates the area of the trangle subtended by 3 [x, y] points.</summary>
		static double Area_Triangle(array<double>^ x, array<double>^ y);

		/// <summary>Returns an interpolation of the specified data at the given interpolation points.</summary>
		/// <param name="xdata">The x-positions of the ydata points to interpolate.</param>
		/// <param name="ydata">The y-values of the data to interpolate.</param>
		/// <param name="xinterp">The x-positions at which to interpolate y-values.</param>
		/// <param name="style">The type of interpolation to compute:
		/// <para> &quot;linear&quot; - linear interpolation</para>
		/// <para> &quot;cubic&quot;  - cubic spline interpolation</para>
		/// <para> &quot;mono&quot;   - monotone cubic spline which preserves monoticity of the data</para>
		/// <para> &quot;catmullrom&quot;   - default Catmull-Rom spline</para>
		/// <para> &quot;akima&quot;  - Akima is a cubic spline which is stable to the outliers, avoiding the oscillations of a cubic spline</para></param>
		static array<double>^ Interpolate1d(array<double>^ xdata, array<double>^ ydata, array<double>^ xinterp, String^ style, bool do_parallel);

		/// <summary>Returns an interpolation of the specified surface data at the given interpolation points with bicubic spline.</summary>
		/// <param name="xdata">The x-positions of the surface points to interpolate. If nullptr is passed a vector will automatically be created of appropriate length.</param>
		/// <param name="ydata">The y-positions of the surface points to interpolate. If nullptr is passed a vector will automatically be created of appropriate length.</param>
		/// <param name="surfdata">The surface data to interpolate.</param>
		/// <param name="xinterpdelta_inv">The inverse of the interpolation delta.  That is, &quot;10&quot; means the grid will interpolated at 1/10th grid scale.</param>
		/// <param name="yinterpdelta_inv">The inverse of the interpolation delta.  That is, &quot;10&quot; means the grid will interpolated at 1/10th grid scale.</param>
		/// <param name="xinterp">The returned interpolated xdata vector. Pass nullptr if not required. If required, must be initialized as an xdata->Length*xinterpdelta_inv length vector.</param>
		/// <param name="yinterp">The returned interpolated ydata vector. Pass nullptr if not required. If required, must be initialized as an ydata->Length*yinterpdelta_inv length vector.</param>
		static array<double, 2>^ Interpolate2d(array<double>^ xdata, array<double>^ ydata, array<double, 2>^ surfdata, int xinterpdelta_inv, int yinterpdelta_inv, array<double>^ xinterp, array<double>^ yinterp, bool do_parallel);

		/// <summary>Returns the 2-D array with gradients removed from a specified dimension.</summary>
		/// <param name="data">The data array to degradient.</param>
		/// <param name="dim">The dimension to degradient: 0 = x, 1 = y.</param>
		static array<double, 2>^ DeGradient(array<double, 2>^ data, int dim, bool do_parallel);

		/// <summary>Returns a circular Gaussian centered on a central pixel in a 2-D array with a given amplitude and Full Width Half Maximum.</summary>
		/// <param name="Amplitude">The amplitude of the Gaussian.</param>
		/// <param name="FWHM">The Full Width Half Maximum of the Gaussian.</param>
		/// <param name="HalfWidth">The half-width or square-radius of the return array at which the Gaussian is calculated.</param>
		static array<double, 2>^ Gaussian(double Amplitude, double FWHM, int HalfWidth, bool do_parallel);

		static array<double, 2>^ Sersic(double Effective_Radius_Re_PIXELS, double Effective_Radius_Io, double SersicFactor_n, int CalculationRadius_NRe);

		static array<double, 2>^ Histogram_IntegerStep(array<double>^ values, int step);
		static array<double, 2>^ Histogram_IntegerDivisions(array<double>^ values, int NDivs);

		/// <summary>Returns true if an integer is even, false if odd.</summary>
		static bool IsEven(int x);

		/// <summary>Returns true if a String can convert to a number, false if it can not.</summary>
		static bool IsNumeric(String^ x);

		/// <summary>Returns true if a number is an integer, false if it is not.</summary>
		static bool IsInteger(double x);

		/// <summary>Returns the dot-product of two equal-length vectors.</summary>
		static double VectorDotProdVector(array<double, 1>^ vector1, array<double, 1>^ vector2, bool do_parallel);

		/// <summary>Returns the sum of a vector with a scalar.</summary>
		static array<double>^ VectorAddScalar(array<double>^ vector, double scalar, bool do_parallel);

		/// <summary>Returns the sum of a vector with a scalar </summary>
		static array<double>^ VectorAddVector(array<double>^ vector1, array<double>^ vector2, bool do_parallel);

		static array<double, 2>^ MatrixAddScalar(array<double, 2>^ matrix, double scalar, bool do_parallel);

		static array<double, 2>^ MatrixAddMatrix(array<double, 2>^ matrix1, array<double, 2>^ matrix2, bool do_parallel);

		static array<double>^ VectorSubScalar(array<double>^ vector, double scalar, bool do_parallel);

		static array<double>^ VectorSubVector(array<double>^ vector1, array<double>^ vector2, bool do_parallel);

		static array<double, 2>^ MatrixSubScalar(array<double, 2>^ matrix, double scalar, bool do_parallel);

		static array<double, 2>^ MatrixSubMatrix(array<double, 2>^ matrix1, array<double, 2>^ matrix2, bool do_parallel);

		static array<double>^ VectorMultScalar(array<double>^ vector, double scalar, bool do_parallel);

		static array<double>^ VectorMultVector(array<double>^ vector1, array<double>^ vector2, bool do_parallel);

		static array<double, 2>^ MatrixMultScalar(array<double, 2>^ matrix, double scalar, bool do_parallel);

		static array<double, 2>^ MatrixMultMatrix(array<double, 2>^ matrix1, array<double, 2>^ matrix2, bool do_parallel);

		static array<double>^ VectorDivScalar(array<double>^ vector, double scalar, bool do_parallel);

		static array<double>^ VectorDivVector(array<double>^ vector1, array<double>^ vector2, bool do_parallel);

		static array<double, 2>^ MatrixDivScalar(array<double, 2>^ matrix, double scalar, bool do_parallel);

		static array<double, 2>^ MatrixDivMatrix(array<double, 2>^ matrix1, array<double, 2>^ matrix2, bool do_parallel);

		//static array<double, 2>^ MatrixMedianFilter(array<double, 2>^ data, int filter_width);

		/// <summary>Filters an 2D array with a median kernel.</summary>
		static array<double, 2>^ MedianFilter(array<double, 2>^ data, int kernelHalfWidth, bool do_parallel);

		/// <summary>Convolves a kernel array into a primary array.  The kernel must have a an odd-numbered with and height.</summary>
		/// <param name="primary">The primary array into which the kernel is convolved.</param>
		/// <param name="kernel">The kernel with which the primary array is convolved.</param>
		static array<double, 2>^ MatrixConvolveMatrix(array<double, 2>^ primary, array<double, 2>^ kernel, bool do_parallel);

		/// <summary>Determines the non-linear least-squares fit parameters for a positively-oriented Gaussian curve G(x|p)
		/// <para>G(x|p) = p(0) * exp( -(x - p(1))^2 / (2*p(2)^2) ) + p(3)</para></summary>
		/// <param name="xdata">The x-data grid positions of the Gaussian data. If nullptr is passed a vector will automatically be created of appropriate size, centered on zero.</param>
		/// <param name="Gdata">The values of the data to be fitted to the Gaussian.</param>
		/// <param name="p">The initial and return parameters of the Gaussian. If p is only initialized and input with all zeros, initial estimates will automatically be computed.
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = sigma; p[3] = bias.</para></param>
		/// <param name="p_lbnd">The lower-bound on the fit parameters. If nullptr is passed they will automatically be set by the Gdata dimensions with allowance.</param>
		/// <param name="p_ubnd">The upper-bound on the fit parameters. If nullptr is passed they will automatically be set by the Gdata dimensions with allowance.</param>
		/// <param name="p_err">The returned errors on the fitted parameters. Pass nullptr if not required.</param>
		/// <param name="fit_residuals">The returned residuals of the fit: Gdata[x] - fit[x].  Pass nullptr if not required.</param>
		static void Fit_Gaussian1d(array<double>^ xdata, array<double>^ Gdata, array<double>^ &p, array<double>^ p_lbnd, array<double>^ p_ubnd, array<double>^ p_err, array<double>^ fit_residuals);

		/// <summary>Computes the elements for a Gaussian curve G(x|p)
		/// <para>G(x|p) = p(0) * exp( -((x - p(1))^2) / (2*p(2)^2) ) + p(3)</para></summary>
		/// <param name="xdata">The x-data grid positions of the Gaussian data. If nullptr is passed a vector will be created of appropriate size, centered on zero.</param>
		/// <param name="G">The values of the data to be computed for the Gaussian.</param>
		/// <param name="p">The parameters of the Gaussian.
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = sigma; p[3] = bias</para></param>
		static void Gaussian1d(array<double>^ xdata, array<double>^ &G, array<double>^ p);

		/// <summary>Determines the non-linear least-squares fit parameters for a positively-oriented 2-d Gaussian surface G(x,y|p)
		/// <para>G(x,y|p) = p(0) * exp(-((x - p(1))^2 + (y - p(2))^2) / (2*p(3)^2)) + p(4)</para>
		/// <para>or</para>
		/// <para>G(x,y|p) = p(0) * exp(-((x - p(1))*cosd(p(3)) + (y - p(2))*sind(p(3)))^2 / (2*p(4)^2) - (-(x - p(1))*sind(p(3)) + (y - p(2))*cosd(p(3))).^2 / (2*p(5)^2) ) + p(6)</para>
		/// <para>The form of G(x,y|p) used is determined by the length of the parameter vector p</para></summary>
		/// <param name="xdata">The x-data grid positions of the Gaussian data.</param>
		/// <param name="ydata">The y-data grid positions of the Gaussian data.</param>
		/// <param name="Gdata">The values of the data to be fitted to the Gaussian.</param>
		/// <param name="p">The initial and return parameters of the Gaussian. If p is only initialized and input with all zeros, initial estimates will automatically be computed. Options are:
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = sigma; p[4] = bias</para>
		/// <para>or</para>
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = phi; p[4] = x-sigma; p[5] = y-sigma; p[6] = bias.</para></param>
		/// <param name="p_LB">The lower bound contraints on the fit parameters. Pass nullptr or an array of length 0 if not required.</param>
		/// <param name="p_UB">The upper bound contraints on the fit parameters. Pass nullptr or an array of length 0 if not required.</param>
		/// <param name="p_err">The return errors on the fitted parameters. Pass an array of length 0 if not required.</param>
		/// <param name="fit_residuals">The return residuals of the fit: Gdata[x, y] - fit[x, y].  Pass an array of length 0 if not required.</param>
		static void Fit_Gaussian2d(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ Gdata, array<double>^ &p, array<double>^ p_LB, array<double>^ p_UB, array<double>^ &p_err, array<double, 2>^ &fit_residuals);

		static void Fit_PointSource(String^ model_name, String^ minimization_type, array<int>^ xdata, array<int>^ ydata, array<double, 2>^ source, array<double>^ &params, array<double>^ params_LB, array<double>^ params_UB, array<double>^ &p_err, array<double, 2>^ &fit_residuals, String^ &termination_msg);

		static void Fit_PointSource_Compound(String^ model_name, String^ minimization_type, array<int>^ xdata, array<int>^ ydata, array<double, 2>^ source, array<double>^ xpositions, array<double>^ ypositions, double position_radius, array<double, 2>^ &params, array<double, 2>^ &p_err, array<double, 2>^ &fit_residuals, String^ &termination_msg);




		/// <summary>Determines the non-linear least-squares fit parameters for a field of n positively-oriented 2-d Gaussian surfaces G(x,y|p_n)
		/// <para>G(x,y|p_n) = Sum[p_n(0) * exp(-((x - p_n(1))^2 + (y - p_n(2))^2) / (2*p_n(3)^2))] + p(4)</para>
		/// <para>or</para>
		/// <para>G(x,y|p_n) =  Sum[p_n(0) * exp(-((x - p_n(1))*cosd(p_n(3)) + (y - p_n(2))*sind(p_n(3)))^2 / (2*p_n(4)^2) - (-(x - p_n(1))*sind(p_n(3)) + (y - p_n(2))*cosd(p_n(3))).^2 / (2*p_n(5)^2) )] + p(6)</para>
		/// <para>The form of G(x,y|p_n) used is determined by the horizontal length of the parameter vector p</para></summary>
		/// <param name="xdata">The x-data grid positions of the Gaussian data.</param>
		/// <param name="ydata">The y-data grid positions of the Gaussian data.</param>
		/// <param name="Gdata">The values of the data to be fitted to the Gaussian.</param>
		/// <param name="p">The initial and return parameters of the Gaussian. If p is only initialized and input with all zeros, initial estimates will automatically be computed. Options are:
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = sigma; p[4] = bias</para>
		/// <para>or</para>
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = phi; p[4] = x-sigma; p[5] = y-sigma; p[6] = bias.</para></param>
		/// <param name="p_LB">The lower bound contraints on the fit parameters.</param>
		/// <param name="p_UB">The upper bound contraints on the fit parameters.</param>
		/// <param name="p_err">The return errors on the fitted parameters. Pass an array of length 0 if not required.</param>
		/// <param name="fit_residuals">The return residuals of the fit: Gdata[x, y] - fit[x, y].  Pass an array of length 0 if not required.</param>
		static void Fit_Gaussian2d_Compound(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ Gdata, array<double, 2>^ &p, array<double, 2>^ p_LB, array<double, 2>^ p_UB, array<double, 2>^ &p_err, array<double, 2>^ &fit_residuals);

		/// <summary>Computes the elements for a 2-d Gaussian surface G(x,y|p)
		/// <para>G(x,y|p) = p(0)*exp( -((x - p(1))^2 + (y - p(2))^2) / (2*p(3)^2) ) + p(4)</para>
		/// <para>or</para>
		/// <para>G(x,y|p) = p(0)*exp( -((x - p(1))*cosd(p(3)) + (y - p(2))*sind(p(3)))^2 / (2*p(4)^2) - (-(x - p(1))*sind(p(3)) + (y - p(2))*cosd(p(3))).^2 / (2*p(5)^2) ) + p(6)</para>
		/// <para>The form of G(x,y|p) used is determined by the length of the parameter vector p</para></summary>
		/// <param name="xdata">The x-data grid positions of the Gaussian data.</param>
		/// <param name="ydata">The y-data grid positions of the Gaussian data.</param>
		/// <param name="G">The values of the data to be computed for the Gaussian.</param>
		/// <param name="p">The parameters of the Gaussian.  Options are:
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = sigma; p[4] = bias</para>
		/// <para>or</para>
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = phi; p[4] = x-sigma; p[5] = y-sigma; p[6] = bias</para></param>
		static void Gaussian2d(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ &G, array<double>^ p, bool do_parallel);

		/// <summary>Determines the non-linear least-squares fit parameters for a 1D Moffat curve M(x|p)
		/// <para>M(x|p) = p(0) * ( 1 + (x - p(1))^2 / p(2)^2 )^(-p(3)) + p(4)</para></summary>
		/// <param name="xdata">The x data grid positions of the Moffat data. If nullptr is passed a vector will automatically be created of appropriate size, centered on zero.</param>
		/// <param name="Mdata">The values of the data to be fitted to the Moffat.</param>
		/// <param name="p">The initial and return parameters of the Moffat. If p is only initialized and input with all zeros, initial estimates will automatically be computed.
		/// <para>p[0] = amplitude, p[1] = x-center, p[2] = theta, p[3] = beta, p[4] = bias.</para></param>
		/// <param name="p_err">The errors on the fitted parameters. Pass nullptr if not required.</param>
		/// <param name="fit_residuals">The residuals of the fit: Mdata[x] - M(x|p). Pass nullptr if not required.</param>
		static void Fit_Moffat1d(array<double>^ xdata, array<double>^ Mdata, array<double>^ &p, array<double>^ p_lbnd, array<double>^ p_ubnd, array<double>^ p_err, array<double>^ fit_residuals);

		/// <summary>Computes the elements for a Moffat curve M(x|p)
		/// <para>M(x|p) = p(0) * ( 1 + (x - p(1))^2 / p(2)^2 )^(-p(3)) + p(4)</para></summary>
		/// <param name="xdata">The x-data grid positions of the Moffat data. If nullptr is passed a vector will be created of appropriate size, centered on zero.</param>
		/// <param name="M">The values of the data to be computed for the Moffat.</param>
		/// <param name="p">The parameters of the Moffat: p[0] = amplitude; p[1] = x-center; p[2] = theta; p[3] = beta; p[4] = bias</param>
		static void Moffat1d(array<double>^ xdata, array<double>^ &M, array<double>^ p);

		/// <summary>Computes the elements for a radial fit Moffat curve to normalized 2D source data of the form M(x|params) = [1 + (x/alpha)^2]^(-beta)</summary>
		/// <param name="radial_x">The radial profile abscissa.</param>
		/// <param name="radial_y">The radial profile ordinate.</param>
		/// <param name="params">The return parameters of the radial Moffat fit. Pass params initialized and input with all zeros, initial estimates will automatically be computed.</param>
		/// <param name="FWHM">User likely wants to know the FWHM of the fit (pixel units).</param>
		/// <param name="interp_radial_x">User likely wants an interpolated (1/10th scale) fit vector result.</param>
		/// <param name="interp_radial_y">User likely wants an interpolated (1/10th scale) fit vector result.</param>
		static void Radial_Profile_Normalized_Fit_Moffat(array<double>^ radial_x, array<double>^ radial_y, array<double>^ &params, double &FWHM, array<double>^& interp_radial_x, array<double>^& interp_radial_y);

		/// <summary>Computes the normalized radial profile.</summary>
		/// <param name="Mdata">The 2D profile to create the radial plot from. Maximum value must be the center pixel, and Mdata array [x, y] size must be odd and square. If the size Mdata is less than 16 elements, the profile is spline-interpolated.</param>
		/// <param name="xdata">The abscissa values for the Mdata array.</param>
		/// <param name="ydata">The ordinate values for the Mdata array.</param>
		/// <param name="axisscale">The unit scale per pixel of the axes, assuming both axes are equal. Pass 0 or 1 for no scale, or any other value greater than zero for scaling.</param>
		/// <param name="radial_x">The radial profile abscissa values (returned).</param>
		/// <param name="radial_y">The radial profile abscissa values (returned).</param>
		static void Radial_Profile_Normalized(array<double, 2>^ Mdata, array<int>^ xdata, array<int>^ ydata, double axisscale, array<double>^& radial_x, array<double>^& radial_y);

		/// <summary>Determines the non-linear least-squares fit parameters for a 2-d Moffat surface M(x,y|p)
		/// <para>M(x,y|p) = p(0) * ( 1 + { (x - p(1))^2 + (y - p(2))^2 } / p(3)^2 )^(-p(4)) + p(5)</para>
		/// <para>or</para>
		/// <para>M(x,y|p) = p(0) * ( 1 + { ((x - p(1))*cosd(p(3)) + (y - p(2))*sind(p(3)))^2 } / p(4)^2 + { (-(x - p(1))*sind(p(3)) + (y - p(2))*cosd(p(3)))^2 } / p(5)^2 )^(-p(6)) + p(7)</para>
		/// <para>The form of M(x,y|p) used is determined by the length of the parameter vector p</para></summary>
		/// <param name="xdata">The x-data grid positions of the Moffat data.</param>
		/// <param name="ydata">The y-data grid positions of the Moffat data.</param>
		/// <param name="Mdata">The values of the data to be fitted to the Moffat.</param>
		/// <param name="p">The initial and return parameters of the Moffat. If p is only initialized and input with all zeros, initial estimates will automatically be computed. Options are:
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = theta; p[4] = beta; p[5] = bias</para>
		/// <para>or</para>
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = phi; p[4] = x-theta; p[5] = y-theta; p[6] = beta; p[7] = bias</para></param>
		/// <param name="p_LB">The lower bound contraints on the fit parameters. Pass nullptr or an array of length 0 if not required.</param>
		/// <param name="p_UB">The upper bound contraints on the fit parameters. Pass nullptr or an array of length 0 if not required.</param>
		/// <param name="p_err">The return errors on the fitted parameters. Pass an array of length 0 if not required.</param>
		/// <param name="fit_residuals">The return residuals of the fit: Mdata[x, y] - fit[x, y].  Pass an array of length 0 if not required.</param>
		static void Fit_Moffat2d(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ Mdata, array<double>^ &p, array<double>^ p_LB, array<double>^ p_UB, array<double>^ &p_err, array<double, 2>^ &fit_residuals);

		/// <summary>Determines the non-linear least-squares fit parameters for a 2-d Moffat surface M(x,y|p)
		/// <para>M(x,y|p_n) = sum[p_n(0) * ( 1 + { (x - p_n(1))^2 + (y - p_n(2))^2 } / p_n(3)^2 )^(-p_n(4))] + p(5)</para>
		/// <para>or</para>
		/// <para>M(x,y|p_n) = sum[p_n(0) * ( 1 + { ((x - p_n(1))*cosd(p_n(3)) + (y - p_n(2))*sind(p_n(3)))^2 } / p_n(4)^2 + { (-(x - p_n(1))*sind(p_n(3)) + (y - p_n(2))*cosd(p_n(3)))^2 } / p_n(5)^2 )^(-p_n(6))] + p_n(7)</para>
		/// <para>The form of M(x,y|p_n) used is determined by the length of the parameter vector p</para></summary>
		/// <param name="xdata">The x-data grid positions of the Moffat data.</param>
		/// <param name="ydata">The y-data grid positions of the Moffat data.</param>
		/// <param name="Mdata">The values of the data to be fitted to the Moffat.</param>
		/// <param name="p">The initial and return parameters of the Moffat. If p is only initialized and input with all zeros, initial estimates will automatically be computed. Options are:
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = theta; p[4] = beta; p[5] = bias</para>
		/// <para>or</para>
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = phi; p[4] = x-theta; p[5] = y-theta; p[6] = beta; p[7] = bias</para></param>
		/// <param name="p_LB">The lower bound contraints on the fit parameters.</param>
		/// <param name="p_UB">The upper bound contraints on the fit parameters.</param>
		/// <param name="p_err">The return errors on the fitted parameters. Pass an array of length 0 if not required.</param>
		/// <param name="fit_residuals">The return residuals of the fit: Mdata[x, y] - fit[x, y].  Pass an array of length 0 if not required.</param>
		static void Fit_Moffat2d_Compound(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ Mdata, array<double, 2>^ &p, array<double, 2>^ p_LB, array<double, 2>^ p_UB, array<double, 2>^ &p_err, array<double, 2>^ &fit_residuals);

		/// <summary>Computes the elements for a 2-d Moffat surface M(x,y|p)
		/// <para>M(x,y|p) = p(0) * ( 1 + { (x-p(1))^2 + (y-p(2))^2 } / p(3)^2 ) ^ (-p(4)) + p(5)</para>
		/// <para>or</para>
		/// <para>M(x,y|p) = p(0) * ( 1 + { ((x-p(1))*cosd(p(3)) + (y-p(2))*sind(p(3)))^2 } / p(4)^2 + { (-(x-p(1))*sind(p(3)) + (y-p(2))*cosd(p(3)))^2 } / p(5)^2 ) ^ (-p(6)) + p(7)</para>
		/// <para>The form of M(x,y|p) used is determined by the length of the parameter vector p</para></summary>
		/// <param name="xdata">The x-data grid positions of the Moffat data.</param>
		/// <param name="ydata">The y-data grid positions of the Moffat data.</param>
		/// <param name="M">The values of the data to be computed for the Moffat.</param>
		/// <param name="p">The parameters of the Moffat. Options are:
		/// <para>p[0] = amplitude, p[1] = x-center, p[2] = y-center, p[3] = theta, p[4] = beta, p[5] = bias</para>
		/// <para>or</para>
		/// <para>p[0] = amplitude, p[1] = x-center, p[2] = y-center, p[3] = phi, p[4] = x-theta, p[5] = y-theta, p[6] = beta, p[7] = bias</para></param>
		static void Moffat2d(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ &M, array<double>^ p, bool do_parallel);

		/// <summary>Computes the 2-D transformation elements between intermediate catalogue coordinates and image pixel coordinates.</summary>
		/// <param name="x_intrmdt">The x-axis reference points of which to determine the transformation to.</param>
		/// <param name="y_intrmdt">The y-axis reference points of which to determine the transformation to.</param>
		/// <param name="x_pix">The x-axis points for which to determine the transformation of.</param>
		/// <param name="y_pix">The y-axis points for which to determine the transformation of.</param>
		/// <param name="p">The initial and return parameters of the tranformation. Options are
		/// <para>p[0] = scale, p[1] = phi (radians), p[2] = x-axis pixel coordinate reference, p[3] = y-axis pixel coordinate reference,</para>
		///<para>or</para>
		///<para>p[0] = Matrix coeff [0, 0], p[1] = Matrix coeff [1, 0], p[2] = Matrix coeff [0, 1], p[3] = Matrix coeff [1, 1], p[4] = x-axis pixel coordinate reference, p[5] = y-axis pixel coordinate reference,</para></param>
		/// <param name="p_lbnd">The lower-bound on the fit parameters.</param>
		/// <param name="p_ubnd">The upper-bound on the fit parameters.</param>
		/// <param name="p_scale">The order of magnitude scale (positive) of the fit parameters.</param>
		static void Fit_WCSTransform2d(array<double>^ x_intrmdt, array<double>^ y_intrmdt, array<double>^ x_pix, array<double>^ y_pix, array<double>^ &p, array<double>^ p_lbnd, array<double>^ p_ubnd, array<double>^ p_scale);

		/// <summary>Computes the 2-D transformation elements between two sets of coordinates.</summary>
		/// <param name="x_ref">The x-axis reference points of which to determine the transformation to.</param>
		/// <param name="y_ref">The y-axis reference points of which to determine the transformation to.</param>
		/// <param name="x_tran">The x-axis points for which to determine the transformation of.</param>
		/// <param name="y_tran">The y-axis points for which to determine the transformation of.</param>
		/// <param name="p">The initial and return parameters of the tranformation. Options are
		/// <para>p[0] = scale, p[1] = phi (radians), p[2] = x-tran pixel coordinate rotation reference, p[3] = y-tran pixel coordinate rotation reference, p[4] = x-tran pixel coordinate shift, p[5] = x-tran pixel coordinate shift</para>
		///<para>or</para>
		///<para>p[0] = Matrix coeff [0, 0], p[1] = Matrix coeff [1, 0], p[2] = Matrix coeff [0, 1], p[3] = Matrix coeff [1, 1], p[4] = x-tran pixel coordinate rotation reference, p[5] = y-tran pixel coordinate rotation reference, p[6] = x-tran pixel coordinate shift, p[7] = x-tran pixel coordinate shift</para></param>
		/// <param name="p_lbnd">The lower-bound on the fit parameters.</param>
		/// <param name="p_ubnd">The upper-bound on the fit parameters.</param>
		/// <param name="p_scale">The order of magnitude scale (positive) of the fit parameters.</param>
		static void Fit_GeneralTransform2d(array<double>^ x_ref, array<double>^ y_ref, array<double>^ x_tran, array<double>^ y_tran, array<double>^ &p, array<double>^ p_lbnd, array<double>^ p_ubnd, array<double>^ p_scale);

		/// <summary>Fits a polynomial to x, y data.</summary>
		/// <param name="xdata">The x-axis data points.</param>
		/// <param name="ydata">The y-axis data points.</param>
		/// <param name="poly_degree">The degree of polynomial to fit: 1 = linear, 2 = quadratic, etc.</param>
		/// <param name="robust">If true, weights will automatically be determined which supress outliers.</param>
		/// <param name="poly_coeffs">The coefficients of the polynomial ordered by increasing power.</param>
		static void Fit_Poly1d(array<double>^ xdata, array<double>^ ydata, int poly_degree, bool robust, array<double>^ &poly_coeffs);

	private:
		static double MedianSTD(double *np, int len);

		static void alglib_Gauss_1d(array<double>^ p, array<double>^ x, double %val, Object^ obj);
		delegate void function_Gauss_1d_delegate(array<double>^ p, array<double>^ x, double %val, Object^ obj);

		static void alglib_Gauss_1d_grad(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj);
		delegate void alglib_Gauss_1d_grad_delegate(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj);

		static void alglib_Moffat_1d(array<double>^ p, array<double>^ x, double %val, Object^ obj);
		delegate void alglib_Moffat_1d_delegate(array<double>^ p, array<double>^ x, double %val, Object^ obj);

		static void alglib_Moffat_1d_grad(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj);
		delegate void alglib_Moffat_1d_grad_delegate(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj);

		///// <summary>Calculates a single point of a 2-d Gaussian surface G(x,y|p)
		///// <para>G(x,y|p) = p(0)*exp(-((x-p(1))^2 + (y - p(2))^2)/(2*p(3)^2)) + p(4)</para>
		///// <para>or</para>
		///// <para>G(x,y|p) = p(0)*exp(-((x-p(1))*cos(p(3)) + (y-p(2))*sin(p(3)))^2 / (2*p(4)^2) - (-(x-p(1))*sin(p(3)) + (y-p(2))*cos(p(3))).^2 / (2*p(5)^2) ) + p(6)</para>
		///// <para>where x[0] is a position on X-axis x, and x[1] is a position on Y-axis y.</para>
		///// <para>The form of G(x,y|p) used is determined by the length of the parmater vector p</para></summary>
		///// <param name="p">The initial parameters of the Gaussian fit.  Options are:
		///// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = sigma; p[4] = bias</para>
		///// <para>or</para>
		///// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = theta; p[4] = x-sigma; p[5] = y-sigma; p[6] = bias</para></param>
		///// <param name="x">The x,y position to calculate the value val of the Gaussian G(x,y|p): x[0] = x, x[1] = y</param>
		///// <param name="val">The calculated value of the Gaussian.</param>
		///// <param name="obj">obj.</param>
		static void alglib_Gauss_2d(array<double>^ p, array<double>^ x, double %val, Object^ obj);
		delegate void function_Gauss_2d_delegate(array<double>^ p, array<double>^ x, double %val, Object^ obj);

		static void alglib_Gauss_2d_grad(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj);
		delegate void alglib_Gauss_2d_grad_delegate(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj);







		static void alglib_Gauss_2d_LM_LS_grad(array<double>^ p, double %f, array<double>^ grad, Object^ obj);
		delegate void function_Gauss_2d_LM_LS_grad_delegate(array<double>^ p, double %f, array<double>^ grad, Object^ obj);

		static void alglib_Gauss_2d_LM_LS_grad_compound(array<double>^ p, double %f, array<double>^ grad, Object^ obj);
		delegate void function_Gauss_2d_LM_LS_grad_compound_delegate(array<double>^ p, double %f, array<double>^ grad, Object^ obj);

		static void alglib_Gauss_2d_LM_LS_CHISQ_grad(array<double>^ p, double% f, array<double>^ grad, Object^ obj);
		delegate void function_Gauss_2d_LM_LS_CHISQ_grad_delegate(array<double>^ p, double% f, array<double>^ grad, Object^ obj);

		static void alglib_Gauss_2d_LM_LS_CHISQ_grad_compound(array<double>^ p, double% f, array<double>^ grad, Object^ obj);
		delegate void function_Gauss_2d_LM_LS_CHISQ_grad__compounddelegate(array<double>^ p, double% f, array<double>^ grad, Object^ obj);

		static void alglib_Gauss_2d_LM_LS_ROBUST_grad(array<double>^ p, double% f, array<double>^ grad, Object^ obj);
		delegate void function_Gauss_2d_LM_LS_ROBUST_grad_delegate(array<double>^ p, double% f, array<double>^ grad, Object^ obj);

		static void alglib_Gauss_2d_LM_LS_ROBUST_grad_compound(array<double>^ p, double% f, array<double>^ grad, Object^ obj);
		delegate void function_Gauss_2d_LM_LS_ROBUST_grad_compound_delegate(array<double>^ p, double% f, array<double>^ grad, Object^ obj);

		static void alglib_Gauss_2d_LM_LS_CSTAT_grad(array<double>^ p, double% f, array<double>^ grad, Object^ obj);
		delegate void function_Gauss_2d_LM_LS_CSTATS_grad_delegate(array<double>^ p, double% f, array<double>^ grad, Object^ obj);

		static void alglib_Gauss_2d_LM_LS_CSTAT_grad_compound(array<double>^ p, double% f, array<double>^ grad, Object^ obj);
		delegate void function_Gauss_2d_LM_LS_CSTATS_grad_compound_delegate(array<double>^ p, double% f, array<double>^ grad, Object^ obj);

		




		

		

		static void alglib_Moffat_2d_LM_LS_grad(array<double>^ p, double %f, array<double>^ grad, Object^ obj);
		delegate void function_Moffat_2d_LM_LS_grad_delegate(array<double>^ p, double %f, array<double>^ grad, Object^ obj);

		static void alglib_Moffat_2d_LM_LS_CHISQ_grad(array<double>^ p, double %f, array<double>^ grad, Object^ obj);
		delegate void function_Moffat_2d_LM_LS_CHISQ_grad_delegate(array<double>^ p, double %f, array<double>^ grad, Object^ obj);

		static void alglib_Moffat_2d_LM_LS_ROBUST_grad(array<double>^ p, double %f, array<double>^ grad, Object^ obj);
		delegate void function_Moffat_2d_LM_LS_ROBUST_grad_delegate(array<double>^ p, double %f, array<double>^ grad, Object^ obj);

		static void alglib_Moffat_2d_LM_LS_CSTAT_grad(array<double>^ p, double %f, array<double>^ grad, Object^ obj);
		delegate void function_Moffat_2d_LM_LS_CSTATS_grad_delegate(array<double>^ p, double f, array<double>^ grad, Object^ obj);












		static array<double>^ Gauss_2D_param_err(array<double>^ p, array<double>^ x, array<double>^ y, array<double, 2>^ Z);
		static array<double>^ Moffat_2D_param_err(array<double>^ p, array<double>^ x, array<double>^ y, array<double, 2>^ Z);

		static void alglib_Gauss_2d_compound(array<double>^ p, array<double>^ x, double %val, Object^ obj);
		delegate void function_Gauss_2d_compound_delegate(array<double>^ p, array<double>^ x, double %val, Object^ obj);

		static void alglib_Gauss_2d_compound_grad(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj);
		delegate void alglib_Gauss_2d_compound_grad_delegate(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj);

		/// <summary>Calculates a single point of a 2-d Moffat surface M(x,y|p)
		/// <para>M(x,y|p) = p(0) * ( 1 + { (x-p(1))^2 + (y-p(2))^2 } / p(3)^2 ) ^ (-p(4)) + p(5)</para>
		/// <para>or</para>
		/// <para>M(x,y|p) = p(0) * ( 1 + { ((x-p(1))*cos(p(3)) + (y-p(2))*sin(p(3)))^2 } / p(4)^2 + { (-(x-p(1))*sin(p(3)) + (y-p(2))*cos(p(3)))^2 } / p(5)^2 ) ^ (-p(6)) + p(7)</para>
		/// <para>where x[0] is a position on X-axis x, and x[1] is a position on Y-axis y.</para>
		/// <para>The form of M(x,y|p) used is determined by the length of the parmater vector "p"</para></summary>
		/// <param name="p">The initial parameters of the Moffat fit.  Options are:
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = theta; p[4] = beta; p[5] = bias</para>
		/// <para>or</para>
		/// <para>p[0] = amplitude; p[1] = x-center; p[2] = y-center; p[3] = phi; p[4] = x-theta; p[5] = y-theta; p[6] = beta; p[7] = bias</para></param>
		/// <param name="x">The x,y position to calculate the value val of the Moffat M(x,y|p): x[0] = x, x[1] = y</param>
		/// <param name="val">The calculated value of the Moffat.</param>
		/// <param name="obj">obj.</param>
		static void alglib_Moffat_2d(array<double>^ p, array<double>^ x, double %val, Object^ obj);
		delegate void alglib_Moffat_2d_delegate(array<double>^ p, array<double>^ x, double %val, Object^ obj);

		static void alglib_Moffat_2d_grad(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj);
		delegate void alglib_Moffat_2d_grad_delegate(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj);

		static void alglib_Moffat_2d_compound(array<double>^ p, array<double>^ x, double %val, Object^ obj);
		delegate void alglib_Moffat_2d_compound_delegate(array<double>^ p, array<double>^ x, double %val, Object^ obj);

		static void alglib_Moffat_2d_compound_grad(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj);
		delegate void alglib_Moffat_2d_grad_compound_delegate(array<double>^ p, array<double>^ x, double %val, array<double>^ grad, Object^ obj);

		static void alglib_WCSTransform2d_fvec(array<double>^ p, array<double>^ f, Object^ obj);
		delegate void alglib_WCSTransform2d_fvec_delegate(array<double>^ p, array<double>^ f, Object^ obj);

		static void alglib_WCSTransform2d_jac(array<double>^ p, array<double>^ f, array<double, 2>^ jac, Object^ obj);
		delegate void alglib_WCSTransform2d_jac_delegate(array<double>^ p, array<double>^ f, array<double, 2>^ jac, Object^ obj);

		static void alglib_GeneralTransform2d_fvec(array<double>^ p, array<double>^ f, Object^ obj);
		delegate void alglib_GeneralTransform2d_fvec_delegate(array<double>^ p, array<double>^ f, Object^ obj);
	};

	/// <summary>SourceExtractor class provides functionality for extracting sources from image arrays.</summary>
	public ref class SourceExtractor
	{
		public:
		SourceExtractor();
		SourceExtractor(array<double>^ XCoords, array<double>^ YCoords);
		SourceExtractor(JPFITS::FITSBinTable^ BinTablePSE);
		~SourceExtractor();

		/// <summary>Gets a metadata table of the extracted sources.</summary>
		property array<String^, 2>^ Source_Table
		{
			array<String^, 2>^ get()
			{
				this->GENERATEPSETABLE();
				return this->PSE_TABLE;
			}
		}

		/// <summary>Gets or Sets the x-axis centroids of extracted sources.</summary>
		property array<double>^ Centroids_X
		{
			array<double>^ get() { return CENTROIDS_X; }
			void set(array<double>^ xcentroids) 
			{
				CENTROIDS_X = xcentroids;
				N_SRC = CENTROIDS_X->Length;
				CENTROID_POINTS = gcnew array<JPMath::PointD^>(N_SRC);
				for (int i = 0; i < N_SRC; i++)
					CENTROID_POINTS[i] = gcnew JPMath::PointD(CENTROIDS_X[i], CENTROIDS_Y[i], CENTROIDS_AMPLITUDE[i]);
			}
		}

		/// <summary>Gets or Sets the y-axis centroids of extracted sources.</summary>
		property array<double>^ Centroids_Y
		{
			array<double>^ get() { return CENTROIDS_Y; }
			void set(array<double>^ ycentroids)
			{
				CENTROIDS_Y = ycentroids;
				N_SRC = CENTROIDS_Y->Length;
				CENTROID_POINTS = gcnew array<JPMath::PointD^>(N_SRC);
				for (int i = 0; i < N_SRC; i++)
					CENTROID_POINTS[i] = gcnew JPMath::PointD(CENTROIDS_X[i], CENTROIDS_Y[i], CENTROIDS_AMPLITUDE[i]);
			}
		}

		property array<double>^ Centroids_Volume
		{
			array<double>^ get() { return CENTROIDS_VOLUME; }
		}

		/// <summary>Gets the total number of extracted sources.</summary>
		property int N_Sources
		{
			int get() { return N_SRC; }
		}

		/// <summary>Gets the number of saturated sources.</summary>
		property int N_SaturatedSources
		{
			int get() { return N_SATURATED; }
		}

		/// <summary>Gets a list of the fitted parameters for all sources.</summary>
		property array<double, 2>^ Fitted_Parameter_List
		{
			array<double, 2>^ get() { return FITS_PARAMS; }
		}

		/// <summary>Gets a String^ of the equation used for least-squares fitting.</summary>
		property String^ LSFit_Equation
		{
			String^ get() { return FIT_EQUATION; }
		}

		/// <summary>Gets a boolean to indicate whether least-squares fits have been performed.</summary>
		property bool Fitted
		{
			bool get() { return FITTED; }
		}

		/// <summary>Returns the boolean source map.</summary>
		property array<bool, 2>^ SourceBooleanMap
		{
			array<bool, 2>^ get() { return SOURCE_BOOLEAN_MAP; }
			void set(array<bool, 2>^ boolmap) { SOURCE_BOOLEAN_MAP = boolmap; }
		}

		/// <summary>Returns the integer index source map.</summary>
		property array<int, 2>^ SourceIndexMap
		{
			array<int, 2>^ get() { return SOURCE_INDEX_MAP; }
			void set(array<int, 2>^ indexmap) { SOURCE_INDEX_MAP = indexmap; }
		}

		property bool IsBusy
		{
			bool get() { return BGWRKR->IsBusy; }
		}

		property double PixelSaturation
		{
			double get() { return PIX_SAT; }
		}

		property double KernelRadius
		{
			double get() { return KERNEL_RADIUS; }
		}

		property double SourceSeparation
		{
			double get() { return SOURCE_SEPARATION; }
		}

		property double PixelMaximum
		{
			double get() { return PIX_MAX; }
		}

		property double PixelMinimum
		{
			double get() { return PIX_MIN; }
		}

		property double KernelMaximum
		{
			double get() { return KERNEL_MAX; }
		}

		property double KernelMinimum
		{
			double get() { return KERNEL_MIN; }
		}

		property bool AutoBackground
		{
			bool get() { return AUTO_BG; }
		}

		property bool SavePointSources
		{
			bool get() { return SAVE_PS; }
		}

		/*property String^ SavePointSourcesFileNameTemplate
		{
			String^ get() { return SAVE_PS_FILENAME; }
		}*/

		property bool SearchROI
		{
			bool get() { return SEARCH_ROI; }
		}

		property bool PSEParametersSet
		{
			bool get() { return PSEPARAMSSET; }
		}

		property int NGroups
		{
			int get() { return NGROUPS; }
		}

		/// <summary>Searches for sources withn a 2D image array.</summary>
		/// <param name="image">The 2D image array to find sources in.</param>
		/// <param name="pix_saturation">The saturation threshold of of the image pixels, for finding saturation islands. Set equal to zero (0) if not needed.</param>
		/// <param name="pix_min">The minimum pixel threshold value (or SN) to consider a potential source.</param>
		/// <param name="pix_max">The maximum pixel threshold value (or SN) to consider a potential source.</param>
		/// <param name="kernel_min">The minimum kernel pixel sum threshold value (or SN) to consider a potential source.</param>
		/// <param name="kernel_max">The maximum kernel pixel sum threshold value (or SN) to consider a potential source.</param>
		/// <param name="threshholds_as_SN">Treat the thresholds as Signal to Noise instead of pixel values.</param>
		/// <param name="kernel_radius">The radius (pixels) of the kernel to find sources within. Secondary sources within the radius will be ignored.</param>
		/// <param name="source_separation">The separation (pixels) between sources. Only the brightest source within the separation radius is kept.</param>
		/// <param name="auto_background">Automatically determine the local background for potential sources.  Not required if background is known to be zeroed, but should have no effect if used in this case.</param>
		/// <param name="kernel_filename_template">The template full file name for the kernels to be saved. Sources will be numbered sequentially. Pass nullptr or empty string for no saving.</param>
		/// <param name="ROI_region">A boolean array of valid area to examine. Pass nullptr or array of equal dimension to source image all true for entire image search.</param>
		/// <param name="show_waitbar">Show a cancellable wait bar.</param>
		void Extract_Sources(array<double, 2>^ image, double pix_saturation, double pix_min, double pix_max, double kernel_min, double kernel_max, bool threshholds_as_SN, int kernel_radius, int source_separation, bool auto_background, String^ kernel_filename_template, array<bool, 2>^ ROI_region, bool show_waitbar);

		/// <summary>Determines centroids and other kernel information for known sources at given coordinates.</summary>
		/// <param name="image">The 2D image array containing the known sources to extract.</param>
		/// <param name="XCoords">The x-axis coordinates of the sources.</param>
		/// <param name="YCoords">The y-axis coordinates of the sources.</param>
		/// <param name="kernel_radius">The radius (pixels) of the kernel to centroid.</param>
		/// <param name="auto_background">Automatically determine the local background for potential sources.  Not required if background is known to be zeroed, but should have no effect if used in this case.</param>
		/// <param name="kernel_filename_template">The template full file name for the kernels to be saved. Sources will be numbered sequentially. Pass nullptr or empty string for no saving.</param>
		void Extract_Sources(array<double, 2>^ image, array<double>^ XCoords, array<double>^ YCoords, int kernel_radius, bool auto_background, String^ kernel_filename_template);

		void Extract_Attempt_N_Sources(int N, array<double, 2>^ image, double pix_saturation, double pix_min, double pix_max, double kernel_min, double kernel_max, bool threshholds_as_SN, int kernel_radius, int source_separation, bool auto_background, String^ kernel_filename_template, array<bool, 2>^ ROI_region, bool show_waitbar);

		/// <summary>Performs a least-squares fit on all sources of the form:
		/// <para>G(x,y|P) = P(0) * exp( -((x - P(1)).^2 + (y - P(2)).^2 ) / (2*P(3)^2)) + P(4).</para></summary>
		void Extract_Source_LSFits_Gaussian_Circular(array<double>^ Pinit, array<double>^ LBnds, array<double>^ UBnds/*bool view,*/);//2-D Circular Gaussian

		/// <summary>Performs a least-squares fit on all sources of the form:
		/// <para>G(x,y|P) = P(0) * exp( -((x - P(1))*cosd(P(3)) + (y - P(2))*sind(P(3))).^2 / (2*P(4)^2) - ( -(x - P(1))*sind(P(3)) + (y - P(2))*cosd(P(3))).^2 / (2*P(5)^2) ) + P(6).</para></summary>
		void Extract_Source_LSFits_Gaussian_Elliptical(array<double>^ Pinit, array<double>^ LBnds, array<double>^ UBnds/*bool view,*/);// 2-D Elliptical Gaussian

		/// <summary>Performs a least-squares fit on all sources of the form:
		/// <para>M(x,y|P) = P(0) * ( 1 + { (x - P(1))^2 + (y - P(2))^2 } / P(3)^2 ) ^ (-P(4)) + P(5).</para></summary>
		void Extract_Source_LSFits_Moffat_Circular(array<double>^ Pinit, array<double>^ LBnds, array<double>^ UBnds/*bool view,*/);// 2-D Circular Moffat

		/// <summary>Performs a least-squares fit on all sources of the form:
		/// <para>M(x,y|P) = P(0) * (1 + { ((x - P(1))*cosd(P(3)) + (y - P(2))*sind(P(3))) ^ 2 } / P(4) ^ 2 + { (-(x - P(1))*sind(P(3)) + (y - P(2))*cosd(P(3))) ^ 2 } / P(5) ^ 2) ^ (-P(6)) + P(7).</para></summary>
		void Extract_Source_LSFits_Moffat_Elliptical(array<double>^ Pinit, array<double>^ LBnds, array<double>^ UBnds/*bool view,*/);// 2-D Elliptical Moffat

		/// <summary>Saves the metadata table of the extracted sources as a delimited text file.</summary>
		/// <param name="delimit">The delimit argument string: &quot;tab&quot; specifies a tab-delimit, otherwise provide a character (such as the comma &quot;,&quot; etc).</param>
		void Save_Source_Table(String^ delimit);

		/// <summary>Generates RA and Dec coordinates for the sources in this instance, using the supplied World Coordinate System instance.</summary>
		/// <param name="wcs">The world coordinate system to use for converting image pixel locations to world coordinates.</param>
		void Generate_Source_RADec_Coords(JPFITS::WorldCoordinateSolution^ wcs);

		/// <summary>Gets a sub-array kernel from a primary image given a center position and square half-width radius.</summary>
		static array<double, 2>^ GetKernel(array<double, 2>^ image, int x0, int y0, int radius);

		/// <summary>Determines the [x, y] centroid location of a given kernel.</summary>
		static void Centroid(array<int>^ xdata, array<int>^ ydata, array<double, 2>^ kernel, double &x_centroid, double& y_centroid);

		void ClipToNBrightest(int NBright);

		void GroupizePSE(double groupRadius);

		ref class Group
		{
			public:
			Group(int id)
			{
				ID = id;
			}

			int ID = -1;			
			int NElements = 0;
			array<int>^ ElementIndices;
			Region^ REGION;
			Color COLOR;
		};

		property array<Group^>^ Groups
		{
			array<Group^>^ get() { return GROUPS; }
		}

		property array<int>^ GroupIDs
		{
			array<int>^ get() { return GROUPIDS; }
		}

		property array<int, 2>^ SourceGroupMap
		{
			array<int, 2>^ get() { return SOURCE_GROUP_MAP; }
		}
		
		private:

		JPWaitBar::WaitBar^ WAITBAR;
		BackgroundWorker^ BGWRKR;
		Object^ BGWRKR_RESULT;
		void BGWRKR_DoWork(System::Object^  sender, System::ComponentModel::DoWorkEventArgs^  e);
		void BGWRKR_ProgressChanged(System::Object^  sender, System::ComponentModel::ProgressChangedEventArgs^  e);
		void BGWRKR_RunWorkerCompleted(System::Object^  sender, System::ComponentModel::RunWorkerCompletedEventArgs^  e);
		void GENERATEPSETABLE();
		double ESTIMATELOCALBACKGROUND(int x, int y, int HW);
		void MAPSATURATIONISLAND(int X, int Y, int sourceindex, int &xmin, int &xmax, int &ymin, int &ymax);
		bool SAFETOMAPSATURATION(int x, int y);
		void DEMAP(int X, int Y, int sourceindex);
		bool SAFETODEMAP(int x, int y, int sourceindex);
		void REMAP(int X, int Y, int oldindex, int newindex);
		bool SAFETOREMAP(int x, int y, int currentindex);
		void INITARRAYS();
		void INITTHIS();

		void RECURSGROUP(int I, int J, double groupRadius, int groupid, ArrayList^ groups);

		bool FITTED = false;
		String^ FIT_EQUATION;
		bool WCS_GENERATED = false;
		bool VIEWFITS = false;
		bool PSEPARAMSSET = false;

		int IMAGEWIDTH, IMAGEHEIGHT;
		int N_SRC = 0;							// number of sources found
		int N_SATURATED = 0;
		int KERNEL_RADIUS;						// source widths for centroiding
		int KERNEL_WIDTH;
		int SOURCE_SEPARATION;
		double PIX_SAT;							//pixel saturation
		double PIX_MIN;							// source min pixel thresh
		double PIX_MAX;							// source max pixel thresh
		double KERNEL_MIN;						// total source min count thresh
		double KERNEL_MAX;						// total source max count thresh
		bool AUTO_BG;							//automatic background determination (corner min method)
		bool SAVE_PS;
		bool THRESHHOLDS_AS_SN;					//interpret pixel value and total count thresholds as SN
		bool SEARCH_ROI;
		array<bool, 2>^ ROI_REGION;
		bool SHOWWAITBAR;
		int NGROUPS = 0;
		
		String^ SAVE_PS_FILENAME;
		array<double, 2>^ IMAGE;
		array<bool, 2>^ SOURCE_BOOLEAN_MAP;
		array<int, 2>^ SOURCE_INDEX_MAP;
		array<int, 2>^ SOURCE_GROUP_MAP;

		array<int>^ GROUPIDS;
		array<Group^>^ GROUPS;
		array<JPMath::PointD^>^ CENTROID_POINTS;
		array<double>^ CENTROIDS_X;				// x centroid positions of sources
		array<double>^ CENTROIDS_Y;				// y centroid positions of sources
		array<double>^ CENTROIDS_RA_DEG;		// right ascension centroid positions of sources - if available
		array<String^>^ CENTROIDS_RA_HMS;		// right ascension centroid positions of sources - if available
		array<double>^ CENTROIDS_DEC_DEG;		// declination centroid positions of sources - if available
		array<String^>^ CENTROIDS_DEC_DMS;		// declination centroid positions of sources - if available
		array<double>^ CENTROIDS_AMPLITUDE;		// sources values (above fcmin)
		array<double>^ CENTROIDS_VOLUME;		// sources energies (above fcmin)
		array<double>^ CENTROIDS_BGESTIMATE;	// corner minimum - estimate of background

		array<double>^ FITS_X;					// x fitted positions of sources
		array<double>^ FITS_Y;					// y fitted positions of sources
		array<double>^ FITS_FWHM_X;				// FWHM of sources
		array<double>^ FITS_FWHM_Y;				// FWHM of sources
		array<double>^ FITS_PHI;				// rotation theta of elliptical fits
		array<double>^ FITS_RA_DEG;				// right ascension centroid positions of sources - if available
		array<String^>^ FITS_RA_HMS;			// right ascension centroid positions of sources - if available
		array<double>^ FITS_DEC_DEG;			// declination centroid positions of sources - if available
		array<String^>^ FITS_DEC_DMS;			// declination centroid positions of sources - if available
		array<double, 2>^ FITS_PARAMS;			// fitted paramaters of sources - 2d because multiple parameters per source
		array<double>^ FITS_AMPLITUDE;			// 
		array<double>^ FITS_VOLUME;				// 
		array<double>^ FITS_BGESTIMATE;			// 
		array<double>^ FITS_CHISQNORM;			//

		array<String^, 2>^ PSE_TABLE;

		array<double>^ LBND;
		array<double>^ UBND;
		array<double>^ PINI;
	};

	/// <summary>WCS_AutoSolver class provides functionality for automatically solving astrometric solutions for FITS image data.</summary>
	public ref class WCS_AutoSolver
	{
	public:

		ref class WCSAutoSolverReportingForm : public System::Windows::Forms::Form
		{
		public:
			/*WCSAutoSolverReportingForm(void)
			{
				InitializeComponent();
				//
				//TODO: Add the constructor code here
				//
			}*/

			WCSAutoSolverReportingForm(JPFITS::WCS_AutoSolver^ solver)
			{
				InitializeComponent();

				WCSAS = solver;
			}

		protected:
			/// <summary>
			/// Clean up any resources being used.
			/// </summary>
			~WCSAutoSolverReportingForm()
			{
				if (components)
				{
					delete components;
				}
			}
		public: System::Windows::Forms::TextBox^ MsgTxt;
		protected:
		public: System::Windows::Forms::Button^ CancelBtn;
		public: System::Windows::Forms::Timer^ WCSAutoReportingTimer;
		public: System::ComponentModel::IContainer^ components;
		public:

			JPFITS::WCS_AutoSolver^ WCSAS;




			#pragma region Windows Form Designer generated code
			/// <summary>
			/// Required method for Designer support - do not modify
			/// the contents of this method with the code editor.
			/// </summary>
			void InitializeComponent(void)
			{
				this->components = (gcnew System::ComponentModel::Container());
				this->MsgTxt = (gcnew System::Windows::Forms::TextBox());
				this->CancelBtn = (gcnew System::Windows::Forms::Button());
				this->WCSAutoReportingTimer = (gcnew System::Windows::Forms::Timer(this->components));
				if (File::Exists("C:\\Program Files\\Astrowerks\\CCDLAB\\trigonometric-algorithm-icon.ico"))
					this->Icon = gcnew System::Drawing::Icon("C:\\Program Files\\Astrowerks\\CCDLAB\\trigonometric-algorithm-icon.ico");
				else if (File::Exists("D:\\Documents\\Visual Studio 2019\\Projects 2019\\Icon Resources\\trigonometric-algorithm-icon.ico"))
					this->Icon = gcnew System::Drawing::Icon("D:\\Documents\\Visual Studio 2019\\Projects 2019\\Icon Resources\\trigonometric-algorithm-icon.ico");
				this->ShowIcon = true;
				this->SuspendLayout();
				// 
				// MsgTxt
				// 
				this->MsgTxt->AcceptsReturn = true;
				this->MsgTxt->Location = System::Drawing::Point(11, 11);
				this->MsgTxt->Margin = System::Windows::Forms::Padding(2);
				this->MsgTxt->Multiline = true;
				this->MsgTxt->Name = L"MsgTxt";
				this->MsgTxt->ReadOnly = true;
				this->MsgTxt->ScrollBars = System::Windows::Forms::ScrollBars::Vertical;
				this->MsgTxt->Size = System::Drawing::Size(423, 471);
				this->MsgTxt->TabIndex = 4;
				// 
				// CancelBtn
				// 
				this->CancelBtn->DialogResult = System::Windows::Forms::DialogResult::Cancel;
				this->CancelBtn->Location = System::Drawing::Point(439, 11);
				this->CancelBtn->Margin = System::Windows::Forms::Padding(2);
				this->CancelBtn->Name = L"CancelBtn";
				this->CancelBtn->Size = System::Drawing::Size(59, 47);
				this->CancelBtn->TabIndex = 3;
				this->CancelBtn->Text = L"Cancel";
				this->CancelBtn->UseVisualStyleBackColor = true;
				// 
				// WCSAutoReportingTimer
				// 
				this->WCSAutoReportingTimer->Tick += gcnew System::EventHandler(this, &WCSAutoSolverReportingForm::WCSAutoReportingTimer_Tick);
				// 
				// WCSAutoSolverReportingForm
				// 
				this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
				this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
				this->ClientSize = System::Drawing::Size(503, 496);
				this->Controls->Add(this->MsgTxt);
				this->Controls->Add(this->CancelBtn);
				this->Name = L"WCSAutoSolverReportingForm";
				this->Text = L"Auto Solver Report...";
				this->TopMost = true;
				this->ResumeLayout(false);
				this->PerformLayout();

			}
			#pragma endregion
		public: System::Void WCSAutoReportingTimer_Tick(System::Object^ sender, System::EventArgs^ e);
		};

		~WCS_AutoSolver();

		/// <summary>Initializes the WCS_AutoSolver class including performing source extraction on a given FITS image.</summary>
		/// <param name="WCS_type">The WCS tranformation type. Solution only uses TAN at this time.</param>
		/// <param name="Number_of_Points">The number of points N to use to compare image coordinates to catalogue coordinates. Suggest N equals 25 for good correspondence, N equals 50 for poor, N equals 100 for very poor.</param>
		/// <param name="Fits_Img">The JPFITS::FITSImage containing the primary image data.</param>
		/// <param name="Image_ROI">The region of interest of the FITS image to search for point sources, of identical size to the FITS image. Pass nullptr or all true for entire image.</param>
		/// <param name="Image_Saturation">The saturation level of the source image for mapping saturated sources. Pass zero if no saturated sources exist.</param>
		/// <param name="auto_background">Automatically determine local background for each centroiding kernel.</param>
		/// <param name="PSE_kernel_radius">The radius of the point-source-extraction kernel, in pixels. PSEkernel_radius greater than or equal to 1.</param>
		/// <param name="PSE_separation_radius">The minimum separation of point sources, in pixels. PSESeparation_radius greater than or equal to PSEkernel_radius.</param>
		/// <param name="Fits_Catalogue_BinTable_File">The full path file name of the FITS binary table containing the catalogue data.</param>
		/// <param name="Catalogue_Extension_Name">The extension name of the FITS binary table which contains the catalogue data. If empty string is passed then the first binary table extension is assumed.</param>
		/// <param name="Catalogue_CVAL1_Name">The name of the entry inside the binary table which lists the CVAL1 (i.e. right ascension) coordinates.</param>
		/// <param name="Catalogue_CVAL2_Name">The name of the entry inside the binary table which lists the CVAL2 (i.e. declination) coordinates.</param>
		/// <param name="Catalogue_Magnitude_Name">The name of the entry inside the binary table which lists the source magnitudes.</param>
		/// <param name="Refine">Option to refine the solution further with additional points after the initial solution is found.</param>
		WCS_AutoSolver(String^ WCS_type, int Number_of_Points, JPFITS::FITSImage^ Fits_Img, array<bool, 2>^ Image_ROI, double Image_Saturation, bool auto_background, int PSE_kernel_radius, int PSE_separation_radius, String^ Fits_Catalogue_BinTable_File, String^ Catalogue_Extension_Name, String^ Catalogue_CVAL1_Name, String^ Catalogue_CVAL2_Name, String^ Catalogue_Magnitude_Name, bool Refine);

		/// <summary>Initializes the WCS_AutoSolver class for a given pair of pixel source and catalogue coordinates.</summary>
		/// <param name="WCS_type">The WCS tranformation type. Solution only uses TAN at this time.</param>
		/// <param name="pixels">The source pixel positions in computer graphics coordinates, i.e., origin top left of screen.</param>
		/// <param name="zero_based_pixels">If the source pixel positions are zero-based.</param>
		/// <param name="pixels_tolerance_radius">The tolerance of the source positions, identical to usage as the PSE_kernel_radius in the other contructor. Typically 2 (pixels).</param>
		/// <param name="image_width">The 1-based width of the source image from where the source pixels points originate.</param>
		/// <param name="image_height">The 1-based height of the source image from where the source pixels points originate.</param>
		/// <param name="wcspoints">The catalogue values corresponding to the region in the image of the source pixel positions.</param>
		WCS_AutoSolver(String^ WCS_type, array<JPMath::PointD^>^ pixels, bool zero_based_pixels, int pixels_tolerance_radius, int image_width, int image_height, array<JPMath::PointD^>^ wcspoints);

		/// <summary>Executes the auto-solver algorithm.</summary>
		/// <param name="scale_init">The initial scale guess, in arcseconds per pixel.</param>
		/// <param name="scale_lb">The lower bound of the scale range, in arcseconds per pixel.</param>
		/// <param name="scale_ub">The upper bound of the scale range, in arcseconds per pixel.</param>
		/// <param name="rotation_init">The initial field rotation guess, in degrees.</param>
		/// <param name="rotation_lb">The lower bound of the field rotation range, in degrees, greater than or equal to -180</param>
		/// <param name="rotation_ub">The upper bound of the field rotation range, in degrees, less than or equal to 180</param>
		/// <param name="vertex_tolerance">The tolerance of the vertex angles when comparing triangles, in degrees. Suggest 0.25.</param>
		/// <param name="N_matches_stop">Stop and solve solution when N matches are found between image and catalogue coordinates. N_matches_stop greater than or equal to 3. Suggest 6. Solution likely requires confirmation at 3 or 4.</param>
		/// <param name="Percentage_matches_stop">Stop and solve solution when Percentage matches are found between image and catalogue coordinates. Suggest 25.</param>
		/// <param name="condition_arrays">Optionally condition the triangle arrays. Suggest true.</param>
		/// <param name="show_report_form">Optionally shows a cancellable Form which displays the solution progress.</param>
		void SolveAsync(double scale_init, double scale_lb, double scale_ub, double rotation_init, double rotation_lb, double rotation_ub, double vertex_tolerance, int N_matches_stop, int Percentage_matches_stop, bool condition_arrays, bool show_report_form);

		/// <summary>Conditions the traingle array so that all threads begin with the brightest triangles.</summary>
		/// <param name="triarray">An array of triangles.</param>
		/// <param name="Nthreads">The number of threads to condition the array for.</param>
		/// <param name="ascending">Brightness is ascending values (i.e. magnitudes) = true, otherwise brightness is descending values (i.e. counts) = false.</param>
		static array<JPMath::Triangle^>^ ConditionTriangleArrayBrightnessThreads(array<JPMath::Triangle^>^ triarray, int Nthreads, bool ascending);

		/// <summary>Queries the Gaia catalogue for entries within a specified region</summary>
		/// <param name="catalogue">A string for the catalogue to query. Options are (case insensitive): "Gaia"</param>
		/// <param name="ra_deg">A string of the right ascension in degrees.</param>
		/// <param name="dec_deg">A string of the declination in degrees.</param>
		/// <param name="result_savepathfilename">The filename to save the catalogue. If saving is not required, pass an empty string.</param>
		/// <param name="radius">A string of the region radius in arcminutes.</param>
		/// <param name="square">Pass 1 if the region is square, 0 for circle.</param>
		static int AstroQuery(String^ catalogue, String^ ra_deg, String^ dec_deg, String^ &result_savepathfilename, String^ radius, String^ square);

		/// <summary>Returns the World Coordinate Solution</summary>
		property WorldCoordinateSolution^ WCS_Solution
		{
			WorldCoordinateSolution^ get() { return this->WCS; }
		}

		/// <summary>Returns the most recent Point Source Extraction</summary>
		property SourceExtractor^ PSE_Extraction
		{
			SourceExtractor^ get() { return this->PSE; }
		}

		/// <summary>Gets or Sets the Status Log</summary>
		property String^ Status_Log
		{
			String^ get() { return this->STATUS_LOG; }
			void set(String^ status) { STATUS_LOG += "\r" + status; }
		}

		/// <summary>Clears the Status Log</summary>
		void Status_Log_Clear()
		{
			STATUS_LOG = "";
		}

		/// <summary>Gets or Sets the Cancel State of the Solver</summary>
		property bool Cancelled
		{
			bool get() { return CANCELLED; }
			void set(bool cancel) { CANCELLED = cancel; }
		}

		/// <summary>Gets the Solving State of the Solver</summary>
		property bool Solving
		{
			bool get() { return SOLVING; }
		}

		/// <summary>Gets the Solving State of the Solver</summary>
		property bool Solved
		{
			bool get() { return SOLVED; }
		}

		/// <summary>Gets Progress of the Solver</summary>
		property int Progress
		{
			int get() { return PROGRESS; }
		}

		/// <summary>Gets or Sets the Solver to run parallelized (default is true)</summary>
		property bool Solver_Parallelized
		{
			bool get() { return DO_PARALLEL; }
			void set(bool parallel) { DO_PARALLEL = parallel; }
		}

	private:
		bool SOLVED = false, CANCELLED = true, DO_PSE = false, DO_PARALLEL = true, SOLVING = false, CONDITION_ARRAYS = false, SHOW_REPORT_FORM = false;
		String^ STATUS_LOG = "";
		array<bool, 2>^ IMAGE_ROI;
		String^ CAT_FILENAME;
		String^ CAT_EXTNAME;
		String^ CAT_CVAL1NAME;
		String^ CAT_CVAL2NAME;
		String^ CAT_MAGNAME;
		FITSImage^ FITS_IMG;
		FITSImage^ FITS_CAT;
		array<double>^ CAT_CVAL1s;
		array<double>^ CAT_CVAL2s;
		array<double>^ CAT_MAGs;
		array<JPMath::PointD^>^ PIX_PTS;
		array<JPMath::PointD^>^ CAT_PTS;
		int PROGRESS, IMAGE_WIDTH, IMAGE_HEIGHT, N_POINTS, PSE_KERNEL_RADIUS, PSE_SEP_RADIUS, N_MATCHES_STOP, PERC_MATCHES_STOP;
		double SCALE_INIT, SCALE_LB, SCALE_UB, ROTATION_INIT, ROTATION_LB, ROTATION_UB, WCS_VERTEX_TOL, PIX_SAT;
		unsigned __int64 NCOMPARES;
		WorldCoordinateSolution^ WCS;
		SourceExtractor^ PSE;
		DateTime DATE;
		String^ WCS_TYPE;
		bool ZERO_BASED_PIX;
		bool AUTO_BACKGROUND;
		bool REFINE = false;
		JPWaitBar::WaitBar^ WAITBAR;
		BackgroundWorker^ BGWRKR;
		void BGWRKR_DoWork(System::Object^  sender, System::ComponentModel::DoWorkEventArgs^  e);
		void BGWRKR_ProgressChanged(System::Object^  sender, System::ComponentModel::ProgressChangedEventArgs^  e);
		void BGWRKR_RunWorkerCompleted(System::Object^  sender, System::ComponentModel::RunWorkerCompletedEventArgs^  e);
		static void MAKEASTROQUERYSCRIPT(String^ script_filename, String^ catalogue);

		WCSAutoSolverReportingForm^ WCSARF;
	};

	/// <summary>JPBitMap class provides functionality for converting 2D arrays to color-mapped Bitmaps.</summary>
	public ref class JPBitMap
	{
		public:

		static Bitmap^ ArrayToBmp(array<double, 2>^ image, int scaling, int colour, bool invert, array<double, 1>^ DImCLim, int WinWidth, int WinHeight, bool invertYaxis);

		//static Bitmap^ ArrayTo16bppGSBmp(array<double, 2>^ image, int scaling, int colour, bool invert, array<double, 1>^ DImCLim, int WinWidth, int WinHeight, bool invertYaxis);
		
		static Bitmap^ RGBBitMap(array<double, 2>^ R, array<double, 2>^ G, array<double, 2>^ B);

		static array<double, 2>^ Bin(array<double, 2>^ data, int Nx, int Ny);

		static double JetR(double val);
		static double JetG(double val);
		static double JetB(double val);
		static double WinterR(double val);
		static double WinterG(double val);
		static double WinterB(double val);
		static double LinesR(double val);
		static double LinesG(double val);
		static double LinesB(double val);
	};

}

