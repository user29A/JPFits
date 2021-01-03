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


#include "StdAfx.h"
#include "FitsFinder.h"

void JPFITS::FitsFinder::CancelBtn_Click(System::Object^  sender, System::EventArgs^  e)
{
	FitsFinder::DialogResult = ::DialogResult::Cancel;
	FitsFinder::Close();
}

void JPFITS::FitsFinder::DirectoryTxt_Click(System::Object^  sender, System::EventArgs^  e)
{
	::FolderBrowserDialog^ fb = gcnew ::FolderBrowserDialog();
	fb->SelectedPath = (String^)GetReg("CCDLAB", "OpenFilesPath");
	fb->ShowDialog();
	String^ dir = fb->SelectedPath;
	DirectoryTxt->Text = dir;
	SetReg("CCDLAB", "OpenFilesPath",dir);
}

void JPFITS::FitsFinder::FitsFinder_Load(System::Object^  sender, System::EventArgs^  e)
{
	DirectoryTxt->Text = (String^)GetReg("CCDLAB", "OpenFilesPath");
	String^ tmplt = (String^)GetReg("CCDLAB", "FindFilesTemplate");
	FileTemplateTxt->Text = tmplt->Substring(0,tmplt->LastIndexOf("."));
	int NumKeyValPairs = System::Convert::ToInt32(GetReg("CCDLAB", "FindFilesNumKeyValPairs"));
	if (NumKeyValPairs > 0)
	{
		String^ key;
		String^ keyval;
		array<System::Windows::Forms::TextBox^>^ k = gcnew array<System::Windows::Forms::TextBox^>{Key1,Key2,Key3,Key4};
		array<System::Windows::Forms::RichTextBox^>^ kv = gcnew array<System::Windows::Forms::RichTextBox^>{Key1Value,Key2Value,Key3Value,Key4Value};
		for (int i = 0; i < NumKeyValPairs; i++)
		{
			key = (String^)GetReg("CCDLAB", String::Concat("FindFilesKey",i));
			keyval = (String^)GetReg("CCDLAB", String::Concat("FindFilesKeyVal",i));
			k[i]->Text = key;
			kv[i]->Text = keyval;
		}
	}
	ExtensionDrop->SelectedIndex = ::Convert::ToInt32(GetReg("CCDLAB", "FindFilesExtIndex"));
	FitsFinder::Tag = ::DialogResult::None;

	SubFoldersChck->Checked = ::Convert::ToBoolean(GetReg("CCDLAB", "SubFoldersChck"));

	CustomExtensionChck->Checked = ::Convert::ToBoolean(GetReg("CCDLAB", "CustomExtChck"));
	CustomExtensionTxtBox->Text = (String^)GetReg("CCDLAB", "CustomExtTxt");

	FOUNDFILES = gcnew array<String^>(0);
}

void JPFITS::FitsFinder::FindBtn_Click(System::Object^  sender, System::EventArgs^  e)
{
	String^ dir = DirectoryTxt->Text;
	if (!::Directory::Exists(dir))
	{
		::MessageBox::Show("Directory doesn't exist...","Error");
		return;
	}

	bool subdirs = SubFoldersChck->Checked;
	SetReg("CCDLAB", "SubFoldersChck",subdirs);

	//need to get search params and write them to reg for later FindFiles()
	String^ extension;
	if (CustomExtensionChck->Checked)
		extension = CustomExtensionTxtBox->Text;
	else
		extension = ExtensionDrop->Items[ExtensionDrop->SelectedIndex]->ToString();

	String^ filetemplate = String::Concat(FileTemplateTxt->Text,extension);//file template for cursory directory search, which we'll start with
	int count = 0;
	array<System::Windows::Forms::TextBox^>^ k = gcnew array<System::Windows::Forms::TextBox^>{Key1,Key2,Key3,Key4};
	array<System::Windows::Forms::RichTextBox^>^ kv = gcnew array<System::Windows::Forms::RichTextBox^>{Key1Value,Key2Value,Key3Value,Key4Value};
	for (int i = 0; i < k->Length; i++)
	{
		if (k[i]->Text->Length == 0)
			continue;
		SetReg("CCDLAB", String::Concat("FindFilesKey",count),k[i]->Text);
		SetReg("CCDLAB", String::Concat("FindFilesKeyVal",count),kv[i]->Text);
		count++;
	}
	SetReg("CCDLAB", "FindFilesNumKeyValPairs",count.ToString());
	SetReg("CCDLAB", "FindFilesExtIndex",ExtensionDrop->SelectedIndex);
	SetReg("CCDLAB", "FindFilesTemplate",filetemplate);
	SetReg("CCDLAB", "CustomExtChck", CustomExtensionChck->Checked);
	SetReg("CCDLAB", "CustomExtTxt", CustomExtensionTxtBox->Text);

	array<String^>^ fullfilesinit;
	if (!subdirs)
		fullfilesinit = ::IO::Directory::GetFiles(dir,filetemplate,System::IO::SearchOption::TopDirectoryOnly);//cursory search
	else
		fullfilesinit = ::IO::Directory::GetFiles(dir,filetemplate,System::IO::SearchOption::AllDirectories);//cursory search

	if (count > 0)//then we're doing more than just a cursory file template search
	{
		this->WAITBAR = gcnew JPWaitBar::WaitBar();
		this->WAITBAR->ProgressBar->Maximum = fullfilesinit->Length;
		this->WAITBAR->Text = "Searching files...";
		FitsFinderWrkr->RunWorkerAsync(fullfilesinit);
		this->WAITBAR->ShowDialog();
	}
	else
	{
		FOUNDFILES = fullfilesinit;
		FitsFinder::Tag = ::DialogResult::OK;
		FitsFinder::DialogResult = ::DialogResult::OK;
	}
}

void JPFITS::FitsFinder::FitsFinderWrkr_DoWork(System::Object^  sender, System::ComponentModel::DoWorkEventArgs^  e)
{
	int numparams = System::Convert::ToInt32(GetReg("CCDLAB", "FindFilesNumKeyValPairs"));
	array<String^>^ fullfilesinit = (array<String^>^)e->Argument;
	ArrayList^ filelist = gcnew ArrayList();
	array<String^>^ KeyParams = gcnew array<String^>(numparams);
	array<String^>^ KeyValParams = gcnew array<String^>(numparams);
	int match = 0;
	for (int i = 0; i < numparams; i++)//get the key/keyvalue pairs
	{
		KeyParams[i] = (String^)GetReg("CCDLAB", String::Concat("FindFilesKey",i));
		KeyValParams[i] = (String^)GetReg("CCDLAB", String::Concat("FindFilesKeyVal",i));
	}
	//done filling the param pairs...now need to do the work

	for (int ii = 0; ii < fullfilesinit->Length; ii++)
	{
		if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
		{
			FitsFinder::DialogResult = ::DialogResult::Cancel;
			FitsFinder::Tag = ::DialogResult::Cancel;
			return;
		}
		FitsFinderWrkr->ReportProgress(ii+1,filelist->Count);

		FITSImage^ f1 = gcnew FITSImage(fullfilesinit[ii], nullptr, true, false, false, false);
		array<String^>^ HeadKeys = f1->HeaderKeys;			//get just header from file in JPFITSImage
		array<String^>^ HeadKeyVals = f1->HeaderKeyValues;	//get just header key values from file in JPFITSImage
		match = 0;
		for (int j = 0; j < HeadKeys->Length; j++)
		{
			String^ key = HeadKeys[j];
			for (int k = 0; k < numparams; k++)
				if (KeyParams[k]->Compare(key,KeyParams[k],false) == 0 && KeyValParams[k]->Compare(KeyValParams[k],HeadKeyVals[j],true) == 0)
					match++;
		}
		if (match == numparams)
			filelist->Add(fullfilesinit[ii]);
	}

	array<String^>^ matchedfiles = gcnew array<String^>(filelist->Count);
	for (int h = 0; h < filelist->Count; h++)
		matchedfiles[h] = (String^)filelist[h];

	e->Result = matchedfiles;

	FitsFinder::DialogResult = ::DialogResult::OK;
}

void JPFITS::FitsFinder::FitsFinderWrkr_ProgressChanged(System::Object^  sender, System::ComponentModel::ProgressChangedEventArgs^  e)
{
	this->WAITBAR->Text = "Searching files. Found " + ::Convert::ToInt32(e->UserState).ToString() + " matches...";
	this->WAITBAR->ProgressBar->Value = e->ProgressPercentage;
	this->WAITBAR->TextMsg->Text = "Examining File: " + e->ProgressPercentage.ToString() + " of " + this->WAITBAR->ProgressBar->Maximum.ToString();
	this->WAITBAR->Refresh();
}

void JPFITS::FitsFinder::FitsFinderWrkr_RunWorkerCompleted(System::Object^  sender, System::ComponentModel::RunWorkerCompletedEventArgs^  e)
{
	if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
	{
		FOUNDFILES = gcnew array<String^>(0);
		FitsFinder::DialogResult = ::DialogResult::Cancel;
		FitsFinder::Tag = ::DialogResult::Cancel;
		FitsFinder::Close();
		this->WAITBAR->Close();
	}
	else
	{
		FOUNDFILES = (array<String^>^)e->Result;
		FitsFinder::Tag = ::DialogResult::OK;
		FitsFinder::Close();
		this->WAITBAR->DialogResult = ::DialogResult::OK;
		this->WAITBAR->Close();
	}
}

void JPFITS::FitsFinder::CustomExtensionTxtBox_TextChanged(System::Object^  sender, System::EventArgs^  e)
{
	/*if (CustomExtensionTxtBox->Text->Substring(0, 1) != ".")
		CustomExtensionTxtBox->Text = "." + CustomExtensionTxtBox->Text;

	CustomExtensionTxtBox->Cur*/
}

