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


#include "StdAfx.h"
#include "FitsFound.h"

void JPFITS::FitsFound::SetFoundFiles(array<String^>^ FoundFiles)
{
	FOUNDFILES = gcnew array<String^>(FoundFiles->Length);
	FOUNDFILES = FoundFiles;
	this->Text = "Found " + FOUNDFILES->Length.ToString() + " files...";
	NumFilesTxt->Text = "Please Select File(s)...";

	FileListTxt->BeginUpdate();
	FileListTxt->Items->AddRange(FOUNDFILES);
	FileListTxt->EndUpdate();
}

void JPFITS::FitsFound::SetFoundFiles(String^ FullFileName)
{
	FileStream^ fs2 = gcnew FileStream(FullFileName, System::IO::FileMode::Open, FileAccess::Read);
	StreamReader^ sr2 = gcnew StreamReader(fs2);
	int numlines = System::Convert::ToInt32(sr2->ReadLine());

	FOUNDFILES = gcnew array<String^>(numlines);
	for (int i = 0; i < numlines; i++)
		FOUNDFILES[i] = sr2->ReadLine();

	sr2->Close();
	fs2->Close();
}

void JPFITS::FitsFound::FileListTxt_SelectedIndexChanged(System::Object^  sender, System::EventArgs^  e)
{
	NumSelectTxt->Text = "(" + FileListTxt->SelectedItems->Count + " selected)";
}

void JPFITS::FitsFound::SelectAllBtn_Click(System::Object^  sender, System::EventArgs^  e) 
{
	FileListTxt->BeginUpdate();
	for (int i = 0; i < FileListTxt->Items->Count; i++)
		FileListTxt->SetSelected(i,true);
	FileListTxt->EndUpdate();
}

void JPFITS::FitsFound::ClearAllBtn_Click(System::Object^  sender, System::EventArgs^  e) 
{
	FileListTxt->BeginUpdate();
	for (int i = 0; i < FileListTxt->Items->Count; i++)
		FileListTxt->SetSelected(i,false);
	FileListTxt->EndUpdate();
}

void JPFITS::FitsFound::SaveListBtn_Click(System::Object^  sender, System::EventArgs^  e)
{
	int Ninds = FileListTxt->SelectedIndices->Count;
	if (Ninds == 0)//no files selected but asked to save files
	{
		::MessageBox::Show("No Files Selected!...","Error");
		return;
	}

	array<String^>^ selectfiles = gcnew array<String^>(Ninds);
	for (int j = 0; j < Ninds; j++)
		selectfiles[j] = (String^)FileListTxt->Items[FileListTxt->SelectedIndices[j]];

	String^ dir = (String^)GetReg("CCDLAB", "OpenFilesPath");
	SaveFileDialog^ dlg = gcnew SaveFileDialog();
	dlg->InitialDirectory = dir;
	dlg->Filter = "CCDLAB File List (*.CFL)|*.CFL";
	::DialogResult res = dlg->ShowDialog();

	if (res == ::DialogResult::OK)
	{
		String^ file = dlg->FileName;
		FileStream^ fs = gcnew FileStream(file,System::IO::FileMode::Create,FileAccess::Write);
		StreamWriter^ sw = gcnew StreamWriter(fs);
		sw->WriteLine(selectfiles->Length);
		for (int u = 0; u < selectfiles->Length; u++)
			 sw->WriteLine(selectfiles[u]);
		sw->Flush();
		fs->Flush();
		sw->Close();
		fs->Close();
		FitsFound::DialogResult = ::DialogResult::Yes;
		FitsFound::Close();

		SetReg("CCDLAB", "FoundFileList",file);
	}
}

void JPFITS::FitsFound::CancelBtn_Click(System::Object^  sender, System::EventArgs^  e) 
{
	FitsFound::DialogResult = ::DialogResult::Cancel;
}

void JPFITS::FitsFound::MoveListBtn_Click(System::Object^  sender, System::EventArgs^  e)
{
	int Ninds = FileListTxt->SelectedIndices->Count;
	if (Ninds == 0)//no files selected but asked to copy files
	{
		::MessageBox::Show("No Files Selected!...","Error");
		return;
	}

	array<String^>^ selectfiles = gcnew array<String^>(Ninds);
	for (int j = 0; j < Ninds; j++)
		selectfiles[j] = (String^)FileListTxt->Items[FileListTxt->SelectedIndices[j]];

	::FolderBrowserDialog^ fdlg = gcnew ::FolderBrowserDialog();
	fdlg->SelectedPath = selectfiles[0]->Substring(0,selectfiles[0]->LastIndexOf("\\")+1);

	if (fdlg->ShowDialog() == ::DialogResult::Cancel)
		return;

	array<Object^>^ arg = gcnew array<Object^>(3){selectfiles,fdlg->SelectedPath,"move"};

	WAITBAR = gcnew JPWaitBar::WaitBar();
	WAITBAR->Text = "Moving Files...";
	WAITBAR->ProgressBar->Maximum = selectfiles->Length;
	FileCopyBGWrkr->RunWorkerAsync(arg);
	WAITBAR->ShowDialog();
}

void JPFITS::FitsFound::CopyListBtn_Click(System::Object^  sender, System::EventArgs^  e)
{
	int Ninds = FileListTxt->SelectedIndices->Count;
	if (Ninds == 0)//no files selected but asked to copy files
	{
		::MessageBox::Show("No Files Selected!...","Error");
		return;
	}

	array<String^>^ selectfiles = gcnew array<String^>(Ninds);
	for (int j = 0; j < Ninds; j++)
		selectfiles[j] = (String^)FileListTxt->Items[FileListTxt->SelectedIndices[j]];

	::FolderBrowserDialog^ fdlg = gcnew ::FolderBrowserDialog();
	fdlg->SelectedPath = selectfiles[0]->Substring(0,selectfiles[0]->LastIndexOf("\\")+1);

	if (fdlg->ShowDialog() == ::DialogResult::Cancel)
		return;

	array<Object^>^ arg = gcnew array<Object^>(3){selectfiles,fdlg->SelectedPath,"copy"};

	WAITBAR = gcnew JPWaitBar::WaitBar();
	WAITBAR->Text = "Copying Files...";
	WAITBAR->ProgressBar->Maximum = selectfiles->Length;
	FileCopyBGWrkr->RunWorkerAsync(arg);
	WAITBAR->ShowDialog();
}

void JPFITS::FitsFound::FileCopyBGWrkr_DoWork(System::Object^  sender, System::ComponentModel::DoWorkEventArgs^  e)
{
	array<Object^>^ arg = (array<Object^>^)e->Argument;
	array<String^>^ selectfiles = (array<String^>^)arg[0];
	String^ selectedpath = (String^)arg[1];
	String^ style = (String^)arg[2];
	bool move = false;
	if (style == "move")
		move = true;

	int Ninds = selectfiles->Length;
	String^ newfile;

	for (int j = 0; j < Ninds; j++)
	{
		if (WAITBAR->DialogResult == ::DialogResult::Cancel)
			return;
		FileCopyBGWrkr->ReportProgress(j+1,move);

		newfile = selectedpath + "\\" + selectfiles[j]->Substring(selectfiles[j]->LastIndexOf("\\"));

		while(::File::Exists(newfile))//then need to add some appendage
		{
			int ind = newfile->LastIndexOf(".");
			if (newfile->Substring(ind-1,1) == ")")
			{
				int num = ::Convert::ToInt32(newfile->Substring(newfile->LastIndexOf("(") + 1,newfile->LastIndexOf(")") - 1 - newfile->LastIndexOf("(")));
				newfile = newfile->Replace("(" + num.ToString() + ").","(" + (num+1).ToString() + ").");
			}
			else
			{
				newfile = newfile->Insert(ind," (1)");
				if (::File::Exists(newfile))
				{
					int num = ::Convert::ToInt32(newfile->Substring(newfile->LastIndexOf("(") + 1,newfile->LastIndexOf(")") - 1 - newfile->LastIndexOf("(")));
					newfile = newfile->Replace("(" + num.ToString() + ").","(" + (num+1).ToString() + ").");
				}
			}
		}
		if (move)
			::File::Move(selectfiles[j], newfile);
		else
			::File::Copy(selectfiles[j], newfile);
	}
}

void JPFITS::FitsFound::FileCopyBGWrkr_ProgressChanged(System::Object^  sender, System::ComponentModel::ProgressChangedEventArgs^  e)
{
	bool move = ::Convert::ToBoolean(e->UserState);
	if (move)
		WAITBAR->TextMsg->Text = "Moving file " + e->ProgressPercentage.ToString() + " of " + WAITBAR->ProgressBar->Maximum.ToString();
	else
		WAITBAR->TextMsg->Text = "Copying file " + e->ProgressPercentage.ToString() + " of " + WAITBAR->ProgressBar->Maximum.ToString();
	WAITBAR->ProgressBar->Value = e->ProgressPercentage;
	WAITBAR->Refresh();
}

void JPFITS::FitsFound::FileCopyBGWrkr_RunWorkerCompleted(System::Object^  sender, System::ComponentModel::RunWorkerCompletedEventArgs^  e)
{
	WAITBAR->Close();
	if (WAITBAR->DialogResult != ::DialogResult::Cancel)
	{
		if (::Convert::ToBoolean(e->UserState))
			::MessageBox::Show("Finished moving files...");
		else
			::MessageBox::Show("Finished copying files...");
	}
	CancelBtn->PerformClick();
}

void JPFITS::FitsFound::FileListTxt_MouseClick(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
}

void JPFITS::FitsFound::FileListTxt_MouseUp(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	int Ninds = FileListTxt->SelectedIndices->Count;
	if (Ninds == 1)
	{
		FoundListContextMenu->Enabled = true;
		OpenFolderContextItem->Enabled = true;
	}
	else
	{
		FoundListContextMenu->Enabled = false;
		OpenFolderContextItem->Enabled = false;
	}
}

void JPFITS::FitsFound::OpenFolderContextItem_Click(System::Object^  sender, System::EventArgs^  e)
{
	String^ selectfile = (String^)FileListTxt->Items[FileListTxt->SelectedIndices[0]];
	selectfile = selectfile->Substring(0,selectfile->LastIndexOf("\\")+1);
	::System::Diagnostics::Process::Start("Explorer.exe",selectfile);
}

void JPFITS::FitsFound::FitsFound_FormClosing(System::Object^  sender, System::Windows::Forms::FormClosingEventArgs^  e)
{
	SetReg("JPFITS","FitsFoundPOSX",this->Location.X);
	SetReg("JPFITS","FitsFoundPOSY",this->Location.Y);
	SetReg("JPFITS","FitsFoundWIDTH",this->Size.Width);
	SetReg("JPFITS","FitsFoundHEIGHT",this->Size.Height);
}

void JPFITS::FitsFound::FitsFound_Load(System::Object^  sender, System::EventArgs^  e)
{
	this->Left = (int)GetReg("JPFITS","FitsFoundPOSX");
	this->Top = (int)GetReg("JPFITS","FitsFoundPOSY");
	this->Width = (int)GetReg("JPFITS","FitsFoundWIDTH");
	this->Height = (int)GetReg("JPFITS","FitsFoundHEIGHT");
}

