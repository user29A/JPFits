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
using namespace System::ComponentModel;
using namespace System::Collections;
using namespace System::Windows::Forms;
using namespace System::Data;
using namespace System::Drawing;
#include "FitsFinder.h"


namespace JPFITS {

	/// <summary>A class to interact with found FITS files from the FitsFinder.</summary>
	public ref class FitsFound : public System::Windows::Forms::Form
	{
	public:
		FitsFound(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~FitsFound()
		{
			if (components)
			{
				delete components;
			}
		}
	public: System::Windows::Forms::Label^  NumSelectTxt;
	protected: 
	private: System::Windows::Forms::Button^  ClearAllBtn;
	public: 
	private: System::Windows::Forms::Button^  SaveListBtn;
	private: System::Windows::Forms::Button^  SelectAllBtn;
	private: System::Windows::Forms::Button^  CancelBtn;


	public: System::Windows::Forms::Label^  NumFilesTxt;
	private: 
	public: System::Windows::Forms::ListBox^  FileListTxt;
	public: System::Windows::Forms::Button^  AddImageSetBtn;
	public: System::Windows::Forms::Button^  LoadImageSetBtn;
	private: System::Windows::Forms::Button^  MoveListBtn;
	private: System::Windows::Forms::Button^  CopyListBtn;

	private: System::Windows::Forms::ContextMenuStrip^  FoundListContextMenu;
	private: System::Windows::Forms::ToolStripMenuItem^  OpenFolderContextItem;
	private: System::ComponentModel::BackgroundWorker^  FileCopyBGWrkr;
	private: System::ComponentModel::IContainer^  components;
	public: 


	public: 


	public: 


	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>


#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(FitsFound::typeid));
			this->NumSelectTxt = (gcnew System::Windows::Forms::Label());
			this->ClearAllBtn = (gcnew System::Windows::Forms::Button());
			this->SaveListBtn = (gcnew System::Windows::Forms::Button());
			this->SelectAllBtn = (gcnew System::Windows::Forms::Button());
			this->CancelBtn = (gcnew System::Windows::Forms::Button());
			this->NumFilesTxt = (gcnew System::Windows::Forms::Label());
			this->FileListTxt = (gcnew System::Windows::Forms::ListBox());
			this->FoundListContextMenu = (gcnew System::Windows::Forms::ContextMenuStrip(this->components));
			this->OpenFolderContextItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->AddImageSetBtn = (gcnew System::Windows::Forms::Button());
			this->LoadImageSetBtn = (gcnew System::Windows::Forms::Button());
			this->MoveListBtn = (gcnew System::Windows::Forms::Button());
			this->CopyListBtn = (gcnew System::Windows::Forms::Button());
			this->FileCopyBGWrkr = (gcnew System::ComponentModel::BackgroundWorker());
			this->FoundListContextMenu->SuspendLayout();
			this->SuspendLayout();
			// 
			// NumSelectTxt
			// 
			this->NumSelectTxt->AutoSize = true;
			this->NumSelectTxt->Location = System::Drawing::Point(26, 25);
			this->NumSelectTxt->Name = L"NumSelectTxt";
			this->NumSelectTxt->Size = System::Drawing::Size(62, 13);
			this->NumSelectTxt->TabIndex = 17;
			this->NumSelectTxt->Text = L"(0 selected)";
			// 
			// ClearAllBtn
			// 
			this->ClearAllBtn->Anchor = System::Windows::Forms::AnchorStyles::None;
			this->ClearAllBtn->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Center;
			this->ClearAllBtn->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->ClearAllBtn->Location = System::Drawing::Point(259, 539);
			this->ClearAllBtn->Name = L"ClearAllBtn";
			this->ClearAllBtn->Size = System::Drawing::Size(75, 23);
			this->ClearAllBtn->TabIndex = 16;
			this->ClearAllBtn->Text = L"Cl&ear All";
			this->ClearAllBtn->UseVisualStyleBackColor = true;
			this->ClearAllBtn->Click += gcnew System::EventHandler(this, &FitsFound::ClearAllBtn_Click);
			// 
			// SaveListBtn
			// 
			this->SaveListBtn->Anchor = System::Windows::Forms::AnchorStyles::None;
			this->SaveListBtn->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Center;
			this->SaveListBtn->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->SaveListBtn->ForeColor = System::Drawing::SystemColors::ControlText;
			this->SaveListBtn->Location = System::Drawing::Point(171, 646);
			this->SaveListBtn->Name = L"SaveListBtn";
			this->SaveListBtn->Size = System::Drawing::Size(75, 23);
			this->SaveListBtn->TabIndex = 15;
			this->SaveListBtn->Text = L"Sa&ve List";
			this->SaveListBtn->UseVisualStyleBackColor = true;
			this->SaveListBtn->Click += gcnew System::EventHandler(this, &FitsFound::SaveListBtn_Click);
			// 
			// SelectAllBtn
			// 
			this->SelectAllBtn->Anchor = System::Windows::Forms::AnchorStyles::None;
			this->SelectAllBtn->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->SelectAllBtn->Location = System::Drawing::Point(171, 539);
			this->SelectAllBtn->Name = L"SelectAllBtn";
			this->SelectAllBtn->Size = System::Drawing::Size(75, 23);
			this->SelectAllBtn->TabIndex = 14;
			this->SelectAllBtn->Text = L"&Select All";
			this->SelectAllBtn->UseVisualStyleBackColor = true;
			this->SelectAllBtn->Click += gcnew System::EventHandler(this, &FitsFound::SelectAllBtn_Click);
			// 
			// CancelBtn
			// 
			this->CancelBtn->Anchor = System::Windows::Forms::AnchorStyles::None;
			this->CancelBtn->DialogResult = System::Windows::Forms::DialogResult::Cancel;
			this->CancelBtn->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->CancelBtn->Location = System::Drawing::Point(259, 646);
			this->CancelBtn->Name = L"CancelBtn";
			this->CancelBtn->Size = System::Drawing::Size(75, 23);
			this->CancelBtn->TabIndex = 13;
			this->CancelBtn->Text = L"&Cancel";
			this->CancelBtn->UseVisualStyleBackColor = true;
			this->CancelBtn->Click += gcnew System::EventHandler(this, &FitsFound::CancelBtn_Click);
			// 
			// NumFilesTxt
			// 
			this->NumFilesTxt->AutoSize = true;
			this->NumFilesTxt->Location = System::Drawing::Point(26, 9);
			this->NumFilesTxt->Name = L"NumFilesTxt";
			this->NumFilesTxt->Size = System::Drawing::Size(71, 13);
			this->NumFilesTxt->TabIndex = 10;
			this->NumFilesTxt->Text = L"Found # Files";
			// 
			// FileListTxt
			// 
			this->FileListTxt->Anchor = static_cast<System::Windows::Forms::AnchorStyles>(((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Left) 
				| System::Windows::Forms::AnchorStyles::Right));
			this->FileListTxt->ContextMenuStrip = this->FoundListContextMenu;
			this->FileListTxt->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->FileListTxt->FormattingEnabled = true;
			this->FileListTxt->HorizontalScrollbar = true;
			this->FileListTxt->Location = System::Drawing::Point(29, 44);
			this->FileListTxt->Name = L"FileListTxt";
			this->FileListTxt->SelectionMode = System::Windows::Forms::SelectionMode::MultiExtended;
			this->FileListTxt->Size = System::Drawing::Size(450, 485);
			this->FileListTxt->Sorted = true;
			this->FileListTxt->TabIndex = 9;
			this->FileListTxt->MouseUp += gcnew System::Windows::Forms::MouseEventHandler(this, &FitsFound::FileListTxt_MouseUp);
			this->FileListTxt->MouseClick += gcnew System::Windows::Forms::MouseEventHandler(this, &FitsFound::FileListTxt_MouseClick);
			this->FileListTxt->SelectedIndexChanged += gcnew System::EventHandler(this, &FitsFound::FileListTxt_SelectedIndexChanged);
			// 
			// FoundListContextMenu
			// 
			this->FoundListContextMenu->Enabled = false;
			this->FoundListContextMenu->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) {this->OpenFolderContextItem});
			this->FoundListContextMenu->Name = L"FoundListContextMenu";
			this->FoundListContextMenu->Size = System::Drawing::Size(140, 26);
			// 
			// OpenFolderContextItem
			// 
			this->OpenFolderContextItem->Enabled = false;
			this->OpenFolderContextItem->Name = L"OpenFolderContextItem";
			this->OpenFolderContextItem->Size = System::Drawing::Size(139, 22);
			this->OpenFolderContextItem->Text = L"Open Folder";
			this->OpenFolderContextItem->Click += gcnew System::EventHandler(this, &FitsFound::OpenFolderContextItem_Click);
			// 
			// AddImageSetBtn
			// 
			this->AddImageSetBtn->Anchor = System::Windows::Forms::AnchorStyles::None;
			this->AddImageSetBtn->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Center;
			this->AddImageSetBtn->DialogResult = System::Windows::Forms::DialogResult::Ignore;
			this->AddImageSetBtn->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->AddImageSetBtn->Location = System::Drawing::Point(259, 568);
			this->AddImageSetBtn->Name = L"AddImageSetBtn";
			this->AddImageSetBtn->Size = System::Drawing::Size(84, 43);
			this->AddImageSetBtn->TabIndex = 19;
			this->AddImageSetBtn->Text = L"&Add Selected to Image Set";
			this->AddImageSetBtn->UseVisualStyleBackColor = true;
			// 
			// LoadImageSetBtn
			// 
			this->LoadImageSetBtn->Anchor = System::Windows::Forms::AnchorStyles::None;
			this->LoadImageSetBtn->DialogResult = System::Windows::Forms::DialogResult::OK;
			this->LoadImageSetBtn->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->LoadImageSetBtn->Location = System::Drawing::Point(162, 568);
			this->LoadImageSetBtn->Name = L"LoadImageSetBtn";
			this->LoadImageSetBtn->Size = System::Drawing::Size(84, 43);
			this->LoadImageSetBtn->TabIndex = 18;
			this->LoadImageSetBtn->Text = L"&Load Selected as Image Set";
			this->LoadImageSetBtn->UseVisualStyleBackColor = true;
			// 
			// MoveListBtn
			// 
			this->MoveListBtn->Anchor = System::Windows::Forms::AnchorStyles::None;
			this->MoveListBtn->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Center;
			this->MoveListBtn->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->MoveListBtn->ForeColor = System::Drawing::SystemColors::ControlText;
			this->MoveListBtn->Location = System::Drawing::Point(171, 617);
			this->MoveListBtn->Name = L"MoveListBtn";
			this->MoveListBtn->Size = System::Drawing::Size(75, 23);
			this->MoveListBtn->TabIndex = 21;
			this->MoveListBtn->Text = L"Move";
			this->MoveListBtn->UseVisualStyleBackColor = true;
			this->MoveListBtn->Click += gcnew System::EventHandler(this, &FitsFound::MoveListBtn_Click);
			// 
			// CopyListBtn
			// 
			this->CopyListBtn->Anchor = System::Windows::Forms::AnchorStyles::None;
			this->CopyListBtn->BackColor = System::Drawing::SystemColors::Control;
			this->CopyListBtn->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Center;
			this->CopyListBtn->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->CopyListBtn->Location = System::Drawing::Point(259, 617);
			this->CopyListBtn->Name = L"CopyListBtn";
			this->CopyListBtn->Size = System::Drawing::Size(75, 23);
			this->CopyListBtn->TabIndex = 20;
			this->CopyListBtn->Text = L"Copy";
			this->CopyListBtn->UseVisualStyleBackColor = true;
			this->CopyListBtn->Click += gcnew System::EventHandler(this, &FitsFound::CopyListBtn_Click);
			// 
			// FileCopyBGWrkr
			// 
			this->FileCopyBGWrkr->WorkerReportsProgress = true;
			this->FileCopyBGWrkr->WorkerSupportsCancellation = true;
			this->FileCopyBGWrkr->DoWork += gcnew System::ComponentModel::DoWorkEventHandler(this, &FitsFound::FileCopyBGWrkr_DoWork);
			this->FileCopyBGWrkr->RunWorkerCompleted += gcnew System::ComponentModel::RunWorkerCompletedEventHandler(this, &FitsFound::FileCopyBGWrkr_RunWorkerCompleted);
			this->FileCopyBGWrkr->ProgressChanged += gcnew System::ComponentModel::ProgressChangedEventHandler(this, &FitsFound::FileCopyBGWrkr_ProgressChanged);
			// 
			// FitsFound
			// 
			this->AcceptButton = this->LoadImageSetBtn;
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::None;
			this->CancelButton = this->CancelBtn;
			this->ClientSize = System::Drawing::Size(504, 691);
			this->Controls->Add(this->MoveListBtn);
			this->Controls->Add(this->CopyListBtn);
			this->Controls->Add(this->AddImageSetBtn);
			this->Controls->Add(this->LoadImageSetBtn);
			this->Controls->Add(this->NumSelectTxt);
			this->Controls->Add(this->ClearAllBtn);
			this->Controls->Add(this->SaveListBtn);
			this->Controls->Add(this->SelectAllBtn);
			this->Controls->Add(this->CancelBtn);
			this->Controls->Add(this->NumFilesTxt);
			this->Controls->Add(this->FileListTxt);
			this->Icon = (cli::safe_cast<System::Drawing::Icon^  >(resources->GetObject(L"$this.Icon")));
			this->MaximizeBox = false;
			this->MinimumSize = System::Drawing::Size(226, 662);
			this->Name = L"FitsFound";
			this->SizeGripStyle = System::Windows::Forms::SizeGripStyle::Show;
			this->StartPosition = System::Windows::Forms::FormStartPosition::CenterScreen;
			this->Text = L"Found FITS Files";
			this->Load += gcnew System::EventHandler(this, &FitsFound::FitsFound_Load);
			this->FormClosing += gcnew System::Windows::Forms::FormClosingEventHandler(this, &FitsFound::FitsFound_FormClosing);
			this->FoundListContextMenu->ResumeLayout(false);
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion

public:

JPWaitBar::WaitBar^ WAITBAR;
void SetFoundFiles(array<String^>^ FoundFiles);
void SetFoundFiles(String^ FullFileName);

private:
array<String^>^ FOUNDFILES;


private: System::Void FileListTxt_SelectedIndexChanged(System::Object^  sender, System::EventArgs^  e);
private: System::Void SelectAllBtn_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void ClearAllBtn_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void SaveListBtn_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void CancelBtn_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void CopyListBtn_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void MoveListBtn_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void FileListTxt_MouseClick(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e);
private: System::Void FileListTxt_MouseUp(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e);
private: System::Void OpenFolderContextItem_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void FileCopyBGWrkr_DoWork(System::Object^  sender, System::ComponentModel::DoWorkEventArgs^  e);
private: System::Void FileCopyBGWrkr_ProgressChanged(System::Object^  sender, System::ComponentModel::ProgressChangedEventArgs^  e);
private: System::Void FileCopyBGWrkr_RunWorkerCompleted(System::Object^  sender, System::ComponentModel::RunWorkerCompletedEventArgs^  e);
private: System::Void FitsFound_FormClosing(System::Object^  sender, System::Windows::Forms::FormClosingEventArgs^  e);
private: System::Void FitsFound_Load(System::Object^  sender, System::EventArgs^  e);
};
}
