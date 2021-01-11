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
#include "JPFITS.h"

using namespace System;
using namespace System::ComponentModel;
using namespace System::Collections;
using namespace System::Windows::Forms;
using namespace System::Data;
using namespace System::Drawing;
using namespace Microsoft::Win32;


namespace JPFITS {

	/// <summary>A class to search for FITS files specifying the file name search pattern and FITS keywords and keyvalues in the primary header.</summary>
	public ref class FitsFinder : public System::Windows::Forms::Form
	{
	public:
		FitsFinder(void)
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
		~FitsFinder()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Label^  label11;
	protected: 
	public: System::Windows::Forms::ComboBox^  ExtensionDrop;
	private: 
	public: System::Windows::Forms::CheckBox^  SubFoldersChck;
	public: System::Windows::Forms::TextBox^  FileTemplateTxt;
	private: System::Windows::Forms::Label^  label10;
	public: 
	private: System::Windows::Forms::Button^  CancelBtn;
	private: System::Windows::Forms::Button^  FindBtn;
	public: System::Windows::Forms::RichTextBox^  Key4Value;
	private: 
	private: System::Windows::Forms::Label^  label8;
	public: 
	private: System::Windows::Forms::Label^  label9;

	private: 
	public: System::Windows::Forms::RichTextBox^  Key3Value;
	private: System::Windows::Forms::Label^  label6;
	public: 
	private: System::Windows::Forms::Label^  label7;

	private: 
	public: System::Windows::Forms::RichTextBox^  Key2Value;
	private: System::Windows::Forms::Label^  label4;
	public: 
	private: System::Windows::Forms::Label^  label5;

	private: 
	public: System::Windows::Forms::RichTextBox^  Key1Value;
	private: System::Windows::Forms::Label^  label3;
	public: 
	private: System::Windows::Forms::Label^  label2;

	private: 
	private: System::Windows::Forms::Label^  label1;
	public: 
	public: System::Windows::Forms::RichTextBox^  DirectoryTxt;
	private: System::ComponentModel::BackgroundWorker^  FitsFinderWrkr;
	public: System::Windows::Forms::TextBox^  Key1;
	public: System::Windows::Forms::TextBox^  Key2;
	public: System::Windows::Forms::TextBox^  Key3;
	public: System::Windows::Forms::TextBox^  Key4;
	public: System::Windows::Forms::CheckBox^  CustomExtensionChck;
	private: System::Windows::Forms::TextBox^  CustomExtensionTxtBox;
	public: 

	public: 

	private: 

	private: 
	public: 
	private: 

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(FitsFinder::typeid));
			this->label11 = (gcnew System::Windows::Forms::Label());
			this->ExtensionDrop = (gcnew System::Windows::Forms::ComboBox());
			this->SubFoldersChck = (gcnew System::Windows::Forms::CheckBox());
			this->FileTemplateTxt = (gcnew System::Windows::Forms::TextBox());
			this->label10 = (gcnew System::Windows::Forms::Label());
			this->CancelBtn = (gcnew System::Windows::Forms::Button());
			this->FindBtn = (gcnew System::Windows::Forms::Button());
			this->Key4Value = (gcnew System::Windows::Forms::RichTextBox());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->label9 = (gcnew System::Windows::Forms::Label());
			this->Key3Value = (gcnew System::Windows::Forms::RichTextBox());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->Key2Value = (gcnew System::Windows::Forms::RichTextBox());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->Key1Value = (gcnew System::Windows::Forms::RichTextBox());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->DirectoryTxt = (gcnew System::Windows::Forms::RichTextBox());
			this->FitsFinderWrkr = (gcnew System::ComponentModel::BackgroundWorker());
			this->Key1 = (gcnew System::Windows::Forms::TextBox());
			this->Key2 = (gcnew System::Windows::Forms::TextBox());
			this->Key3 = (gcnew System::Windows::Forms::TextBox());
			this->Key4 = (gcnew System::Windows::Forms::TextBox());
			this->CustomExtensionChck = (gcnew System::Windows::Forms::CheckBox());
			this->CustomExtensionTxtBox = (gcnew System::Windows::Forms::TextBox());
			this->SuspendLayout();
			// 
			// label11
			// 
			this->label11->AutoSize = true;
			this->label11->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Underline, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label11->Location = System::Drawing::Point(100, 99);
			this->label11->Name = L"label11";
			this->label11->Size = System::Drawing::Size(56, 13);
			this->label11->TabIndex = 45;
			this->label11->Text = L"Extension:";
			// 
			// ExtensionDrop
			// 
			this->ExtensionDrop->DropDownStyle = System::Windows::Forms::ComboBoxStyle::DropDownList;
			this->ExtensionDrop->FormattingEnabled = true;
			this->ExtensionDrop->Items->AddRange(gcnew cli::array< System::Object^  >(4) { L".fts", L".fit", L".fits", L".raw" });
			this->ExtensionDrop->Location = System::Drawing::Point(103, 115);
			this->ExtensionDrop->Name = L"ExtensionDrop";
			this->ExtensionDrop->Size = System::Drawing::Size(45, 21);
			this->ExtensionDrop->TabIndex = 32;
			// 
			// SubFoldersChck
			// 
			this->SubFoldersChck->AutoSize = true;
			this->SubFoldersChck->Location = System::Drawing::Point(162, 72);
			this->SubFoldersChck->Name = L"SubFoldersChck";
			this->SubFoldersChck->Size = System::Drawing::Size(136, 17);
			this->SubFoldersChck->TabIndex = 30;
			this->SubFoldersChck->Text = L"Include Sub Directories";
			this->SubFoldersChck->UseVisualStyleBackColor = true;
			// 
			// FileTemplateTxt
			// 
			this->FileTemplateTxt->Location = System::Drawing::Point(14, 116);
			this->FileTemplateTxt->MaxLength = 100;
			this->FileTemplateTxt->Name = L"FileTemplateTxt";
			this->FileTemplateTxt->Size = System::Drawing::Size(71, 20);
			this->FileTemplateTxt->TabIndex = 31;
			// 
			// label10
			// 
			this->label10->AutoSize = true;
			this->label10->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Underline, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label10->Location = System::Drawing::Point(11, 100);
			this->label10->Name = L"label10";
			this->label10->Size = System::Drawing::Size(73, 13);
			this->label10->TabIndex = 29;
			this->label10->Text = L"File Template:";
			// 
			// CancelBtn
			// 
			this->CancelBtn->DialogResult = System::Windows::Forms::DialogResult::Cancel;
			this->CancelBtn->Location = System::Drawing::Point(162, 305);
			this->CancelBtn->Name = L"CancelBtn";
			this->CancelBtn->Size = System::Drawing::Size(75, 23);
			this->CancelBtn->TabIndex = 42;
			this->CancelBtn->Text = L"&Cancel";
			this->CancelBtn->UseVisualStyleBackColor = true;
			this->CancelBtn->Click += gcnew System::EventHandler(this, &FitsFinder::CancelBtn_Click);
			// 
			// FindBtn
			// 
			this->FindBtn->Location = System::Drawing::Point(73, 305);
			this->FindBtn->Name = L"FindBtn";
			this->FindBtn->Size = System::Drawing::Size(71, 23);
			this->FindBtn->TabIndex = 41;
			this->FindBtn->Text = L"&Find";
			this->FindBtn->UseVisualStyleBackColor = true;
			this->FindBtn->Click += gcnew System::EventHandler(this, &FitsFinder::FindBtn_Click);
			// 
			// Key4Value
			// 
			this->Key4Value->Location = System::Drawing::Point(103, 271);
			this->Key4Value->MaxLength = 30;
			this->Key4Value->Multiline = false;
			this->Key4Value->Name = L"Key4Value";
			this->Key4Value->ScrollBars = System::Windows::Forms::RichTextBoxScrollBars::None;
			this->Key4Value->Size = System::Drawing::Size(197, 18);
			this->Key4Value->TabIndex = 40;
			this->Key4Value->Text = L"";
			// 
			// label8
			// 
			this->label8->AutoSize = true;
			this->label8->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Underline, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label8->Location = System::Drawing::Point(100, 255);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(108, 13);
			this->label8->TabIndex = 25;
			this->label8->Text = L"Keyname Four Value:";
			// 
			// label9
			// 
			this->label9->AutoSize = true;
			this->label9->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Underline, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label9->Location = System::Drawing::Point(11, 255);
			this->label9->Name = L"label9";
			this->label9->Size = System::Drawing::Size(78, 13);
			this->label9->TabIndex = 24;
			this->label9->Text = L"Keyname Four:";
			// 
			// Key3Value
			// 
			this->Key3Value->Location = System::Drawing::Point(103, 234);
			this->Key3Value->MaxLength = 30;
			this->Key3Value->Multiline = false;
			this->Key3Value->Name = L"Key3Value";
			this->Key3Value->ScrollBars = System::Windows::Forms::RichTextBoxScrollBars::None;
			this->Key3Value->Size = System::Drawing::Size(197, 18);
			this->Key3Value->TabIndex = 38;
			this->Key3Value->Text = L"";
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Underline, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label6->Location = System::Drawing::Point(100, 218);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(115, 13);
			this->label6->TabIndex = 27;
			this->label6->Text = L"Keyname Three Value:";
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Underline, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label7->Location = System::Drawing::Point(11, 218);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(85, 13);
			this->label7->TabIndex = 23;
			this->label7->Text = L"Keyname Three:";
			// 
			// Key2Value
			// 
			this->Key2Value->Location = System::Drawing::Point(103, 197);
			this->Key2Value->MaxLength = 30;
			this->Key2Value->Multiline = false;
			this->Key2Value->Name = L"Key2Value";
			this->Key2Value->ScrollBars = System::Windows::Forms::RichTextBoxScrollBars::None;
			this->Key2Value->Size = System::Drawing::Size(197, 18);
			this->Key2Value->TabIndex = 36;
			this->Key2Value->Text = L"";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Underline, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label4->Location = System::Drawing::Point(100, 181);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(108, 13);
			this->label4->TabIndex = 21;
			this->label4->Text = L"Keyname Two Value:";
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Underline, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label5->Location = System::Drawing::Point(11, 181);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(78, 13);
			this->label5->TabIndex = 22;
			this->label5->Text = L"Keyname Two:";
			// 
			// Key1Value
			// 
			this->Key1Value->Location = System::Drawing::Point(103, 160);
			this->Key1Value->MaxLength = 30;
			this->Key1Value->Multiline = false;
			this->Key1Value->Name = L"Key1Value";
			this->Key1Value->ScrollBars = System::Windows::Forms::RichTextBoxScrollBars::None;
			this->Key1Value->Size = System::Drawing::Size(197, 18);
			this->Key1Value->TabIndex = 34;
			this->Key1Value->Text = L"";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Underline, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label3->Location = System::Drawing::Point(100, 144);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(107, 13);
			this->label3->TabIndex = 30;
			this->label3->Text = L"Keyname One Value:";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Underline, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label2->Location = System::Drawing::Point(11, 144);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(77, 13);
			this->label2->TabIndex = 28;
			this->label2->Text = L"Keyname One:";
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Underline, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label1->Location = System::Drawing::Point(11, 27);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(52, 13);
			this->label1->TabIndex = 26;
			this->label1->Text = L"Directory:";
			// 
			// DirectoryTxt
			// 
			this->DirectoryTxt->BackColor = System::Drawing::SystemColors::Control;
			this->DirectoryTxt->DetectUrls = false;
			this->DirectoryTxt->Location = System::Drawing::Point(69, 15);
			this->DirectoryTxt->MaxLength = 500;
			this->DirectoryTxt->Name = L"DirectoryTxt";
			this->DirectoryTxt->ReadOnly = true;
			this->DirectoryTxt->Size = System::Drawing::Size(233, 51);
			this->DirectoryTxt->TabIndex = 29;
			this->DirectoryTxt->Text = L"";
			this->DirectoryTxt->Click += gcnew System::EventHandler(this, &FitsFinder::DirectoryTxt_Click);
			// 
			// FitsFinderWrkr
			// 
			this->FitsFinderWrkr->WorkerReportsProgress = true;
			this->FitsFinderWrkr->WorkerSupportsCancellation = true;
			this->FitsFinderWrkr->DoWork += gcnew System::ComponentModel::DoWorkEventHandler(this, &FitsFinder::FitsFinderWrkr_DoWork);
			this->FitsFinderWrkr->ProgressChanged += gcnew System::ComponentModel::ProgressChangedEventHandler(this, &FitsFinder::FitsFinderWrkr_ProgressChanged);
			this->FitsFinderWrkr->RunWorkerCompleted += gcnew System::ComponentModel::RunWorkerCompletedEventHandler(this, &FitsFinder::FitsFinderWrkr_RunWorkerCompleted);
			// 
			// Key1
			// 
			this->Key1->CharacterCasing = System::Windows::Forms::CharacterCasing::Upper;
			this->Key1->Location = System::Drawing::Point(14, 160);
			this->Key1->MaxLength = 8;
			this->Key1->Name = L"Key1";
			this->Key1->Size = System::Drawing::Size(71, 20);
			this->Key1->TabIndex = 33;
			// 
			// Key2
			// 
			this->Key2->CharacterCasing = System::Windows::Forms::CharacterCasing::Upper;
			this->Key2->Location = System::Drawing::Point(13, 197);
			this->Key2->MaxLength = 8;
			this->Key2->Name = L"Key2";
			this->Key2->Size = System::Drawing::Size(71, 20);
			this->Key2->TabIndex = 35;
			// 
			// Key3
			// 
			this->Key3->CharacterCasing = System::Windows::Forms::CharacterCasing::Upper;
			this->Key3->Location = System::Drawing::Point(14, 234);
			this->Key3->MaxLength = 8;
			this->Key3->Name = L"Key3";
			this->Key3->Size = System::Drawing::Size(71, 20);
			this->Key3->TabIndex = 37;
			// 
			// Key4
			// 
			this->Key4->CharacterCasing = System::Windows::Forms::CharacterCasing::Upper;
			this->Key4->Location = System::Drawing::Point(14, 271);
			this->Key4->MaxLength = 8;
			this->Key4->Name = L"Key4";
			this->Key4->Size = System::Drawing::Size(71, 20);
			this->Key4->TabIndex = 39;
			// 
			// CustomExtensionChck
			// 
			this->CustomExtensionChck->AutoSize = true;
			this->CustomExtensionChck->Location = System::Drawing::Point(162, 99);
			this->CustomExtensionChck->Name = L"CustomExtensionChck";
			this->CustomExtensionChck->Size = System::Drawing::Size(61, 17);
			this->CustomExtensionChck->TabIndex = 46;
			this->CustomExtensionChck->Text = L"Custom";
			this->CustomExtensionChck->UseVisualStyleBackColor = true;
			this->CustomExtensionChck->CheckedChanged += gcnew System::EventHandler(this, &FitsFinder::CustomExtensionChck_CheckedChanged);
			// 
			// CustomExtensionTxtBox
			// 
			this->CustomExtensionTxtBox->Enabled = false;
			this->CustomExtensionTxtBox->Location = System::Drawing::Point(162, 116);
			this->CustomExtensionTxtBox->Name = L"CustomExtensionTxtBox";
			this->CustomExtensionTxtBox->Size = System::Drawing::Size(100, 20);
			this->CustomExtensionTxtBox->TabIndex = 47;
			this->CustomExtensionTxtBox->Text = L".events";
			this->CustomExtensionTxtBox->TextChanged += gcnew System::EventHandler(this, &FitsFinder::CustomExtensionTxtBox_TextChanged);
			// 
			// FitsFinder
			// 
			this->AcceptButton = this->FindBtn;
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->CancelButton = this->CancelBtn;
			this->ClientSize = System::Drawing::Size(313, 343);
			this->ControlBox = false;
			this->Controls->Add(this->CustomExtensionTxtBox);
			this->Controls->Add(this->CustomExtensionChck);
			this->Controls->Add(this->Key4);
			this->Controls->Add(this->Key3);
			this->Controls->Add(this->Key2);
			this->Controls->Add(this->Key1);
			this->Controls->Add(this->label11);
			this->Controls->Add(this->ExtensionDrop);
			this->Controls->Add(this->SubFoldersChck);
			this->Controls->Add(this->FileTemplateTxt);
			this->Controls->Add(this->label10);
			this->Controls->Add(this->CancelBtn);
			this->Controls->Add(this->FindBtn);
			this->Controls->Add(this->Key4Value);
			this->Controls->Add(this->label8);
			this->Controls->Add(this->label9);
			this->Controls->Add(this->Key3Value);
			this->Controls->Add(this->label6);
			this->Controls->Add(this->label7);
			this->Controls->Add(this->Key2Value);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->label5);
			this->Controls->Add(this->Key1Value);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->label1);
			this->Controls->Add(this->DirectoryTxt);
			this->FormBorderStyle = System::Windows::Forms::FormBorderStyle::FixedToolWindow;
			this->Icon = (cli::safe_cast<System::Drawing::Icon^>(resources->GetObject(L"$this.Icon")));
			this->Name = L"FitsFinder";
			this->SizeGripStyle = System::Windows::Forms::SizeGripStyle::Hide;
			this->StartPosition = System::Windows::Forms::FormStartPosition::CenterParent;
			this->Text = L"Search for FITS File(s)...";
			this->Load += gcnew System::EventHandler(this, &FitsFinder::FitsFinder_Load);
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion




public:

JPWaitBar::WaitBar^ WAITBAR;
array<String^>^ FOUNDFILES;


private: System::Void DirectoryTxt_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void FitsFinder_Load(System::Object^  sender, System::EventArgs^  e);
private: System::Void FindBtn_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void FitsFinderWrkr_DoWork(System::Object^  sender, System::ComponentModel::DoWorkEventArgs^  e);
private: System::Void FitsFinderWrkr_ProgressChanged(System::Object^  sender, System::ComponentModel::ProgressChangedEventArgs^  e);
private: System::Void FitsFinderWrkr_RunWorkerCompleted(System::Object^  sender, System::ComponentModel::RunWorkerCompletedEventArgs^  e);
private: System::Void CancelBtn_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void CustomExtensionChck_CheckedChanged(System::Object^  sender, System::EventArgs^  e)
{
	if (CustomExtensionChck->Checked)
	{
		CustomExtensionTxtBox->Enabled = true;
		ExtensionDrop->Enabled = false;
	}
	else
	{
		CustomExtensionTxtBox->Enabled = false;
		ExtensionDrop->Enabled = true;
	}
}
private: System::Void CustomExtensionTxtBox_TextChanged(System::Object^  sender, System::EventArgs^  e);


};
}
