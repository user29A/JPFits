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
#include "FitsExtensionTablePlotter.h"

namespace JPFITS {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for FitsExtensionTableViewer
	/// </summary>
	public ref class FitsExtensionTableViewer : public System::Windows::Forms::Form
	{
	public:
		FitsExtensionTableViewer(String^ FileName)
		{
			InitializeComponent();
			OpenFITSImage(FileName);
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~FitsExtensionTableViewer()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::DataGridView^  ExtensionTableGrid;
	private: System::Windows::Forms::MenuStrip^  menuStrip1;
	private: System::Windows::Forms::ToolStripMenuItem^  MenuChooseTableEntries;

	private: System::Windows::Forms::ToolStripMenuItem^  ViewAllChck;
	private: System::Windows::Forms::ToolStripSeparator^  toolStripSeparator1;
	private: System::Windows::Forms::ToolStripMenuItem^  PlotEntryMenu;


	private: System::Windows::Forms::ToolStripComboBox^  XDrop;
	private: System::Windows::Forms::ToolStripSeparator^  toolStripSeparator2;

	private: System::Windows::Forms::ToolStripComboBox^  YDrop;
	private: System::Windows::Forms::ToolStripMenuItem^  PlotXChck;

	private: System::Windows::Forms::ToolStripMenuItem^  toolStripMenuItem2;
	private: System::Windows::Forms::ToolStripSeparator^  toolStripSeparator3;
	private: System::Windows::Forms::ToolStripMenuItem^  PlotMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  MenuChooseTable;

	private: System::Windows::Forms::ToolStripMenuItem^  MenuFile;
	private: System::Windows::Forms::ToolStripMenuItem^  FileOpenMenu;
	private: System::Windows::Forms::Button^  button1;
	private: System::Windows::Forms::ToolStripMenuItem^  ViewHeaderMenu;
	private: System::Windows::Forms::ListBox^  HeaderListBox;



	protected:

	protected:

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
			System::Windows::Forms::DataGridViewCellStyle^  dataGridViewCellStyle1 = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(FitsExtensionTableViewer::typeid));
			this->ExtensionTableGrid = (gcnew System::Windows::Forms::DataGridView());
			this->menuStrip1 = (gcnew System::Windows::Forms::MenuStrip());
			this->MenuFile = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->FileOpenMenu = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->MenuChooseTable = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->MenuChooseTableEntries = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->ViewAllChck = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->toolStripSeparator1 = (gcnew System::Windows::Forms::ToolStripSeparator());
			this->PlotEntryMenu = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->PlotXChck = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->XDrop = (gcnew System::Windows::Forms::ToolStripComboBox());
			this->toolStripSeparator2 = (gcnew System::Windows::Forms::ToolStripSeparator());
			this->toolStripMenuItem2 = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->YDrop = (gcnew System::Windows::Forms::ToolStripComboBox());
			this->toolStripSeparator3 = (gcnew System::Windows::Forms::ToolStripSeparator());
			this->PlotMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->ViewHeaderMenu = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->HeaderListBox = (gcnew System::Windows::Forms::ListBox());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->ExtensionTableGrid))->BeginInit();
			this->menuStrip1->SuspendLayout();
			this->SuspendLayout();
			// 
			// ExtensionTableGrid
			// 
			this->ExtensionTableGrid->AllowUserToAddRows = false;
			this->ExtensionTableGrid->AllowUserToDeleteRows = false;
			this->ExtensionTableGrid->ClipboardCopyMode = System::Windows::Forms::DataGridViewClipboardCopyMode::EnableWithoutHeaderText;
			this->ExtensionTableGrid->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->ExtensionTableGrid->Dock = System::Windows::Forms::DockStyle::Fill;
			this->ExtensionTableGrid->Location = System::Drawing::Point(0, 24);
			this->ExtensionTableGrid->Name = L"ExtensionTableGrid";
			this->ExtensionTableGrid->ReadOnly = true;
			dataGridViewCellStyle1->Alignment = System::Windows::Forms::DataGridViewContentAlignment::MiddleLeft;
			dataGridViewCellStyle1->BackColor = System::Drawing::SystemColors::Control;
			dataGridViewCellStyle1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular,
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			dataGridViewCellStyle1->ForeColor = System::Drawing::SystemColors::WindowText;
			dataGridViewCellStyle1->Format = L"N0";
			dataGridViewCellStyle1->NullValue = nullptr;
			dataGridViewCellStyle1->SelectionBackColor = System::Drawing::SystemColors::Highlight;
			dataGridViewCellStyle1->SelectionForeColor = System::Drawing::SystemColors::HighlightText;
			dataGridViewCellStyle1->WrapMode = System::Windows::Forms::DataGridViewTriState::True;
			this->ExtensionTableGrid->RowHeadersDefaultCellStyle = dataGridViewCellStyle1;
			this->ExtensionTableGrid->RowHeadersWidth = 75;
			this->ExtensionTableGrid->Size = System::Drawing::Size(583, 329);
			this->ExtensionTableGrid->TabIndex = 0;
			this->ExtensionTableGrid->VirtualMode = true;
			this->ExtensionTableGrid->CellValueNeeded += gcnew System::Windows::Forms::DataGridViewCellValueEventHandler(this, &FitsExtensionTableViewer::ExtensionTableGrid_CellValueNeeded);
			this->ExtensionTableGrid->NewRowNeeded += gcnew System::Windows::Forms::DataGridViewRowEventHandler(this, &FitsExtensionTableViewer::ExtensionTableGrid_NewRowNeeded);
			this->ExtensionTableGrid->Scroll += gcnew System::Windows::Forms::ScrollEventHandler(this, &FitsExtensionTableViewer::ExtensionTableGrid_Scroll);
			// 
			// menuStrip1
			// 
			this->menuStrip1->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(5) {
				this->MenuFile, this->MenuChooseTable,
					this->MenuChooseTableEntries, this->PlotEntryMenu, this->ViewHeaderMenu
			});
			this->menuStrip1->Location = System::Drawing::Point(0, 0);
			this->menuStrip1->Name = L"menuStrip1";
			this->menuStrip1->Size = System::Drawing::Size(583, 24);
			this->menuStrip1->TabIndex = 1;
			this->menuStrip1->Text = L"menuStrip1";
			// 
			// MenuFile
			// 
			this->MenuFile->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) { this->FileOpenMenu });
			this->MenuFile->Name = L"MenuFile";
			this->MenuFile->Size = System::Drawing::Size(37, 20);
			this->MenuFile->Text = L"File";
			// 
			// FileOpenMenu
			// 
			this->FileOpenMenu->Name = L"FileOpenMenu";
			this->FileOpenMenu->Size = System::Drawing::Size(103, 22);
			this->FileOpenMenu->Text = L"Open";
			this->FileOpenMenu->Click += gcnew System::EventHandler(this, &FitsExtensionTableViewer::FileOpenMenu_Click);
			// 
			// MenuChooseTable
			// 
			this->MenuChooseTable->Name = L"MenuChooseTable";
			this->MenuChooseTable->Size = System::Drawing::Size(51, 20);
			this->MenuChooseTable->Text = L"Tables";
			// 
			// MenuChooseTableEntries
			// 
			this->MenuChooseTableEntries->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(2) {
				this->ViewAllChck,
					this->toolStripSeparator1
			});
			this->MenuChooseTableEntries->Name = L"MenuChooseTableEntries";
			this->MenuChooseTableEntries->Size = System::Drawing::Size(84, 20);
			this->MenuChooseTableEntries->Text = L"Table Entries";
			this->MenuChooseTableEntries->DropDownItemClicked += gcnew System::Windows::Forms::ToolStripItemClickedEventHandler(this, &FitsExtensionTableViewer::MenuChooseTableEntries_DropDownItemClicked);
			this->MenuChooseTableEntries->Click += gcnew System::EventHandler(this, &FitsExtensionTableViewer::MenuChooseTableEntries_Click);
			this->MenuChooseTableEntries->MouseLeave += gcnew System::EventHandler(this, &FitsExtensionTableViewer::MenuChooseTableEntries_MouseLeave);
			// 
			// ViewAllChck
			// 
			this->ViewAllChck->Name = L"ViewAllChck";
			this->ViewAllChck->Size = System::Drawing::Size(131, 22);
			this->ViewAllChck->Text = L"View None";
			this->ViewAllChck->CheckedChanged += gcnew System::EventHandler(this, &FitsExtensionTableViewer::ViewAllChck_CheckedChanged);
			this->ViewAllChck->Click += gcnew System::EventHandler(this, &FitsExtensionTableViewer::ViewAllChck_Click);
			// 
			// toolStripSeparator1
			// 
			this->toolStripSeparator1->Name = L"toolStripSeparator1";
			this->toolStripSeparator1->Size = System::Drawing::Size(128, 6);
			// 
			// PlotEntryMenu
			// 
			this->PlotEntryMenu->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(7) {
				this->PlotXChck,
					this->XDrop, this->toolStripSeparator2, this->toolStripMenuItem2, this->YDrop, this->toolStripSeparator3, this->PlotMenuItem
			});
			this->PlotEntryMenu->Name = L"PlotEntryMenu";
			this->PlotEntryMenu->Size = System::Drawing::Size(70, 20);
			this->PlotEntryMenu->Text = L"Plot Entry";
			// 
			// PlotXChck
			// 
			this->PlotXChck->Checked = true;
			this->PlotXChck->CheckOnClick = true;
			this->PlotXChck->CheckState = System::Windows::Forms::CheckState::Checked;
			this->PlotXChck->Name = L"PlotXChck";
			this->PlotXChck->Size = System::Drawing::Size(181, 22);
			this->PlotXChck->Text = L"X";
			this->PlotXChck->CheckedChanged += gcnew System::EventHandler(this, &FitsExtensionTableViewer::toolStripMenuItem1_CheckedChanged);
			// 
			// XDrop
			// 
			this->XDrop->BackColor = System::Drawing::Color::Gainsboro;
			this->XDrop->DropDownStyle = System::Windows::Forms::ComboBoxStyle::DropDownList;
			this->XDrop->Name = L"XDrop";
			this->XDrop->Size = System::Drawing::Size(121, 23);
			// 
			// toolStripSeparator2
			// 
			this->toolStripSeparator2->Name = L"toolStripSeparator2";
			this->toolStripSeparator2->Size = System::Drawing::Size(178, 6);
			// 
			// toolStripMenuItem2
			// 
			this->toolStripMenuItem2->Name = L"toolStripMenuItem2";
			this->toolStripMenuItem2->Size = System::Drawing::Size(181, 22);
			this->toolStripMenuItem2->Text = L"Y";
			// 
			// YDrop
			// 
			this->YDrop->BackColor = System::Drawing::Color::Gainsboro;
			this->YDrop->DropDownStyle = System::Windows::Forms::ComboBoxStyle::DropDownList;
			this->YDrop->Name = L"YDrop";
			this->YDrop->Size = System::Drawing::Size(121, 23);
			// 
			// toolStripSeparator3
			// 
			this->toolStripSeparator3->Name = L"toolStripSeparator3";
			this->toolStripSeparator3->Size = System::Drawing::Size(178, 6);
			// 
			// PlotMenuItem
			// 
			this->PlotMenuItem->Name = L"PlotMenuItem";
			this->PlotMenuItem->Size = System::Drawing::Size(181, 22);
			this->PlotMenuItem->Text = L"Plot Selection";
			this->PlotMenuItem->Click += gcnew System::EventHandler(this, &FitsExtensionTableViewer::PlotMenuItem_Click);
			// 
			// ViewHeaderMenu
			// 
			this->ViewHeaderMenu->Name = L"ViewHeaderMenu";
			this->ViewHeaderMenu->Size = System::Drawing::Size(85, 20);
			this->ViewHeaderMenu->Text = L"View Header";
			this->ViewHeaderMenu->Click += gcnew System::EventHandler(this, &FitsExtensionTableViewer::ViewHeaderMenu_Click);
			// 
			// button1
			// 
			this->button1->DialogResult = System::Windows::Forms::DialogResult::Cancel;
			this->button1->Location = System::Drawing::Point(239, 133);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(75, 23);
			this->button1->TabIndex = 2;
			this->button1->Text = L"button1";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &FitsExtensionTableViewer::button1_Click);
			// 
			// HeaderListBox
			// 
			this->HeaderListBox->Dock = System::Windows::Forms::DockStyle::Fill;
			this->HeaderListBox->FormattingEnabled = true;
			this->HeaderListBox->Location = System::Drawing::Point(0, 0);
			this->HeaderListBox->Name = L"HeaderListBox";
			this->HeaderListBox->Size = System::Drawing::Size(583, 353);
			this->HeaderListBox->TabIndex = 3;
			// 
			// FitsExtensionTableViewer
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->CancelButton = this->button1;
			this->ClientSize = System::Drawing::Size(583, 353);
			this->Controls->Add(this->ExtensionTableGrid);
			this->Controls->Add(this->menuStrip1);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->HeaderListBox);
			this->Icon = (cli::safe_cast<System::Drawing::Icon^>(resources->GetObject(L"$this.Icon")));
			this->MainMenuStrip = this->menuStrip1;
			this->Name = L"FitsExtensionTableViewer";
			this->StartPosition = System::Windows::Forms::FormStartPosition::CenterParent;
			this->Text = L"FitsExtensionTableViewer";
			this->FormClosing += gcnew System::Windows::Forms::FormClosingEventHandler(this, &FitsExtensionTableViewer::FitsExtensionTableViewer_FormClosing);
			this->Load += gcnew System::EventHandler(this, &FitsExtensionTableViewer::FitsExtensionTableViewer_Load);
			this->Shown += gcnew System::EventHandler(this, &FitsExtensionTableViewer::FitsExtensionTableViewer_Shown);
			this->ResizeBegin += gcnew System::EventHandler(this, &FitsExtensionTableViewer::FitsExtensionTableViewer_ResizeBegin);
			this->ResizeEnd += gcnew System::EventHandler(this, &FitsExtensionTableViewer::FitsExtensionTableViewer_ResizeEnd);
			this->LocationChanged += gcnew System::EventHandler(this, &FitsExtensionTableViewer::FitsExtensionTableViewer_LocationChanged);
			this->SizeChanged += gcnew System::EventHandler(this, &FitsExtensionTableViewer::FitsExtensionTableViewer_SizeChanged);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->ExtensionTableGrid))->EndInit();
			this->menuStrip1->ResumeLayout(false);
			this->menuStrip1->PerformLayout();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion

String^ FILENAME;
String^ EXTENSIONNAME;
array<double, 2>^ DATATABLE;

void OpenFITSImage(String^ FileName);
void PopulateTable(String^ ExtensionName);
bool MenuChooseTable_OPENED = false;
bool headerfront = false;

private: System::Void FitsExtensionTableViewer_ResizeBegin(System::Object^  sender, System::EventArgs^  e);
private: System::Void FitsExtensionTableViewer_ResizeEnd(System::Object^  sender, System::EventArgs^  e);
private: System::Void ViewAllChck_CheckedChanged(System::Object^  sender, System::EventArgs^  e);
private: System::Void MenuChooseTableEntries_DropDownItemClicked(System::Object^  sender, System::Windows::Forms::ToolStripItemClickedEventArgs^  e);
private: System::Void MenuChooseTableEntries_MouseLeave(System::Object^  sender, System::EventArgs^  e) {}
private: System::Void MenuChooseTableEntries_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void MenuChooseTable_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void MenuChooseTable_CheckedChanged(System::Object^  sender, System::EventArgs^  e);
private: System::Void ViewAllChck_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void FitsExtensionTableViewer_FormClosing(System::Object^  sender, System::Windows::Forms::FormClosingEventArgs^  e);
private: System::Void FitsExtensionTableViewer_Load(System::Object^  sender, System::EventArgs^  e);
private: System::Void toolStripMenuItem1_CheckedChanged(System::Object^  sender, System::EventArgs^  e);
private: System::Void PlotMenuItem_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void FileOpenMenu_Click(System::Object^  sender, System::EventArgs^  e);
private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) {
		this->Close();
	}
private: System::Void ViewHeaderMenu_Click(System::Object^  sender, System::EventArgs^  e);


private: System::Void FitsExtensionTableViewer_SizeChanged(System::Object^  sender, System::EventArgs^  e) 
{
	SetReg("JPChart", this->Text + "FitsTableLeft", this->Left);
	SetReg("JPChart", this->Text + "FitsTableTop", this->Top);
	SetReg("JPChart", this->Text + "FitsTableWidth", this->Width);
	SetReg("JPChart", this->Text + "FitsTableHeight", this->Height);
}
private: System::Void FitsExtensionTableViewer_LocationChanged(System::Object^  sender, System::EventArgs^  e) 
{
	SetReg("JPChart", this->Text + "FitsTableLeft", this->Left);
	SetReg("JPChart", this->Text + "FitsTableTop", this->Top);
	SetReg("JPChart", this->Text + "FitsTableWidth", this->Width);
	SetReg("JPChart", this->Text + "FitsTableHeight", this->Height);
}
private: System::Void FitsExtensionTableViewer_Shown(System::Object^  sender, System::EventArgs^  e);
private: System::Void ExtensionTableGrid_Scroll(System::Object^  sender, System::Windows::Forms::ScrollEventArgs^  e);
private: System::Void ExtensionTableGrid_NewRowNeeded(System::Object^  sender, System::Windows::Forms::DataGridViewRowEventArgs^  e);
private: System::Void ExtensionTableGrid_CellValueNeeded(System::Object^  sender, System::Windows::Forms::DataGridViewCellValueEventArgs^  e);
};
}
