#include "stdafx.h"
#include "FitsExtensionTableViewer.h"
#include "JPFITS.h"

void JPFITS::FitsExtensionTableViewer::OpenFITSImage(String^ FileName)
{
	try
	{
		FILENAME = FileName;
		String^ file = FILENAME->Substring(FILENAME->LastIndexOf("\\"));
		this->Text = file->Substring(1);
		array<String^>^ list = JPFITS::FITSBinTable::GetAllExtensionNames(FileName);
		SetReg("JPFITS", "BinTableOpenFilesPath", FileName->Substring(0, FileName->LastIndexOf("\\")));

		for (int i = 0; i < list->Length; i++)
		{
			if (list[i] == "")
				MenuChooseTable->DropDownItems->Add("UNNAMED" + i);
			else
				MenuChooseTable->DropDownItems->Add(list[i]);

			((ToolStripMenuItem^)MenuChooseTable->DropDownItems[MenuChooseTable->DropDownItems->Count - 1])->CheckOnClick = true;
			this->MenuChooseTable->DropDownItems[MenuChooseTable->DropDownItems->Count - 1]->Click += gcnew System::EventHandler(this, &FitsExtensionTableViewer::MenuChooseTable_Click);
			((ToolStripMenuItem^)this->MenuChooseTable->DropDownItems[MenuChooseTable->DropDownItems->Count - 1])->CheckedChanged += gcnew System::EventHandler(this, &FitsExtensionTableViewer::MenuChooseTable_CheckedChanged);
		}

		this->Show();

		if (list->Length == 1)
			PopulateTable(list[0]);
		else
			MenuChooseTable->ShowDropDown();
	}
	catch (Exception^ e)
	{
		MessageBox::Show(e->Data + "	" + e->InnerException + "	" + e->Message + "	" + e->Source + "	" + e->StackTrace + "	" + e->TargetSite);
	}
}

void JPFITS::FitsExtensionTableViewer::PopulateTable(String^ ExtensionName)
{
	try
	{
		this->ExtensionTableGrid->SuspendLayout();

		EXTENSIONNAME = ExtensionName;

		FITSBinTable^ bt = gcnew FITSBinTable(FILENAME, EXTENSIONNAME);
		array<String^>^ labels = bt->TableDataLabelsTTYPE;

		//array<String^>^ labels = JPFITS::FITSBinTable::GetExtensionEntryLabels(FILENAME, EXTENSIONNAME);

		for (int i = 0; i < labels->Length; i++)
			labels[i] = labels[i]->Trim();

		for (int i = 0; i < labels->Length; i++)
		{
			MenuChooseTableEntries->DropDownItems->Add(labels[i]);
			((ToolStripMenuItem^)MenuChooseTableEntries->DropDownItems[MenuChooseTableEntries->DropDownItems->Count - 1])->CheckOnClick = true;
			((ToolStripMenuItem^)MenuChooseTableEntries->DropDownItems[MenuChooseTableEntries->DropDownItems->Count - 1])->Checked = true;
			this->MenuChooseTableEntries->DropDownItems[MenuChooseTableEntries->DropDownItems->Count - 1]->Click += gcnew System::EventHandler(this, &FitsExtensionTableViewer::MenuChooseTableEntries_Click);
			((ToolStripMenuItem^)this->MenuChooseTableEntries->DropDownItems[MenuChooseTableEntries->DropDownItems->Count - 1])->CheckedChanged += gcnew System::EventHandler(this, &FitsExtensionTableViewer::ViewAllChck_CheckedChanged);
		}

		XDrop->Items->AddRange(labels);
		YDrop->Items->AddRange(labels);

		ExtensionTableGrid->ColumnCount = labels->Length;

		for (int i = 0; i < labels->Length; i++)
			ExtensionTableGrid->Columns[i]->HeaderText = labels[i];

		//DATATABLE = FITSImage::ReadBinaryTableExtensionEntries(FILENAME, EXTENSIONNAME, labels);
		int width = 0, height = 0;
		array<double>^ entry = bt->GetTTYPEEntry(labels[0], width, height);
		ExtensionTableGrid->RowCount = height;

		/*#pragma omp parallel for
		for (int i = 0; i < labels->Length; i++)
			for (int j = 0; j < DATATABLE->GetLength(1); j++)
			{
				#pragma omp critical
				{
					ExtensionTableGrid[i, j]->Value = DATATABLE[i, j];
				}
			}*/

		DATATABLE = gcnew array<double, 2>(labels->Length, height);
		for (int i = 0; i < labels->Length; i++)
		{
			entry = bt->GetTTYPEEntry(labels[i], width, height);
			if (width != 1)
				for (int j = 0; j < height; j++)
				{
					ExtensionTableGrid[i, j]->Value = "multi";
					DATATABLE[i, j] = Double::NaN;
				}
			else
				for (int j = 0; j < height; j++)
				{
					ExtensionTableGrid[i, j]->Value = entry[j];
					DATATABLE[i, j] = entry[j];
				}
		}

		this->ExtensionTableGrid->ResumeLayout();
	}
	catch (Exception^ e)
	{
		MessageBox::Show(e->Data + "	" + e->InnerException + "	" + e->Message + "	" + e->Source + "	" + e->StackTrace + "	" + e->TargetSite);
	}
}

void JPFITS::FitsExtensionTableViewer::MenuChooseTable_Click(System::Object^  sender, System::EventArgs^  e)
{
	MenuChooseTable->ShowDropDown();

	for (int i = 0; i < MenuChooseTable->DropDownItems->Count; i++)
	{
		if (MenuChooseTable->DropDownItems[i]->Text != ((ToolStripMenuItem^)sender)->Text)
			((ToolStripMenuItem^)MenuChooseTable->DropDownItems[i])->Checked = false;
	}

	if (((ToolStripMenuItem^)sender)->Checked == false)
		((ToolStripMenuItem^)sender)->Checked = true;
}

void JPFITS::FitsExtensionTableViewer::MenuChooseTable_CheckedChanged(System::Object^  sender, System::EventArgs^  e)
{
	if (((ToolStripMenuItem^)sender)->Checked == true)
	{
		String^ text = ((ToolStripMenuItem^)sender)->Text;
		if (text->Contains("UNNAMED"))
			text = "";

		XDrop->Items->Clear();
		YDrop->Items->Clear();
		int c = MenuChooseTableEntries->DropDownItems->Count;
		for (int i = 2; i < c; i++)
			MenuChooseTableEntries->DropDownItems->RemoveAt(2);

		PopulateTable(text);
	}
}

void JPFITS::FitsExtensionTableViewer::FitsExtensionTableViewer_Load(System::Object^  sender, System::EventArgs^  e)
{
	//try
	{
		this->Left = (int)GetReg("JPChart", /*this->Text + */"FitsTableLeft");
		this->Top = (int)GetReg("JPChart", /*this->Text + */"FitsTableTop");
		this->Width = (int)GetReg("JPChart", /*this->Text + */"FitsTableWidth");
		this->Height = (int)GetReg("JPChart", /*this->Text + */"FitsTableHeight");

		//MessageBox::Show(((int)GetReg("JPChart", this->Text + "FitsTableWidth")).ToString());
		//MessageBox::Show(this->Width.ToString());
	}
	//catch (...) {}
}

void JPFITS::FitsExtensionTableViewer::FitsExtensionTableViewer_Shown(System::Object^  sender, System::EventArgs^  e)
{
	//this->Left = (int)GetReg("JPChart", /*this->Text + */"FitsTableLeft");
	//this->Top = (int)GetReg("JPChart", /*this->Text + */"FitsTableTop");
	//this->Width = (int)GetReg("JPChart", /*this->Text + */"FitsTableWidth");
	//this->Height = (int)GetReg("JPChart", /*this->Text + */"FitsTableHeight");
}

void JPFITS::FitsExtensionTableViewer::FileOpenMenu_Click(System::Object^  sender, System::EventArgs^  e)
{
	OpenFileDialog^ ofd = gcnew OpenFileDialog();
	ofd->Filter = "FITS|*.fits;*.fit;*.fts|All Files|*.*";
	ofd->InitialDirectory = (String^)GetReg("JPFITS", "BinTableOpenFilesPath");

	if (ofd->ShowDialog() == ::DialogResult::Cancel)
		return;

	XDrop->Items->Clear();
	YDrop->Items->Clear();
	int c = MenuChooseTableEntries->DropDownItems->Count;
	for (int i = 2; i < c; i++)
		MenuChooseTableEntries->DropDownItems->RemoveAt(2);
	MenuChooseTable->DropDownItems->Clear();

	OpenFITSImage(ofd->FileName);
}

void JPFITS::FitsExtensionTableViewer::FitsExtensionTableViewer_FormClosing(System::Object^  sender, System::Windows::Forms::FormClosingEventArgs^  e)
{
	SetReg("JPChart", /*this->Text + */"FitsTableLeft", this->Left);
	SetReg("JPChart", /*this->Text + */"FitsTableTop", this->Top);
	SetReg("JPChart", /*this->Text + */"FitsTableWidth", this->Width);
	SetReg("JPChart", /*this->Text + */"FitsTableHeight", this->Height);

	//MessageBox::Show(this->Text + "FitsTableWidth" + this->Width.ToString());
}

void JPFITS::FitsExtensionTableViewer::toolStripMenuItem1_CheckedChanged(System::Object^  sender, System::EventArgs^  e)
{
	if (PlotXChck->Checked == false)
		XDrop->Enabled = false;
	else
		XDrop->Enabled = true;

	PlotEntryMenu->ShowDropDown();
}

void JPFITS::FitsExtensionTableViewer::FitsExtensionTableViewer_ResizeBegin(System::Object^  sender, System::EventArgs^  e)
{
	this->SuspendLayout();
}

void JPFITS::FitsExtensionTableViewer::FitsExtensionTableViewer_ResizeEnd(System::Object^  sender, System::EventArgs^  e)
{
	this->ResumeLayout();
}

void JPFITS::FitsExtensionTableViewer::ViewAllChck_Click(System::Object^  sender, System::EventArgs^  e)
{
	this->ExtensionTableGrid->SuspendLayout();

	if (ViewAllChck->Text == "View None")
	{
		ViewAllChck->Text = "View All";
		for (int i = 2; i < MenuChooseTableEntries->DropDownItems->Count; i++)
		{
			MenuChooseTableEntries->DropDownItems[i]->Tag = "ViewAll";
			((ToolStripMenuItem^)MenuChooseTableEntries->DropDownItems[i])->Checked = false;
			MenuChooseTableEntries->DropDownItems[i]->Tag = nullptr;
		}
		ExtensionTableGrid->Columns->Clear();

		ExtensionTableGrid->ColumnCount = MenuChooseTableEntries->DropDownItems->Count - 2;

		for (int i = 0; i < ExtensionTableGrid->ColumnCount; i++)
			ExtensionTableGrid->Columns[i]->HeaderText = MenuChooseTableEntries->DropDownItems[i + 2]->Text;
	}
	else
	{
		ViewAllChck->Text = "View None";
		for (int i = 2; i < MenuChooseTableEntries->DropDownItems->Count; i++)
		{
			MenuChooseTableEntries->DropDownItems[i]->Tag = "ViewAll";
			((ToolStripMenuItem^)MenuChooseTableEntries->DropDownItems[i])->Checked = true;
			MenuChooseTableEntries->DropDownItems[i]->Tag = nullptr;
		}

		ExtensionTableGrid->RowCount = DATATABLE->GetLength(1);

		for (int i = 0; i < ExtensionTableGrid->ColumnCount; i++)
			for (int j = 0; j < DATATABLE->GetLength(1); j++)
				ExtensionTableGrid[i, j]->Value = DATATABLE[i, j];
	}

	this->ExtensionTableGrid->ResumeLayout();
}

void JPFITS::FitsExtensionTableViewer::MenuChooseTableEntries_Click(System::Object^  sender, System::EventArgs^  e)
{
	MenuChooseTableEntries->ShowDropDown();
}

void JPFITS::FitsExtensionTableViewer::ViewAllChck_CheckedChanged(System::Object^  sender, System::EventArgs^  e)
{
	if (((ToolStripMenuItem^)sender)->Text->Contains("View"))
		return;

	if ((String^)(((ToolStripMenuItem^)sender)->Tag) == "ViewAll")
		return;	

	this->ExtensionTableGrid->SuspendLayout();

	ArrayList^ checked = gcnew ArrayList();
	for (int i = 2; i < MenuChooseTableEntries->DropDownItems->Count; i++)
		if (((ToolStripMenuItem^)MenuChooseTableEntries->DropDownItems[i])->Checked)
			checked->Add(i);

	ExtensionTableGrid->Columns->Clear();
	ExtensionTableGrid->ColumnCount = checked->Count;
	ExtensionTableGrid->RowCount = DATATABLE->GetLength(1);

	for (int i = 0; i < checked->Count; i++)
		ExtensionTableGrid->Columns[i]->HeaderText = MenuChooseTableEntries->DropDownItems[(int)checked[i]]->Text;

	for (int i = 0; i < ExtensionTableGrid->ColumnCount; i++)
		for (int j = 0; j < DATATABLE->GetLength(1); j++)
			ExtensionTableGrid[i, j]->Value = DATATABLE[(int)checked[i] - 2, j];

	this->ExtensionTableGrid->ResumeLayout();
}

void JPFITS::FitsExtensionTableViewer::MenuChooseTableEntries_DropDownItemClicked(System::Object^  sender, System::Windows::Forms::ToolStripItemClickedEventArgs^  e)
{

}

void JPFITS::FitsExtensionTableViewer::PlotMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
	FitsExtensionTablePlotter^ plot = gcnew FitsExtensionTablePlotter();
	array<double>^ x = gcnew array<double>(DATATABLE->GetLength(1));
	array<double>^ y = gcnew array<double>(DATATABLE->GetLength(1));

	int xind = -1;
	if (PlotXChck->Checked)
		xind = XDrop->SelectedIndex;
	int yind = YDrop->SelectedIndex;

	for (int i = 0; i < DATATABLE->GetLength(1); i++)
	{
		if (xind == -1)
			x[i] = i;
		else
			x[i] = (double)ExtensionTableGrid[xind, i]->Value;// DATATABLE[xind, i];
		y[i] = (double)ExtensionTableGrid[yind, i]->Value; //DATATABLE[yind, i];
	}

	String^ xlabel;
	if (xind == -1)
		xlabel = "index";
	else
		xlabel = ExtensionTableGrid->Columns[xind]->HeaderText;
	String^ ylabel = ExtensionTableGrid->Columns[yind]->HeaderText;
	String^ title = ylabel;
	if (xind != -1)
		title += " vs. " + xlabel;

	plot->Chart1->PlotXYData(x, y, title, xlabel, ylabel, ::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line, title->Replace(" ", ""));
	plot->Text = title;
	plot->Show();
}

void JPFITS::FitsExtensionTableViewer::ViewHeaderMenu_Click(System::Object^  sender, System::EventArgs^  e)
{
	HeaderListBox->Items->Clear();

	FITSBinTable^ bt = gcnew FITSBinTable(FILENAME, EXTENSIONNAME);
	//array<String^>^ header = JPFITS::FITSBinTable::GetExtensionHeader(FILENAME, EXTENSIONNAME);
	
	for (int i = 0; i < bt->Header->Length; i++)
		HeaderListBox->Items->Add(bt->Header[i]);

	if (!headerfront)
		HeaderListBox->BringToFront();
	else
		HeaderListBox->SendToBack();
	headerfront = !headerfront;
}
