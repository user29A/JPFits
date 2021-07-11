#pragma once

namespace JPFits {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for FITSHeaderViewer
	/// </summary>
	public ref class FITSHeaderViewer : public System::Windows::Forms::Form
	{
	public:
		FITSHeaderViewer(void)
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
		~FITSHeaderViewer()
		{
			if (components)
			{
				delete components;
			}
		}
		public: System::Windows::Forms::ListBox^ HeaderKeysListBox;
		protected:

		protected:

		protected:

		protected:
		private: System::Windows::Forms::MenuStrip^ MenuStrip;
		private: System::Windows::Forms::ToolStripMenuItem^ EditMenuItem;

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
			this->HeaderKeysListBox = (gcnew System::Windows::Forms::ListBox());
			this->MenuStrip = (gcnew System::Windows::Forms::MenuStrip());
			this->EditMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->MenuStrip->SuspendLayout();
			this->SuspendLayout();
			// 
			// HeaderKeysListBox
			// 
			this->HeaderKeysListBox->BackColor = System::Drawing::Color::Gainsboro;
			this->HeaderKeysListBox->Dock = System::Windows::Forms::DockStyle::Fill;
			this->HeaderKeysListBox->Font = (gcnew System::Drawing::Font(L"Courier New", 9, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->HeaderKeysListBox->FormattingEnabled = true;
			this->HeaderKeysListBox->ItemHeight = 15;
			this->HeaderKeysListBox->Location = System::Drawing::Point(0, 24);
			this->HeaderKeysListBox->Name = L"HeaderKeysListBox";
			this->HeaderKeysListBox->SelectionMode = System::Windows::Forms::SelectionMode::MultiExtended;
			this->HeaderKeysListBox->Size = System::Drawing::Size(654, 672);
			this->HeaderKeysListBox->TabIndex = 0;
			// 
			// MenuStrip
			// 
			this->MenuStrip->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) { this->EditMenuItem });
			this->MenuStrip->Location = System::Drawing::Point(0, 0);
			this->MenuStrip->Name = L"MenuStrip";
			this->MenuStrip->Size = System::Drawing::Size(654, 24);
			this->MenuStrip->TabIndex = 1;
			this->MenuStrip->Text = L"menuStrip1";
			// 
			// EditMenuItem
			// 
			this->EditMenuItem->Name = L"EditMenuItem";
			this->EditMenuItem->Size = System::Drawing::Size(39, 20);
			this->EditMenuItem->Text = L"Edit";
			// 
			// FITSHeaderViewer
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(654, 696);
			this->Controls->Add(this->HeaderKeysListBox);
			this->Controls->Add(this->MenuStrip);
			this->MainMenuStrip = this->MenuStrip;
			this->Name = L"FITSHeaderViewer";
			this->ShowIcon = false;
			this->Text = L"Header";
			this->MenuStrip->ResumeLayout(false);
			this->MenuStrip->PerformLayout();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	};
}
