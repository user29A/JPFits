#pragma once

namespace JPFITS {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for FitsExtensionTablePlotter
	/// </summary>
	public ref class FitsExtensionTablePlotter : public System::Windows::Forms::Form
	{
	public:
		FitsExtensionTablePlotter(void)
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
		~FitsExtensionTablePlotter()
		{
			if (components)
			{
				delete components;
			}
		}
	public: CustomChart::MyChart^  Chart1;
	protected:

	protected:

	protected:
	private: System::ComponentModel::IContainer^  components;

	protected:

	protected:

	protected:

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
			System::Windows::Forms::DataVisualization::Charting::ChartArea^  chartArea2 = (gcnew System::Windows::Forms::DataVisualization::Charting::ChartArea());
			System::Windows::Forms::DataVisualization::Charting::Series^  series2 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			System::Windows::Forms::DataVisualization::Charting::Title^  title2 = (gcnew System::Windows::Forms::DataVisualization::Charting::Title());
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(FitsExtensionTablePlotter::typeid));
			this->Chart1 = (gcnew CustomChart::MyChart());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->Chart1))->BeginInit();
			this->SuspendLayout();
			// 
			// Chart1
			// 
			this->Chart1->BackColor = System::Drawing::Color::DarkGray;
			this->Chart1->BackGradientStyle = System::Windows::Forms::DataVisualization::Charting::GradientStyle::TopBottom;
			this->Chart1->BackSecondaryColor = System::Drawing::Color::White;
			chartArea2->AxisX->Enabled = System::Windows::Forms::DataVisualization::Charting::AxisEnabled::True;
			chartArea2->AxisX->IntervalAutoMode = System::Windows::Forms::DataVisualization::Charting::IntervalAutoMode::VariableCount;
			chartArea2->AxisX->IsStartedFromZero = false;
			chartArea2->AxisX->LabelAutoFitStyle = System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::DecreaseFont;
			chartArea2->AxisX->LabelStyle->Format = L"G6";
			chartArea2->AxisX->LabelStyle->TruncatedLabels = true;
			chartArea2->AxisX->MajorGrid->Enabled = false;
			chartArea2->AxisX->MajorGrid->LineColor = System::Drawing::Color::Gray;
			chartArea2->AxisX->MajorGrid->LineDashStyle = System::Windows::Forms::DataVisualization::Charting::ChartDashStyle::Dash;
			chartArea2->AxisX->MajorTickMark->Interval = 0;
			chartArea2->AxisX->MajorTickMark->IntervalOffset = 0;
			chartArea2->AxisX->MajorTickMark->IntervalOffsetType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Auto;
			chartArea2->AxisX->MajorTickMark->IntervalType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Auto;
			chartArea2->AxisX->MajorTickMark->Size = 2;
			chartArea2->AxisX->MinorTickMark->Enabled = true;
			chartArea2->AxisX->MinorTickMark->IntervalType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Number;
			chartArea2->AxisX->ScaleView->Zoomable = false;
			chartArea2->AxisX->ScrollBar->Enabled = false;
			chartArea2->AxisX->Title = L"X Axis";
			chartArea2->AxisX2->Enabled = System::Windows::Forms::DataVisualization::Charting::AxisEnabled::False;
			chartArea2->AxisX2->IsStartedFromZero = false;
			chartArea2->AxisX2->MajorTickMark->Interval = 0;
			chartArea2->AxisX2->ScrollBar->Enabled = false;
			chartArea2->AxisY->Enabled = System::Windows::Forms::DataVisualization::Charting::AxisEnabled::True;
			chartArea2->AxisY->IntervalAutoMode = System::Windows::Forms::DataVisualization::Charting::IntervalAutoMode::VariableCount;
			chartArea2->AxisY->IsStartedFromZero = false;
			chartArea2->AxisY->LabelAutoFitStyle = System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::DecreaseFont;
			chartArea2->AxisY->LabelStyle->Format = L"G5";
			chartArea2->AxisY->LabelStyle->TruncatedLabels = true;
			chartArea2->AxisY->MajorGrid->Enabled = false;
			chartArea2->AxisY->MajorGrid->LineColor = System::Drawing::Color::Gray;
			chartArea2->AxisY->MajorGrid->LineDashStyle = System::Windows::Forms::DataVisualization::Charting::ChartDashStyle::Dash;
			chartArea2->AxisY->MajorTickMark->Interval = 0;
			chartArea2->AxisY->MajorTickMark->IntervalOffset = 0;
			chartArea2->AxisY->MajorTickMark->IntervalOffsetType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Auto;
			chartArea2->AxisY->MajorTickMark->IntervalType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Auto;
			chartArea2->AxisY->MajorTickMark->Size = 2;
			chartArea2->AxisY->MinorTickMark->Enabled = true;
			chartArea2->AxisY->ScaleView->Zoomable = false;
			chartArea2->AxisY->ScrollBar->Enabled = false;
			chartArea2->AxisY->Title = L"Y Axis";
			chartArea2->AxisY2->Enabled = System::Windows::Forms::DataVisualization::Charting::AxisEnabled::False;
			chartArea2->AxisY2->IsStartedFromZero = false;
			chartArea2->AxisY2->ScrollBar->Enabled = false;
			chartArea2->BackColor = System::Drawing::Color::White;
			chartArea2->BorderDashStyle = System::Windows::Forms::DataVisualization::Charting::ChartDashStyle::Solid;
			chartArea2->CursorX->AutoScroll = false;
			chartArea2->CursorX->Interval = 0.001;
			chartArea2->CursorX->IntervalType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Number;
			chartArea2->CursorX->SelectionColor = System::Drawing::Color::LightGreen;
			chartArea2->CursorY->AutoScroll = false;
			chartArea2->CursorY->Interval = 0.001;
			chartArea2->CursorY->IntervalType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Number;
			chartArea2->CursorY->SelectionColor = System::Drawing::Color::LightGreen;
			chartArea2->Name = L"ChartArea1";
			chartArea2->Position->Auto = false;
			chartArea2->Position->Height = 85;
			chartArea2->Position->Width = 85;
			chartArea2->Position->X = 2;
			chartArea2->Position->Y = 12.64803F;
			//this->Chart1->ChartAreas->Add(chartArea2);
			this->Chart1->Dock = System::Windows::Forms::DockStyle::Fill;
			this->Chart1->IsSoftShadows = false;
			this->Chart1->Location = System::Drawing::Point(0, 0);
			this->Chart1->Name = L"Chart1";
			series2->BorderColor = System::Drawing::Color::White;
			series2->ChartArea = L"ChartArea1";
			series2->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			series2->Color = System::Drawing::Color::Teal;
			series2->CustomProperties = L"DrawingStyle=Emboss, PointWidth=1";
			series2->MarkerBorderWidth = 0;
			series2->MarkerColor = System::Drawing::Color::Teal;
			series2->Name = L"Series1";
			//this->Chart1->Series->Add(series2);
			this->Chart1->Size = System::Drawing::Size(564, 471);
			this->Chart1->TabIndex = 0;
			this->Chart1->Text = L"myChart1";
			title2->DockedToChartArea = L"ChartArea1";
			title2->Font = (gcnew System::Drawing::Font(L"Times New Roman", 16, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			title2->IsDockedInsideChartArea = false;
			title2->Name = L"Title1";
			title2->Position->Auto = false;
			title2->Position->Height = 6.828034F;
			title2->Position->Width = 79.9F;
			title2->Position->X = 10;
			title2->Position->Y = 5.819999F;
			title2->Text = L"Title";
			//this->Chart1->Titles->Add(title2);
			// 
			// FitsExtensionTablePlotter
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(564, 471);
			this->Controls->Add(this->Chart1);
			this->Icon = (cli::safe_cast<System::Drawing::Icon^>(resources->GetObject(L"$this.Icon")));
			this->Name = L"FitsExtensionTablePlotter";
			this->Text = L"FitsExtensionTablePlotter";
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->Chart1))->EndInit();
			this->ResumeLayout(false);

		}
#pragma endregion
	};
}
