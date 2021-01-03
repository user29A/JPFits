#include "stdafx.h"
#include "JPFITS.h"

using namespace JPFITS;

JPFITS::WCS_AutoSolver::~WCS_AutoSolver() {}

JPFITS::WCS_AutoSolver::WCS_AutoSolver(String^ WCS_type, array<JPMath::PointD^>^ pixels, bool zero_based_pixels, int pixels_tolerance_radius, int image_width, int image_height, array<JPMath::PointD^>^ wcscoords)
{
	this->BGWRKR = gcnew BackgroundWorker();
	this->BGWRKR->WorkerReportsProgress = true;
	this->BGWRKR->WorkerSupportsCancellation = true;
	this->BGWRKR->DoWork += gcnew System::ComponentModel::DoWorkEventHandler(this, &WCS_AutoSolver::BGWRKR_DoWork);
	this->BGWRKR->ProgressChanged += gcnew System::ComponentModel::ProgressChangedEventHandler(this, &WCS_AutoSolver::BGWRKR_ProgressChanged);
	this->BGWRKR->RunWorkerCompleted += gcnew System::ComponentModel::RunWorkerCompletedEventHandler(this, &WCS_AutoSolver::BGWRKR_RunWorkerCompleted);

	WCS_TYPE = WCS_type;
	PIX_PTS = pixels;
	ZERO_BASED_PIX = zero_based_pixels;
	PSE_KERNEL_RADIUS = pixels_tolerance_radius;
	IMAGE_WIDTH = image_width;
	IMAGE_HEIGHT = image_height;
	CAT_PTS = wcscoords;
	CANCELLED = false;
	PROGRESS = 0;
	DO_PSE = false;
	N_POINTS = pixels->Length;//????????????????????????????????????????????????
	SOLVED = false;
	PSE = nullptr;
	WCS = gcnew JPFITS::WorldCoordinateSolution();
}

JPFITS::WCS_AutoSolver::WCS_AutoSolver(String^ WCS_type, int Number_of_Points, JPFITS::FITSImage^ Fits_Img, array<bool, 2>^ Image_ROI, double Image_Saturation, bool auto_background, int PSE_kernel_radius, int PSE_separation_radius, String^ Fits_Catalogue_BinTable_File, String^ Catalogue_Extension_Name, String^ Catalogue_CVAL1_Name, String^ Catalogue_CVAL2_Name, String^ Catalogue_Magnitude_Name, bool Refine)
{
	this->BGWRKR = gcnew BackgroundWorker();
	this->BGWRKR->WorkerReportsProgress = true;
	this->BGWRKR->WorkerSupportsCancellation = true;
	this->BGWRKR->DoWork += gcnew System::ComponentModel::DoWorkEventHandler(this, &WCS_AutoSolver::BGWRKR_DoWork);
	this->BGWRKR->ProgressChanged += gcnew System::ComponentModel::ProgressChangedEventHandler(this, &WCS_AutoSolver::BGWRKR_ProgressChanged);
	this->BGWRKR->RunWorkerCompleted += gcnew System::ComponentModel::RunWorkerCompletedEventHandler(this, &WCS_AutoSolver::BGWRKR_RunWorkerCompleted);

	WCS_TYPE = WCS_type;
	N_POINTS = Number_of_Points;
	FITS_IMG = Fits_Img;
	IMAGE_ROI = Image_ROI;
	PIX_SAT = Image_Saturation;
	AUTO_BACKGROUND = auto_background;
	PSE_KERNEL_RADIUS = PSE_kernel_radius;
	PSE_SEP_RADIUS = PSE_separation_radius;
	CAT_FILENAME = Fits_Catalogue_BinTable_File;
	CAT_EXTNAME = Catalogue_Extension_Name;
	CAT_CVAL1NAME = Catalogue_CVAL1_Name;
	CAT_CVAL2NAME = Catalogue_CVAL2_Name;
	CAT_MAGNAME = Catalogue_Magnitude_Name;
	REFINE = Refine;
	DO_PSE = true;
	CANCELLED = false;
	PROGRESS = 0;
	SOLVED = false;
	PSE = gcnew JPFITS::SourceExtractor();
	WCS = gcnew JPFITS::WorldCoordinateSolution();
}

void JPFITS::WCS_AutoSolver::SolveAsync(double scale_init, double scale_lb, double scale_ub, double rotation_init, double rotation_lb, double rotation_ub, double vertex_tolerance, int N_matches_stop, int Percentage_matches_stop, bool condition_arrays, bool show_report_form)
{
	SCALE_INIT = scale_init;
	SCALE_LB = scale_lb;
	SCALE_UB = scale_ub;
	ROTATION_INIT = rotation_init;
	ROTATION_LB = rotation_lb;
	ROTATION_UB = rotation_ub;
	WCS_VERTEX_TOL = vertex_tolerance;
	N_MATCHES_STOP = N_matches_stop;
	PERC_MATCHES_STOP = Percentage_matches_stop;
	CONDITION_ARRAYS = condition_arrays;
	SHOW_REPORT_FORM = show_report_form;
	SOLVING = true;

	if (SHOW_REPORT_FORM)
	{
		SHOW_REPORT_FORM = false;
		WCSARF = gcnew JPFITS::WCS_AutoSolver::WCSAutoSolverReportingForm(this);
		WCSARF->WCSAutoReportingTimer->Enabled = true;
		BGWRKR->RunWorkerAsync();
		::DialogResult res = WCSARF->ShowDialog();
		if (res == ::DialogResult::Cancel)
			CANCELLED = true;
	}
	else
		BGWRKR->RunWorkerAsync();
}

void JPFITS::WCS_AutoSolver::WCSAutoSolverReportingForm::WCSAutoReportingTimer_Tick(System::Object^ sender, System::EventArgs^ e)
{
	this->MsgTxt->Text = this->WCSAS->Status_Log;
	this->MsgTxt->SelectionStart = this->MsgTxt->TextLength;
	this->MsgTxt->ScrollToCaret();

	if (this->DialogResult == ::DialogResult::Cancel)
		this->WCSAS->Cancelled = true;

	if (!this->WCSAS->Solving)
		this->WCSAutoReportingTimer->Enabled = false;
}

void JPFITS::WCS_AutoSolver::BGWRKR_DoWork(System::Object^ sender, System::ComponentModel::DoWorkEventArgs^ e)
{
	if (DO_PSE)
	{
		//get catalogue RA, Dec, and mag's
		BGWRKR->ReportProgress(0, "Reading the Catalogue FITS binary tables...");
		CAT_CVAL1s = JPFITS::FITSBinTable::GetExtensionEntry(CAT_FILENAME, CAT_EXTNAME, CAT_CVAL1NAME);
		CAT_CVAL2s = JPFITS::FITSBinTable::GetExtensionEntry(CAT_FILENAME, CAT_EXTNAME, CAT_CVAL2NAME);
		CAT_MAGs = JPFITS::FITSBinTable::GetExtensionEntry(CAT_FILENAME, CAT_EXTNAME, CAT_MAGNAME);

		//need to check mag for NaN's and re-form ra dec mag
		BGWRKR->ReportProgress(0, "Formatting the Catalogue FITS binary tables...");
		ArrayList^ ralist = gcnew ArrayList(CAT_CVAL1s->Length);
		ArrayList^ declist = gcnew ArrayList(CAT_CVAL1s->Length);
		ArrayList^ maglist = gcnew ArrayList(CAT_CVAL1s->Length);
		#pragma omp parallel for
		for (int i = 0; i < CAT_CVAL1s->Length; i++)
		{
			if (Double::IsNaN(CAT_MAGs[i]))
				continue;

			#pragma omp critical
			{
				ralist->Add(CAT_CVAL1s[i]);
				declist->Add(CAT_CVAL2s[i]);
				maglist->Add(CAT_MAGs[i]);
			}
		}
		#pragma omp parallel for
		for (int i = 0; i < ralist->Count; i++)
		{
			CAT_CVAL1s[i] = Convert::ToDouble(ralist[i]);
			CAT_CVAL2s[i] = Convert::ToDouble(declist[i]);
			CAT_MAGs[i] = Convert::ToDouble(maglist[i]);
		}
		Array::Resize(CAT_CVAL1s, ralist->Count);
		Array::Resize(CAT_CVAL2s, ralist->Count);
		Array::Resize(CAT_MAGs, ralist->Count);

		//sort the catalogue list by magnitude
		BGWRKR->ReportProgress(0, "Sorting the Catalogue FITS binary tables...");
		array<double>^ keysref = gcnew array<double>(CAT_MAGs->Length);
		Array::Copy(CAT_MAGs, keysref, CAT_MAGs->Length);
		Array::Sort(CAT_MAGs, CAT_CVAL1s);
		Array::Copy(keysref, CAT_MAGs, CAT_MAGs->Length);
		Array::Sort(CAT_MAGs, CAT_CVAL2s);

		//get the brightest few catlaogue points
		BGWRKR->ReportProgress(0, "Making Catalogue points...");
		CAT_PTS = gcnew array<JPMath::PointD^>(N_POINTS);
		for (int i = 0; i < CAT_PTS->Length; i++)
			CAT_PTS[i] = gcnew JPMath::PointD(CAT_CVAL1s[i], CAT_CVAL2s[i], CAT_MAGs[i]);

		//process the fits image
		BGWRKR->ReportProgress(0, Environment::NewLine + "Processing the FITS Image...");
		IMAGE_WIDTH = FITS_IMG->Width;
		IMAGE_HEIGHT = FITS_IMG->Height;

		BGWRKR->ReportProgress(0, Environment::NewLine + "Seraching '" + FITS_IMG->FileName + "' for " + N_POINTS + " point sources...");
		double immax = FITS_IMG->Max;
		double pixthresh = immax / 256;
		double div = 2;
		double amp = pixthresh;
		int PSEiters = 0;
		int maxPSEiters = 15;

		while (PSE->N_Sources < N_POINTS || PSE->N_Sources > N_POINTS + 1)
		{
			if (CANCELLED)
				return;

			PSEiters++;
			if (PSEiters > maxPSEiters)
				break;

			if (PSE->N_SaturatedSources >= N_POINTS || PSE->N_Sources >= N_POINTS)
				break;

			PSE->Extract_Sources(FITS_IMG->Image, PIX_SAT, pixthresh, Double::MaxValue, 0, Double::MaxValue, false, PSE_KERNEL_RADIUS, PSE_SEP_RADIUS, AUTO_BACKGROUND, "", IMAGE_ROI, false);

			if (PSE->N_Sources < N_POINTS)
				pixthresh -= amp / div;
			if (PSE->N_Sources > N_POINTS + 1)
				pixthresh += amp / div;
			div *= 2;

			BGWRKR->ReportProgress(0, "Found " + PSE->N_Sources + " point sources on iteration " + PSEiters);
		}
		if (PSE->N_Sources > N_POINTS)
			PSE->ClipToNBrightest(N_POINTS);

		BGWRKR->ReportProgress(0, Environment::NewLine + "Stopped searching on iteration " + PSEiters + " with " + PSE->N_Sources + " point sources");

		//turn the PSE results into points
		BGWRKR->ReportProgress(0, Environment::NewLine + "Making point source points...");
		PIX_PTS = gcnew array<JPMath::PointD^>(PSE->N_Sources);
		for (int i = 0; i < PIX_PTS->Length; i++)
			PIX_PTS[i] = gcnew JPMath::PointD(IMAGE_WIDTH - 1 - PSE->Centroids_X[i], IMAGE_HEIGHT - 1 - PSE->Centroids_Y[i], PSE->Centroids_Volume[i]);

		//now run the auto solver by returning from run worker compeleted
		return;
	}

	if (CANCELLED)
		return;

	//convert field and rotation parameters to radians here:
	SCALE_INIT *= Math::PI / 180 / 3600;
	SCALE_LB *= Math::PI / 180 / 3600;
	SCALE_UB *= Math::PI / 180 / 3600;
	ROTATION_INIT *= Math::PI / 180;
	ROTATION_LB *= Math::PI / 180;
	ROTATION_UB *= Math::PI / 180;
	WCS_VERTEX_TOL *= Math::PI / 180;

	//get pixel initial conditions and bounadries, and CRVAL values
	BGWRKR->ReportProgress(0, "Determining pixel initial conditions and boundaries, and sky coordinate references...");
	double crpix1_init = 0, crpix2_init = 0, crpix1_lb = Double::MaxValue, crpix1_ub = Double::MinValue, crpix2_lb = Double::MaxValue, crpix2_ub = Double::MinValue, crval1 = 0, crval2 = 0;
	for (int i = 0; i < PIX_PTS->Length; i++)
	{
		PIX_PTS[i] = gcnew JPMath::PointD(IMAGE_WIDTH - 1 - PSE->Centroids_X[i], IMAGE_HEIGHT - 1 - PSE->Centroids_Y[i], PSE->Centroids_Volume[i]);
		crpix1_init += PIX_PTS[i]->X;
		crpix2_init += PIX_PTS[i]->Y;
		if (crpix1_ub < PIX_PTS[i]->X)
			crpix1_ub = PIX_PTS[i]->X;
		if (crpix1_lb > PIX_PTS[i]->X)
			crpix1_lb = PIX_PTS[i]->X;
		if (crpix2_ub < PIX_PTS[i]->Y)
			crpix2_ub = PIX_PTS[i]->Y;
		if (crpix2_lb > PIX_PTS[i]->Y)
			crpix2_lb = PIX_PTS[i]->Y;
	}
	crpix1_init /= (double)PIX_PTS->Length;//the reference value initial guesses can be the means
	crpix2_init /= (double)PIX_PTS->Length;

	for (int i = 0; i < CAT_PTS->Length; i++)
	{
		CAT_PTS[i] = gcnew JPMath::PointD(CAT_CVAL1s[i], CAT_CVAL2s[i], CAT_MAGs[i]);
		crval1 += CAT_PTS[i]->X;
		crval2 += CAT_PTS[i]->Y;
	}
	crval1 /= (double)CAT_PTS->Length;//the reference value can be the mean
	crval2 /= (double)CAT_PTS->Length;//the reference value can be the mean

	if (CANCELLED)
		return;

	//make PSE triangles
	BGWRKR->ReportProgress(0, "Making point source triangles...");
	int nPSEtriangles = PIX_PTS->Length * (PIX_PTS->Length - 1) * (PIX_PTS->Length - 2) / 6;
	array<JPMath::Triangle^>^ PSEtriangles = gcnew array<JPMath::Triangle^>(nPSEtriangles);
	int c = 0;
	for (int i = 0; i < PIX_PTS->Length - 2; i++)
		for (int j = i + 1; j < PIX_PTS->Length - 1; j++)
			for (int k = j + 1; k < PIX_PTS->Length; k++)
			{
				PSEtriangles[c] = gcnew JPMath::Triangle(PIX_PTS[i], PIX_PTS[j], PIX_PTS[k]/*, true*/);
				c++;
			}

	//convert the catalogue points to intermediate points
	BGWRKR->ReportProgress(0, "Making catalogue intermediate points...");
	array<JPMath::PointD^>^ CATpts_intrmdt = gcnew array<JPMath::PointD^>(N_POINTS);
	double a0 = crval1 * Math::PI / 180, d0 = crval2 * Math::PI / 180;
	for (int i = 0; i < CATpts_intrmdt->Length; i++)
	{
		double a = CAT_PTS[i]->X * Math::PI / 180;//radians
		double d = CAT_PTS[i]->Y * Math::PI / 180;//radians

		//for tangent plane Gnomic
		double Xint = Math::Cos(d) * Math::Sin(a - a0) / (Math::Cos(d0) * Math::Cos(d) * Math::Cos(a - a0) + Math::Sin(d0) * Math::Sin(d));
		double Yint = (Math::Cos(d0) * Math::Sin(d) - Math::Cos(d) * Math::Sin(d0) * Math::Cos(a - a0)) / (Math::Sin(d0) * Math::Sin(d) + Math::Cos(d0) * Math::Cos(d) * Math::Cos(a - a0));

		CATpts_intrmdt[i] = gcnew JPMath::PointD(Xint, Yint, CAT_PTS[i]->Value);
	}

	if (CANCELLED)
		return;

	//make intermediate coordinate triangles
	BGWRKR->ReportProgress(0, "Making catalogue intermediate triangles...");
	int nCATtriangles = CATpts_intrmdt->Length * (CATpts_intrmdt->Length - 1) * (CATpts_intrmdt->Length - 2) / 6;
	array<JPMath::Triangle^>^ CATtriangles_intrmdt = gcnew array<JPMath::Triangle^>(nCATtriangles);
	c = 0;
	for (int i = 0; i < CATpts_intrmdt->Length - 2; i++)
		for (int j = i + 1; j < CATpts_intrmdt->Length - 1; j++)
			for (int k = j + 1; k < CATpts_intrmdt->Length; k++)
			{
				CATtriangles_intrmdt[c] = gcnew JPMath::Triangle(CATpts_intrmdt[i], CATpts_intrmdt[j], CATpts_intrmdt[k]/*, true*/);
				c++;
			}

	if (CANCELLED)
		return;

	if (DO_PARALLEL)
		if (CONDITION_ARRAYS)
		{
			BGWRKR->ReportProgress(0, Environment::NewLine + "Conditioning the triangle arrays...");
			PSEtriangles = JPFITS::WCS_AutoSolver::ConditionTriangleArrayBrightnessThreads(PSEtriangles, omp_get_max_threads(), false);
			CATtriangles_intrmdt = JPFITS::WCS_AutoSolver::ConditionTriangleArrayBrightnessThreads(CATtriangles_intrmdt, 1, true);
		}

	if (CANCELLED)
		return;

	//for each PSE triangle, fit it to a CAT intermediate triangle, and then check if this fit satisfies the other CAT points to the PSE points
	//rotation transformation p[0] = scale; p[1] = phi (radians); p[2] = x-axis coordinate reference; p[3] = y-axis coordinate reference;
	array<double>^ plb = gcnew array<double>(4) { SCALE_LB, ROTATION_LB, crpix1_lb, crpix2_lb };
	array<double>^ pub = gcnew array<double>(4) { SCALE_UB, ROTATION_UB, crpix1_ub, crpix2_ub };
	array<double>^ psc = gcnew array<double>(4) { SCALE_INIT, 1, Math::Abs(crpix1_init), Math::Abs(crpix2_init) };
	double kern_diam = double(2 * PSE_KERNEL_RADIUS) + 1;
	double p00, p01, p02, p03;
	int total_pt_matches = 0;
	DATE = DateTime::Now;
	TimeSpan ts;
	int prog = 0, threadnum = 0;
	unsigned __int64 ncompares = 0, nfalse_sols = 0;

	bool compare_fieldvectors = ROTATION_LB != -Math::PI && ROTATION_UB != Math::PI;

	BGWRKR->ReportProgress(0, Environment::NewLine + "Starting search for matching triangles among " + (PSEtriangles->Length * CATtriangles_intrmdt->Length).ToString("0.00##e+00") + " possible comparisons...");

	#pragma omp parallel for reduction(+:ncompares, nfalse_sols) if(DO_PARALLEL)
	for (int i = 0; i < PSEtriangles->Length; i++)
	{
		if (SOLVED)
			break;
		if (CANCELLED)
			break;

		if (i < PSEtriangles->Length / omp_get_num_threads())
			if (omp_get_num_threads() * i * 100 / PSEtriangles->Length > prog)
			{
				prog++;
				BGWRKR->ReportProgress(prog);
			}

		//create these here so that each thread when parallel has own copy
		array<double>^ xpix_triplet = gcnew array<double>(3);
		array<double>^ ypix_triplet = gcnew array<double>(3);
		xpix_triplet[0] = PSEtriangles[i]->Vertex[0]->X;
		ypix_triplet[0] = PSEtriangles[i]->Vertex[0]->Y;
		xpix_triplet[1] = PSEtriangles[i]->Vertex[1]->X;
		ypix_triplet[1] = PSEtriangles[i]->Vertex[1]->Y;
		xpix_triplet[2] = PSEtriangles[i]->Vertex[2]->X;
		ypix_triplet[2] = PSEtriangles[i]->Vertex[2]->Y;

		//create these here so that each thread when parallel has own copy
		array<double>^ Xintrmdt_triplet = gcnew array<double>(3);
		array<double>^ Yintrmdt_triplet = gcnew array<double>(3);
		array<double>^ P0 = gcnew array<double>(4);
		//double minlength0 = SCALE_LB * (PSEtriangles[i]->SideLength[0] - kern_diam);//redundant as per below
		//double maxlength0 = SCALE_UB * (PSEtriangles[i]->SideLength[0] + kern_diam);//redundant as per below
		//double minlength1 = SCALE_LB * (PSEtriangles[i]->SideLength[1] - kern_diam);//redundant as per below
		//double maxlength1 = SCALE_UB * (PSEtriangles[i]->SideLength[1] + kern_diam);//redundant as per below
		double minlength2 = SCALE_LB * (PSEtriangles[i]->SideLength[2] - kern_diam);
		double maxlength2 = SCALE_UB * (PSEtriangles[i]->SideLength[2] + kern_diam);

		for (int j = 0; j < CATtriangles_intrmdt->Length; j++)
		{
			if (SOLVED)
				break;
			if (CANCELLED)
				break;

			ncompares++;

			if (Math::Abs(PSEtriangles[i]->VertexAngle[0] - CATtriangles_intrmdt[j]->VertexAngle[0]) > WCS_VERTEX_TOL)
				continue;
			if (Math::Abs(PSEtriangles[i]->VertexAngle[1] - CATtriangles_intrmdt[j]->VertexAngle[1]) > WCS_VERTEX_TOL)
				continue;
			/*if (Math::Abs(PSEtriangles[i]->VertexAngle[2] - CATtriangles_intrmdt[j]->VertexAngle[2]) > WCS_VERTEX_TOL)
				continue;*/

				//these are unneccessary because they're redundant

			/*if (CATtriangles_intrmdt[j]->SideLength[0] <  minlength0 || CATtriangles_intrmdt[j]->SideLength[0] > maxlength0)
				continue;*/
			/*if (CATtriangles_intrmdt[j]->SideLength[1] < minlength1 || CATtriangles_intrmdt[j]->SideLength[1] > maxlength1)
				continue;*/
			if (CATtriangles_intrmdt[j]->SideLength[2] < minlength2 || CATtriangles_intrmdt[j]->SideLength[2] > maxlength2)
				continue;

			if (compare_fieldvectors)
			{
				double theta = Math::Atan2(PSEtriangles[i]->FieldVector->X * CATtriangles_intrmdt[j]->FieldVector->Y - PSEtriangles[i]->FieldVector->Y * CATtriangles_intrmdt[j]->FieldVector->X, PSEtriangles[i]->FieldVector->X * CATtriangles_intrmdt[j]->FieldVector->X + PSEtriangles[i]->FieldVector->Y * CATtriangles_intrmdt[j]->FieldVector->Y);
				//double theta = CATtriangles_intrmdt[j]->FieldVectorRadAngle - PSEtriangles[i]->FieldVectorRadAngle;
				
				if (theta > ROTATION_UB || theta < ROTATION_LB)
					continue;
			}

			Xintrmdt_triplet[0] = CATtriangles_intrmdt[j]->Vertex[0]->X;
			Yintrmdt_triplet[0] = CATtriangles_intrmdt[j]->Vertex[0]->Y;
			Xintrmdt_triplet[1] = CATtriangles_intrmdt[j]->Vertex[1]->X;
			Yintrmdt_triplet[1] = CATtriangles_intrmdt[j]->Vertex[1]->Y;
			Xintrmdt_triplet[2] = CATtriangles_intrmdt[j]->Vertex[2]->X;
			Yintrmdt_triplet[2] = CATtriangles_intrmdt[j]->Vertex[2]->Y;

			//reset P0 for j'th iteration
			P0[0] = SCALE_INIT;
			P0[1] = ROTATION_INIT;
			P0[2] = crpix1_init;
			P0[3] = crpix2_init;

			//try a fit
			JPMath::Fit_WCSTransform2d(Xintrmdt_triplet, Yintrmdt_triplet, xpix_triplet, ypix_triplet, P0, plb, pub, psc);

			int N_pt_matches = 0;
			for (int k = 0; k < 3; k++)
			{
				int x = (int)Math::Round((double)IMAGE_WIDTH - 1 - (1 / P0[0] * (Math::Cos(-P0[1]) * Xintrmdt_triplet[k] - Math::Sin(-P0[1]) * Yintrmdt_triplet[k]) + P0[2]));
				int y = (int)Math::Round((double)IMAGE_HEIGHT - 1 - (1 / P0[0] * (Math::Sin(-P0[1]) * Xintrmdt_triplet[k] + Math::Cos(-P0[1]) * Yintrmdt_triplet[k]) + P0[3]));

				if (x > 0 && y > 0 && x < IMAGE_WIDTH && y < IMAGE_HEIGHT && PSE->SourceIndexMap[x, y] == PSE->SourceIndexMap[IMAGE_WIDTH - 1 - (int)Math::Round(xpix_triplet[k]), IMAGE_HEIGHT - 1 - (int)Math::Round(ypix_triplet[k])])
					N_pt_matches++;
			}

			if (N_pt_matches != 3)//not a possible solution
				continue;

			if (!SOLVED)
			{
				#pragma omp critical
				{
					if (!SOLVED)
					{
						//need to check if the other CAT points match the PSE pts
						N_pt_matches = 0;
						for (int k = 0; k < CATpts_intrmdt->Length; k++)
						{
							double X_int = CATpts_intrmdt[k]->X;
							double Y_int = CATpts_intrmdt[k]->Y;

							int x = (int)Math::Round((double)IMAGE_WIDTH - 1 - (1 / P0[0] * (Math::Cos(-P0[1]) * X_int - Math::Sin(-P0[1]) * Y_int) + P0[2]));
							int y = (int)Math::Round((double)IMAGE_HEIGHT - 1 - (1 / P0[0] * (Math::Sin(-P0[1]) * X_int + Math::Cos(-P0[1]) * Y_int) + P0[3]));

							if (x > 0 && y > 0 && x < IMAGE_WIDTH && y < IMAGE_HEIGHT && PSE->SourceBooleanMap[x, y])
								N_pt_matches++;
						}

						if (N_pt_matches >= N_MATCHES_STOP || N_pt_matches * 100 / CATpts_intrmdt->Length >= PERC_MATCHES_STOP)
						{
							SOLVED = true;
							ts = DateTime::Now - DATE;
							total_pt_matches = N_pt_matches;
							p00 = P0[0];
							p01 = P0[1];
							p02 = P0[2];
							p03 = P0[3];
							threadnum = omp_get_thread_num();
						}
						else
							nfalse_sols++;
					}
				}
			}
		}
	}

	if (!SOLVED)
		BGWRKR->ReportProgress(0, "No solution...");
	if (CANCELLED)
		BGWRKR->ReportProgress(0, "Cancelled...");
	if (!SOLVED || CANCELLED)
		return;

	//solve on all matching coordinates
	BGWRKR->ReportProgress(0, "Found an initial solution given the stopping criteria of either " + N_MATCHES_STOP + " matching points or " + PERC_MATCHES_STOP + "% matching points" + Environment::NewLine);
	BGWRKR->ReportProgress(0, "Field Scale: " + Math::Round(p00 * 180 / Math::PI * 3600, 4));
	BGWRKR->ReportProgress(0, "Field Rotation: " + Math::Round(p01 * 180 / 3.14159265, 3));
	BGWRKR->ReportProgress(0, "N Pt. Matches: " + total_pt_matches + " (" + (total_pt_matches * 100 / CATpts_intrmdt->Length).ToString("00.0") + "%)");
	BGWRKR->ReportProgress(0, "N Comparisons: " + ncompares.ToString("0.00e00") + " (" + Math::Round(double(ncompares * 100) / double(PSEtriangles->Length) / double(CATtriangles_intrmdt->Length), 1) + "%)");
	BGWRKR->ReportProgress(0, "N False Postives: " + nfalse_sols);
	BGWRKR->ReportProgress(0, "Thread: " + threadnum);
	BGWRKR->ReportProgress(0, "Completed in: " + ts.Minutes.ToString() + "m" + (double(ts.Seconds) + (double)ts.Milliseconds / 1000).ToString() + "s");
	BGWRKR->ReportProgress(0, "Comparisons per Second: " + (ncompares / ts.TotalSeconds).ToString("0.00e00") + Environment::NewLine);

	if (CANCELLED)
		return;

	BGWRKR->ReportProgress(0, "Matching all the points that I can given the " + N_POINTS + " input point pairs...");
	array<double>^ cval1 = gcnew array<double>(total_pt_matches);
	array<double>^ cval2 = gcnew array<double>(total_pt_matches);
	array<double>^ cpix1 = gcnew array<double>(total_pt_matches);
	array<double>^ cpix2 = gcnew array<double>(total_pt_matches);
	c = 0;
	for (int k = 0; k < CATpts_intrmdt->Length; k++)
	{
		double X_int = CATpts_intrmdt[k]->X;
		double Y_int = CATpts_intrmdt[k]->Y;

		int x = (int)Math::Round((double)IMAGE_WIDTH - 1 - (1 / p00 * (Math::Cos(-p01) * X_int - Math::Sin(-p01) * Y_int) + p02));
		int y = (int)Math::Round((double)IMAGE_HEIGHT - 1 - (1 / p00 * (Math::Sin(-p01) * X_int + Math::Cos(-p01) * Y_int) + p03));

		if (x > 0 && y > 0 && x < IMAGE_WIDTH && y < IMAGE_HEIGHT && PSE->SourceBooleanMap[x, y])
		{
			int index = PSE->SourceIndexMap[x, y];
			cpix1[c] = PSE->Centroids_X[index];
			cpix2[c] = PSE->Centroids_Y[index];
			cval1[c] = CAT_PTS[k]->X;
			cval2[c] = CAT_PTS[k]->Y;
			c++;
		}
	}

	if (CANCELLED)
		return;

	BGWRKR->ReportProgress(0, "Solving for " + c + " point pair matches out of a possible " + N_POINTS);
	WCS->Solve_WCS("TAN", cpix1, cpix2, true, cval1, cval2, FITS_IMG);
	BGWRKR->ReportProgress(0, "Solution:" + Environment::NewLine);
	BGWRKR->ReportProgress(0, "CRPIX1 = " + WCS->CRPIXn[1]);
	BGWRKR->ReportProgress(0, "CRPIX2 = " + WCS->CRPIXn[2]);
	BGWRKR->ReportProgress(0, "CRVAL1 = " + WCS->CRVALn[1]);
	BGWRKR->ReportProgress(0, "CRVAL2 = " + WCS->CRVALn[2]);
	BGWRKR->ReportProgress(0, "CD_1_1 = " + WCS->CDi_j[1, 1]);
	BGWRKR->ReportProgress(0, "CD_1_2 = " + WCS->CDi_j[1, 2]);
	BGWRKR->ReportProgress(0, "CD_2_1 = " + WCS->CDi_j[2, 1]);
	BGWRKR->ReportProgress(0, "CD_2_2 = " + WCS->CDi_j[2, 2]);
	BGWRKR->ReportProgress(0, "CDELT1 = " + WCS->CDELTn[1]);
	BGWRKR->ReportProgress(0, "CDELT2 = " + WCS->CDELTn[2]);
	BGWRKR->ReportProgress(0, "CROTA1 = " + WCS->CROTAn[1]);
	BGWRKR->ReportProgress(0, "CROTA2 = " + WCS->CROTAn[2]);
	BGWRKR->ReportProgress(0, "WCS Fit Residual Mean (pixels) = " + WCS->WCSFitResidual_MeanPix);
	BGWRKR->ReportProgress(0, "WCS Fit Residual Stdv (pixels) = " + WCS->WCSFitResidual_StdvPix);
	BGWRKR->ReportProgress(0, "WCS Fit Residual Mean (arcsec) = " + WCS->WCSFitResidual_MeanSky);
	BGWRKR->ReportProgress(0, "WCS Fit Residual Stdv (arcsec) = " + WCS->WCSFitResidual_StdvSky + Environment::NewLine);

	if (!REFINE || FITS_IMG == nullptr)
	{
		BGWRKR->ReportProgress(0, "Finished...");
		return;
	}

	if (CANCELLED)
		return;

	BGWRKR->ReportProgress(0, "Refining solution...");
	PSE = gcnew JPFITS::SourceExtractor();
	double immax = FITS_IMG->Max;
	double pixthresh = immax / 256;
	double div = 2;
	double amp = pixthresh;
	int PSEiters = 0;
	int maxPSEiters = 15;
	N_POINTS *= 3;
	BGWRKR->ReportProgress(0, "Searching for " + N_POINTS + " point sources..." + Environment::NewLine);

	while (PSE->N_Sources < N_POINTS || PSE->N_Sources > N_POINTS + 1)
	{
		if (CANCELLED)
			return;

		PSEiters++;
		if (PSEiters > maxPSEiters)
			break;

		if (PSE->N_SaturatedSources >= N_POINTS || PSE->N_Sources >= N_POINTS)
			break;

		PSE->Extract_Sources(FITS_IMG->Image, PIX_SAT, pixthresh, Double::MaxValue, 0, Double::MaxValue, false, PSE_KERNEL_RADIUS, PSE_SEP_RADIUS, AUTO_BACKGROUND, "", IMAGE_ROI, false);

		if (PSE->N_Sources < N_POINTS)
			pixthresh -= amp / div;
		if (PSE->N_Sources > N_POINTS + 1)
			pixthresh += amp / div;
		div *= 2;

		BGWRKR->ReportProgress(0, "Found " + PSE->N_Sources + " point sources on iteration " + PSEiters);
	}
	if (PSE->N_Sources > N_POINTS)
		PSE->ClipToNBrightest(N_POINTS);
	N_POINTS = PSE->N_Sources;
	BGWRKR->ReportProgress(0, Environment::NewLine + "Stopped searching on iteration " + PSEiters + " with " + N_POINTS + " point sources");

	if (CANCELLED)
		return;

	if (N_POINTS > CAT_CVAL1s->Length)
		N_POINTS = CAT_CVAL1s->Length;

	//get the brightest catlaogue points
	cval1 = gcnew array<double>(N_POINTS);
	cval2 = gcnew array<double>(N_POINTS);
	for (int i = 0; i < N_POINTS; i++)
	{
		cval1[i] = CAT_CVAL1s[i];
		cval2[i] = CAT_CVAL2s[i];
	}

	//get the wcs pixel locations
	cpix1 = gcnew array<double>(N_POINTS);
	cpix2 = gcnew array<double>(N_POINTS);
	WCS->Get_Pixels(cval1, cval2, "TAN", cpix1, cpix2, true);

	//check for WCS pixels which fall onto PSE pixels
	int nmatches = 0;
	array<bool>^ match = gcnew array<bool>(N_POINTS);
	array<int>^ matchinds = gcnew array<int>(N_POINTS);
	for (int i = 0; i < N_POINTS; i++)
	{
		int x = (int)Math::Round(cpix1[i]);
		int y = (int)Math::Round(cpix2[i]);
		if (x > 0 && x < IMAGE_WIDTH && y > 0 && y < IMAGE_HEIGHT)
			if (PSE->SourceBooleanMap[x, y])
			{
				nmatches++;
				match[i] = true;
				matchinds[i] = PSE->SourceIndexMap[x, y];
			}
	}

	int n = 0;
	for (int i = 0; i < N_POINTS; i++)
	{
		if (!match[i])
			continue;

		cval1[n] = CAT_CVAL1s[i];
		cval2[n] = CAT_CVAL2s[i];
		cpix1[n] = PSE->Centroids_X[matchinds[i]];
		cpix2[n] = PSE->Centroids_Y[matchinds[i]];
		n++;
	}
	Array::Resize(cval1, nmatches);
	Array::Resize(cval2, nmatches);
	Array::Resize(cpix1, nmatches);
	Array::Resize(cpix2, nmatches);

	if (CANCELLED)
		return;

	WCS->Clear(FITS_IMG);
	WCS->Solve_WCS("TAN", cpix1, cpix2, true, cval1, cval2, FITS_IMG);
	BGWRKR->ReportProgress(0, Environment::NewLine + nmatches + " sources of " + N_POINTS + " were able to be used for WCS refinement.");
	BGWRKR->ReportProgress(0, Environment::NewLine + "Refined solution:" + Environment::NewLine);
	BGWRKR->ReportProgress(0, "CRPIX1 = " + WCS->CRPIXn[1]);
	BGWRKR->ReportProgress(0, "CRPIX2 = " + WCS->CRPIXn[2]);
	BGWRKR->ReportProgress(0, "CRVAL1 = " + WCS->CRVALn[1]);
	BGWRKR->ReportProgress(0, "CRVAL2 = " + WCS->CRVALn[2]);
	BGWRKR->ReportProgress(0, "CD_1_1 = " + WCS->CDi_j[1, 1]);
	BGWRKR->ReportProgress(0, "CD_1_2 = " + WCS->CDi_j[1, 2]);
	BGWRKR->ReportProgress(0, "CD_2_1 = " + WCS->CDi_j[2, 1]);
	BGWRKR->ReportProgress(0, "CD_2_2 = " + WCS->CDi_j[2, 2]);
	BGWRKR->ReportProgress(0, "CDELT1 = " + WCS->CDELTn[1]);
	BGWRKR->ReportProgress(0, "CDELT2 = " + WCS->CDELTn[2]);
	BGWRKR->ReportProgress(0, "CROTA1 = " + WCS->CROTAn[1]);
	BGWRKR->ReportProgress(0, "CROTA2 = " + WCS->CROTAn[2]);
	BGWRKR->ReportProgress(0, "WCS Fit Residual Mean (pixels) = " + WCS->WCSFitResidual_MeanPix);
	BGWRKR->ReportProgress(0, "WCS Fit Residual Stdv (pixels) = " + WCS->WCSFitResidual_StdvPix);
	BGWRKR->ReportProgress(0, "WCS Fit Residual Mean (arcsec) = " + WCS->WCSFitResidual_MeanSky);
	BGWRKR->ReportProgress(0, "WCS Fit Residual Stdv (arcsec) = " + WCS->WCSFitResidual_StdvSky + Environment::NewLine);
	BGWRKR->ReportProgress(0, "Finished..." + Environment::NewLine);
}

void JPFITS::WCS_AutoSolver::BGWRKR_ProgressChanged(System::Object^ sender, System::ComponentModel::ProgressChangedEventArgs^ e)
{
	PROGRESS = e->ProgressPercentage;

	if (PROGRESS == 0)
		STATUS_LOG += Environment::NewLine + (String^)e->UserState;
	else
	{
		STATUS_LOG = STATUS_LOG->Remove(STATUS_LOG->LastIndexOf(Environment::NewLine));
		STATUS_LOG += Environment::NewLine + "Approximate progress: " + (PROGRESS + 1).ToString() + "%";
	}

	
}

void JPFITS::WCS_AutoSolver::BGWRKR_RunWorkerCompleted(System::Object^ sender, System::ComponentModel::RunWorkerCompletedEventArgs^ e)
{
	if (CANCELLED)
	{
		WCS = nullptr;
		SOLVED = false;
		SOLVING = false;
		return;
	}

	if (DO_PSE)
	{
		DO_PSE = false;
		SolveAsync(SCALE_INIT, SCALE_LB, SCALE_UB, ROTATION_INIT, ROTATION_LB, ROTATION_UB, WCS_VERTEX_TOL, N_MATCHES_STOP, PERC_MATCHES_STOP, CONDITION_ARRAYS, SHOW_REPORT_FORM);
		return;
	}

	if (!SOLVED)
	{
		SOLVING = false;
		if (WCSARF != nullptr)
			WCSARF->CancelBtn->DialogResult = ::DialogResult::No;
		return;
	}

	SOLVING = false;
	if (WCSARF != nullptr)
	{
		WCSARF->CancelBtn->Text = "OK";
		WCSARF->CancelBtn->DialogResult = ::DialogResult::OK;
	}
}

array<JPMath::Triangle^>^ JPFITS::WCS_AutoSolver::ConditionTriangleArrayBrightnessThreads(array<JPMath::Triangle^>^ triarray, int Nthreads, bool ascending)
{
	//reformat traingle arrays for threading
	int Ntris = triarray->Length;
	array<double>^ trivals = gcnew array<double>(Ntris);

	#pragma omp parallel for
	for (int i = 0; i < Ntris; i++)
		trivals[i] = triarray[i]->VertexPointSum;

	Array::Sort(trivals, triarray);

	if (!ascending)
		Array::Reverse(triarray);

	if (Nthreads == 1)
	{
		array<JPMath::Triangle^>^ retarray = gcnew array<JPMath::Triangle^>(Ntris);
		Array::Copy(triarray, retarray, Ntris);
		return retarray;
	}

	Ntris = (Ntris / Nthreads) * Nthreads;
	Array::Resize(triarray, Ntris);

	int dim1 = triarray->Length / Nthreads;
	array<JPMath::Triangle^, 2>^ temptris = gcnew array<JPMath::Triangle^, 2>(Nthreads, dim1);

	#pragma omp parallel for
	for (int i = 0; i < dim1; i++)
		for (int j = 0; j < Nthreads; j++)
			temptris[j, i] = triarray[i * Nthreads + j];

	#pragma omp parallel for
	for (int i = 0; i < dim1; i++)
		for (int j = 0; j < Nthreads; j++)
			triarray[i + j * dim1] = temptris[j, i];

	array<JPMath::Triangle^>^ retarray = gcnew array<JPMath::Triangle^>(Ntris);
	Array::Copy(triarray, retarray, Ntris);
	return retarray;
}

int JPFITS::WCS_AutoSolver::AstroQuery(String^ catalogue, String^ ra_deg, String^ dec_deg, String^ & result_savepathfilename, String^ radius, String^ square)
{
	String^ pypath = (String^)GetReg("CCDLAB", "PythonExePath");
	if (pypath == nullptr || !File::Exists(pypath))
	{
		array<String^>^ dirs = ::Directory::GetDirectories("C:\\Program Files\\", "*Python*", ::SearchOption::TopDirectoryOnly);
		if (dirs->Length == 1)
		{
			::DialogResult res = MessageBox::Show("Use " + dirs[0] + "\\python.exe as your Python environment?", "Finding python.exe...", ::MessageBoxButtons::YesNoCancel);
			if (res == ::DialogResult::Cancel)
				return -2;
			else if (res == ::DialogResult::Yes)
				pypath = dirs[0] + "\\python.exe";
			else if (res == ::DialogResult::No)
			{
				if (MessageBox::Show("Do you want to show me where your Python installation is located, then?", "Find *python.exe*...", ::MessageBoxButtons::OKCancel) == ::DialogResult::Cancel)
					return -2;

				::OpenFileDialog^ ofd = gcnew ::OpenFileDialog();
				ofd->InitialDirectory = "C:\\Program Files\\";
				ofd->Filter = "exe|*.exe;*";
				if (ofd->ShowDialog() == ::DialogResult::Cancel)
					return -2;

				pypath = ofd->FileName;
			}
		}
		else if (dirs->Length == 0 || dirs->Length > 1)
		{
			if (MessageBox::Show("Please show me where your Python installation is located, OK?", "I cannot find *python.exe*...", ::MessageBoxButtons::OKCancel) == ::DialogResult::Cancel)
				return -2;

			::OpenFileDialog^ ofd = gcnew ::OpenFileDialog();
			ofd->InitialDirectory = "C:\\Program Files\\";
			ofd->Filter = "exe|*.exe;*";
			if (ofd->ShowDialog() == ::DialogResult::Cancel)
				return -2;

			pypath = ofd->FileName;
		}

		SetReg("CCDLAB", "PythonExePath", pypath);
	}

	catalogue = catalogue->ToLower();
	String^ script = "C:\\ProgramData\\Astrowerks\\CCDLABx64\\astro_query.py";
	MAKEASTROQUERYSCRIPT(script, catalogue);

	if (result_savepathfilename == "")
	{
		if (!::Directory::Exists("C:\\ProgramData\\Astrowerks\\CCDLABx64\\"))
			::Directory::CreateDirectory("C:\\ProgramData\\Astrowerks\\CCDLABx64\\");
		result_savepathfilename = "C:\\ProgramData\\Astrowerks\\CCDLABx64\\queryCatalog.fit";
	}

	::Diagnostics::ProcessStartInfo^ psi = gcnew ::Diagnostics::ProcessStartInfo();
	psi->FileName = pypath;
	psi->Arguments = String::Format("\"" + script + "\"" + " {0} {1} {2} {3} {4}", ra_deg, dec_deg, "\"" + result_savepathfilename + "\"", radius, square);

	/*psi->UseShellExecute = false;//??????
	psi->CreateNoWindow = true;//????
	psi->RedirectStandardError = true;
	psi->RedirectStandardOutput = true;
	String^ errs = "";
	String^ res = "";*/

	::Diagnostics::Process^ proc = ::Diagnostics::Process::Start(psi);
	proc->WaitForExit();
	int res = proc->ExitCode;
	if (res != 0)
		return res;

	/*errs = proc->StandardError->ReadToEnd();
	res = proc->StandardOutput->ReadToEnd();
	MessageBox::Show(errs + "\r\n" + res);*/

	array<String^>^ ExtensionEntryLabels = FITSBinTable::GetExtensionEntryLabels(result_savepathfilename, "");
	array<TypeCode>^ ExtensionEntryDataTypes = FITSBinTable::GetExtensionEntryDataTypes(result_savepathfilename, "");
	array<String^>^ ExtensionEntryDataUnits = FITSBinTable::GetExtensionEntryUnits(result_savepathfilename, "");
	array<double, 2>^ table = FITSBinTable::GetExtensionEntries(result_savepathfilename, "", ExtensionEntryLabels);
	FITSImage^ fits = gcnew FITSImage(result_savepathfilename, nullptr, true, true, false, false);
	//fits->AddKey("RA", ra_deg, "Right Ascension of query field center, degrees", -1);
	//fits->AddKey("DEC", dec_deg, "Declination of query field center, degrees", -1);
	fits->WriteFile(TypeCode::Double, false);
	array<String^>^ exkeys = gcnew array<String^>(2) { "RA", "DEC" };
	array<String^>^ exvals = gcnew array<String^>(2) { ra_deg->ToString(), dec_deg->ToString() };
	array<String^>^ excoms = gcnew array<String^>(2) { "Right Ascension of query field center, degrees", "Declination of query field center, degrees" };
	FITSBinTable::WriteExtension(result_savepathfilename, "", true, ExtensionEntryLabels, ExtensionEntryDataTypes, ExtensionEntryDataUnits, exkeys, exvals, excoms, table);

	return res;
}

void JPFITS::WCS_AutoSolver::MAKEASTROQUERYSCRIPT(String^ script_filename, String^ catalogue)
{
	String^ script = "";
	script += "import argparse" + Environment::NewLine;
	script += "import sys" + Environment::NewLine;
	script += "from astroquery.simbad import Simbad" + Environment::NewLine;
	script += "from astropy.coordinates import SkyCoord" + Environment::NewLine;
	script += "import astropy.units as u" + Environment::NewLine;
	if (catalogue == "gaia")
		script += "from astroquery.gaia import Gaia" + Environment::NewLine;

	script += "ra = float(sys.argv[1])" + Environment::NewLine;
	script += "dec = float(sys.argv[2])" + Environment::NewLine;
	script += "filename = str(sys.argv[3])" + Environment::NewLine;
	script += "radius = float(sys.argv[4])" + Environment::NewLine;
	script += "square = float(sys.argv[5])" + Environment::NewLine;

	script += "radvals = radius * u.arcmin" + Environment::NewLine;
	script += "if square == 1:" + Environment::NewLine;
	script += "    radvals = radvals * 2;" + Environment::NewLine;

	script += "coords = SkyCoord(ra = ra, dec = dec, unit = (u.deg, u.deg), frame = 'fk5')" + Environment::NewLine;

	if (catalogue == "gaia")
	{
		script += "jobstr = \"SELECT * FROM gaiadr2.gaia_source\"" + Environment::NewLine;
		script += "jobstr += \" WHERE CONTAINS(POINT('ICRS', gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),\"" + Environment::NewLine;
	}
	script += "if square == 1:" + Environment::NewLine;
	script += "    jobstr += \"BOX('ICRS',{0},{1},{2},{2}))=1;\".format(coords.ra.deg, coords.dec.deg, radvals.to(u.deg).value)" + Environment::NewLine;
	script += "else:" + Environment::NewLine;
	script += "    jobstr += \"CIRCLE('ICRS',{0},{1},{2}))=1;\".format(coords.ra.deg, coords.dec.deg, radvals.to(u.deg).value)" + Environment::NewLine;

	if (catalogue == "gaia")
		script += "print(\"Launching job query to Gaia archive\")" + Environment::NewLine;
	script += "print(jobstr)" + Environment::NewLine;
	script += "print(\" \")" + Environment::NewLine;
	script += "print(\"Waiting for query results...\")" + Environment::NewLine;
	if (catalogue == "gaia")
		script += "job = Gaia.launch_job_async(jobstr, dump_to_file = False)" + Environment::NewLine;
	script += "print(job)" + Environment::NewLine;
	script += "results = job.get_results()" + Environment::NewLine;
	script += "removelist = []" + Environment::NewLine;

	//Strip object columns from FITS table
	script += "for col in results.columns:" + Environment::NewLine;
	script += "    if results[col].dtype == 'object' :" + Environment::NewLine;
	script += "        removelist += [col]" + Environment::NewLine;
	script += "results.remove_columns(removelist)" + Environment::NewLine;
	script += "results.write(filename, overwrite = True, format = 'fits')";

	StreamWriter^ sw = gcnew StreamWriter(script_filename);
	sw->Write(script);
	sw->Close();
}

