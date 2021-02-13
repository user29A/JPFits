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


#include "stdafx.h"
#include "JPFITS.h"
using namespace JPFITS;


JPFITS::FITSImageSet::FITSImageSet()
{
	FITSLIST = gcnew ArrayList();
	CODIMENSIONAL = true;
	this->BGWRKR = gcnew BackgroundWorker();
	this->BGWRKR->WorkerReportsProgress = true;
	this->BGWRKR->WorkerSupportsCancellation = true;
	this->BGWRKR->DoWork += gcnew System::ComponentModel::DoWorkEventHandler(this, &JPFITS::FITSImageSet::BGWRKR_DoWork);
	this->BGWRKR->ProgressChanged += gcnew System::ComponentModel::ProgressChangedEventHandler(this, &JPFITS::FITSImageSet::BGWRKR_ProgressChanged);
	this->BGWRKR->RunWorkerCompleted += gcnew System::ComponentModel::RunWorkerCompletedEventHandler(this, &JPFITS::FITSImageSet::BGWRKR_RunWorkerCompleted);
}

bool JPFITS::FITSImageSet::Write(TypeCode precision, bool do_parallel, bool show_waitbar, String^ waitbar_message)
{
	if (!show_waitbar)
	{
		#pragma omp parallel for if (do_parallel)
		for (int i = 0; i < FITSLIST->Count; i++)
			((FITSImage^)FITSLIST[i])->WriteImage(precision, !do_parallel);
		return true;
	}
	else
	{
		this->WAITBAR = gcnew JPWaitBar::WaitBar();
		this->WAITBAR->ProgressBar->Maximum = FITSLIST->Count;
		this->WAITBAR->Text = "Saving Image Set: " + FITSLIST->Count + " files...";
		//this->WAITBAR->StartPosition = Windows::Forms::FormStartPosition::;
		array<Object^>^ arg = gcnew array<Object^>(5);
		arg[0] = "";
		arg[1] = "save";
		arg[2] = do_parallel;
		arg[3] = precision;
		arg[4] = waitbar_message;
		BGWRKR->RunWorkerAsync(arg);
		this->WAITBAR->ShowDialog();
		if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
			return false;
		else
			return true;
	}
}

bool JPFITS::FITSImageSet::Load(array<String^>^ files, array<int>^ range, bool do_stats, bool do_parallel, bool show_waitbar, String^ waitbar_message)
{
	if (!show_waitbar)
	{
		array<JPFITS::FITSImage^>^ set = gcnew array<JPFITS::FITSImage^>(files->Length);
		#pragma omp parallel for if (do_parallel)
		for (int i = 0; i < files->Length; i++)
			set[i] = gcnew FITSImage(files[i], range, true, true, do_stats, !do_parallel);
		for (int i = 0; i < files->Length; i++)
			FITSLIST->Add(set[i]);
		return true;
	}
	else
	{
		this->WAITBAR = gcnew JPWaitBar::WaitBar();
		this->WAITBAR->ProgressBar->Maximum = files->Length;
		this->WAITBAR->Text = "Loading Image Set: " + files->Length + " files...";
		//this->WAITBAR->StartPosition = Windows::Forms::FormStartPosition::CenterScreen;
		array<Object^>^ arg = gcnew array<Object^>(6);
		arg[0] = files;
		arg[1] = "load";
		arg[2] = do_stats;
		arg[3] = do_parallel;
		arg[4] = range;
		arg[5] = waitbar_message;
		BGWRKR->RunWorkerAsync(arg);
		this->WAITBAR->ShowDialog();
		if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
			return false;
		else
			return true;
	}
}

void JPFITS::FITSImageSet::BGWRKR_DoWork(System::Object^  sender, System::ComponentModel::DoWorkEventArgs^  e)
{
	array<Object^>^ arg = (array<Object^>^)e->Argument;
	String^ op = arg[1]->ToString();

	if (op == "save")
	{
		bool do_parallel = bool(arg[2]);
		TypeCode precision = (TypeCode)arg[3];
		String^ waitbar_message = (String^)arg[4];

		int n = 1;
		if (do_parallel)
			n = omp_get_max_threads();
		array<int>^ count = gcnew array<int>(n);
		int prog = 0;

		#pragma omp parallel for if (do_parallel)
		for (int i = 0; i < FITSLIST->Count; i++)
		{
			if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
				continue;

			count[omp_get_thread_num()]++;
			#pragma omp critical
			{
				int sum = 0;
				for (int si = 0; si < count->Length; si++)
					sum += count[si];
				if (sum * 100 / FITSLIST->Count > prog)
				{
					BGWRKR->ReportProgress(sum, "save " + waitbar_message);
					prog = sum * 100 / FITSLIST->Count;
				}
			}

			((FITSImage^)FITSLIST[i])->WriteImage(precision, !do_parallel);
		}
		return;
	}

	if (op == "load")
	{
		array<String^>^ files = (array<String^>^)arg[0];
		bool do_stats = bool(arg[2]);
		bool do_parallel = bool(arg[3]);
		array<int>^ range = (array<int>^)arg[4];
		String^ waitbar_message = (String^)arg[5];

		array<JPFITS::FITSImage^>^ set = gcnew array<JPFITS::FITSImage^>(files->Length);

		int n = 1;
		if (do_parallel)
			n = omp_get_max_threads();
		array<int>^ count = gcnew array<int>(n);
		int prog = 0;

		#pragma omp parallel for if (do_parallel)
		for (int i = 0; i < files->Length; i++)
		{
			if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
				continue;

			count[omp_get_thread_num()]++;
			#pragma omp critical
			{
				int sum = 0;
				for (int si = 0; si < count->Length; si++)
					sum += count[si];
				if (sum * 100 / files->Length > prog)
				{
					BGWRKR->ReportProgress(sum, "load " + waitbar_message);
					prog = sum * 100 / files->Length;
				}
			}

			set[i] = gcnew FITSImage(files[i], range, true, true, do_stats, !do_parallel);
		}
		if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
			return;

		for (int i = 0; i < files->Length; i++)
			FITSLIST->Add(set[i]);
		return;
	}
	
	if (op->Equals("Min"))
	{
		FITSImageSet^ ImageSet = (FITSImageSet^)arg[0];
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;
		double N = (double)L;
		int n = 1;//def num threads
		int intprog = 0;
		int prog = -1;

		#pragma omp parallel for
		for (int i=0; i < width; i++)
		{
			if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
			{
				e->Result = nullptr;
				continue;
			}
			n = omp_get_num_threads();
			intprog = i*100*n/width;
			if (i < (int)(width/n) && intprog-prog > 1)//keep the update of progress bar to only one thread of the team...avoids locks
			{
				prog++;
				BGWRKR->ReportProgress(intprog);
			}

			for (int j=0; j < height; j++)
			{
				double min = Double::MaxValue;
				for (int k = 0; k < L; k++)
					if (ImageSet[k]->Image[i, j] < min)
						min = ImageSet[k]->Image[i, j];
				img[i,j] = min;
			}
		}
		e->Result = img;
		return;
	}

	if (op->Equals("Max"))
	{
		FITSImageSet^ ImageSet = (FITSImageSet^)arg[0];
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;
		double N = (double)L;
		int n = 1;//def num threads
		int intprog = 0;
		int prog = -1;

		#pragma omp parallel for
		for (int i=0; i < width; i++)
		{
			if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
			{
				e->Result = nullptr;
				continue;
			}
			n = omp_get_num_threads();
			intprog = i*100*n/width;
			if (i < (int)(width/n) && intprog-prog > 1)//keep the update of progress bar to only one thread of the team...avoids locks
			{
				prog++;
				BGWRKR->ReportProgress(intprog);
			}

			for (int j=0; j < height; j++)
			{
				double max = Double::MinValue;
				for (int k = 0; k < L; k++)
					if (ImageSet[k]->Image[i, j] > max)
						max = ImageSet[k]->Image[i, j];
				img[i,j] = max;
			}
		}
		e->Result = img;
		return;
	}

	if (op->Equals("Mean"))
	{
		FITSImageSet^ ImageSet = (FITSImageSet^)arg[0];
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;
		double N = (double)L;
		int n = 1;//def num threads
		int intprog = 0;
		int prog = -1;

		#pragma omp parallel for
		for (int i=0; i < width; i++)
		{
			if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
			{
				e->Result = nullptr;
				continue;
			}
			n = omp_get_num_threads();
			intprog = i*100*n/width;
			if (i < (int)(width/n) && intprog-prog > 1)//keep the update of progress bar to only one thread of the team...avoids locks
			{
				prog++;
				BGWRKR->ReportProgress(intprog);
			}

			for (int j=0; j < height; j++)
			{
				double mean = 0;
				for (int k = 0; k < L; k++)
					mean = mean + ImageSet[k]->Image[i,j];
				img[i,j] = mean/N;
			}
		}
		e->Result = img;
		return;
	}

	if (op->Equals("Sum"))
	{
		FITSImageSet^ ImageSet = (FITSImageSet^)arg[0];
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;
		double N = (double)L;
		int n = 1;//def num threads
		int intprog = 0;
		int prog = -1;

		#pragma omp parallel for
		for (int i=0; i < width; i++)
		{
			if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
			{
				e->Result = nullptr;
				continue;
			}
			n = omp_get_num_threads();
			intprog = i*100*n/width;
			if (i < (int)(width/n) && intprog-prog > 1)//keep the update of progress bar to only one thread of the team...avoids locks
			{
				prog++;
				BGWRKR->ReportProgress(intprog);
			}

			for (int j=0; j < height; j++)
			{
				double sum = 0;
				for (int k = 0; k < L; k++)
					sum = sum + ImageSet[k]->Image[i,j];
				img[i,j] = sum;
			}
		}
		e->Result = img;
		return;
	}

	if (op->Equals("Quadrature"))
	{
		FITSImageSet^ ImageSet = (FITSImageSet^)arg[0];
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;
		double N = (double)L;
		int n = 1;//def num threads
		int intprog = 0;
		int prog = -1;

		#pragma omp parallel for
		for (int i = 0; i < width; i++)
		{
			if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
			{
				e->Result = nullptr;
				continue;
			}
			n = omp_get_num_threads();
			intprog = i * 100 * n / width;
			if (i < (int)(width / n) && intprog - prog > 1)//keep the update of progress bar to only one thread of the team...avoids locks
			{
				prog++;
				BGWRKR->ReportProgress(intprog);
			}

			for (int j = 0; j < height; j++)
			{
				double quadsum = 0;
				for (int k = 0; k < L; k++)
					quadsum += (ImageSet[k]->Image[i, j] * ImageSet[k]->Image[i, j]);
				img[i, j] = Math::Sqrt(quadsum);
			}
		}
		e->Result = img;
		return;
	}

	if (op->Equals("Median"))
	{
		FITSImageSet^ ImageSet = (FITSImageSet^)arg[0];
		String^ waitbar_message = (String^)(arg[2]);
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;
		double N = (double)L;
		int n = 1;//def num threads
		int intprog = 0;
		int prog = -1;

		#pragma omp parallel for
		for (int i = 0; i < width; i++)
		{
			if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
			{
				e->Result = nullptr;
				continue;
			}
			n = omp_get_num_threads();
			intprog = i * 100 * n / width;
			if (i < (int)(width/n) && intprog-prog > 1)//keep the update of progress bar to only one thread of the team...avoids locks
			{
				prog++;
				BGWRKR->ReportProgress(intprog, waitbar_message);
			}

			for (int j = 0; j < height; j++)
			{
				array<double,1>^ medarray = gcnew array<double,1>(L);
				for (int k = 0; k < L; k++)
					medarray[k] = ImageSet[k]->Image[i, j];
				img[i, j] = JPMath::Median(medarray);
			}
		}
		e->Result = img;
		return;
	}

	if (op->Equals("Stdv"))
	{
		FITSImageSet^ ImageSet = (FITSImageSet^)arg[0];
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;
		double N = (double)L;
		int n = 1;//def num threads
		int intprog = 0;
		int prog = -1;

		#pragma omp parallel for
		for (int i=0; i < width; i++)
		{
			if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
			{
				e->Result = nullptr;
				continue;
			}
			n = omp_get_num_threads();
			intprog = i*100*n/width;
			if (i < (int)(width/n) && intprog-prog > 1)//keep the update of progress bar to only one thread of the team...avoids locks
			{
				prog++;
				BGWRKR->ReportProgress(intprog);
			}

			for (int j=0; j < height; j++)
			{
				double std = 0;
				double mean = 0;
				for (int k = 0; k < L; k++)
					mean = mean + ImageSet[k]->Image[i,j];
				mean = mean/N;
				for (int q = 0; q < L; q++)
					std = std + (ImageSet[q]->Image[i,j]-mean)*(ImageSet[q]->Image[i,j]-mean);
				std = Math::Sqrt(std/(N-1.0));
				img[i, j] = std;
			}
		}
		e->Result = img;
		return;
	}

	if (op->Equals("AutoReg"))
	{
		FITSImageSet^ ImageSet = (FITSImageSet^)arg[0];
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;
		double N = (double)L;
		int n = 1;//def num threads
		int intprog = 0;
		int prog = -1;

		int RefImgIndex = (int)arg[2];
		bool Do_Stats = (bool)arg[3];

		array<double, 2>^ ref = gcnew array<double, 2>(ImageSet[RefImgIndex]->Width, ImageSet[RefImgIndex]->Height);
		Array::Copy(ImageSet[RefImgIndex]->Image, ref, ImageSet[RefImgIndex]->Image->LongLength);

		ref = JPMath::DeGradient(ref, 0, true);
		ref = JPMath::DeGradient(ref, 1, true);
		ref = JPMath::Hanning(ref, true);
		array<double, 1>^ Href = JPMath::Sum(ref, 1, true);
		array<double, 1>^ Vref = JPMath::Sum(ref, 0, true);
		Href = JPMath::VectorSubScalar(Href, JPMath::Mean(Href, true), true);
		Vref = JPMath::VectorSubScalar(Vref, JPMath::Mean(Vref, true), true);

		for (int i = 0; i < ImageSet->Count; i++)//create the array with the missing reference index
		{
			BGWRKR->ReportProgress(i * 100 / (ImageSet->Count - 1));

			if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
				continue;
			if (i == RefImgIndex)
				continue;//don't register to one's self

			double xshift = 0, yshift = 0;
			JPMath::XCorrImageLagShifts(Href, Vref, ImageSet[i]->Image, true, true, true, xshift, yshift, true);

			int intxshift = (int)Math::Round(xshift);
			int intyshift = (int)Math::Round(yshift);

			array<double, 2>^ COM = gcnew array<double, 2>(ImageSet[i]->Width, ImageSet[i]->Height);
			Array::Copy(ImageSet[i]->Image, COM, ImageSet[i]->Image->LongLength);

			double mean = JPMath::Median(COM);

			if (intxshift > 0)//shift REL in horizontal dim (0) to the left
			{
				#pragma omp parallel for
				for (int j = 0; j < COM->GetLength(1); j++)
				{
					for (int i = 0; i < COM->GetLength(0) - intxshift; i++)
						COM[i, j] = COM[i + intxshift, j];
					for (int i = COM->GetLength(0) - intxshift; i < COM->GetLength(0); i++)
						COM[i, j] = mean;
				}
			}

			if (intxshift < 0)//shift REL in horizontal dim (0) to the right
			{
				intxshift = -intxshift;
				#pragma omp parallel for
				for (int j = 0; j < COM->GetLength(1); j++)
				{
					for (int i = 0; i < COM->GetLength(0) - intxshift; i++)
						COM[COM->GetLength(0) - i - 1, j] = COM[COM->GetLength(0) - intxshift - i - 1, j];
					for (int i = 0; i < intxshift; i++)
						COM[i, j] = mean;
				}
			}

			if (intyshift > 0)//shift REL in horizontal dim (0) to the left
			{
				#pragma omp parallel for
				for (int i = 0; i < COM->GetLength(0); i++)
				{
					for (int j = 0; j < COM->GetLength(1) - intyshift; j++)
						COM[i, j] = COM[i, j + intyshift];
					for (int j = COM->GetLength(1) - intyshift; j < COM->GetLength(1); j++)
						COM[i, j] = mean;
				}
			}

			if (intyshift < 0)//shift REL in horizontal dim (0) to the right
			{
				intyshift = -intyshift;
				#pragma omp parallel for
				for (int i = 0; i < COM->GetLength(0); i++)
				{
					for (int j = 0; j < COM->GetLength(1) - intyshift; j++)
						COM[i, COM->GetLength(1) - j - 1] = COM[i, COM->GetLength(1) - intyshift - j - 1];
					for (int j = 0; j < intyshift; j++)
						COM[i, j] = mean;
				}
			}

			ImageSet[i]->SetImage(COM, Do_Stats, true);
		}
		return;
	}

	if (op->Equals("SCM"))
	{
		FITSImageSet^ ImageSet = (FITSImageSet^)arg[0];
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;
		double N = (double)L;
		int n = 1;//def num threads
		int intprog = 0;
		int prog = -1;

		double sigma = (double)arg[2];
		array<double,2>^ simg = (array<double,2>^)arg[3];
		img  = (array<double,2>^)arg[4];

		bool global_need = true;
		int global_count = 0;
		int psuedo_global_count = 0;

		int nptsrepeat1 = 0;
		int nptsrepeat2 = 0;

		while (global_need == true)
		{
			if (this->WAITBAR->DialogResult == ::DialogResult::Cancel)
			{
				//still return image because we are just stopping the iteration
				global_need = false;
				goto stop_scm;
			}

			double std = JPMath::Stdv(simg, true);
			double mean = JPMath::Mean(simg, true);
			array<int,2>^ pts = JPMath::Find(simg,mean + sigma*std,">", true);//find pts > clip range

			if (pts == nullptr)
				goto stop_scm;

			BGWRKR->ReportProgress(psuedo_global_count + 1,String::Concat("Iteration: ",global_count + 1,". # of Offending Points: ",pts->GetLength(0)));

			if (pts != nullptr)//then do some clippin!
			{
				if (nptsrepeat1 != pts->GetLength(0))
					nptsrepeat1 = pts->GetLength(0);
				else
					nptsrepeat2++;

				double med;
				bool need = true;
				bool first = true;
				int count = 0;//breakout counter if too high
				array<double,1> ^clipvec = gcnew array<double,1>(L);
				double max;
				int index;

				for (int c = 0; c < pts->GetLength(0); c++)
				{
					need = true;
					first = true;
					while (need == true)
					{
						if (first == true)
						{
							for (int i = 0; i < L; i++)
							{
								clipvec[i] = ImageSet[i]->Image[pts[c, 0], pts[c, 1]];
							}
							first = false;
						}
						med = JPMath::Median(clipvec);
						max = JPMath::Max(JPMath::Abs(JPMath::VectorSubScalar(clipvec, med, false), false), index, false);
						if (max > sigma*std)
						{
							clipvec[index] = med;
							img[pts[c,0],pts[c,1]] = JPMath::Mean(clipvec, true);
							simg[pts[c,0],pts[c,1]] = JPMath::Stdv(clipvec, true);
						}
						else
						{
							need = false;
						}
						count++;
						if (count > 3*pts->GetLength(0))//?arbitrary? number of iterations limit
							need = false;//and go to next point
					}
				}
			}
			else
			{
				global_need = false;
			}
			global_count++;
			if (global_count > 2000 || nptsrepeat2 > 5)
				global_need = false;
			psuedo_global_count++;
			if (psuedo_global_count >= 100)
				psuedo_global_count = 0;
		}
		stop_scm:;
		return;
	}
}

void JPFITS::FITSImageSet::BGWRKR_ProgressChanged(System::Object^  sender, System::ComponentModel::ProgressChangedEventArgs^  e)
{
	this->WAITBAR->ProgressBar->Value = e->ProgressPercentage;

	if (e->UserState == nullptr)
	{
		this->WAITBAR->TextMsg->Text = e->ProgressPercentage + "%";
	}
	else if (((String^)e->UserState)->Substring(0, 4) == "load")
	{
		this->WAITBAR->TextMsg->Text = e->ProgressPercentage + " " + ((String^)e->UserState)->Substring(4);
	}
	else if (((String^)e->UserState)->Substring(0, 4) == "save")
	{
		this->WAITBAR->TextMsg->Text = e->ProgressPercentage + " " + ((String^)e->UserState)->Substring(4);
	}
	else
	{
		this->WAITBAR->TextMsg->Text = e->ProgressPercentage + "% " + (String^)e->UserState;
	}
	this->WAITBAR->Refresh();
}

void JPFITS::FITSImageSet::BGWRKR_RunWorkerCompleted(System::Object^  sender, System::ComponentModel::RunWorkerCompletedEventArgs^  e)
{
	this->BGWRKR_RESULT = e->Result;
	this->WAITBAR->DialogResult = ::DialogResult::OK;
	this->WAITBAR->Close();
}

FITSImage^ JPFITS::FITSImageSet::Mean(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar)
{
	if (ImageSet->CoDimensional == false)
	{
		MessageBox::Show("Can Not Perform Mean: Data Stack not Co-Dimensional", "Invalid Operation...");
		return nullptr;
	}
	if (ImageSet->Count <= 1)
	{
		MessageBox::Show("Can Not Perform Mean: Data Stack Contains One or Fewer Images", "Invalid Operation...");
		return nullptr;
	}

	if (Show_Waitbar)
	{
		WAITBAR = gcnew JPWaitBar::WaitBar();
		WAITBAR->ProgressBar->Maximum = 100;
		WAITBAR->Text = "Computing Mean Data Stack Image";
		array<Object^>^ arg = gcnew array<Object^>{ImageSet, "Mean"};
		BGWRKR->RunWorkerAsync(arg);
		WAITBAR->ShowDialog();
		if (WAITBAR->DialogResult == ::DialogResult::Cancel)
		{
			WAITBAR->Close();
			return nullptr;
		}
		
		return gcnew FITSImage("Mean", (array<double, 2>^)BGWRKR_RESULT, Do_Stats, true);
	}
	else
	{
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;
		double N = (double)L;

		#pragma omp parallel for
		for (int i = 0; i < width; i++)
			for (int j = 0; j < height; j++)
			{
				double sum = 0;
				for (int k = 0; k < L; k++)
				{
					sum = sum + ImageSet[k]->Image[i, j];
				}
				img[i, j] = sum / N;
			}
		return  gcnew FITSImage("c:\\Mean.fits", img, Do_Stats, true);
	}
}

FITSImage^ JPFITS::FITSImageSet::MeanClipped(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, double sigma)
{
	if (ImageSet->CoDimensional == false)
	{
		MessageBox::Show("Can Not Perform Sigma Clipped Mean: Data Stack not Co-Dimensional","Invalid Operation...");
		return nullptr;
	}
	if (ImageSet->Count <= 1)
	{
		MessageBox::Show("Can Not Perform Sigma Clipped Mean: Data Stack Contains One or Fewer Images","Invalid Operation...");
		return nullptr;
	}

	FITSImage^ f = FITSImageSet::Stdv(ImageSet, false, true);
	if (f == nullptr)
		return nullptr;
	array<double, 2>^ simg = f->Image;
	f = FITSImageSet::Mean(ImageSet, false, true);
	if (f == nullptr)
		return nullptr;
	array<double, 2>^ img = f->Image;

	WAITBAR = gcnew JPWaitBar::WaitBar();
	WAITBAR->ProgressBar->Maximum = 100;
	WAITBAR->Text = "Computing Sigma Clipped Mean Data Stack Image";
	WAITBAR->CancelBtn->Text = "Stop Iterating (Max = 2000)";
	array<Object^>^ arg = gcnew array<Object^>{ImageSet,"SCM",sigma,simg,img};
	BGWRKR->RunWorkerAsync(arg);
	WAITBAR->ShowDialog();
	return gcnew FITSImage("c:\\ClippedMean.fits", (array<double, 2>^)BGWRKR_RESULT, Do_Stats, true);
}

FITSImage^ JPFITS::FITSImageSet::Max(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar)
{
	if (ImageSet->CoDimensional == false)
	{
		MessageBox::Show("Can Not Perform Mean: Data Stack not Co-Dimensional","Invalid Operation...");
		return nullptr;
	}
	if (ImageSet->Count <= 1)
	{
		MessageBox::Show("Can Not Perform Mean: Data Stack Contains One or Fewer Images","Invalid Operation...");
		return nullptr;
	}

	if (Show_Waitbar)
	{
		WAITBAR = gcnew JPWaitBar::WaitBar();
		WAITBAR->ProgressBar->Maximum = 100;
		WAITBAR->Text = "Computing Max Data Stack Image";
		array<Object^>^ arg = gcnew array<Object^>{ImageSet, "Max"};
		BGWRKR->RunWorkerAsync(arg);
		WAITBAR->ShowDialog();
		if (WAITBAR->DialogResult == ::DialogResult::OK)
		{
			return gcnew FITSImage("Max", (array<double, 2>^)BGWRKR_RESULT, Do_Stats, true);
		}
		else
		{
			WAITBAR->Close();
			return nullptr;
		}
	}
	else
	{
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;

		#pragma omp parallel for
		for (int i = 0; i < width; i++)
			for (int j = 0; j < height; j++)
			{
				double max = Double::MinValue;
				for (int k = 0; k < L; k++)
					if (ImageSet[k]->Image[i, j] > max)
						max = ImageSet[k]->Image[i, j];
				img[i, j] = max;
			}
		return gcnew FITSImage("c:\\Max.fits", img, Do_Stats, true);
	}
}

FITSImage^ JPFITS::FITSImageSet::Min(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar)
{
	if (ImageSet->CoDimensional == false)
	{
		MessageBox::Show("Can Not Perform Mean: Data Stack not Co-Dimensional","Invalid Operation...");
		return nullptr;
	}
	if (ImageSet->Count <= 1)
	{
		MessageBox::Show("Can Not Perform Mean: Data Stack Contains One or Fewer Images","Invalid Operation...");
		return nullptr;
	}


	if (Show_Waitbar)
	{
		WAITBAR = gcnew JPWaitBar::WaitBar();
		WAITBAR->ProgressBar->Maximum = 100;
		WAITBAR->Text = "Computing Min Data Stack Image";
		array<Object^>^ arg = gcnew array<Object^>{ImageSet, "Min"};
		BGWRKR->RunWorkerAsync(arg);
		WAITBAR->ShowDialog();
		if (WAITBAR->DialogResult == ::DialogResult::OK)
		{
			return gcnew FITSImage("Min", (array<double, 2>^)BGWRKR_RESULT, Do_Stats, true);
		}
		else
		{
			WAITBAR->Close();
			return nullptr;
		}
	}
	else
	{
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;

		#pragma omp parallel for
		for (int i = 0; i < width; i++)
			for (int j = 0; j < height; j++)
			{
				double min = Double::MaxValue;
				for (int k = 0; k < L; k++)
					if (ImageSet[k]->Image[i, j] < min)
						min = ImageSet[k]->Image[i, j];
				img[i, j] = min;
			}
		return gcnew FITSImage("c:\\Min.fits", img, Do_Stats, true);
	}
}

FITSImage^ JPFITS::FITSImageSet::Sum(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar)
{
	if (ImageSet->CoDimensional == false)
	{
		MessageBox::Show("Can Not Perform Sum: Data Stack not Co-Dimensional","Invalid Operation...");
		return nullptr;
	}

	if (ImageSet->Count <= 1)
	{
		MessageBox::Show("Can Not Perform Sum: Data Stack Contains One or Fewer Images","Invalid Operation...");
		return nullptr;
	}

	if (Show_Waitbar)
	{
		WAITBAR = gcnew JPWaitBar::WaitBar();
		WAITBAR->ProgressBar->Maximum = 100;
		WAITBAR->Text = "Computing Summed Data Stack Image";
		array<Object^>^ arg = gcnew array<Object^>{ImageSet, "Sum"};
		BGWRKR->RunWorkerAsync(arg);
		WAITBAR->ShowDialog();
		if (WAITBAR->DialogResult == ::DialogResult::OK)
		{
			return gcnew FITSImage("Sum", (array<double, 2>^)BGWRKR_RESULT, Do_Stats, true);
		}
		else
			return nullptr;
	}
	else
	{
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;

		#pragma omp parallel for
		for (int i = 0; i < width; i++)
			for (int j = 0; j < height; j++)
			{
				double sum = 0;
				for (int k = 0; k < L; k++)
					sum = sum + ImageSet[k]->Image[i, j];
				img[i, j] = sum;
			}
		return gcnew FITSImage("c:\\Sum.fits", img, Do_Stats, true);
	}
}

FITSImage^ JPFITS::FITSImageSet::Quadrature(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar)
{
	if (ImageSet->CoDimensional == false)
	{
		MessageBox::Show("Can Not Perform Quadrature Sum: Data Stack not Co-Dimensional", "Invalid Operation...");
		return nullptr;
	}

	if (ImageSet->Count <= 1)
	{
		MessageBox::Show("Can Not Perform Quadrature Sum: Data Stack Contains One or Fewer Images", "Invalid Operation...");
		return nullptr;
	}

	if (Show_Waitbar)
	{
		WAITBAR = gcnew JPWaitBar::WaitBar();
		WAITBAR->ProgressBar->Maximum = 100;
		WAITBAR->Text = "Computing Quadrature Summed Data Stack Image";
		array<Object^>^ arg = gcnew array<Object^>{ImageSet, "Quadrature"};
		BGWRKR->RunWorkerAsync(arg);
		WAITBAR->ShowDialog();
		if (WAITBAR->DialogResult == ::DialogResult::OK)
		{
			return gcnew FITSImage("Quadrature", (array<double, 2>^)BGWRKR_RESULT, Do_Stats, true);
		}
		else
			return nullptr;
	}
	else
	{
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;

		#pragma omp parallel for
		for (int i = 0; i < width; i++)
			for (int j = 0; j < height; j++)
			{
				double quadsum = 0;
				for (int k = 0; k < L; k++)
					quadsum += (ImageSet[k]->Image[i, j] * ImageSet[k]->Image[i, j]);
				img[i, j] = Math::Sqrt(quadsum);
			}
		return gcnew FITSImage("c:\\Quadrature.fits", img, Do_Stats, true);
	}
}

FITSImage^ JPFITS::FITSImageSet::Median(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar, String^ waitbar_message)
{
	if (ImageSet->CoDimensional == false)
	{
		MessageBox::Show("Can Not Perform Median: Data Stack not Co-Dimensional","Invalid Operation...");
		return nullptr;
	}
	if (ImageSet->Count <= 1)
	{
		MessageBox::Show("Can Not Perform Median: Data Stack Contains One or Fewer Images","Invalid Operation...");
		return nullptr;
	}

	if (Show_Waitbar)
	{
		WAITBAR = gcnew JPWaitBar::WaitBar();
		WAITBAR->ProgressBar->Maximum = 100;
		WAITBAR->Text = "Computing Median Data Stack Image";
		array<Object^>^ arg = gcnew array<Object^>{ImageSet, "Median", waitbar_message};
		BGWRKR->RunWorkerAsync(arg);
		WAITBAR->ShowDialog();
		if (WAITBAR->DialogResult == ::DialogResult::OK)
		{
			return gcnew FITSImage("Median", (array<double, 2>^)BGWRKR_RESULT, Do_Stats, true);
		}
		else
			return nullptr;
	}
	else
	{
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;

		#pragma omp parallel for
		for (int i = 0; i < width; i++)
			for (int j = 0; j < height; j++)
			{
				array<double, 1>^ medarray = gcnew array<double, 1>(L);
				for (int k = 0; k < L; k++)
					medarray[k] = ImageSet[k]->Image[i, j];
				img[i, j] = JPMath::Median(medarray);
			}
		return gcnew FITSImage("c:\\Median.fits", img, Do_Stats, true);
	}
}

FITSImage^ JPFITS::FITSImageSet::Stdv(JPFITS::FITSImageSet^ ImageSet, bool Do_Stats, bool Show_Waitbar)
{
	if (ImageSet->CoDimensional == false)
	{
		MessageBox::Show("Can Not Perform Standard Deviation: Data Stack not Co-Dimensional","Invalid Operation...");
		return nullptr;
	}
	if (ImageSet->Count <= 1)
	{
		MessageBox::Show("Can Not Perform Standard Deviation: Data Stack Contains One or Fewer Images","Invalid Operation...");
		return nullptr;
	}

	if (Show_Waitbar)
	{
		WAITBAR = gcnew JPWaitBar::WaitBar();
		WAITBAR->ProgressBar->Maximum = 100;
		WAITBAR->Text = "Computing Stdv Data Stack Image";
		array<Object^>^ arg = gcnew array<Object^>{ImageSet, "Stdv"};
		BGWRKR->RunWorkerAsync(arg);
		WAITBAR->ShowDialog();
		if (WAITBAR->DialogResult == ::DialogResult::OK)
		{
			return gcnew FITSImage("Stdv", (array<double, 2>^)BGWRKR_RESULT, Do_Stats, true);
		}
		else
			return nullptr;
	}
	else
	{
		array<double, 2>^ img = gcnew array<double, 2>(ImageSet[0]->Width, ImageSet[0]->Height);
		int width = img->GetLength(0);
		int height = img->GetLength(1);
		int L = ImageSet->Count;
		double N = double(L);

		#pragma omp parallel for
		for (int i = 0; i < width; i++)
			for (int j = 0; j < height; j++)
			{
				double std = 0;
				double mean = 0;
				for (int k = 0; k < L; k++)
					mean = mean + ImageSet[k]->Image[i, j];
				mean = mean / N;
				for (int q = 0; q < L; q++)
					std = std + (ImageSet[q]->Image[i, j] - mean)*(ImageSet[q]->Image[i, j] - mean);
				std = Math::Sqrt(std / (N - 1.0));
				img[i, j] = std;
			}
		return gcnew FITSImage("c:\\Stdv.fits", img, Do_Stats, true);
	}
}

void JPFITS::FITSImageSet::Register(JPFITS::FITSImageSet ^ImageSet, int RefImgIndex, bool Do_Stats)
{
	WAITBAR = gcnew JPWaitBar::WaitBar();
	WAITBAR->ProgressBar->Maximum = 100;
	WAITBAR->Text = "Auto-Registering Images";
	array<Object^>^ arg = gcnew array<Object^>{ImageSet,"AutoReg",RefImgIndex, Do_Stats};
	BGWRKR->RunWorkerAsync(arg);
	WAITBAR->ShowDialog();
}

void JPFITS::FITSImageSet::CHECK_CODIMENSIONAL()
{
	CODIMENSIONAL = true;
	for (int i = 1; i < FITSLIST->Count; i++)
		if (((FITSImage^)(FITSLIST[0]))->Width != ((FITSImage^)(FITSLIST[i]))->Width || ((FITSImage^)(FITSLIST[0]))->Height != ((FITSImage^)(FITSLIST[i]))->Height)
			CODIMENSIONAL = false;
}

void JPFITS::FITSImageSet::GatherHeaders(JPFITS::FITSImageSet^ FITS_Set, FITSImage^ FITS_destination)
{
	FITS_destination->Header->RemoveAllKeys(FITS_destination->Image);

	bool skip = false;
	for (int j = 0; j < FITS_Set[0]->Header->HeaderKeys->Length; j++)
	{
		skip = false;
		String^ key = FITS_Set[0]->Header->GetKeyName(j);
		String^ val = FITS_Set[0]->Header->GetKeyValue(j);
		String^ com = FITS_Set[0]->Header->GetKeyComment(j);

		if (!JPFITS::FITSImageHeader::VALIDKEYEDIT(key)/*FITS_Set[0]->Header->VALIDKEYEDIT(key)*/)
			continue;

		for (int i = 1; i < FITS_Set->Count; i++)
			if (FITS_Set[i]->Header->GetKeyIndex(key, val, com) == -1)
			{
				skip = true;
				break;
			}
		
		if (!skip)
			if (FITS_Set[0]->Header->HeaderLineIsComment[j])
				FITS_destination->Header->AddCommentKeyLine(com, -1);
			else
				FITS_destination->Header->AddKey(key, val, com, -1);
	}
}

void JPFITS::FITSImageSet::GatherHeaders(array<String^>^ filenames, JPFITS::FITSImage^ FITS_destination)
{
	JPFITS::FITSImageSet^ mergeset = gcnew JPFITS::FITSImageSet();
	for (int i = 0; i < filenames->Length; i++)
		mergeset->Add(gcnew JPFITS::FITSImage(filenames[i], nullptr, true, false, false, false));

	JPFITS::FITSImageSet::GatherHeaders(mergeset, FITS_destination);
}

int JPFITS::FITSImageSet::Sort(String^ key)
{
	if (key == "filename")//filenames are nice because they are always unique
	{
		array<String^>^ keys = gcnew array<String^>(FITSLIST->Count);

		for (int i = 0; i < FITSLIST->Count; i++)
			keys[i] = ((FITSImage^)(FITSLIST[i]))->FullFileName;

		Array::Sort(keys);

		for (int i = 0; i < FITSLIST->Count; i++)
		{
			if (keys[i] == ((FITSImage^)(FITSLIST[i]))->FullFileName)
				continue;

			for (int j = i+1; j < FITSLIST->Count; j++)
				if (keys[i] == ((FITSImage^)(FITSLIST[j]))->FullFileName)
				{
					FITSImage^ tempfits = (FITSImage^)FITSLIST[i];
					FITSLIST[i] = FITSLIST[j];
					FITSLIST[j] = tempfits;
				}
		}
		return 0;
	}

	//else use header key value, returned already if for filename
	//key values aren't as nice because they might not always be unique

	//check for either numeric or alphabetical case
	bool numeric = true;
	try
	{
		double d = Convert::ToDouble(((FITSImage^)(FITSLIST[0]))->Header->GetKeyValue(key));
	}
	catch (...)
	{
		numeric = false;
	}

	if (!numeric)
	{
		array<String^>^ keys = gcnew array<String^>(FITSLIST->Count);

		bool keycheck = true;

		for (int i = 0; i < FITSLIST->Count; i++)
		{
			keys[i] = ((FITSImage^)(FITSLIST[i]))->Header->GetKeyValue(key);

			if (keys[i] == "" && keycheck)
			{
				if (MessageBox::Show("Key not found in at least one FITS header: continue?", "Warning", MessageBoxButtons::YesNo) == DialogResult::No)
					return -1;
				keycheck = false;
			}
		}

		Array::Sort(keys);

		for (int i = 0; i < FITSLIST->Count; i++)
		{
			if (keys[i] == ((FITSImage^)(FITSLIST[i]))->Header->GetKeyValue(key))
				continue;

			for (int j = i + 1; j < FITSLIST->Count; j++)
				if (keys[i] == ((FITSImage^)(FITSLIST[j]))->Header->GetKeyValue(key))
				{
					FITSImage^ tempfits = (FITSImage^)FITSLIST[i];
					FITSLIST[i] = FITSLIST[j];
					FITSLIST[j] = tempfits;
				}
		}
		return 0;
	}
	else//numeric
	{
		array<double>^ keys = gcnew array<double>(FITSLIST->Count);

		bool keycheck = true;

		for (int i = 0; i < FITSLIST->Count; i++)
		{
			String^ k = ((FITSImage^)(FITSLIST[i]))->Header->GetKeyValue(key);

			if (k == "" && keycheck)
			{
				if (MessageBox::Show("Key not found in at least one FITS header: continue?", "Warning", MessageBoxButtons::YesNo) == DialogResult::No)
					return -1;
				keycheck = false;
			}

			try
			{
				keys[i] = Convert::ToDouble(((FITSImage^)(FITSLIST[i]))->Header->GetKeyValue(key));
			}
			catch (...)
			{
				MessageBox::Show("Tried to convert what had been numeric key values for the sorting keys.  Check FITSImageSet index (1-based) " + (i+1).ToString() + ".", "Sorting Error");
				return -1;
			}
		}

		Array::Sort(keys);

		for (int i = 0; i < FITSLIST->Count; i++)
		{
			if (keys[i] == Convert::ToDouble(((FITSImage^)(FITSLIST[i]))->Header->GetKeyValue(key)))
				continue;

			for (int j = i + 1; j < FITSLIST->Count; j++)
				if (keys[i] == Convert::ToDouble(((FITSImage^)(FITSLIST[j]))->Header->GetKeyValue(key)))
				{
					FITSImage^ tempfits = (FITSImage^)FITSLIST[i];
					FITSLIST[i] = FITSLIST[j];
					FITSLIST[j] = tempfits;
				}
		}
		return 0;
	}
}

String^ JPFITS::FITSImageSet::GetCommonDirectory()
{
	return FITSImageSet::GetCommonDirectory(this->FilePaths);
}

String^ JPFITS::FITSImageSet::GetCommonDirectory(array<String^>^ filelist)
{
	String^ first = Path::GetDirectoryName(filelist[0]) + "\\";
	for (int i = 1; i < filelist->Length; i++)
	{
		String^ second = filelist[i];
		int N = first->Length;
		for (int j = 0; j < second->Length; j++)
		{
			if (j == N)
				break;
			if (first[j] != second[j])
			{
				first = first->Substring(0, j);
				first = first->Substring(0, first->LastIndexOf("\\"));
				break;
			}
		}
	}
	return first;
}

