
#include "LambdaMoments.h"

// ROOT headers
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "Math/IFunction.h"
#include "TGraphErrors.h"

// CPP headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

//function definitions

void Emptying1DArrayInt(Int_t *array, Int_t raws)
{
	for(Int_t i = 0; i < raws; i++) array[i] = 0;
}

void Emptying1DArray(Double_t *array, Int_t raws)
{
	for(Int_t i = 0; i < raws; i++) array[i] = 0;
}

void Emptying2DArray(Double_t **array, Int_t raws, Int_t columns)
{
 	for(Int_t i = 0; i < raws; i++)
 		for(Int_t j = 0; j < columns; j++) array[i][j] = 0.;
}

void Emptying3DArray(Double_t ***array, Int_t dim1, Int_t dim2, Int_t dim3)
{
	for(Int_t i = 0; i < dim1; i++)
 		for(Int_t j = 0; j < dim2; j++)
 			for(Int_t k = 0; k < dim3; k++) array[i][j][k] = 0.;
}

void DeleteDynamic1DArrayInt(Int_t *Array, TString ArrayName)
{
  delete[] Array;
}

void DeleteDynamic1DArray(Double_t *Array, TString ArrayName)
{
  delete[] Array;
}


void DeleteDynamic2DArray(Double_t **Array, Int_t rows, TString ArrayName)
{
 	{
 		for(Int_t i = 0; i < rows; i++) delete[] Array[i];
 		delete[] Array;
 	}
}

void DeleteDynamic3DArray(Double_t ***Array, Int_t rows, Int_t columns, TString ArrayName)
{
  for(Int_t i = 0; i < rows; i++)
  {
    for(Int_t j = 0; j < columns; j++)
    {
      delete[] Array[i][j];
    }
    delete[] Array[i];
  }
  delete[] Array;
}

void DeleteDynamic5DArray(Double_t *****Array, Int_t Dim_1, Int_t Dim_2, Int_t Dim_3, Int_t Dim_4, TString ArrayName)
{
  for(Int_t i = 0; i < Dim_1; i++)
  {
    for(Int_t j = 0; j < Dim_2; j++)
    {
      for(Int_t k = 0; k < Dim_3; k++)
      {
        for(Int_t l = 0; l < Dim_4; l++)
        {
          delete[] Array[i][j][k][l];
        }
        delete[] Array[i][j][k];
      }
      delete[] Array[i][j];
    }
    delete[] Array[i];
  }
  delete[] Array;
}


Int_t S_1(Int_t n, Int_t i)
{
  if      (n < i)           {s_1_n_i = 0;}
  else if (n == i)          {s_1_n_i = 1;}
  else if (n > i && i == 0) {s_1_n_i = 0;}
  else                      {s_1_n_i = S_1(n-1, i-1) - (n-1)*S_1(n-1, i);}

  return s_1_n_i;
}


Int_t S_2(Int_t n, Int_t i)
{
  if      (n < i)           {s_2_n_i = 0;}
  else if (n == i)          {s_2_n_i = 1;}
  else if (n > i && i == 0) {s_2_n_i = 0;}
  else                      {s_2_n_i = S_2(n-1, i-1) + i*S_2(n-1, i);}

  return s_2_n_i;
}


long double FACT(Int_t Num)
{
  long double FacNum = 1;
  for (Int_t i = 0; i < Num; i++) FacNum *= (Num - i);
  return FacNum;
}


void CentralMomCalcRefMult(Int_t *BinEntriesRefMult, Double_t **CentralMomRefMult)
{
	/*======================================= Central moments calculation =============================*/

  TFile *LambdaMom_file = new TFile(InRootFile, "read");

  TProfile *NetL_moment_1 = (TProfile*) LambdaMom_file->Get(hh_M_1);
  TProfile *NetL_moment_2 = (TProfile*) LambdaMom_file->Get(hh_M_2);
  TProfile *NetL_moment_3 = (TProfile*) LambdaMom_file->Get(hh_M_3);
  TProfile *NetL_moment_4 = (TProfile*) LambdaMom_file->Get(hh_M_4);
  TProfile *NetL_moment_5 = (TProfile*) LambdaMom_file->Get(hh_M_5);
  TProfile *NetL_moment_6 = (TProfile*) LambdaMom_file->Get(hh_M_6);
  TProfile *NetL_moment_7 = (TProfile*) LambdaMom_file->Get(hh_M_7);
  TProfile *NetL_moment_8 = (TProfile*) LambdaMom_file->Get(hh_M_8);

  Double_t M_1 = 0.0, M_2 = 0.0, M_3 = 0.0, M_4 = 0.0, M_6 = 0.0, M_5 = 0.0, M_7 = 0.0, M_8 = 0.0;
  Double_t m1 =0.0, m2 = 0.0, m3 = 0.0, m4 = 0.0, m5 = 0.0, m6 = 0.0, m7 = 0.0, m8 = 0.0;

  Int_t BinEntries = 0;

  for(Int_t i = 0; i < 1000; i++)
  {
  	BinEntriesRefMult[i] = 0;
  	for(Int_t j = 0; j < 8; j++)
  	{
  		CentralMomRefMult[i][j] = 0;
  	}
  }

  //cout << " RefMult3 " << "   \tCm_1 " << "   \tCm_2 " << "   \tCm_3 " << "   \tCm_4 " << "   \tCm_5 " << "   \tCm_6 " << "   \tCm_7 " << "   \tCm_8 " << "\n " << endl;
  for(Int_t RefMultBin = 0; RefMultBin < 1000; RefMultBin++)
  {

  	M_1        = NetL_moment_1 -> GetBinContent(RefMultBin);
  	BinEntries = NetL_moment_1 -> GetBinEntries(RefMultBin);
  	BinEntriesRefMult[RefMultBin] = BinEntries;

  	M_2 = NetL_moment_2 -> GetBinContent(RefMultBin);
  	M_3 = NetL_moment_3 -> GetBinContent(RefMultBin);
  	M_4 = NetL_moment_4 -> GetBinContent(RefMultBin);
  	M_5 = NetL_moment_5 -> GetBinContent(RefMultBin);
  	M_6 = NetL_moment_6 -> GetBinContent(RefMultBin);
  	M_7 = NetL_moment_7 -> GetBinContent(RefMultBin);
  	M_8 = NetL_moment_8 -> GetBinContent(RefMultBin);

    // Calculation of central moments.
	m1 = M_1;
  	m2 = M_2  -  TMath::Power(M_1, 2);
  	m3 = M_3  -  3 * M_1 * m2  -   TMath::Power(M_1, 3);
  	m4 = M_4  -  4 * M_1 * m3  -   6 * TMath::Power(M_1, 2) * m2   -  TMath::Power(M_1, 4);
  	m5 = M_5  -  5 * M_1 * m4  -  10 * TMath::Power(M_1, 2) * m3   -  10 * TMath::Power(M_1, 3) * m2  -  TMath::Power(M_1, 5);
  	m6 = M_6  -  6 * M_1 * m5  -  15 * TMath::Power(M_1, 2) * m4   -  20 * TMath::Power(M_1, 3) * m3  -  15 * TMath::Power(M_1, 4) * m2  -  TMath::Power(M_1, 6);
  	m7 = M_7  -  7 * M_1 * m6  -  21 * TMath::Power(M_1, 2) * m5   -  35 * TMath::Power(M_1, 3) * m4  -  35 * TMath::Power(M_1, 4) * m3  -  21 * TMath::Power(M_1, 5) * m2  -  TMath::Power(M_1, 7);
        m8 = M_8  -  8 * M_1 * m7  -  28 * TMath::Power(M_1, 2) * m6   -  56 * TMath::Power(M_1, 3) * m5  -  70 * TMath::Power(M_1, 4) * m4  -  56 * TMath::Power(M_1, 5) * m3  -  28 * TMath::Power(M_1, 6) * m2 - TMath::Power(M_1, 8);

    //cout << " " << RefMultBin << " \t"<< m1 <<"     \t"<< m2 <<"     \t"<< m3 <<"     \t"<< m4 <<"     \t"<< m5 <<"     \t"<< m6 <<"     \t"<< m7 <<"     \t"<< m8 << endl;

    CentralMomRefMult[RefMultBin][0] = m1;
    CentralMomRefMult[RefMultBin][1] = m2;
    CentralMomRefMult[RefMultBin][2] = m3;
    CentralMomRefMult[RefMultBin][3] = m4;
    CentralMomRefMult[RefMultBin][4] = m5;
    CentralMomRefMult[RefMultBin][5] = m6;
    CentralMomRefMult[RefMultBin][6] = m7;
    CentralMomRefMult[RefMultBin][7] = m8;
  }

	/*======================================= Central Moment calculation ends   =============================*/
	LambdaMom_file -> Close();
}


void CBWC_Central_mom(Double_t **mu, Double_t **cbwc_mu, Int_t *events, Int_t CentBins, Int_t RefMultLow, Int_t RefMultHigh, Int_t Type)
{
	Int_t RefMultBinWidth = 0;
	if(Type == 0)
	{
		RefMultBinWidth = (Int_t) ((RefMultHigh - RefMultLow)/CentBins);
		cout << CentBins << " " << RefMultLow << " " << RefMultHigh << " " << RefMultBinWidth << endl;
	}

	Int_t binLow  = 0;
	Int_t binHigh = 0;

	if(Type == 0){binLow  = RefMultLow; binHigh = RefMultLow + RefMultBinWidth;}
	else{binLow  = RefMult3[0]; binHigh = RefMult3[1];}

	Double_t Sum_MU_1 = 0., Sum_MU_2 = 0., Sum_MU_3 = 0., Sum_MU_4 = 0., Sum_MU_5 = 0., Sum_MU_6 = 0., Sum_MU_7 = 0., Sum_MU_8 = 0.;
	Int_t eventsSum = 0;

	for(Int_t i = 0; i < CentBins; i++)
	{
		for(Int_t j = binLow; j < binHigh; j++)
		{
			Sum_MU_1 += mu[j][0]*events[j];
			Sum_MU_2 += mu[j][1]*events[j];
			Sum_MU_3 += mu[j][2]*events[j];
			Sum_MU_4 += mu[j][3]*events[j];
			Sum_MU_5 += mu[j][4]*events[j];
			Sum_MU_6 += mu[j][5]*events[j];
			Sum_MU_7 += mu[j][6]*events[j];
			Sum_MU_8 += mu[j][7]*events[j];
			eventsSum += events[j];

		}
		cbwc_mu[i][0] =  Sum_MU_1/eventsSum;
		cbwc_mu[i][1] =  Sum_MU_2/eventsSum;
		cbwc_mu[i][2] =  Sum_MU_3/eventsSum;
		cbwc_mu[i][3] =  Sum_MU_4/eventsSum;
		cbwc_mu[i][4] =  Sum_MU_5/eventsSum;
		cbwc_mu[i][5] =  Sum_MU_6/eventsSum;
		cbwc_mu[i][6] =  Sum_MU_7/eventsSum;
		cbwc_mu[i][7] =  Sum_MU_8/eventsSum;

		Sum_MU_1 = 0.; Sum_MU_2 = 0.; Sum_MU_3 = 0.; Sum_MU_4 = 0.; Sum_MU_5 = 0.; Sum_MU_6 = 0.; Sum_MU_7 = 0.; Sum_MU_8 = 0.;
		eventsSum = 0;

		if(Type == 0){binLow  += RefMultBinWidth; binHigh += RefMultBinWidth;}
		else{binLow  = RefMult3[i+1]; binHigh = RefMult3[i+2];}
	}
}


void Moments(Double_t **mu, Double_t **moments, Int_t RefMultLow, Int_t RefMultHigh)
{
	for(Int_t i = RefMultLow; i < RefMultHigh; i++)
	{
		if(mu[i][1] > 0)
		{
			moments[i][0] = mu[i][0];					//Mean.
			moments[i][1] = TMath::Sqrt(mu[i][1]);			        //Standared deviation.
			moments[i][2] = mu[i][2]/TMath::Power(mu[i][1], 1.5);		//Skewness.
			moments[i][3] = (mu[i][3]/TMath::Power(mu[i][1], 2)) - 3;	//Kurtosis.
		}
	}
}


void CBWC_Moments(Double_t **moments, Double_t **cbwc_moments, Int_t *events, Int_t CentBins, Int_t RefMultLow, Int_t RefMultHigh, Int_t Type)
{
	Int_t RefMultBinWidth = 0;
	if(Type == 0)
	{
		RefMultBinWidth = (Int_t) ((RefMultHigh - RefMultLow)/CentBins);
		cout << CentBins << " " << RefMultLow << " " << RefMultHigh << " " << RefMultBinWidth << endl;
	}

	Int_t binLow  = 0;
	Int_t binHigh = 0;

	if(Type == 0){binLow  = RefMultLow; binHigh = RefMultLow + RefMultBinWidth;}
	else{binLow  = RefMult3[0]; binHigh = RefMult3[1];}

	Double_t Sum_Mom_1 = 0., Sum_Mom_2 = 0., Sum_Mom_3 = 0., Sum_Mom_4 = 0.;
	Int_t eventsSum = 0;

	for(Int_t i = 0; i < CentBins; i++)
	{
		for(Int_t j = binLow; j < binHigh; j++)
		{
			Sum_Mom_1 += moments[j][0]*events[j];
			Sum_Mom_2 += moments[j][1]*events[j];
			Sum_Mom_3 += moments[j][2]*events[j];
			Sum_Mom_4 += moments[j][3]*events[j];

			eventsSum += events[j];
		}
		cbwc_moments[i][0] =  Sum_Mom_1/eventsSum;
		cbwc_moments[i][1] =  Sum_Mom_2/eventsSum;
		cbwc_moments[i][2] =  Sum_Mom_3/eventsSum;
		cbwc_moments[i][3] =  Sum_Mom_4/eventsSum;

		Sum_Mom_1 = 0.; Sum_Mom_2 = 0.; Sum_Mom_3 = 0.; Sum_Mom_4 = 0.;
		eventsSum = 0;

		if(Type == 0){binLow  += RefMultBinWidth; binHigh += RefMultBinWidth;}
		else{binLow  = RefMult3[i+1]; binHigh = RefMult3[i+2];}
	}
}

void Cumulants(Double_t **mu, Double_t **cumulants, Int_t RefMultLow, Int_t RefMultHigh)
{
	for(Int_t i = RefMultLow; i < RefMultHigh; i++)
	{
		//if(mu[i][1] > 0)
		if(mu[i][0] > 0 && mu[i][1] > 0)
		{
			cumulants[i][0] = mu[i][0];					//C_1.
			cumulants[i][1] = mu[i][1];				        //C_2.
			cumulants[i][2] = mu[i][2];					//C_3.
			cumulants[i][3] = mu[i][3] - 3*TMath::Power(mu[i][1], 2);	//C_4.
			cumulants[i][4] = cumulants[i][1]/cumulants[i][0];		//C_2/C_1
			cumulants[i][5] = cumulants[i][2]/cumulants[i][1];		//C_3/C_2
			cumulants[i][6] = cumulants[i][3]/cumulants[i][1];		//C_4/C_2
		}
	}
}


void CBWC_Cumulants(Double_t **cumulants, Double_t **cbwc_cumulants, Int_t *events, Int_t CentBins, Int_t RefMultLow, Int_t RefMultHigh, Int_t Type)
{
	Int_t RefMultBinWidth = 0;
	if(Type == 0)
	{
		RefMultBinWidth = (Int_t) ((RefMultHigh - RefMultLow)/CentBins);
		cout << CentBins << " " << RefMultLow << " " << RefMultHigh << " " << RefMultBinWidth << endl;
	}

	Int_t binLow  = 0;
	Int_t binHigh = 0;

	if(Type == 0){binLow  = RefMultLow; binHigh = RefMultLow + RefMultBinWidth;}
	else{binLow  = RefMult3[0]; binHigh = RefMult3[1];}

	Double_t Sum_Cum_1 = 0., Sum_Cum_2 = 0., Sum_Cum_3 = 0., Sum_Cum_4 = 0., Sum_Cum_5 = 0., Sum_Cum_6 = 0., Sum_Cum_7 = 0.;
	Int_t eventsSum = 0;

	for(Int_t i = 0; i < CentBins; i++)
	{
		for(Int_t j = binLow; j < binHigh; j++)
		{
			Sum_Cum_1 += cumulants[j][0]*events[j];
			Sum_Cum_2 += cumulants[j][1]*events[j];
			Sum_Cum_3 += cumulants[j][2]*events[j];
			Sum_Cum_4 += cumulants[j][3]*events[j];
			Sum_Cum_5 += cumulants[j][4]*events[j];
			Sum_Cum_6 += cumulants[j][5]*events[j];
			Sum_Cum_7 += cumulants[j][6]*events[j];

			eventsSum += events[j];
		}
		cbwc_cumulants[i][0] =  Sum_Cum_1/eventsSum;
		cbwc_cumulants[i][1] =  Sum_Cum_2/eventsSum;
		cbwc_cumulants[i][2] =  Sum_Cum_3/eventsSum;
		cbwc_cumulants[i][3] =  Sum_Cum_4/eventsSum;
		cbwc_cumulants[i][4] =  Sum_Cum_5/eventsSum;
		cbwc_cumulants[i][5] =  Sum_Cum_6/eventsSum;
		cbwc_cumulants[i][6] =  Sum_Cum_7/eventsSum;

		Sum_Cum_1 = 0.; Sum_Cum_2 = 0.; Sum_Cum_3 = 0.; Sum_Cum_4 = 0.; Sum_Cum_5 = 0.; Sum_Cum_6 = 0.; Sum_Cum_7 = 0.;
		eventsSum = 0;

		if(Type == 0){binLow  += RefMultBinWidth; binHigh += RefMultBinWidth;}
		else{binLow  = RefMult3[i+1]; binHigh = RefMult3[i+2];}
	}
}


void FactorialCalc(Double_t ***f, Int_t Dim_1, Int_t Dim_2, Int_t Dim_3)
{
  TFile *LambdaCent_file      = new TFile(InRootFile, "read");

  TH3D *h_La_vs_antiLa_Cent   = (TH3D*) LambdaCent_file->Get(hh_La_vs_antiLa_Cent);
  TH1D *h_Cent                = (TH1D*) LambdaCent_file->Get(hh_Cent);

  h_La_vs_antiLa_Cent ->GetXaxis()->SetRangeUser(0., 10.0);
  h_La_vs_antiLa_Cent ->GetYaxis()->SetRangeUser(0., 10.0);

  Int_t Tot_N_of_Entries = h_La_vs_antiLa_Cent ->GetEntries();
  Double_t P_Array[9][9][9];

  for(Int_t CentBin = 0; CentBin < Dim_1; CentBin++)
  {
    Int_t Events_in_Cent = h_Cent ->GetBinContent(CentBin+1);

    for(Int_t bin_Lambda = 1; bin_Lambda < 10; bin_Lambda++)
    {
      for(Int_t bin_aLambda = 1; bin_aLambda < 10; bin_aLambda++)
      {
        P_Array[CentBin][bin_Lambda -1][bin_aLambda -1] = (h_La_vs_antiLa_Cent -> GetBinContent(bin_Lambda,bin_aLambda, CentBin+2))/Events_in_Cent;
      }
    }

    for(Int_t i = 0; i < Dim_2; i++)
    {
      for(Int_t j = 0; j < Dim_3; j++)
      {
        f[CentBin][i][j] = 0;
        for(Int_t nl = i; nl < 9; nl++)
        {
          for(Int_t nal = j; nal < 9; nal++)
          {
            f[CentBin][i][j] += P_Array[CentBin][nl][nal]*(FACT(nl)/FACT(nl-i))*(FACT(nal)/FACT(nal-j));
          }
        }
      }
    }
  }
  LambdaCent_file -> Close();
}

void EffCorrectFactorial(Double_t ***fac, Double_t ***effCorFac, Double_t *eff, Int_t Dim_1, Int_t Dim_2, Int_t Dim_3)
{
  for(Int_t i = 0; i < Dim_1; i++)
  {
    for(Int_t j = 0; j < Dim_2; j++)
    {
      for(Int_t k = 0; k < Dim_3; k++)
      {
      	effCorFac[i][j][k] = fac[i][j][k]/TMath::Power(eff[i], (j+k));
      }
    }
  }
}

/*_________________ Test for the Fac moments with different eff. for Lambda and anti-Lambda _________________*/

void EffCorrectFactorial_DiffEff(Double_t ***fac, Double_t ***effCorFac_DiffEff, Int_t Dim_1, Int_t Dim_2, Int_t Dim_3)
{
  for(Int_t i = 0; i < Dim_1; i++)
  {
    for(Int_t j = 0; j < Dim_2; j++)
    {
      for(Int_t k = 0; k < Dim_3; k++)
      {
      	effCorFac_DiffEff[i][j][k] = fac[i][j][k]/(TMath::Power(Eff_la[i], j) * TMath::Power(Eff_ala[i], k));
      }
    }
  }
}

void EC_Cumulants_DiffEff(Double_t **ec_cumulants_diffeff)
{
	Int_t OD1 = 9; 				//for the fac_mom - order_1.
	Int_t OD2 = 9; 				//for the fac_mom - order_2.

	Double_t ***FF_DF;    //Pointer to store the factorial moments array.

	//Malloc for ff.
	FF_DF = new Double_t**[9];
	for(Int_t i = 0; i < 9; i++)
	{
		FF_DF[i] = new Double_t*[OD1];
		for(Int_t j = 0; j < OD1; j++) FF_DF[i][j] = new Double_t[OD2];
	}

	Double_t ***ff_DF;    //Pointer to store the factorial moments array.

	//Malloc for ff.
	ff_DF = new Double_t**[9];
	for(Int_t i = 0; i < 9; i++)
	{
		ff_DF[i] = new Double_t*[OD1];
		for(Int_t j = 0; j < OD1; j++) ff_DF[i][j] = new Double_t[OD2];
	}

	FactorialCalc(ff_DF, 9, OD1, OD2);
	EffCorrectFactorial_DiffEff(ff_DF, FF_DF, 9, OD1, OD2);

	Double_t N[9] = {};

	for(Int_t i = 0; i < 9; i++)
	{
		N[i] = FF_DF[i][1][0] + FF_DF[i][0][1];
		ec_cumulants_diffeff[i][0] = FF_DF[i][1][0] - FF_DF[i][0][1];
		ec_cumulants_diffeff[i][1] = N[i] - TMath::Power(ec_cumulants_diffeff[i][0], 2) + FF_DF[i][0][2] - 2*FF_DF[i][1][1] + FF_DF[i][2][0];
                ec_cumulants_diffeff[i][2] = ec_cumulants_diffeff[i][0] + 2*TMath::Power(ec_cumulants_diffeff[i][0], 3) - FF_DF[i][0][3] - 3*FF_DF[i][0][2]
																 + 3*FF_DF[i][1][2] + 3*FF_DF[i][2][0] - 3*FF_DF[i][2][1] + FF_DF[i][3][0]
																 - 3*ec_cumulants_diffeff[i][0]*(N[i] + FF_DF[i][0][2] - 2*FF_DF[i][1][1] + FF_DF[i][2][0]);
  }
}

void Avg_Tot_Part_counting(Double_t *AvgTotPartArray)
{
	/*____________________Average number of total(<n+> + <n->) particles caculation for efficiancy correction.________________*/

  TFile *LambdaCent_file 			= new TFile(InRootFile, "read");

  TH2D *h_Lambdas_Cent 	   	  = (TH2D*) LambdaCent_file->Get(hh_Lambdas_Cent);
  TH2D *h_AntiLambdas_Cent 	  = (TH2D*) LambdaCent_file->Get(hh_AntiLambdas_Cent);
  TH1D *h_Cent            	  = (TH1D*) LambdaCent_file->Get(hh_Cent);

  Double_t LambdaCount 		  = 0;
  Double_t AntiLambdaCount 	  = 0;
  Double_t TotalLambdas 	  = 0;

  for(int j = 2; j < 11; j++)
  {
    for(int i = 1; i < 15; i++) //21
    {
      Double_t LambdaBinContent        = h_Lambdas_Cent     -> GetBinContent(i,j);
      Double_t AntiLambdaBinContent    = h_AntiLambdas_Cent -> GetBinContent(i,j);
      if(i > 1)
      {
        LambdaCount += (i-1)*LambdaBinContent;         // Count lambdas.
        AntiLambdaCount += (i-1)*AntiLambdaBinContent; // Count Anti-Lambdas.
      }
		}
    TotalLambdas = LambdaCount + AntiLambdaCount;

    Double_t NumberOfEvents    = h_Cent -> GetBinContent(j-1);
    Double_t AvarageParticles  = TotalLambdas/NumberOfEvents;

    AvgTotPartArray[j-2] = AvarageParticles;
    LambdaCount = 0;
    AntiLambdaCount = 0;
    TotalLambdas = 0;
  }

  LambdaCent_file -> Close();
}


void CovarianceCalc(Double_t *****Cov, Int_t Dim_1, Int_t Dim_2, Int_t Dim_3, Int_t Dim_4, Int_t Dim_5)
{
  //cout << "\n Covariance calculation starts here .... \n" << endl;

  Int_t OD = 5;      //for the cov - order_1 = order_2 = OD.
  Int_t fac_OD1 = 9; //for the fac_mom - order_1.
  Int_t fac_OD2 = 9; //for the fac_mom - order_2.

  Double_t ***ff;    //Pointer to store the factorial moments array.

  //Malloc for ff.
  ff = new Double_t**[9];
  for(Int_t i = 0; i < 9; i++)
  {
    ff[i] = new Double_t*[fac_OD1];
    for(Int_t j = 0; j < fac_OD1; j++) ff[i][j] = new Double_t[fac_OD2];
  }

  FactorialCalc(ff, 9, fac_OD1, fac_OD2);

  TFile *LambdaCent_file      = new TFile(InRootFile, "read");
  TH1D  *h_Cent               = (TH1D*) LambdaCent_file->Get(hh_Cent);

  Double_t f[9][fac_OD1][fac_OD2];

  for(Int_t i = 0; i < 9; i++)
    for(Int_t j = 0; j < fac_OD1; j++)
      for(Int_t k = 0; k < fac_OD2; k++)
        f[i][j][k] = ff[i][j][k];

  DeleteDynamic3DArray(ff, 9, fac_OD1, "ff");

  for(Int_t i = 0; i < Dim_1; i++)
    for(Int_t j = 0; j < Dim_2; j++)
      for(Int_t k = 0; k < Dim_3; k++)
        for(Int_t l = 0; l < Dim_4; l++)
          for(Int_t m = 0; m < Dim_5; m++)
            Cov[i][j][k][l][m] = 0.;

  Int_t CONFG  = TMath::Power(OD, 4); //number of rsuv configurations.

  Int_t rsuv[CONFG][4];

  Int_t index1 = 0;
  Int_t index2 = 0;
  Int_t index3 = 0;
  Int_t index4 = 0;

  for(Int_t row = 0; row < CONFG; row++)
  {
    if(row%OD == 0) index1 = 0;
    else index1++;

    if(row%OD == 0) index2++;
    if(row%(OD*OD) == 0) index2 = 0;

    if(row%(OD*OD) == 0) index3++;
    if(row%(OD*OD*OD) == 0) index3 = 0;

    if(row%(OD*OD*OD) == 0) index4++;
    if(row%(OD*OD*OD*OD) == 0) index4 = 0;

    for(Int_t column = 0; column < 4; column++)
    {
      if(column == 3)       {rsuv[row][column] = index1;}
      else if (column == 2) {rsuv[row][column] = index2;}
      else if (column == 1) {rsuv[row][column] = index3;}
      else rsuv[row][column] = index4;
    }
  }

  /*
	for(Int_t row = 0; row < CONFG; row++)
  {
    if(row%(OD*OD*OD) == 0) cout << "\n";
    for(Int_t column = 0; column < 4; column++)
    {
      cout << rsuv[row][column] << ",";
    }
    cout << "\n";
  }

  cout << "\n\nCovariance metrix printing... ==>\n" << endl;
	*/

  Double_t CovVal = 0;
  Int_t r = 0;
  Int_t s = 0;
  Int_t u = 0;
  Int_t v = 0;

  for(Int_t CentBin = 0; CentBin < 9; CentBin++)
  {
    long double Events_in_Cent = h_Cent ->GetBinContent(CentBin+1);
    for(Int_t config = 0; config < CONFG; config++)
    {
      r = rsuv[config][0];
      s = rsuv[config][1];
      u = rsuv[config][2];
      v = rsuv[config][3];

      for(Int_t i = 0; i <= r; i++)
      {
        for(Int_t j = 0; j <= s; j++)
        {
          for(Int_t k = 0; k <= u; k++)
          {
            for(Int_t h = 0; h <= v; h++)
            {
              for(Int_t alfa = 0; alfa <= i+k; alfa++)
              {
                for(Int_t beta = 0; beta <= j+h; beta++)
                {
                  CovVal += S_1(r,i)*S_1(s,j)*S_1(u,k)*S_1(v,h)*S_2(i+k,alfa)*S_2(j+h,beta)*f[CentBin][alfa][beta];
                }
              }
            }
          }
        }
      }
      Cov[CentBin][r][s][u][v] = (1/Events_in_Cent)*(CovVal - (f[CentBin][r][s] * f[CentBin][u][v]));
      CovVal = 0;
      //cout << "  Cov("<< rsuv[config][0] << rsuv[config][1] << rsuv[config][2] << rsuv[config][3] << ") = "
      //<< Cov[CentBin][r][s][u][v] << " \n ";
    }
  }

  /*
	cout << "\n\n\n **** Manual calculation for Cov(f_10, f_10) for the cross-check.. starts here .... \n" << endl;

  Double_t CovCrosschk1010  = 0.;

  for(Int_t i = 0; i < 9; i++)
  {
    cout << "\n__________ " << CentName[8 - i] << " __________\n";
    long double Events_in_Cent_cc = h_Cent ->GetBinContent(i+1);
    CovCrosschk1010 = (1/Events_in_Cent_cc)*(f[i][1][0] - f[i][1][0]*f[i][1][0] + f[i][2][0]);
    cout << "\n CrossCheck1010 = " << CovCrosschk1010 << endl;
  }

  cout << "\n Cross-check for covariance calculation ends here .... \n" << endl;
	*/

  LambdaCent_file -> Close();
}


void EC_Cumulants(Double_t **cumulants, Double_t **ec_cumulants)
{
	Int_t fac_OD1 = 4;
	Int_t fac_OD2 = 4;

  Double_t ***f;
	f = new Double_t**[9];
  for(Int_t i = 0; i < 9; i++)
  {
    f[i] = new Double_t*[fac_OD1];
    for(Int_t j = 0; j < fac_OD1; j++) f[i][j] = new Double_t[fac_OD2];
  }

	Emptying3DArray(f, 9, fac_OD1, fac_OD2);
        FactorialCalc(f, 9, fac_OD1, fac_OD2);

	Double_t *AVG_TOT_PART;
	AVG_TOT_PART = new Double_t[9];
	Emptying1DArray(AVG_TOT_PART, 9);

        Avg_Tot_Part_counting(AVG_TOT_PART);

	for(Int_t i = 0; i < 9; i++)
	{
		ec_cumulants[i][0] = cumulants[i][0]/Eff[i];  //CC_1
		ec_cumulants[i][1] = (cumulants[i][1] - AVG_TOT_PART[i]*(1 - Eff[i]))/TMath::Power(Eff[i], 2);  //CC_2
		ec_cumulants[i][2] = (cumulants[i][2] - cumulants[i][0]*(1 - TMath::Power(Eff[i], 2)) - 3*(1 - Eff[i])*(f[i][2][0] - f[i][0][2] - AVG_TOT_PART[i]*cumulants[i][0]))/TMath::Power(Eff[i], 3);
		ec_cumulants[i][3] = (cumulants[i][3]
		                     - AVG_TOT_PART[i] * TMath::Power(Eff[i], 2) * (1 - Eff[i])
		                     - 3 * TMath::Power(AVG_TOT_PART[i], 2) * TMath::Power((1 - Eff[i]), 2)
		                     - 6 * Eff[i] * (1 - Eff[i]) * (f[i][2][0] + f[i][0][2])
		                     + 12 * cumulants[i][0] * (1 - Eff[i]) * (f[i][2][0] - f[i][0][2])
		                     - (1 - TMath::Power(Eff[i], 2)) * (cumulants[i][1] - 3 * TMath::Power(cumulants[i][0], 2))
		                     - 6 * AVG_TOT_PART[i] * (1 - Eff[i]) * (TMath::Power(cumulants[i][0], 2) - cumulants[i][1])
		                     - 6 * (1 - Eff[i]) * (f[i][0][3] - f[i][1][2] + f[i][0][2] + f[i][2][0] - f[i][2][1] + f[i][3][0]))/TMath::Power(Eff[i], 4);
		ec_cumulants[i][4] = ec_cumulants[i][1]/ec_cumulants[i][0];
		ec_cumulants[i][5] = ec_cumulants[i][2]/ec_cumulants[i][1];
		ec_cumulants[i][6] = ec_cumulants[i][3]/ec_cumulants[i][1];
	}

	DeleteDynamic3DArray(f, 9, 4, "f");
	DeleteDynamic1DArray(AVG_TOT_PART, "AVG_TOT_PART");
}


void DeltaStatErrorForMoments(Double_t **mu, Double_t **error)
{

	TFile *LambdaCent_file      = new TFile(InRootFile, "read");
	TH1D *h_Cent                = (TH1D*) LambdaCent_file -> Get(hh_Cent);

	for(Int_t i = 0; i < 9; i++)
	{
		Int_t nevents = h_Cent ->GetBinContent(i+1);
		Double_t sigma = TMath::Sqrt(mu[i][1]);

		Double_t m_1 = mu[i][0]/sigma;
		Double_t m_2 = mu[i][1]/TMath::Power(sigma,2);
		Double_t m_3 = mu[i][2]/TMath::Power(sigma,3);
		Double_t m_4 = mu[i][3]/TMath::Power(sigma,4);
		Double_t m_5 = mu[i][4]/TMath::Power(sigma,5);
		Double_t m_6 = mu[i][5]/TMath::Power(sigma,6);
		Double_t m_7 = mu[i][6]/TMath::Power(sigma,7);
		Double_t m_8 = mu[i][7]/TMath::Power(sigma,8);

		error[i][0] = sigma/TMath::Sqrt(nevents);
		error[i][1] = TMath::Sqrt((m_4 - 1)*TMath::Power(sigma, 2)/(4*nevents));
		error[i][2] = TMath::Sqrt((9 - 6*m_4 + (m_3*m_3*(35 + 9*m_4)/4) - 3*m_3*m_5 + m_6)/nevents);
		error[i][3] = TMath::Sqrt((-1*m_4*m_4 + 4*TMath::Power(m_4, 3) + 16*TMath::Power(m_3, 2)*(1 + m_4) - 8*m_3*m_5 - 4*m_4*m_6 + m_8)/nevents);

		//cout << "\t " << CentName[8-i] << "\t\t\t " << error[i][0] << "\t\t " << error[i][1] << "\t\t " << error[i][2] << "\t\t " << error[i][3] << endl;
	}
	LambdaCent_file -> Close();
}


void EC_C1_Err(Double_t *C1_Error, Double_t *eff)
{
  TFile *LambdaCent_file      = new TFile(InRootFile, "read");
  TH1D *h_Cent                = (TH1D*) LambdaCent_file->Get(hh_Cent);

  Double_t ***f;

  f = new Double_t**[9];
  for(Int_t i = 0; i < 9; i++)
  {
    f[i] = new Double_t*[9];
    for(Int_t j = 0; j < 9; j++)f[i][j] = new Double_t[9];
  }

  FactorialCalc(f, 9, 9, 9);

  for(Int_t i = 0; i < 9; i++)
	{
    Int_t n = h_Cent ->GetBinContent(i+1);
    Double_t var_M = (1/(n*TMath::Power(eff[i], 2)))*(f[i][1][0] + f[i][0][1] - TMath::Power((f[i][1][0] - f[i][0][1]), 2) + f[i][0][2] - 2*f[i][1][1] + f[i][2][0]);
    C1_Error[i] = TMath::Sqrt(var_M);
    cout << "\n\tC1_[StatError]_from_Cov = " << C1_Error[i] << endl;
  }
  DeleteDynamic3DArray(f, 9, 9, "f");
  LambdaCent_file -> Close();
}


Double_t DC2(Int_t CentBin, Int_t Dim_1, Int_t Dim_2)
{
  Double_t D[9][3][3];
  for(Int_t i = 0; i < 9; i++)
    for(Int_t j = 0; j < 3; j++)
      for(Int_t k = 0; k < 3; k++)
        D[i][j][k] = 0.;

  Double_t ***F;
  Double_t ***f;

  F = new Double_t**[9];
  for(Int_t i = 0; i < 9; i++)
  {
    F[i] = new Double_t*[9];
    for(Int_t j = 0; j < 9; j++)F[i][j] = new Double_t[9];
  }

  f = new Double_t**[9];
  for(Int_t i = 0; i < 9; i++)
  {
  	f[i] = new Double_t*[9];
    for(Int_t j = 0; j < 9; j++)f[i][j] = new Double_t[9];
  }

  FactorialCalc(f, 9, 9, 9);
  EffCorrectFactorial(f, F, Eff, 9, 9, 9);

  for(Int_t i = 0; i < 9; i++)
  {
    Double_t M = F[i][1][0] - F[i][0][1];
    D[i][0][0] = 0.;
    D[i][0][1] = 1 + 2*M;
    D[i][0][2] = 1;
    D[i][1][0] = 1 - 2*M;
    D[i][1][1] = -2;
    D[i][1][2] = 0.;
    D[i][2][0] = 1;
    D[i][2][1] = 0.;
    D[i][2][2] = 0.;
  }

  DeleteDynamic3DArray(f, 9, 9, "f");
  DeleteDynamic3DArray(F, 9, 9, "F");

  Double_t Dout = D[CentBin][Dim_1][Dim_2];
  return Dout;
}

void EC_C2_Err(Double_t *C2_Error, Double_t *eff)
{
  Double_t *****COV;

  COV = new Double_t****[9];
  for(Int_t i = 0; i < 9; i++)
  {
    COV[i] = new Double_t***[5];
    for(Int_t j = 0; j < 5; j++)
    {
      COV[i][j] = new Double_t**[5];
      for(Int_t k = 0; k < 5; k++)
      {
        COV[i][j][k] = new Double_t*[5];
        for(Int_t l = 0; l < 5; l++) COV[i][j][k][l] = new Double_t[5];
      }
    }
  }

  CovarianceCalc(COV, 9, 5, 5, 5, 5);
  Double_t CovSum = 0;

  for(Int_t i = 0; i < 9; i++)
  {
    CovSum = 0;
    for(Int_t j = 0; j < 3; j++)
    {
      for(Int_t k = 0; k < 3; k++)
      {
        for(Int_t l = 0; l < 3; l++)
        {
        	for(Int_t m = 0; m < 3; m++)
        	{
        		CovSum += (1/(TMath::Power(eff[i], (j+k)) * TMath::Power(eff[i], (l+m)))) * DC2(i,j,k) * DC2(i,l,m) * COV[i][j][k][l][m];
        	}
        }
      }
    }
    C2_Error[i] = TMath::Sqrt(CovSum);
    cout << "\n\tC2_[StatError]_from_Cov = " << C2_Error[i] << endl;
  }

  DeleteDynamic5DArray(COV, 9, 5, 5, 5, "COV");
}

void EC_C2_C1_Err(Double_t *C2_C1_Error, Double_t *C_1_Error, Double_t *C_2_Error)
{
	Double_t ***f;
	f = new Double_t**[9];
  for(Int_t i = 0; i < 9; i++)
  {
    f[i] = new Double_t*[9];
    for(Int_t j = 0; j < 9; j++)f[i][j] = new Double_t[9];
  }

  Double_t ***F;
  F = new Double_t**[9];
  for(Int_t i = 0; i < 9; i++)
  {
    F[i] = new Double_t*[9];
    for(Int_t j = 0; j < 9; j++)F[i][j] = new Double_t[9];
  }

	FactorialCalc(f, 9, 9, 9);
  EffCorrectFactorial(f, F, Eff, 9, 9, 9);

	for(Int_t i = 0; i < 9;  i++)
	{
		Double_t C2 		= F[i][1][0] + F[i][0][1] - TMath::Power((F[i][1][0] - F[i][0][1]), 2) + F[i][0][2] - 2*F[i][1][1] + F[i][2][0];
		Double_t C1 		= F[i][1][0] - F[i][0][1];
		Double_t C2_C1 		= C2/C1;
		Double_t err_C2_C1 	= C2_C1*(TMath::Sqrt(TMath::Power(C_2_Error[i]/C2, 2) + TMath::Power(C_1_Error[i]/C1, 2)));
		C2_C1_Error[i]		= err_C2_C1;
		cout << "\n\tC2/C1_[StatError]_from_Cov = " << C2_C1_Error[i] << endl;
	}

	DeleteDynamic3DArray(f, 9, 9, "f");
  DeleteDynamic3DArray(F, 9, 9, "F");
}

Double_t DC3(Int_t CentBin, Int_t Dim_1, Int_t Dim_2)
{
  Double_t D[9][4][4];
  for(Int_t i = 0; i < 9; i++)
    for(Int_t j = 0; j < 4; j++)
      for(Int_t k = 0; k < 4; k++)
        D[i][j][k] = 0.;

  Double_t ***F;
  Double_t ***f;

  F = new Double_t**[9];
  for(Int_t i = 0; i < 9; i++)
  {
    F[i] = new Double_t*[9];
    for(Int_t j = 0; j < 9; j++)F[i][j] = new Double_t[9];
  }

  f = new Double_t**[9];
  for(Int_t i = 0; i < 9; i++)
  {
    f[i] = new Double_t*[9];
    for(Int_t j = 0; j < 9; j++)f[i][j] = new Double_t[9];
  }

  FactorialCalc(f, 9, 9, 9);
  EffCorrectFactorial(f, F, Eff, 9, 9, 9);

  for(Int_t i = 0; i < 9; i++)
  {
    Double_t M = F[i][1][0] - F[i][0][1];
    D[i][0][0] = 0.;
    D[i][0][1] = 3*F[i][0][1] + 3*F[i][0][2] + F[i][1][0] - 6*F[i][1][1] + 3*F[i][2][0] - 6*TMath::Power(M, 2) - 3*M - 1;
    D[i][0][2] = 3 - 3*M;
    D[i][0][3] = -1;
    D[i][1][0] = -3*F[i][0][1] - 3*F[i][0][2] -3*F[i][1][0] + 6*F[i][1][1] - 3*F[i][2][0] + 6*TMath::Power(M, 2) - 3*M + 1;
    D[i][1][1] = 6*M;
    D[i][1][2] = 3;
    D[i][1][3] = 0.;
    D[i][2][0] = 3 - 3*M;
    D[i][2][1] = -3;
    D[i][2][2] = 0.;
    D[i][2][3] = 0.;
    D[i][3][0] = 1;
    D[i][3][1] = 0.;
    D[i][3][2] = 0.;
    D[i][3][3] = 0.;
  }

  DeleteDynamic3DArray(f, 9, 9, "f");
  DeleteDynamic3DArray(F, 9, 9, "F");

  Double_t Dout = D[CentBin][Dim_1][Dim_2];
  return Dout;
}


void EC_C3_Err(Double_t *C3_Error, Double_t *eff)
{
	Double_t *****COV;

  COV = new Double_t****[9];
  for(Int_t i = 0; i < 9; i++)
  {
    COV[i] = new Double_t***[5];
    for(Int_t j = 0; j < 5; j++)
    {
      COV[i][j] = new Double_t**[5];
      for(Int_t k = 0; k < 5; k++)
      {
        COV[i][j][k] = new Double_t*[5];
        for(Int_t l = 0; l < 5; l++) COV[i][j][k][l] = new Double_t[5];
      }
    }
  }

  for(Int_t i = 0; i < 9; i++)
    for(Int_t j = 0; j < 5; j++)
      for(Int_t k = 0; k < 5; k++)
        for(Int_t l = 0; l < 5; l++)
          for(Int_t m = 0; m < 5; m++)
            COV[i][j][k][l][m] = 0.;

  CovarianceCalc(COV, 9, 5, 5, 5, 5);
  Double_t CovSum = 0;

	for(Int_t i = 0; i < 9; i++)
  {
    CovSum = 0;
    for(Int_t j = 0; j < 4; j++)
    {
      for(Int_t k = 0; k < 4; k++)
      {
        for(Int_t l = 0; l < 4; l++)
        {
        	for(Int_t m = 0; m < 4; m++)
        	{
        		CovSum += (1/(TMath::Power(eff[i], (j+k)) * TMath::Power(eff[i], (l+m)))) *
        		DC3(i,j,k) * DC3(i,l,m) * COV[i][j][k][l][m];
        	}
        }
      }
    }
    C3_Error[i] = TMath::Sqrt(CovSum);
    cout << "\n\tC3_[StatError]_from_Cov = " << C3_Error[i] << endl;
  }

  DeleteDynamic5DArray(COV, 9, 5, 5, 5, "COV");
}



Double_t DC3C2(Int_t CentBin, Int_t Dim_1, Int_t Dim_2)
{
  Double_t D[9][4][4];
  for(Int_t i = 0; i < 9; i++)
    for(Int_t j = 0; j < 4; j++)
      for(Int_t k = 0; k < 4; k++)
        D[i][j][k] = 0.;

  Double_t ***F;
  Double_t ***f;

  F = new Double_t**[9];
  for(Int_t i = 0; i < 9; i++)
  {
    F[i] = new Double_t*[9];
    for(Int_t j = 0; j < 9; j++)F[i][j] = new Double_t[9];
  }

  f = new Double_t**[9];
  for(Int_t i = 0; i < 9; i++)
  {
    f[i] = new Double_t*[9];
    for(Int_t j = 0; j < 9; j++)f[i][j] = new Double_t[9];
  }

  FactorialCalc(f, 9, 9, 9);
  EffCorrectFactorial(f, F, Eff, 9, 9, 9);

  Int_t *BIN_ENTRIES_REFMULT;
  BIN_ENTRIES_REFMULT = new Int_t[1000];

  Double_t **CENTRAL_MOM_REFMULT;
  CENTRAL_MOM_REFMULT = new Double_t*[1000];
  for(Int_t i = 0; i < 1000; i++) CENTRAL_MOM_REFMULT[i] = new Double_t[8];

  CentralMomCalcRefMult(BIN_ENTRIES_REFMULT, CENTRAL_MOM_REFMULT);

  Double_t **CBWC_MU;
  CBWC_MU = new Double_t*[9];
  for(Int_t i = 0; i < 9; i++) CBWC_MU[i] = new Double_t[8];

  CBWC_Central_mom(CENTRAL_MOM_REFMULT, CBWC_MU, BIN_ENTRIES_REFMULT, 9, 0, 1000, 1);

  for(Int_t i = 0; i < 9; i++)
  {
    Double_t mu1 = CBWC_MU[i][0];
    Double_t mu2 = CBWC_MU[i][1];
    Double_t mu3 = CBWC_MU[i][2];
    Double_t M = F[i][1][0] - F[i][0][1];

    D[i][0][0] = 0.;
    D[i][0][1] = ((3*(F[i][0][1] + F[i][0][2] + F[i][1][0] - 2*F[i][1][1] + F[i][2][0]) - 6*TMath::Power(M, 2) - 3*M -1)/mu2) - (1 + 2*M)*mu3/TMath::Power(mu2, 2);
    D[i][0][2] = (-3 - 3*M)/mu2 - mu3/TMath::Power(mu2, 2);
    D[i][0][3] = -1/mu2;
    D[i][1][0] = ((-3*(F[i][0][1] + F[i][0][2] + F[i][1][0] - 2*F[i][1][1] + F[i][2][0]) + 6*TMath::Power(M, 2) - 3*M + 1)/mu2) - (1 - 2*M)*mu3/TMath::Power(mu2, 2);
    D[i][1][1] = 6*M/mu2 + 2*mu3/TMath::Power(mu2, 2);
    D[i][1][2] = 3/mu2;
    D[i][1][3] = 0.;
    D[i][2][0] = (3 - 3*M)/mu2;
    D[i][2][1] = -3/mu2;
    D[i][2][2] = 0.;
    D[i][2][3] = 0.;
    D[i][3][0] = 1/mu2;
    D[i][3][1] = 0.;
    D[i][3][2] = 0.;
    D[i][3][3] = 0.;
  }

  DeleteDynamic3DArray(f, 9, 9, "f");
  DeleteDynamic3DArray(F, 9, 9, "F");
  DeleteDynamic1DArrayInt(BIN_ENTRIES_REFMULT, "BIN_ENTRIES_REFMULT");
  DeleteDynamic2DArray(CENTRAL_MOM_REFMULT, 1000, "CENTRAL_MOM_REFMULT");
  DeleteDynamic2DArray(CBWC_MU, 9, "CBWC_MU");

  Double_t Dout = D[CentBin][Dim_1][Dim_2];
  return Dout;
}

void EC_C3_C2_Err(Double_t *C3_C2_Error, Double_t *eff)
{
  Double_t *****COV;

  COV = new Double_t****[9];
  for(Int_t i = 0; i < 9; i++)
  {
    COV[i] = new Double_t***[5];
    for(Int_t j = 0; j < 5; j++)
    {
      COV[i][j] = new Double_t**[5];
      for(Int_t k = 0; k < 5; k++)
      {
        COV[i][j][k] = new Double_t*[5];
        for(Int_t l = 0; l < 5; l++) COV[i][j][k][l] = new Double_t[5];
      }
    }
	}

  for(Int_t i = 0; i < 9; i++)
    for(Int_t j = 0; j < 5; j++)
      for(Int_t k = 0; k < 5; k++)
        for(Int_t l = 0; l < 5; l++)
          for(Int_t m = 0; m < 5; m++)
            COV[i][j][k][l][m] = 0.;

  CovarianceCalc(COV, 9, 5, 5, 5, 5);
  Double_t CovSum = 0;

  for(Int_t i = 0; i < 9; i++)
  {
    CovSum = 0;
    for(Int_t j = 0; j < 4; j++)
    {
      for(Int_t k = 0; k < 4; k++)
      {
        for(Int_t l = 0; l < 4; l++)
        {
          for(Int_t m = 0; m < 4; m++)
          {
            CovSum += (1/(TMath::Power(eff[i], (j+k)) * TMath::Power(eff[i], (l+m)))) *
            DC3C2(i,j,k) * DC3C2(i,l,m) * COV[i][j][k][l][m];
          }
        }
      }
  	}
    C3_C2_Error[i] = TMath::Sqrt(CovSum);
    cout << "\n\tC3/C2_[StatError]_from_Cov = " << C3_C2_Error[i] << endl;
  }

  DeleteDynamic5DArray(COV, 9, 5, 5, 5, "COV");
}

void PlotPoints_CompareCorrection(TString Title_01, TCanvas *CanName, Int_t division, Double_t YrangeLow, Double_t YrangeHigh, Double_t *plot_1, Double_t *plot_2, Double_t *plot_1_Er, Double_t *plot_2_Er)
{

  TString Title_02 ="EC_";
  Title_02 += Title_01;

	TGraphErrors *graph_01 = new TGraphErrors(9, RefMult3BinCntre, plot_1, RefMult3BinCntreEr, plot_1_Er);//,&NpartErrorArray[0],&MeanErrorArray[0]);
	graph_01->SetTitle(Title_01);
	graph_01->GetXaxis()->SetTitle("Refmult3");
	graph_01->GetXaxis()->SetTitleSize(0.04);
	graph_01->GetYaxis()->SetTitle(Title_01);
	graph_01->GetYaxis()->SetTitleSize(0.04);
	graph_01->GetXaxis()->CenterTitle();
	graph_01->GetYaxis()->CenterTitle();

	graph_01->SetMarkerStyle(25);
	graph_01->SetMarkerColor(4);
	graph_01->GetYaxis()->SetRangeUser(YrangeLow, YrangeHigh);

	TGraphErrors *graph_02 = new TGraphErrors(9, RefMult3BinCntre, plot_2, RefMult3BinCntreEr, plot_2_Er);
	graph_02->SetTitle(Title_02);
	graph_02->GetXaxis()->SetTitle("Refmult3");
	graph_02->GetXaxis()->SetTitleSize(0.04);
	graph_02->GetYaxis()->SetTitle(Title_02);
	graph_02->GetYaxis()->SetTitleSize(0.04);
	graph_02->GetXaxis()->CenterTitle();
	graph_02->GetYaxis()->CenterTitle();

	graph_02->SetMarkerStyle(21);
	graph_02->SetMarkerColor(2);
	graph_02->GetYaxis()->SetRangeUser(YrangeLow, YrangeHigh);

	CanName  -> cd(division);
	graph_01 -> Draw("AP");
	graph_02 -> Draw("PSAME");
}

void PlotPoints_SingleMoments(TString Title_01, TCanvas *CanName, Int_t division, Double_t *plot, Double_t *plot_Er)
{
	TGraphErrors *graph_01 = new TGraphErrors(9, RefMult3BinCntre, plot, RefMult3BinCntreEr, plot_Er);
	graph_01->SetTitle(Title_01);
	graph_01->GetXaxis()->SetTitle("Refmult3");
	graph_01->GetXaxis()->SetTitleSize(0.04);
	graph_01->GetYaxis()->SetTitle(Title_01);
	graph_01->GetYaxis()->SetTitleSize(0.04);
	graph_01->GetXaxis()->CenterTitle();
	graph_01->GetYaxis()->CenterTitle();

	graph_01->SetMarkerStyle(25);
	graph_01->SetMarkerColor(4);

	CanName  -> cd(division);
	graph_01 -> Draw("AP");
}

/*________________________________________________________ THE END __________________________________________________*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
