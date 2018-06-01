#include "LambdaMoments.h"

#include <iostream>
#include <string>

#include "TFile.h"

using namespace std;

void MainProgram(TString Append_E, TString Append_cut = "d")
{
	cout << "\n\t\t\t\t\t~~~-<><>- Main program starts here -<><>-~~~ \n" << endl;

	Int_t *BIN_ENTRIES_REFMULT;
	BIN_ENTRIES_REFMULT = new Int_t[1000];
	Emptying1DArrayInt(BIN_ENTRIES_REFMULT, 1000);

	Double_t **CENTRAL_MOM_REFMULT;
	CENTRAL_MOM_REFMULT = new Double_t*[1000];
	for(Int_t i = 0; i < 1000; i++) CENTRAL_MOM_REFMULT[i] = new Double_t[8];
  Emptying2DArray(CENTRAL_MOM_REFMULT, 1000, 8);

	CentralMomCalcRefMult(BIN_ENTRIES_REFMULT, CENTRAL_MOM_REFMULT);

  Double_t **CBWC_MU;
	CBWC_MU = new Double_t*[9];
	for(Int_t i = 0; i < 9; i++) CBWC_MU[i] = new Double_t[8];
	Emptying2DArray(CBWC_MU, 9, 8);

  CBWC_Central_mom(CENTRAL_MOM_REFMULT, CBWC_MU, BIN_ENTRIES_REFMULT, 9, 0, 1000, 1);

	Double_t **MOMENTS_REFMULT;
	MOMENTS_REFMULT = new Double_t*[1000];
	for(Int_t i = 0; i < 1000; i++) MOMENTS_REFMULT[i] = new Double_t[4];
	Emptying2DArray(MOMENTS_REFMULT, 1000, 4);

	Moments(CENTRAL_MOM_REFMULT, MOMENTS_REFMULT, 0, 1000);

	Double_t **CBWC_MOMENTS;
	CBWC_MOMENTS = new Double_t*[9];
	for(Int_t i = 0; i < 9; i++) CBWC_MOMENTS[i] = new Double_t[4];
	Emptying2DArray(CBWC_MOMENTS, 9, 4);

	CBWC_Moments(MOMENTS_REFMULT, CBWC_MOMENTS, BIN_ENTRIES_REFMULT, 9, 0, 1000, 1);

  Double_t **CUMULANTS_REFMULT;
	CUMULANTS_REFMULT = new Double_t*[1000];
	for(Int_t i = 0; i < 1000; i++) CUMULANTS_REFMULT[i] = new Double_t[7];
	Emptying2DArray(CUMULANTS_REFMULT, 1000, 7);

	Cumulants(CENTRAL_MOM_REFMULT, CUMULANTS_REFMULT, 0, 1000);

	Double_t **CBWC_CUMULANTS;
	CBWC_CUMULANTS = new Double_t*[9];
	for(Int_t i = 0; i < 9; i++) CBWC_CUMULANTS[i] = new Double_t[7];
	Emptying2DArray(CBWC_CUMULANTS, 9, 7);

	CBWC_Cumulants(CUMULANTS_REFMULT, CBWC_CUMULANTS, BIN_ENTRIES_REFMULT, 9, 0, 1000, 1);

	Double_t **EC_CUMULANTS;
	EC_CUMULANTS = new Double_t*[9];
	for(Int_t i = 0; i < 9; i++) EC_CUMULANTS[i] = new Double_t[7];
	Emptying2DArray(EC_CUMULANTS, 9, 7);

	EC_Cumulants(CBWC_CUMULANTS, EC_CUMULANTS);

	Double_t **EC_CUMULANTS_DF;
	EC_CUMULANTS_DF = new Double_t*[9];
	for(Int_t i = 0; i < 9; i++) EC_CUMULANTS_DF[i] = new Double_t[7];
	Emptying2DArray(EC_CUMULANTS_DF, 9, 7);

	EC_Cumulants_DiffEff(EC_CUMULANTS_DF);

	Double_t **ERROR_MOMENTS;
	ERROR_MOMENTS = new Double_t*[9];
	for(Int_t i = 0; i < 9; i++) ERROR_MOMENTS[i] = new Double_t[4];
	Emptying2DArray(ERROR_MOMENTS, 9, 4);

	DeltaStatErrorForMoments(CBWC_MU, ERROR_MOMENTS);

  /*______________ EC_Errors -> C1, C2, C3, C2/C1, C3/C2 _________________*/

	Double_t *EC_C1_ERROR 	 = new Double_t[9];
	Double_t *EC_C2_ERROR 	 = new Double_t[9];
	Double_t *EC_C3_ERROR 	 = new Double_t[9];
	Double_t *EC_C2_C1_ERROR = new Double_t[9];
	Double_t *EC_C3_C2_ERROR = new Double_t[9];

	Double_t C1_ERROR[9] 	 = {};
	Double_t C2_ERROR[9] 	 = {};
	Double_t C3_ERROR[9] 	 = {};
	Double_t C2_C1_ERROR[9]  = {};
	Double_t C3_C2_ERROR[9]  = {};

  cout << "\n\t_________C1 Statistical uncertainty_________\n" << endl;

	EC_C1_Err(EC_C1_ERROR, Eff);

	cout << "\n\t_________C2 Statistical uncertainty_________\n" << endl;

	EC_C2_Err(EC_C2_ERROR, Eff);

	cout << "\n\t_________C3 Statistical uncertainty_________\n" << endl;

	EC_C3_Err(EC_C3_ERROR, Eff);

	cout << "\n\t________C2/C1 Statistical uncertainty_________\n" << endl;

	EC_C2_C1_Err(EC_C2_C1_ERROR, EC_C1_ERROR, EC_C2_ERROR);

	cout << "\n\t________C3/C2 Statistical uncertainty_________\n" << endl;

	EC_C3_C2_Err(EC_C3_C2_ERROR, Eff);

  /*________________________________________ Ploting Starts __________________________________________*/

  Double_t plot_c1[9];				Double_t plot_c_c1[9];
  Double_t plot_c2[9];				Double_t plot_c_c2[9];
	Double_t plot_c3[9];				Double_t plot_c_c3[9];
  Double_t plot_c2_c1[9];			Double_t plot_c_c2_c1[9];
	Double_t plot_c3_c2[9];			Double_t plot_c_c3_c2[9];

	Emptying1DArray(plot_c1, 9);		Emptying1DArray(plot_c_c1, 9);
	Emptying1DArray(plot_c2, 9);		Emptying1DArray(plot_c_c2, 9);
	Emptying1DArray(plot_c3, 9);		Emptying1DArray(plot_c_c3, 9);
	Emptying1DArray(plot_c2_c1, 9);	Emptying1DArray(plot_c_c2_c1, 9);
	Emptying1DArray(plot_c3_c2, 9);	Emptying1DArray(plot_c_c3_c2, 9);

	for(Int_t i = 0; i < 9; i++) {plot_c1[i] 			= CBWC_CUMULANTS[i][0]; }
	for(Int_t i = 0; i < 9; i++) {plot_c_c1[i] 		= EC_CUMULANTS[i][0]; }

	for(Int_t i = 0; i < 9; i++) {plot_c2[i] 			= CBWC_CUMULANTS[i][1]; }
	for(Int_t i = 0; i < 9; i++) {plot_c_c2[i] 		= EC_CUMULANTS[i][1]; }

	for(Int_t i = 0; i < 9; i++) {plot_c3[i] 			= CBWC_CUMULANTS[i][2]; }
	for(Int_t i = 0; i < 9; i++) {plot_c_c3[i] 		= EC_CUMULANTS[i][2]; }

	for(Int_t i = 0; i < 9; i++) {plot_c2_c1[i] 	= CBWC_CUMULANTS[i][4]; }
	for(Int_t i = 0; i < 9; i++) {plot_c_c2_c1[i] = EC_CUMULANTS[i][4]; }

	for(Int_t i = 0; i < 9; i++) {plot_c3_c2[i] 	= CBWC_CUMULANTS[i][5]; }
	for(Int_t i = 0; i < 9; i++) {plot_c_c3_c2[i] = EC_CUMULANTS[i][5]; }

	TCanvas *canvas_A = new TCanvas("canvas_A", "canvas_A", 400, 900);
  canvas_A -> Divide(1,3);

	PlotPoints_CompareCorrection("C_{1}", canvas_A, 1, -0.5, 8, plot_c1, plot_c_c1, C1_ERROR, EC_C1_ERROR);
	PlotPoints_CompareCorrection("C_{2}", canvas_A, 2, -0.5, 16, plot_c2, plot_c_c2, C2_ERROR, EC_C2_ERROR);
	PlotPoints_CompareCorrection("C_{3}", canvas_A, 3, -5, 30, plot_c3, plot_c_c3, C3_ERROR, EC_C3_ERROR);

	TCanvas *canvas_B = new TCanvas("canvas_B", "canvas_B", 400, 700);
  canvas_B -> Divide(1,2);

	PlotPoints_CompareCorrection("C_{2}/C_{1}", canvas_B, 1, -0.5, 6, plot_c2_c1, plot_c_c2_c1, C2_C1_ERROR, EC_C2_C1_ERROR);
	PlotPoints_CompareCorrection("C_{3}/C_{2}", canvas_B, 2, -2, 9, plot_c3_c2, plot_c_c3_c2, C3_C2_ERROR, EC_C3_C2_ERROR);

	TFile *OutFile_01     = new TFile("OutFileCumulantsCompareCorrection.root", "RECREATE");
	canvas_A 		-> Write();
	canvas_B 		-> Write();
	OutFile_01  -> Close();

	Double_t plot_M[9];				Double_t plot_M_Er[9];
	Double_t plot_Sigma[9];		Double_t plot_Sigma_Er[9];
	Double_t plot_S[9];				Double_t plot_S_Er[9];
	Double_t plot_K[9];				Double_t plot_K_Er[9];

	Emptying1DArray(plot_M, 9);					Emptying1DArray(plot_M_Er, 9);
	Emptying1DArray(plot_Sigma, 9);			Emptying1DArray(plot_Sigma_Er, 9);
	Emptying1DArray(plot_S, 9);					Emptying1DArray(plot_S_Er, 9);
	Emptying1DArray(plot_K, 9);					Emptying1DArray(plot_K_Er, 9);

	for(Int_t i = 0; i < 9; i++) {plot_M[i] 				= CBWC_MOMENTS[i][0]; }
	for(Int_t i = 0; i < 9; i++) {plot_M_Er[i] 			= ERROR_MOMENTS[i][0]; }

	for(Int_t i = 0; i < 9; i++) {plot_Sigma[i] 		= CBWC_MOMENTS[i][1]; }
	for(Int_t i = 0; i < 9; i++) {plot_Sigma_Er[i] 	= ERROR_MOMENTS[i][1]; }

	for(Int_t i = 0; i < 9; i++) {plot_S[i] 				= CBWC_MOMENTS[i][2]; }
	for(Int_t i = 0; i < 9; i++) {plot_S_Er[i] 			= ERROR_MOMENTS[i][2]; }

	for(Int_t i = 0; i < 9; i++) {plot_K[i] 				= CBWC_MOMENTS[i][3]; }
	for(Int_t i = 0; i < 9; i++) {plot_K_Er[i] 			= ERROR_MOMENTS[i][3]; }

	TCanvas *canvas_C = new TCanvas("canvas_C", "canvas_C", 850, 700);
  canvas_C -> Divide(2,2);

	PlotPoints_SingleMoments("Mean(M)", canvas_C, 1, plot_M, plot_M_Er);
	PlotPoints_SingleMoments("St.Dev(#sigma)", canvas_C, 2, plot_Sigma, plot_Sigma_Er);
	PlotPoints_SingleMoments("Skewness(S)", canvas_C, 3, plot_S, plot_S_Er);
	PlotPoints_SingleMoments("Kurtosis(#kappa)", canvas_C, 4, plot_K, plot_K_Er);

	TFile *OutFile_02     = new TFile("OutFileSingleMoments_01.root", "RECREATE");
	canvas_C 		-> Write();
	OutFile_02  -> Close();

	//TString Append_eta = "_050_"; //Not needed for sys. error calc.
	TString Append = "";
	//Append += Append_E; //Not needed for sys. error calc.
	//Append += Append_eta; //Not needed for sys. error calc.
	Append += Append_cut;

	cout << "\n \t\t\t\t[__________ Cumulants calculated on " << Append_E << "GeV data set. ___________]\n\n" << endl;

	/*________________________________________ Ploting Points Ends __________________________________________*/

	cout << "\n \t=====================================================|> Uncorr. Moments Arrays - " << Append << " - STARTS.\n" << endl;

	cout << "\n\tDouble_t M_" << Append << "[9] \t= {";
	for(Int_t i = 0; i < 9; i++){cout << plot_M[i]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\n\tDouble_t M_Er_" << Append << "[9] \t= {";
	for(Int_t i = 0; i < 9; i++){cout << plot_M_Er[i]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\n\tDouble_t Sigma_" << Append << "[9] \t= {";
	for(Int_t i = 0; i < 9; i++){cout << plot_Sigma[i]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\n\tDouble_t Sigma_Er_" << Append << "[9] \t= {";
	for(Int_t i = 0; i < 9; i++){cout << plot_Sigma_Er[i]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\n\tDouble_t S_" << Append << "[9] \t= {";
	for(Int_t i = 0; i < 9; i++){cout << plot_S[i]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\n\tDouble_t S_Er_" << Append << "[9] \t= {";
	for(Int_t i = 0; i < 9; i++){cout << plot_S_Er[i]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\n\tDouble_t K_" << Append << "[9] \t= {";
	for(Int_t i = 0; i < 9; i++){cout << plot_K[i]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\n\tDouble_t K_Er_" << Append << "[9] \t= {";
	for(Int_t i = 0; i < 9; i++){cout << plot_K_Er[i]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\n \t=====================================================|> Moment Arrays END.\n\n" << endl;

	/*________________________________________ Cumulant Arrays for Sys. Error _______________________________*/

  cout << "\n \t=====================================================|> Ci and Ci/Cj Arrays - " << Append << " - STARTS.\n" << endl;

  cout << "\n\tDouble_t C1_" << Append << "[9] \t= {";
  for(Int_t i = 0; i < 9; i++){cout << plot_c_c1[i]; if(i != 8){cout << ", ";}}
  cout << "}; \n" << endl;

  cout << "\n\tDouble_t C1_Er_" << Append << "[9] \t= {";
  for(Int_t i = 0; i < 9; i++){cout << EC_C1_ERROR[i]; if(i != 8){cout << ", ";}}
  cout << "}; \n" << endl;

  cout << "\n\tDouble_t C2_" << Append << "[9] \t= {";
  for(Int_t i = 0; i < 9; i++){cout << plot_c_c2[i]; if(i != 8){cout << ", ";}}
  cout << "}; \n" << endl;

  cout << "\n\tDouble_t C2_Er_" << Append << "[9] \t= {";
  for(Int_t i = 0; i < 9; i++){cout << EC_C2_ERROR[i]; if(i != 8){cout << ", ";}}
  cout << "}; \n" << endl;

  cout << "\n\tDouble_t C3_" << Append << "[9] \t= {";
  for(Int_t i = 0; i < 9; i++){cout << plot_c_c3[i]; if(i != 8){cout << ", ";}}
  cout << "}; \n" << endl;

  cout << "\n\tDouble_t C3_Er_" << Append << "[9] \t= {";
  for(Int_t i = 0; i < 9; i++){cout << EC_C3_ERROR[i]; if(i != 8){cout << ", ";}}
  cout << "}; \n" << endl;

  cout << "\n\tDouble_t C2C1_" << Append << "[9] \t= {";
  for(Int_t i = 0; i < 9; i++){cout << plot_c_c2_c1[i]; if(i != 8){cout << ", ";}}
  cout << "}; \n" << endl;

  cout << "\n\tDouble_t C2C1_Er_" << Append << "[9] \t= {";
  for(Int_t i = 0; i < 9; i++){cout << EC_C2_C1_ERROR[i]; if(i != 8){cout << ", ";}}
  cout << "}; \n" << endl;

  cout << "\n\tDouble_t C3C2_" << Append << "[9] \t= {";
  for(Int_t i = 0; i < 9; i++){cout << plot_c_c3_c2[i]; if(i != 8){cout << ", ";}}
  cout << "}; \n" << endl;

  cout << "\n\tDouble_t C3C2_Er_" << Append << "[9] \t= {";
  for(Int_t i = 0; i < 9; i++){cout << EC_C3_C2_ERROR[i]; if(i != 8){cout << ", ";}}
  cout << "}; \n" << endl;

  cout << "\n \t=====================================================|> Cumulant Arrays END.\n" << endl;

  /*
  cout << "\n\n =====  Via FacMom_Moments (Eff_La == Eff_aLa is not nessasary) - STARTS  ===== \n" << endl;

	cout << "\nDouble_t C1_Fac_" << Append << "[9] = {";
	for(Int_t i = 0; i < 9; i++){cout << EC_CUMULANTS_DF[i][0]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\nDouble_t C1_Fac_Er_" << Append << "[9] = {";
	for(Int_t i = 0; i < 9; i++){cout << EC_C1_ERROR[i]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\nDouble_t C2_Fac_" << Append << "[9] = {";
	for(Int_t i = 0; i < 9; i++){cout << EC_CUMULANTS_DF[i][1]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\nDouble_t C2_Fac_Er_" << Append << "[9] = {";
	for(Int_t i = 0; i < 9; i++){cout << EC_C2_ERROR[i]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\nDouble_t C3_Fac_" << Append << "[9] = {";
	for(Int_t i = 0; i < 9; i++){cout << EC_CUMULANTS_DF[i][2]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\nDouble_t C3_Fac_Er_" << Append << "[9] = {";
	for(Int_t i = 0; i < 9; i++){cout << EC_C3_ERROR[i]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\nDouble_t C2C1_Fac_" << Append << "[9] = {";
	for(Int_t i = 0; i < 9; i++){cout << EC_CUMULANTS_DF[i][1]/EC_CUMULANTS_DF[i][0]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\nDouble_t C2C1_Fac_Er_" << Append << "[9] = {";
	for(Int_t i = 0; i < 9; i++){cout << EC_C2_C1_ERROR[i]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\nDouble_t C3C2_Fac_" << Append << "[9] = {";
	for(Int_t i = 0; i < 9; i++){cout << EC_CUMULANTS_DF[i][2]/EC_CUMULANTS_DF[i][1]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\nDouble_t C3C2_Fac_Er_" << Append << "[9] = {";
	for(Int_t i = 0; i < 9; i++){cout << EC_C3_C2_ERROR[i]; if(i != 8){cout << ", ";}}
	cout << "}; \n" << endl;

	cout << "\n\n ===  FacMom_Moments arrays (Eff_La == Eff_aLa not nessasary) - ENDS  ===== \n" << endl;
	*/

	/*_______________________________________________________________________________________________________*/

	DeleteDynamic1DArrayInt(BIN_ENTRIES_REFMULT, "BIN_ENTRIES_REFMULT");
	DeleteDynamic2DArray(CENTRAL_MOM_REFMULT, 1000, "CENTRAL_MOM_REFMULT");
	DeleteDynamic2DArray(CBWC_MU, 9, "CBWC_MU");
	DeleteDynamic2DArray(MOMENTS_REFMULT, 1000, "MOMENTS_REFMULT");
	DeleteDynamic2DArray(CBWC_MOMENTS, 9, "CBWC_MOMENTS");
	DeleteDynamic2DArray(CUMULANTS_REFMULT, 1000, "CUMULANTS_REFMULT");
	DeleteDynamic2DArray(CBWC_CUMULANTS, 9, "CBWC_CUMULANTS");
	DeleteDynamic2DArray(EC_CUMULANTS, 9, "EC_CUMULANTS");
	DeleteDynamic2DArray(EC_CUMULANTS_DF, 9, "EC_CUMULANTS_DF");
	DeleteDynamic2DArray(ERROR_MOMENTS, 9, "ERROR_MOMENTS");

	DeleteDynamic1DArray(EC_C1_ERROR, "EC_C1_ERROR");
	DeleteDynamic1DArray(EC_C2_ERROR, "EC_C2_ERROR");
	DeleteDynamic1DArray(EC_C3_ERROR, "EC_C3_ERROR");
	DeleteDynamic1DArray(EC_C2_C1_ERROR, "EC_C2_C1_ERROR");
	DeleteDynamic1DArray(EC_C3_C2_ERROR, "EC_C3_C2_ERROR");

	cout << "\n\t\t\t\t\t________________________________________" << endl;
	cout << "\t\t\t\t\t-<><><>- Main program ends here -<><><>- \n\n" << endl;
}
