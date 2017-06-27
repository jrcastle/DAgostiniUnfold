#include "/home/j550c590/tdrstyle.C"
#include "DAgostiniUnfold.h"
#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/interface/differentailBinning.h"
#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfoldResponse.h"
#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfoldBayes.h"
#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfold.h"
#include "TLine.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TH2D.h" 
#include "TLatex.h"
#include "TFile.h"
#include "TCanvas.h"

void testDAgostiniUnfold(){

  bool validate = 1;
  bool dosys    = 0;
  int ptBin     = 0;
  int centBin   = 1;
  int VN        = 2;

  static const int NITER  = 4;
  static const int iter[] = {1, 2, 4, 8, 16, 32, 64, 128};
  int col[]               = {kOrange-2, kGreen+3, kCyan, kMagenta, kViolet-1, kBlue, kRed, kGray+2};

  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;

  setTDRStyle();
  TLatex latex;
  latex.SetNDC();

  TFile * fMC   = new TFile( Form("/rfs/jcastle/PbPb2015/mimicDataMC/cent%i-%i/pt%.2f-%.2f/MCTree.root", cent_min[centBin], cent_max[centBin], pt_min[ptBin], pt_max[ptBin]) );
  TFile * fData = new TFile( Form("/rfs/jcastle/PbPb2015/mimicDataMC/cent%i-%i/pt%.2f-%.2f/CastleEbyE.root", cent_min[centBin], cent_max[centBin], pt_min[ptBin], pt_max[ptBin]) );

  //-- CastleEbyE.root data members
  static const int NVn = 7;
  TH1D * hVnFull[NVn][NCENT_PBPB];
  TH2D * hResp[NVn][NCENT_PBPB];

  //-- Unfolding data
  TH2D * hresp[NCENT];
  TH1D * h1D[NCENT];
  TH1D * hreco1[NCENT][NITER];
  TH1D * hrefold[NCENT][NITER];
  TH1D * hreco1_RU[NCENT][NITER];
  TH1D * hratio[NCENT][NITER];
  TH1D * hratioErr[NCENT][NITER];
  RooUnfoldResponse * RUresp[NCENT];
  RooUnfoldBayes * RUBayes[NCENT][NITER];

  TLine * line = new TLine(0., 1., v2Max[ptBin][centBin], 1.);
  line->SetLineColor(2);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TLegend * leg[NCENT][NITER];
  TCanvas * cIter[NCENT][NITER];
  TCanvas * cIterErr[NCENT][NITER];

  //-- Get DATA
  for(int ivn = 1; ivn <=4; ivn++){
    if(ivn != VN) continue;
    for(int icent = 0; icent < NCENT_PBPB; icent++){
      if(icent == NCENT40_cutoff) break;
      hVnFull[ivn][icent]   = (TH1D*) fData->Get(Form("qwebye/hVnFull_%i_%i",ivn,icent));
      hResp[ivn][icent]     = (TH2D*) fData->Get(Form("qwebye/hResp_%i_%i",ivn,icent));
    }
  }

  //-- Condense from 2.5 to 10% cent binning
  for(int c = 0; c < NCENT; c++){

    if(c != centBin) continue;
    std::cout<<"!! Processing Cent = "<<c<<std::endl;

    //-- Observed distribution
    h1D[c] = (TH1D*)hVnFull[VN][4*c]->Clone(Form("h1D_%i_%i", VN, c));
    h1D[c]->Add(hVnFull[VN][4*c+1]);
    h1D[c]->Add(hVnFull[VN][4*c+2]);
    h1D[c]->Add(hVnFull[VN][4*c+3]);

    //-- Response matrix
    hresp[c] = (TH2D*) hResp[VN][4*c]->Clone(Form("hresp_%i_%i", VN, c));
    hresp[c]->SetOption("colz");
    hresp[c]->Add(hResp[VN][4*c+1]);
    hresp[c]->Add(hResp[VN][4*c+2]);
    hresp[c]->Add(hResp[VN][4*c+3]);

    RUresp[c] = new RooUnfoldResponse(0, 0, hresp[c], Form("response_%i_%i", VN, c));

    for(int i = 0; i<NITER; i++){

      if( !validate && i > 0 ) break;

      hreco1[c][i]    = 0;
      hreco1_RU[c][i] = 0;

      DAgostiniUnfold unfolder(hresp[c]);
      if(dosys) unfolder.DoSystematics();
      unfolder.Unfold(h1D[c], iter[i]);

      hreco1[c][i] = (TH1D*) unfolder.GetHReco( Form("hreco1_%ci_iter%i", c, iter[i]) );
      hrefold[c][i] = (TH1D*) unfolder.refold(hreco1[c][i], Form("hrefold_%ci_iter%i", c, iter[i]) );
      hreco1[c][i]->GetXaxis()->SetTitle("v_{2}");
      hreco1[c][i]->GetYaxis()->SetTitle("Events");
      hreco1[c][i]->SetLineColor(col[i]);
      hreco1[c][i]->SetMarkerColor(col[i]);
      hreco1[c][i]->SetMarkerStyle(20);

      RUBayes[c][i] = new RooUnfoldBayes(RUresp[c], h1D[c], iter[i]);
      hreco1_RU[c][i] = (TH1D*) RUBayes[c][i]->Hreco(errorTreatment);
      hreco1_RU[c][i]->GetXaxis()->SetTitle("v_{2}");
      hreco1_RU[c][i]->GetYaxis()->SetTitle("Events");
      hreco1_RU[c][i]->SetLineColor(col[i+4]);
      hreco1_RU[c][i]->SetMarkerColor(col[i+4]);
      hreco1_RU[c][i]->SetMarkerStyle(24);

      hratio[c][i] = (TH1D*) hreco1[c][i]->Clone( Form("ratio_%i", i ) );
      hratio[c][i]->Divide( hreco1_RU[c][i] );
      hratio[c][i]->GetXaxis()->SetTitle("v_{2}");
      hratio[c][i]->GetYaxis()->SetTitle("Ratio");

      
      hratioErr[c][i] = (TH1D*) hratio[c][i]->Clone( Form("hratioErr_%ci_iter%i", c, iter[i]) );
      hratioErr[c][i]->Reset();
      hratioErr[c][i]->GetXaxis()->SetTitle("v_{2}");
      hratioErr[c][i]->GetYaxis()->SetTitle("Error Bar Ratio");
      hratioErr[c][i]->SetLineWidth(3);

      for(int j = 1; j <= hratio[c][i]->GetNbinsX(); j++){
	
	double myErr = hreco1[c][i]->GetBinError(j);
	double ruErr = hreco1_RU[c][i]->GetBinError(j);
	if(ruErr == 0) continue;

	double ratio = myErr/ruErr;
	hratioErr[c][i]->SetBinContent(j, ratio);
      }
      

      leg[c][i] = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg[c][i]->SetFillStyle(0);
      leg[c][i]->SetBorderSize(0);
      leg[c][i]->AddEntry(hreco1[c][i], "CastleUnfold", "lp");
      leg[c][i]->AddEntry(hreco1_RU[c][i], "RooUnfold", "lp");
      /*
      cIter[c][i] = new TCanvas( Form("citer_c%i_iter%i", c, iter[i]), Form("citer_c%i_iter%i", c, iter[i]), 500, 1000);
      cIter[c][i]->Divide(1,2);
      cIter[c][i]->cd(1);
      cIter[c][i]->cd(1)->SetLogy();
      hreco1[c][i]->Draw();
      hreco1_RU[c][i]->Draw("same");
      leg[c][i]->Draw("same");
      latex.DrawLatex(0.2, 0.22, Form("N_{iter} = %i", iter[i]) );
      cIter[c][i]->cd(2);
      hratio[c][i]->Draw();

      cIter[c][i]->SaveAs( Form("plots/citer_c%i_iter%i.png", c, iter[i]) );
      */
      cIterErr[c][i] = new TCanvas( Form("citerErr_c%i_iter%i", c, iter[i]), Form("citerErr_c%i_iter%i", c, iter[i]), 500, 500);
      cIterErr[c][i]->cd();
      hratioErr[c][i]->Draw();
      line->Draw("same");

    } //-- End iter loop

  } //-- End cent loop

}
