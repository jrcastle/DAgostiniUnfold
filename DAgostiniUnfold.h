//-- Unfolding procedure follows that proposed by D'Agostini in: 
//--       "A multidimensional unfolding method based on Bayes' theorem", D'Agostini, NIM-A 362, 487
//-- Error propagation based on improved approach to D'Agostini as outlined in:
//--       http://hepunx.rl.ac.uk/~adye/software/unfold/bayes_errors.pdf
//-- Some methods in this class were stolen directly from the RooUnfoldBayes class.  Class reference and credit can be found here: 
//--       http://arxiv.org/pdf/1105.1160.pdf

#include <iostream>
#include <vector>
#include <string>

#include "TH1D.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"

class DAgostiniUnfold{

 public:
  //-- Constructor
  DAgostiniUnfold(TH2D * hresp, bool smoothResp = 0);
  ~DAgostiniUnfold();

  //-- Methods
  TMatrixD ABAT(TMatrixD a, TVectorD b, TMatrixD c);
  int computeObsCovMatrix(int iter);
  int computeObsErrMatrix(int iter);
  int computeRespCovMatrix(int iter);
  int computeRespErrMatrix(int iter);
  void DebugMessages();
  void DoSystematics();
  TMatrixD H2M(TH2D * h);
  TMatrixD H2M(TH2D * h, TMatrixD &Mw2);
  TVectorD H2V(TH1D * h);
  TVectorD H2V(TH1D * h, TVectorD &Vw2);
  double kronDelta(int i, int j);
  TH2D * M2H(TMatrixD M,TH2D * h);
  void makeEfficiencyVector();
  void makeRespCov0();
  TH1D * refold(TH1D * hreco, string name);
  void smooth(TVectorD &V);
  int Unfold(TH1D * hObserved, int niter, TH1D * hInitialPrior = 0);
  TH1D * V2H(TH1D * h, TVectorD v);
  TH1D * V2H(TH1D * h, TVectorD v, TVectorD vw2);

  //-- Getters
  TH1D * GetHReco(string name);
  TH2D * GetHRecoCov(string name);

  //-- Setters
  void SetSmoothIter();

 private:
  TMatrixD mResp_;        //-- Response matrix p( ei | cj ) with dimensions truth VS observed
  TMatrixD mRespw2_;      //-- Variance of each element of the response matrix with dimensions truth VS observed

  TMatrixD Mij_;          //-- The unfolding matrix with dimensions truth vs observed 
  TMatrixD obsErrMat_;    //-- The error matrix for each iteration that arises from the statistical uncertainties in the observed distributions
  TMatrixD respErrMat_;   //-- The error matrix for each iteration that arises from the uncertainties on the respose matrix elements
  TMatrixD obsCovMat_;    //-- The covariance matrix for each iteration that arises from the statistical uncertainties in the observed distributions
  TMatrixD respCovMat_;   //-- The covariance matrix for each iteration that arises from the statistical uncertainties in the response matrix elements 
  TMatrixD recoCovMat_;   //-- The covariance matrix for the final, unfolded distribution.  Equal to the sum of obsCovMat_ + respCovMat_

  TVectorD respCov0_;     //-- Special vector version of mRespw2_. Made specifically for the method computeRespCovMatrix(int iter)

  TVectorD vObs_;         //-- Vector of the observed, smeared histogram
  TVectorD vObsw2_;       //-- Vector of the variances on the observed, smeared histogram
  TVectorD vPrior_;       //-- "Prior" vector, before iterations it will be the guess at the true, underlying vector
  TVectorD vReco_;        //-- Unfolded vector
  TVectorD vRecow2_;      //-- Variances on the unfolded vector
  TVectorD foldPrior_;    //-- The response matrix applied to the prior vector

  TH2D * hResp_;          //-- Response matrix p( ei | cj ) as a histogram with dimensions truth VS observed

  vector<double> effVec_; //-- Detection efficiency vector = Sum[ P( ei | cj ) , { i, 1, Nobs} ]

  int kiter_;
  int Nobs;               //-- Number of bins in the observed histogram
  int Ntru;               //-- Number of bins in the true histogram
  bool dosys_;            //-- Choice to propagate uncertainties on the response matrix. These are regarded as systematic uncertainties
  bool DEBUG_;            //-- If true, adds some cout statements so the user can track where the code is if/when it breaks
  bool smoothResp_;       //-- If true, calls the TH2::Smooth method on the input response matrix
  bool smoothIter_;       //-- If true, the unfolded vector of each iteration will be smoothed using the TH1::SmoothArray method
};


//-- =======================================
//--               CONSTRUCTOR             
//-- =======================================

DAgostiniUnfold::DAgostiniUnfold(TH2D * hresp, bool smoothResp){

  std::cout<<"Constructing a new DAgostiniUnfold object..."<<std::endl;
  hResp_      = 0;
  dosys_      = 0;
  DEBUG_      = 0;
  smoothIter_ = 0;
  smoothResp_ = smoothResp;

  //-- Store the response matrix into a class member
  hResp_ = (TH2D*) hresp->Clone("hResp_");
  if( smoothResp_ ){
    std::cout<<"Response matrix smoothing enabled"<<std::endl;
    hResp_->Smooth(1, "k5b");
  }

  TH1D * hy = (TH1D*) hResp_->ProjectionY();
  Nobs   = hResp_->GetNbinsX();
  Ntru   = hResp_->GetNbinsY();

  //-- Initialize the vector corresponding to the prior distribution
  vPrior_.ResizeTo(Ntru);
  TVectorD vy = H2V(hy);
  vPrior_ = vy;

  //-- Normalize the response matrix by the projection onto the truth axis
  for(int y = 1; y <= Ntru; y++){
    double ntrue = hy->GetBinContent(y);
    if(ntrue==0){
      for(int x = 1; x <= Nobs; x++){
	hResp_->SetBinContent(x, y, 0.);
      }
    }
    else{
      for(int x = 1; x <= Nobs; x++){
	double binCont   = hResp_->GetBinContent(x, y);
	double binErr    = hResp_->GetBinError(x, y);
	hResp_->SetBinContent(x, y, binCont / ntrue);
	hResp_->SetBinError(x, y, binErr / ntrue);
      }
    }
  }
  delete hy;  

  //-- Convert the TH2D hResp_ into a TMatrixD
  mResp_.ResizeTo(Ntru, Nobs);
  mRespw2_.ResizeTo(Ntru, Nobs);
  mResp_ = H2M(hResp_, mRespw2_);
  makeRespCov0();

  //-- Create a std::vector that contains the detection efficiencies for each truth bin
  makeEfficiencyVector();

  //-- Initialize the error propagation matrix for the observed distribution
  //-- dn(cj)/dn(ei)
  //-- x <==> ei
  //-- y <==> cj
  obsErrMat_.ResizeTo(Ntru, Nobs);

  //-- Initialize the error propagation matrix for the response matrix
  //-- dn(cj)/dP( ei | ck )
  //-- Collapse i and k onto the "x" axis of this matrix
  //-- x <==> ei ck
  //-- y <==> cj
  if( dosys_ ) respErrMat_.ResizeTo(Ntru, Ntru*Nobs);

  //-- Initialize the covariance matrices for the respective error propagation matrices, which will build
  //-- the final covariance martrix for the unfolded distribution Cov_unfold = Cov_obs + Cov_resp
  obsCovMat_.ResizeTo(Ntru, Ntru);
  respCovMat_.ResizeTo(Ntru, Ntru);
  recoCovMat_.ResizeTo(Ntru, Ntru);

  //-- Initialize what will be the final unfolded distribution
  vReco_.ResizeTo(Ntru);
  vRecow2_.ResizeTo(Ntru);

  //-- Initialize what will be the unfolding matrix
  //-- x <==> ei
  //-- y <==> cj
  Mij_.ResizeTo(Ntru, Nobs);

  //-- Initialize what will be the folded prior vector
  foldPrior_.ResizeTo(Nobs);

}

//-- =======================================
//--                DESTRUCTOR
//-- =======================================

DAgostiniUnfold::~DAgostiniUnfold(){

  if( DEBUG_ ) std::cout<<"Destroying DAgostiniUnfold object..."<<std::endl;
  if( hResp_ ) delete hResp_;
  
}

//-- =======================================
//--                  ABAT                
//-- =======================================

TMatrixD DAgostiniUnfold::ABAT (TMatrixD a, TVectorD b, TMatrixD c){

  //-- Fills C such that C = A * B * A^T, where B is a diagonal matrix specified by the vector.
  //-- Note that C cannot be the same object as A.
  //-- NB this method was stolen directly from the RooUnfoldBayes Class.  Credit to Tim Ayde
  TMatrixD d (TMatrixD::kTransposed, a);
  d.NormByColumn (b, "M");
  c.Mult (a, d);
  return c;

}

//-- =======================================
//--            computeObsCovMatrix
//-- =======================================

int DAgostiniUnfold::computeObsCovMatrix(int iter){

  if( DEBUG_ ) std::cout<<"Converting the observed error propagation matrix to a covariance matrix for iteration "<<iter<<std::endl;

  TMatrixD dummy(Ntru, Ntru);
  obsCovMat_ = ABAT(obsErrMat_, vObsw2_, dummy);

  return 1;

}

//-- =======================================
//--           computeObsErrMatrix        
//-- =======================================

int DAgostiniUnfold::computeObsErrMatrix(int iter){

  if( DEBUG_ ) std::cout<<"Computing the observed error propagation matrix for iteration "<<iter<<std::endl;

  if(iter == 0){
    obsErrMat_ = Mij_;
    return 1;
  }
  else{

    TVectorD ksum(Nobs);

    for (Int_t j = 0 ; j < Nobs ; j++) {
      for (Int_t k = 0 ; k < Nobs ; k++) {

	Double_t sum = 0.0;
	for (Int_t l = 0 ; l < Ntru ; l++) {
	  if ( vPrior_(l) > 0.0 ) sum += effVec_[l] * Mij_(l, k) * obsErrMat_(l, j) / vPrior_(l);
	}
	ksum(k) = sum;

      }

      for (Int_t i = 0 ; i < Ntru ; i++) {

	Double_t dsum;
	if( vPrior_(i) > 0 ) dsum = obsErrMat_(i,j) * vReco_(i) / vPrior_(i);
	else                 dsum = 0.;

	for (Int_t k = 0 ; k < Nobs ; k++) {
	  dsum -= Mij_(i,k) * vObs_(k) * ksum(k);
	}

	//-- Update obsErrMat_. Note that we can do this in-place due to the ordering of the accesses.
	obsErrMat_(i,j) = Mij_(i,j) + dsum;

      }

    }

    return 1;
  
  } // End else 

  return 0; //-- This return should never be reached, but is placed here to keep the compiler happy

}

//-- =======================================
//--            computeRespCovMatrix
//-- =======================================

int DAgostiniUnfold::computeRespCovMatrix(int iter){

  if( DEBUG_ ) std::cout<<"Converting the response error propagation matrix to a covariance matrix for iteration "<<iter<<std::endl;

  TMatrixD dummy(Ntru, Ntru);
  respCovMat_ = ABAT(respErrMat_, respCov0_, dummy);
  return 1;

}

//-- =======================================
//--            computeRespErrMatrix        
//-- =======================================

int DAgostiniUnfold::computeRespErrMatrix(int iter){

  if( DEBUG_ ) std::cout<<"Computing the response error propagation matrix for iteration "<<iter<<std::endl;

  //-- dn(cj)/dP( ei | ck )
  //-- x <==> ei ck
  //-- y <==> cj
  /*
  if(iter == 0){

    for(int i = 0; i < Nobs; i++){
      for(int j = 0; j < Ntru; j++){
        for(int k = 0; k < Ntru; k++){

          double epsj   = effVec_[j];
          double n0Cj   = vPrior_(j);
          double nEi    = vObs_(i);
          double nHatCj = vReco_(j);
          double n0Ck   = vPrior_(k);
          double Mij    = Mij_(j, i);
          double fi     = foldPrior_(i);
          if(fi == 0) continue;

	  double dncj_dPeick = 0.;
	  if( epsj == 0){
	    dncj_dPeick = -( n0Ck * nEi / fi ) * Mij;
	  }
	  else{
	    dncj_dPeick = (1./epsj) * ( n0Cj * nEi / fi - nHatCj ) * kronDelta(j,k) - ( n0Ck * nEi / fi ) * Mij;
	  }

	  //-- Collapse i and k to the "x" dimension of the matrix
          //--
          //--          __                                                             __
          //--  j=0    |                                                                 |
          //--   .     |                                                                 |
          //--   .     |   k=0           k=1             k=2      ...         k=nC-1     |
          //--   .     |                                                                 |
          //--   .     |                                                                 |
          //--  j=nC-1 |__                                                             __|
          //--          i=1...nE-1   i=nE...2nE-1   i=2nE...3nE-1 ... i=(nC-1)nE...nCnE-1

          int bin = i + k * Nobs;
          //respErrMat_(j, bin) = dncj_dPeick;

        }
      }
    }

    return 1;

  } //-- End if(iter == 0)
  else{
    
    TMatrixD M1( Ntru, Ntru*Nobs );
    TMatrixD suml( Nobs, Ntru*Nobs );

    for(int i = 0; i < Nobs; i++){
      for(int k = 0; k < Ntru; k++){

	int bin = i + k * Nobs;

	for(int l = 0; l < Nobs; l++){

	  double lsum = 0.;
	  for(int r = 0; r < Ntru; r++){

	    double Mlr            = Mij_(l, r);
	    double errMatPrev_irk = respErrMat_(r, bin);

	    lsum += Mlr * errMatPrev_irk;

	  }

	  suml(l, bin) = lsum;

	}

	for(int j = 0; j < Ntru; j++){

	  int bin = i + k * Nobs;

	  //-- M1 part 1/4
	  double epsj   = effVec_[j];
	  double n0Cj   = vPrior_(j);
	  double nhatCj = vReco_(j);
	  double nEi    = vObs_(i);
	  double fi     = foldPrior_(i);

	  if( epsj == 0 || fi == 0 ) M1(j, bin) = 0.;
	  else                       M1(j, bin) = (1./epsj) * ( n0Cj * nEi / fi - nhatCj ) * kronDelta(j, k);

	  //-- M1 part 2/4
	  double n0Ck = vPrior_(k);
	  double Mij  = Mij_(j, i);

	  if( fi == 0 ) M1(j, bin) -= 0.;
	  else          M1(j, bin) -= ( n0Ck * nEi / fi ) * Mij;

	  //-- M1 part 3/4
	  double errMatPrev_ijk = respErrMat_(j, bin);

	  if( n0Cj == 0) M1(j, bin) += 0.;
	  else           M1(j, bin) += ( nhatCj / n0Cj ) * errMatPrev_ijk;

	  //-- M1 part 4/4
	  //-- To save some computer cycles, check to see if n0cj == 0.  
	  //-- If so, there is not need to do part 4/4 as this contrubution will be 0;
	  if( n0Cj == 0 ) continue;
	  double doubSum = 0.;
	  for( int l = 0; l < Nobs; l++){

	    double nEl      = vObs_(l);
	    double Mjl      = Mij_(j, l);
	    double lsum_lik = suml(l, bin);

	    doubSum += nEl * Mjl * lsum_lik;

	  }

	  if( n0Cj == 0 ) M1(j, bin) -= 0.;
	  else            M1(j, bin) -= (epsj / n0Cj) * doubSum;

	}
      }
    }

    //-- Update respErrMat_
    respErrMat_ = M1;

    return 1;

    } //-- End else

    return 0; //-- This return should never be reached, but is placed here to keep the compiler happy

  */

  //-- How it's done in RooUnfold
  if (iter > 0) {
    TVectorD mbyu(Nobs);
    for (Int_t i = 0 ; i < Nobs ; i++) {
      double fi = foldPrior_(i);
      if(fi == 0) continue;
      double nEi = vObs_(i);
      mbyu(i)= nEi /fi;
    }
    TMatrixD A= Mij_;
    A.NormByRow (mbyu, "M");
    TMatrixD B(A, TMatrixD::kMult, mResp_);
    TMatrixD dnCidPjkUpd (B, TMatrixD::kMult, respErrMat_);
    Int_t nec= Nobs*Ntru;
    for (Int_t j = 0 ; j < Ntru ; j++) {
      if (vPrior_(j)<=0.0) continue;  // skip loop: dnCidPjkUpd(i,jk) will also be 0 because _Mij(i,j) will be 0                                                               
      Double_t r= vReco_(j) / vPrior_(j);
      for (Int_t ik= 0; ik<nec; ik++)
	respErrMat_(j,ik)= r*respErrMat_(j,ik) - dnCidPjkUpd(j,ik);
    }
  }

  for (Int_t i = 0 ; i < Nobs; i++) {
    if (foldPrior_(i)==0.0) continue;
    Double_t mbyu= 1./foldPrior_(i)*vObs_(i);
    Int_t i0= i*Ntru;
    for (Int_t j = 0 ; j < Ntru ; j++) {
      Double_t b= -mbyu * Mij_(j,i);
      for (Int_t k = 0 ; k < Ntru ; k++) respErrMat_(j,i0+k) += b*vPrior_(k);
      if (effVec_[j]!=0.0)
	respErrMat_(j,i0+j) += (vPrior_(j)*mbyu - vReco_(j)) / effVec_[j];
    }
  }

  return 1;

}

//-- =======================================
//--              DebugMessages
//-- =======================================

void DAgostiniUnfold::DebugMessages(){

  DEBUG_ = 1;

}

//-- =======================================
//--              DoSystematics            
//-- =======================================

void DAgostiniUnfold::DoSystematics(){

  std::cout<<"dosys_ set to true.  Will propagate uncertainties on the response matrix..."<<std::endl;
  dosys_ = 1;

}

//-- =======================================
//--                   H2M
//-- =======================================

TMatrixD DAgostiniUnfold::H2M(TH2D * h){

  //-- Convert a TH2D into a TMatrixD.

  int nx = h->GetNbinsX();
  int ny = h->GetNbinsY();

  TMatrixD M(ny, nx);
  for(int i = 0; i<nx; i++){
    for(int j = 0; j<ny; j++){
      M(j, i) = h->GetBinContent(i+1, j+1);
    }
  }

  return M;

}

//-- =======================================

TMatrixD DAgostiniUnfold::H2M(TH2D * h, TMatrixD &Mw2){

  //-- Convert a TH2D into a TMatrixD.  Pass another matrix in to build it's respective error matrix

  int nx = h->GetNbinsX();
  int ny = h->GetNbinsY();

  TMatrixD M(ny, nx);
  for(int i = 0; i<nx; i++){
    for(int j = 0; j<ny; j++){
      M(j, i)   = h->GetBinContent(i+1, j+1);
      Mw2(j, i) = pow( h->GetBinError(i+1, j+1), 2 );
    }
  }

  return M;

}

//-- =======================================
//--                   H2V
//-- ======================================= 

TVectorD DAgostiniUnfold::H2V(TH1D * h){

  //-- Converts a TH1D into a TVectorD
  //-- If Vw2 != 0, Vw2 will be filled with the square of the error bars on the TH1D

  int nx = h->GetNbinsX();

  TVectorD V(nx);
  for(int i = 0; i<nx; i++) V(i) = h->GetBinContent(i+1);

  return V;

}

//-- =======================================

TVectorD DAgostiniUnfold::H2V(TH1D * h, TVectorD &Vw2){

  //-- Converts a TH1D into a TVectorD
  //-- If Vw2 != 0, Vw2 will be filled with the square of the error bars on the TH1D

  int nx = h->GetNbinsX();

  TVectorD V(nx);
  for(int i = 0; i<nx; i++){
    V(i)   = h->GetBinContent(i+1);
    Vw2(i) = pow( h->GetBinError(i+1), 2 );
  }

  return V;

}

//-- =======================================
//--                kronDelta
//-- =======================================

double DAgostiniUnfold::kronDelta(int i, int j){

  //-- Kronecker delta

  if(i == j) return 1.;
  return 0.;

}

//-- =======================================
//--                   M2H                  
//-- =======================================

TH2D * DAgostiniUnfold::M2H(TMatrixD M, TH2D * h){

  //-- Converts a TMatrixD into a TH2D without error bars

  //-- Check dimensions
  int hdimx = h->GetNbinsX();
  int hdimy = h->GetNbinsY();
  int mdimx = M.GetNcols();
  int mdimy = M.GetNrows();

  if( hdimx != mdimx){
    std::cout<<"Warning, in H2M the x dimensions of the matrix and TH2D do not match!"<<std::endl;
    std::cout<<"Returning a null pointer for the TH2D."<<std::endl;
    return 0;
  }
  if( hdimy != mdimy){
    std::cout<<"Warning, in H2M the x dimensions of the matrix and TH2D do not match!"<<std::endl;
    std::cout<<"Returning a null pointer for the TH2D."<<std::endl;
    return 0;
  }

  for(int i = 1; i <= hdimx; i++){
    for(int j = 1; j <= hdimy; j++){
      h->SetBinContent(i, j, M(j-1, i-1) );
    }
  }

  return h;

}

//-- =======================================
//--           makeEfficiencyVector         
//-- =======================================

void DAgostiniUnfold::makeEfficiencyVector(){

  //-- Make the detection efficiency vector eff_j = Sum[ P( ei | cj ) , { i, 1, Nobs} ]

  for(int y = 0; y < Ntru; y++){

    double sum = 0.;
    for(int m = 0; m < Nobs; m++) sum += mResp_(y, m);
    effVec_.push_back(sum);

  }

}

//-- =======================================
//--               makeRespCov0            
//-- =======================================

void DAgostiniUnfold::makeRespCov0(){

  //-- Makes a diagonal matrix of the variances of the initial response matrix

  respCov0_.ResizeTo(Nobs*Ntru);

  int binE = 0;
  int binC = 0;
  for(int i = 0; i < Nobs*Ntru; i++){

    if(i != 0) binE++;
    if(binE >= Nobs){
      binC++;
      binE = 0;
    }

    respCov0_(i) = mRespw2_(binE, binC);

  }

}

//-- =======================================
//--                  refold
//-- =======================================

TH1D * DAgostiniUnfold::refold(TH1D * hreco, string name){

  int nt = hreco->GetNbinsX();
  if( nt != Ntru ){
    std::cout<<"Warning!!! Binning in hreco does not match Ntru!"<<std::endl;
    std::cout<<"Returning an null pointer for the refolded histogram"<<std::endl;
    return 0;
  }

  TVectorD vSmear(Ntru);
  TMatrixD mResp(TMatrixD::kTransposed, mResp_);
  vSmear = H2V(hreco);
  vSmear *= mResp;

  double smL = hResp_->ProjectionX()->GetBinLowEdge(1);
  double smH = hResp_->ProjectionX()->GetBinLowEdge(Nobs) + hResp_->ProjectionY()->GetBinWidth(Nobs);

  TH1D * hsmear = new TH1D(name.data(), name.data(), Nobs, smL, smH);
  for(int i = 1; i <= Nobs; i++){
    hsmear->SetBinContent( i, vSmear(i-1) );
  }

  return hsmear;

}

//-- =======================================
//--                  smooth
//-- =======================================

void DAgostiniUnfold::smooth(TVectorD &V){

  std::cout<<"Smoothing..."<<std::endl;
  TH1::SmoothArray(Ntru, V.GetMatrixArray(), 1);

}

//-- =======================================
//--                  Unfold                
//-- =======================================

int DAgostiniUnfold::Unfold(TH1D * hObserved, int niter, TH1D * hInitialPrior){

  kiter_ = niter;

  vObs_.ResizeTo(Nobs);
  vObsw2_.ResizeTo(Nobs);
  vObs_ = H2V(hObserved, vObsw2_);

  //-- If the initial prior is not specified, use the projection onto the truth axis as the starting point
  if(hInitialPrior){
    std::cout<<"Manually setting the starting distribution"<<std::endl;
    TVectorD v = H2V(hInitialPrior);
    vPrior_ = v;
  }
  else std::cout<<"Using the projection of the response onto the truth axis as the starting distribution"<<std::endl;

  for(int k = 0; k < kiter_; k++){

    //-- Update the prior for subsequent iterations
    if(k > 0) vPrior_ = vReco_;

    //-- Fold the prior
    for(int x = 0; x < Nobs; x++){
      double sum = 0;
      for(int L = 0; L < Ntru; L++) sum += mResp_(L, x) * vPrior_(L);
      foldPrior_(x) = sum;
    }

    //-- Construct the unfolding matrix
    std::cout<<"Constructing the unfolding matrix for iteration "<<k<<"..."<<std::endl;
    for(int x = 0; x < Nobs; x++){
      for(int y = 0; y < Ntru; y++){

	double epsj  = effVec_[y];
	if(epsj == 0) continue;

	double PeicLc0L = foldPrior_(x);
	if(PeicLc0L == 0) continue;

	double Peicj = mResp_(y, x);
	double c0j   = vPrior_(y);

	double Mij = Peicj * c0j / PeicLc0L / epsj;
	Mij_(y, x) = Mij;

      }
    } //-- End unfolding matrix construction loop

    //-- Unfold!
    TVectorD vDummy = vObs_;
    vDummy *= Mij_;
    vReco_ = vDummy;
    //-- If enabled, smooth all but the last iteration
    if( smoothIter_ && k < kiter_-1 ) smooth(vReco_);


    //-- Set up the error propagation matrices    
    int obsErrCheck  = computeObsErrMatrix(k);
    int respErrCheck;
    if(dosys_) respErrCheck = computeRespErrMatrix(k);

    if( !obsErrCheck ){
      std::cout<<"WARNING!!! Procedure broke when constructing the observed error propagation matrix"<<std::endl;
      return 0;
    }
    if( dosys_ && !respErrCheck ){
      std::cout<<"WARNING!!! Procedure broke when constructing the response error propagation matrix"<<std::endl;
      return 0;
    }
    
    //-- Set up the covaciance matrices
    int obsCovCheck  = computeObsCovMatrix(k);
    int respCovCheck;
    if(dosys_) respCovCheck = computeRespCovMatrix(k);

    if( !obsCovCheck ){
      std::cout<<"WARNING!!! Procedure broke when constructing the observed covariance matrix"<<std::endl;
      return 0;
    }
    if( dosys_ && !respCovCheck ){
      std::cout<<"WARNING!!! Procedure broke when constructing the response covariance matrix"<<std::endl;
      return 0;
    }

    //-- Make the covariance matrix for the unfolded distribution
    TMatrixD mDummy(Ntru, Ntru);
    mDummy += obsCovMat_;
    if( dosys_ ) mDummy += respCovMat_;
    recoCovMat_ = mDummy;

    //Error bars on hReco_ set as the square root of the diagonal elements of the cov matrix for hReco_
    for(int i = 0; i<Ntru; i++){
      double w2 = recoCovMat_(i, i);
      vRecow2_(i) = w2;
    }

  } //-- End iteration loop

  return 1;
}

//-- =======================================
//--                   V2H                 
//-- =======================================

TH1D * DAgostiniUnfold::V2H(TH1D * h, TVectorD v){

  //-- Converts a TVectorD into a TH1D
  
  //-- Check dimensions
  int hdim = h->GetNbinsX();
  int vdim = v.GetNoElements();
  if(hdim != vdim){
    std::cout<<"Warning! vReco dimension does not match hReco dimension!"<<std::endl;
    std::cout<<"Returning a null pointer for hReco"<<std::endl;
    return 0;
  }

  for(int i = 1; i<= hdim; i++) h->SetBinContent( i, vReco_(i-1) );

  return h;

}

//-- =======================================

TH1D * DAgostiniUnfold::V2H(TH1D * h, TVectorD v, TVectorD vw2){

  //-- Converts a TVectorD into a TH1D
  //-- If vw2 != 0, sets the error bars of the TH1D as the sqrt of the elements of vw2

  //-- Check dimensions
  int hdim = h->GetNbinsX();
  int vdim = v.GetNoElements();
  if(hdim != vdim){
    std::cout<<"Warning! vReco dimension does not match hReco dimension!"<<std::endl;
    std::cout<<"Returning a null pointer for hReco"<<std::endl;
    return 0;
  }

  for(int i = 1; i<= hdim; i++){
    h->SetBinContent( i, vReco_(i-1) );
    h->SetBinError( i, TMath::Sqrt( vRecow2_(i-1) ) );
  }

  return h;

}


//-- =======================================
//--                 GetHReco               
//-- =======================================
TH1D * DAgostiniUnfold::GetHReco(string name){

  double tL = hResp_->ProjectionY()->GetBinLowEdge(1);
  double tH = hResp_->ProjectionY()->GetBinLowEdge(Ntru) + hResp_->ProjectionY()->GetBinWidth(Ntru);

  TH1D * hReco_ = new TH1D(name.data(), name.data(), Ntru, tL, tH);
  V2H(hReco_, vReco_, vRecow2_);

  //-- Returns the unfolded distribution as a TH1D object
  return hReco_;

}

//-- =======================================
//--               GetHRecoCov              
//-- =======================================

TH2D * DAgostiniUnfold::GetHRecoCov(string name){

  double tL = hResp_->ProjectionY()->GetBinLowEdge(1);
  double tH = hResp_->ProjectionY()->GetBinLowEdge(Ntru) + hResp_->ProjectionY()->GetBinWidth(Ntru);

  TH2D * hRecoCov_ = new TH2D(name.data(), name.data(), Ntru, tL, tH, Ntru, tL, tH);
  M2H(recoCovMat_, hRecoCov_);

  //-- Returns the unfolded distribution's covariance matrix as a TH2D object
  return hRecoCov_;

}

//-- =======================================
//--              SetSmoothIter             
//-- =======================================

void DAgostiniUnfold::SetSmoothIter(){

  std::cout<<"Smoothing of each unfolding iteration enabled"<<std::endl;
  smoothIter_ = 1;

}
