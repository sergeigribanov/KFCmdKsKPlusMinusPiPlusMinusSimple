#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <set>

#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>

#include <KFCmd/Hypo3ChPionsKPlus.hpp>
#include <KFCmd/Hypo3ChPionsKMinus.hpp>

#include "TrPh.h"

const double TrPh::_dZ = 30;
const double TrPh::_dRho = 30;
const double TrPh::_mindEdX = 0;
const double TrPh::_maxdEdX = 15000;
const double TrPh::_minTPtot = 5;
const double TrPh::_maxTPtot = 1000;

TrPh::TrPh(TTree *tree) :
  KFCmd::TrPh(tree) {
}

TrPh::~TrPh() {
}

bool TrPh::cutTracks() {
  _trackIndices.clear();
  for (int i = 0; i < nt; i++) {
    bool point = (std::fabs(tz[i]) < _dZ) && (std::fabs(trho[i]) < _dRho);
    bool dedx = (tdedx[i] > _mindEdX) && (tdedx[i] < _maxdEdX);
    bool ptot = (tptot[i] > _minTPtot) && (tptot[i] < _maxTPtot);
    if (point && dedx && ptot) _trackIndices.push_back(i);
  }
  if (_trackIndices.size() == 4) {
    int totalCharge = 0;
    for (int i = 0; i < 4; ++i) totalCharge += tcharge[_trackIndices[i]];
    return (totalCharge == 0);
  }
  return false;
}

bool TrPh::cutPhotons() {
  double phEnergy = 0;
  for (int i = 0; i < nph; ++i) {
    phEnergy += phen[i];
  }
  return phEnergy < 50;
}

Int_t TrPh::Cut(Long64_t) {
  if (nt < 4) return -1;
  if (!cutTracks()) return -1;
  // if (!cutPhotons()) return -1;
  std::vector<Int_t> charges(nt);
  std::copy(tcharge, tcharge + nt, charges.begin());
  std::sort(_trackIndices.begin(), _trackIndices.end(),
            [&charges](int i, int j) { return charges[i] < charges[j]; });
  return 1;
}

void TrPh::Loop(const std::string& outpath, double magneticField) {
  if (fChain == 0) return;
  std::set<std::string> sKs = {"pi+_0", "pi-_0"};
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TH1F h_kf_mks("h_kf_mks", "", 512, 0, 2000);
  TH1F h_in_mks("h_in_mks", "", 512, 0, 2000);
  TH1F h_kf_chi2("h_kf_chi2", "", 1024, 0, 1024);
  TH1F h_vtx0_x("h_vtx0_x", "", 600, -30, 30);
  TH1F h_vtx0_y("h_vtx0_y", "", 600, -30, 30);
  TH1F h_vtx0_z("h_vtx0_z", "", 600, -30, 30);
  TH1F h_vtx1_x("h_vtx1_x", "", 600, -30, 30);
  TH1F h_vtx1_y("h_vtx1_y", "", 600, -30, 30);
  TH1F h_vtx1_z("h_vtx1_z", "", 600, -30, 30);
  TH1F h_vtx0_r("h_vtx0_r", "", 600, 0, 60);
  TH1F h_vtx1_dr("h_vtx1_dr", "", 600, 0, 60);
  TH1F h_ksminv("ksminv", "", 512, 0, 2000);
  TH2F h_kfm_ksminv("kfm_ksminv", "", 512, 0, 2000, 512, 0, 2000);
  TH2F h_inm_ksminv("inm_ksminv", "", 512, 0, 2000, 512, 0, 2000);
  TH2F h_pdedx_k("pdedx_k", "", 1024, 0, 16384, 1024, 0, 1024);
  TH2F h_piv("pdedx_piv", "", 1024, 0, 16384, 1024, 0, 1024);
  TH2F h_pib("pdedx_pib", "", 1024, 0, 16384, 1024, 0, 1024);
  fChain->GetEntry(0);
  KFCmd::Hypo3ChPionsKPlus hypo_plus(2 * emeas, magneticField);
  hypo_plus.setBeamXY(xbeam, ybeam);
  hypo_plus.fixVertexComponent("vtx0", xbeam, KFBase::VERTEX_X);
  hypo_plus.fixVertexComponent("vtx0", ybeam, KFBase::VERTEX_Y);
  KFCmd::Hypo3ChPionsKMinus hypo_minus(2 * emeas, magneticField);
  hypo_minus.setBeamXY(xbeam, ybeam);
  hypo_minus.fixVertexComponent("vtx0", xbeam, KFBase::VERTEX_X);
  hypo_minus.fixVertexComponent("vtx0", ybeam, KFBase::VERTEX_Y);
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  double kf_chi2_plus;
  double tchi2_plus;
  double kf_chi2_minus;
  double tchi2_minus;
  bool flag_plus;
  bool flag_minus;
  int errCode;
  double v_kf_mks_plus = 0;
  double v_in_mks_plus = 0;
  double v_kf_mks_minus = 0;
  double v_in_mks_minus = 0;

  double v_dedx_plus = 0;
  double v_p_plus = 0;
  double v_dedx_minus = 0;
  double v_p_minus = 0;

  double v_dedx_piv_plus[2] = {0., 0.};
  double v_p_piv_plus[2] = {0., 0.};
  double v_dedx_piv_minus[2] = {0., 0.};
  double v_p_piv_minus[2] = {0., 0.};

  double v_dedx_pib_plus = 0;
  double v_p_pib_plus = 0;
  double v_dedx_pib_minus = 0;
  double v_p_pib_minus = 0;
  
  TVector3 v_vtx0_plus;
  TVector3 v_vtx1_plus;
  TVector3 v_vtx0_minus;
  TVector3 v_vtx1_minus;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (Cut(ientry) < 0) continue;
    if (nks != 1) continue;
    hypo_plus.setBeamXY(xbeam, ybeam);
    hypo_plus.fixVertexComponent("vtx0", xbeam, KFBase::VERTEX_X);
    hypo_plus.fixVertexComponent("vtx0", ybeam, KFBase::VERTEX_Y);
    hypo_minus.setBeamXY(xbeam, ybeam);
    hypo_minus.fixVertexComponent("vtx0", xbeam, KFBase::VERTEX_X);
    hypo_minus.fixVertexComponent("vtx0", ybeam, KFBase::VERTEX_Y);
    flag_plus = false;
    flag_minus = false;
    kf_chi2_plus = std::numeric_limits<double>::infinity();
    kf_chi2_minus = std::numeric_limits<double>::infinity();
    for (int im = 0; im < 2; ++im) {
      if (!hypo_plus.fillTrack("pi-_0", _trackIndices[im], *this)) continue;
      if (!hypo_plus.fillTrack("pi-_1", _trackIndices[1 - im], *this)) continue;
      for (int ip = 2; ip < 4; ++ip) {
	if (!hypo_plus.fillTrack("pi+_0", _trackIndices[ip], *this)) continue;
	if (!hypo_plus.fillTrack("k+", _trackIndices[5 - ip], *this)) continue;
	hypo_plus.optimize();
	errCode = hypo_plus.getErrorCode();
	if (errCode != 0) continue;
	tchi2_plus = hypo_plus.getChiSquare();
	if (tchi2_plus < kf_chi2_plus) {
	  flag_plus = true;
	  kf_chi2_plus = tchi2_plus;
	  v_in_mks_plus = hypo_plus.getInitialMomentum(sKs).M();
	  v_kf_mks_plus = hypo_plus.getFinalMomentum(sKs).M();
	  v_vtx0_plus = hypo_plus.getFinalVertex("vtx0");
	  v_vtx1_plus = hypo_plus.getFinalVertex("vtx1");

	  v_dedx_plus = tdedx[_trackIndices[5 - ip]];
	  v_p_plus = hypo_plus.getFinalMomentum("k+").P();
	  v_dedx_piv_plus[0] = tdedx[_trackIndices[im]];
	  v_dedx_piv_plus[1] = tdedx[_trackIndices[ip]];
	  v_p_piv_plus[0] = hypo_plus.getFinalMomentum("pi-_0").P();
	  v_p_piv_plus[1] = hypo_plus.getFinalMomentum("pi+_0").P();
	  v_dedx_pib_plus = tdedx[_trackIndices[1 - im]];
	  v_p_pib_plus = hypo_plus.getFinalMomentum("pi-_1").P();
	}
      }
    }

    for (int im = 0; im < 2; ++im) {
      if (!hypo_minus.fillTrack("pi-_0", _trackIndices[im], *this)) continue;
      if (!hypo_minus.fillTrack("k-", _trackIndices[1 - im], *this)) continue;
      for (int ip = 2; ip < 4; ++ip) {
	if (!hypo_minus.fillTrack("pi+_0", _trackIndices[ip], *this)) continue;
	if (!hypo_minus.fillTrack("pi+_1", _trackIndices[5 - ip], *this)) continue;
	hypo_minus.optimize();
	errCode = hypo_minus.getErrorCode();
	if (errCode != 0) continue;
	tchi2_minus = hypo_minus.getChiSquare();
	if (tchi2_minus < kf_chi2_minus) {
	  flag_minus = true;
	  kf_chi2_minus = tchi2_minus;
	  v_in_mks_minus = hypo_minus.getInitialMomentum(sKs).M();
	  v_kf_mks_minus = hypo_minus.getFinalMomentum(sKs).M();
	  v_vtx0_minus = hypo_minus.getFinalVertex("vtx0");
	  v_vtx1_minus = hypo_minus.getFinalVertex("vtx1");

	  v_dedx_minus = tdedx[_trackIndices[1 - im]];
	  v_p_minus = hypo_minus.getFinalMomentum("k-").P();
	  v_dedx_piv_minus[0] = tdedx[_trackIndices[im]];
	  v_dedx_piv_minus[1] = tdedx[_trackIndices[ip]];
	  v_p_piv_minus[0] = hypo_minus.getFinalMomentum("pi-_0").P();
	  v_p_piv_minus[1] = hypo_minus.getFinalMomentum("pi+_0").P();
	  v_dedx_pib_minus = tdedx[_trackIndices[5 - ip]];
	  v_p_pib_minus = hypo_minus.getFinalMomentum("pi+_1").P();
	}
      }
    }
    
    if (flag_plus && flag_minus) {
      if (kf_chi2_plus > kf_chi2_minus) {
	flag_plus = false;
      } else {
	flag_minus = false;
      }
    }

    if (flag_plus) {
      h_kf_chi2.Fill(kf_chi2_plus);
      if (kf_chi2_plus < 100) {
	h_inm_ksminv.Fill(v_in_mks_plus, ksminv[0]);
	h_kfm_ksminv.Fill(v_kf_mks_plus, ksminv[0]);
	h_ksminv.Fill(ksminv[0]);
	h_in_mks.Fill(v_in_mks_plus);
	h_kf_mks.Fill(v_kf_mks_plus);
	h_vtx0_x.Fill(v_vtx0_plus.X());
	h_vtx0_y.Fill(v_vtx0_plus.Y());
	h_vtx0_z.Fill(v_vtx0_plus.Z());
	h_vtx1_x.Fill(v_vtx1_plus.X());
	h_vtx1_y.Fill(v_vtx1_plus.Y());
	h_vtx1_z.Fill(v_vtx1_plus.Z());
	h_vtx0_r.Fill(v_vtx0_plus.Perp());
	h_vtx1_dr.Fill((v_vtx1_plus - v_vtx0_plus).Perp());
	h_pdedx_k.Fill(v_dedx_plus, v_p_plus);
	h_piv.Fill(v_dedx_piv_plus[0], v_p_piv_plus[0]);
	h_piv.Fill(v_dedx_piv_plus[1], v_p_piv_plus[1]);
	h_pib.Fill(v_dedx_pib_plus, v_p_pib_plus);
      }
    }

    if (flag_minus) {
      h_kf_chi2.Fill(kf_chi2_minus);
      if (kf_chi2_minus < 100) {
	h_inm_ksminv.Fill(v_in_mks_minus, ksminv[0]);
	h_kfm_ksminv.Fill(v_kf_mks_minus, ksminv[0]);
	h_ksminv.Fill(ksminv[0]);
	h_in_mks.Fill(v_in_mks_minus);
	h_kf_mks.Fill(v_kf_mks_minus);
	h_vtx0_x.Fill(v_vtx0_minus.X());
	h_vtx0_y.Fill(v_vtx0_minus.Y());
	h_vtx0_z.Fill(v_vtx0_minus.Z());
	h_vtx1_x.Fill(v_vtx1_minus.X());
	h_vtx1_y.Fill(v_vtx1_minus.Y());
	h_vtx1_z.Fill(v_vtx1_minus.Z());
	h_vtx0_r.Fill(v_vtx0_minus.Perp());
	h_vtx1_dr.Fill((v_vtx1_minus - v_vtx0_minus).Perp());
	h_pdedx_k.Fill(v_p_minus, v_dedx_minus);
	h_piv.Fill(v_dedx_piv_minus[0], v_p_piv_minus[0]);
	h_piv.Fill(v_dedx_piv_minus[1], v_p_piv_minus[1]);
	h_pib.Fill(v_dedx_pib_minus, v_p_pib_minus);
      }
    }
    
  }
  h_pdedx_k.SetMarkerColor(2);
  h_piv.SetMarkerColor(4);
  h_pib.SetMarkerColor(3);
  h_pdedx_k.SetMarkerStyle(6);
  h_piv.SetMarkerStyle(6);
  h_pib.SetMarkerStyle(6);
  outfl->cd();
  h_kf_chi2.Write();
  h_in_mks.Write();
  h_kf_mks.Write();
  h_vtx0_x.Write();
  h_vtx0_y.Write();
  h_vtx0_z.Write();
  h_vtx1_x.Write();
  h_vtx1_y.Write();
  h_vtx1_z.Write();
  h_vtx0_r.Write();
  h_vtx1_dr.Write();
  h_ksminv.Write();
  h_inm_ksminv.Write();
  h_kfm_ksminv.Write();
  h_pdedx_k.Write();
  h_piv.Write();
  h_pib.Write();
  outfl->Close();
}
