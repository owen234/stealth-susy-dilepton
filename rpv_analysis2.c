#define rpv_analysis2_cxx
#include "rpv_analysis2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "histio.c"
#include "analysis_utils.c"

void rpv_analysis2::Loop(int max_events, bool verb, const char* hist_file )
{

   bool skip_sl_top_events(true) ;

   gDirectory -> Delete( "h*" ) ;

   TH1F* h_stop_pt = new TH1F( "h_stop_pt", "Stop pt", 100, 0., 1000. ) ;
   TH1F* h_top_pt = new TH1F( "h_top_pt", "top pt", 100, 0., 1000. ) ;
   TH2F* h_top1_vs_top2_pt = new TH2F( "h_top1_vs_top2_pt", "top1 vs top2 pt", 100, 0., 1000., 100, 0., 1000. ) ;
   TH1F* h_chi_pt = new TH1F( "h_chi_pt", "Chi pt", 100, 0., 1000. ) ;
   TH1F* h_w_pt = new TH1F( "h_w_pt", "w pt", 100, 0., 1000. ) ;
   TH1F* h_b_pt = new TH1F( "h_b_pt", "b pt", 100, 0., 1000. ) ;
   TH1F* h_w_dh_pt = new TH1F( "h_w_dh_pt", "W daughter pt (harder one)", 100, 0., 1000. ) ;
   TH1F* h_w_ds_pt = new TH1F( "h_w_ds_pt", "W daughter pt (softer one)", 100, 0., 1000. ) ;

   TH1F* h_stop_eta = new TH1F( "h_stop_eta", "Stop eta", 100, -6., 6. ) ;
   TH1F* h_top_eta = new TH1F( "h_top_eta", "top eta", 100, -6., 6.  ) ;
   TH1F* h_chi_eta = new TH1F( "h_chi_eta", "Chi eta", 100, -6., 6.  ) ;
   TH1F* h_w_eta = new TH1F( "h_w_eta", "w eta", 100, -6., 6.  ) ;
   TH1F* h_b_eta = new TH1F( "h_b_eta", "b eta", 100, -6., 6.  ) ;

   TH1F* h_twostop_pt = new TH1F( "h_twostop_pt", "Two stop system, pt", 100, 0., 1000. ) ;
   TH1F* h_twostop_pz = new TH1F( "h_twostop_pz", "Two stop system, pz", 100, -5000., 5000. ) ;

   TH1F* h_stopstop_cosdphi = new TH1F( "h_stopstop_cosdphi", "stop1*stop2 cos dphi", 110, -1.1, 1.1 ) ;
   TH2F* h_stopstop_pt_vs_cosdphi = new TH2F( "h_stopstop_pt_vs_cosdphi", "stop1,stop2 pt vs cos dphi", 100, -1.1, 1.1, 100, 0., 1000. ) ;
   TH1F* h_stopstop_cosdphi_two_taggable_tops = new TH1F( "h_stopstop_cosdphi_two_taggable_tops", "stop1*stop2 cos dphi, two taggable tops", 110, -1.1, 1.1 ) ;
   TH1F* h_stopstop_costhetasum_cmf = new TH1F( "h_stopstop_costhetasum_cmf", "stop1,stop2 cos thetasum, CMF", 110, -1.1, 1.1 ) ;
   TH1F* h_stopstop_costhetasum_cmf_two_taggable_tops = new TH1F( "h_stopstop_costhetasum_cmf_two_taggable_tops", "stop1,stop2 cos thetasum, CMF, two taggable tops", 110, -1.1, 1.1 ) ;

   TH1F* h_drmin_top_dau = new TH1F( "h_drmin_top_dau", "dR min, top daughters", 100, 0., 5. ) ;
   TH1F* h_drmax_top_dau = new TH1F( "h_drmax_top_dau", "dR max, top daughters", 100, 0., 5. ) ;
   TH2F* h_drmin_vs_drmax_top_dau = new TH2F( "h_drmin_vs_drmax_top_dau", "dR min vs dR max, top daughters", 100, 0., 5., 100, 0., 5. ) ;
   TH2F* h_top_dau_drmax_vs_top_pt = new TH2F( "h_top_dau_drmax_vs_top_pt", "top, dR max vs pt", 100, 0., 1000., 100, 0., 5. ) ;
   TH1F* h_drmax_wrt_top_axis = new TH1F( "h_drmax_wrt_top_axis", "Max top daughter Delta R w.r.t. top axis", 100, 0., 5. ) ;
   TH2F* h_drmax_wrt_top_axis_vs_top_pt = new TH2F( "h_drmax_wrt_top_axis_vs_top_pt", "top, dR max wrt top axis vs top pt", 100, 0., 1000., 100, 0., 5. ) ;
   TH2F* h_drmax_wrt_top_axis_vs_top_pt_dpt20 = new TH2F( "h_drmax_wrt_top_axis_vs_top_pt_dpt20", "top, dR max wrt top axis vs top pt, daughter pt>20", 100, 0., 1000., 100, 0., 5. ) ;
   TH2F* h_drmax_wrt_top_axis_vs_top_pt_dpt30 = new TH2F( "h_drmax_wrt_top_axis_vs_top_pt_dpt30", "top, dR max wrt top axis vs top pt, daughter pt>30", 100, 0., 1000., 100, 0., 5. ) ;

   TH1F* h_w_gen_mass_check = new TH1F( "h_w_gen_mass_check", "Check on gen W mass", 100, 0., 150. ) ;
   TH1F* h_w_gen_dsum_mass_check = new TH1F( "h_w_gen_dsum_mass_check", "Check on gen W mass, sum of daughters", 100, 0., 150. ) ;

   TH1F* h_chi_gen_mass_check = new TH1F( "h_chi_gen_mass_check", "Check on gen Chi mass", 100, 0., 150. ) ;
   TH1F* h_chi_gen_dsum_mass_check = new TH1F( "h_chi_gen_dsum_mass_check", "Check on gen Chi mass, sum of daughters", 100, 0., 150. ) ;

   TH1F* h_top_gen_mass_check = new TH1F( "h_top_gen_mass_check", "Check on gen top mass", 100, 0., 300. ) ;
   TH1F* h_top_gen_dsum_mass_check = new TH1F( "h_top_gen_dsum_mass_check", "Check on gen top mass, sum of daughters", 100, 0., 300. ) ;

   TH1F* h_toptop_cosdphi = new TH1F( "h_toptop_cosdphi", "top*top cos dphi", 110, -1.1, 1.1 ) ;
   TH1F* h_toptop_cosdphi_two_taggable_tops = new TH1F( "h_toptop_cosdphi_two_taggable_tops", "top*top cos dphi, two taggable tops", 110, -1.1, 1.1 ) ;
   TH1F* h_toptop_costhetasum_cmf = new TH1F( "h_toptop_costhetasum_cmf", "top,top cos theta sum, CMF", 110, -1.1, 1.1 ) ;
   TH1F* h_toptop_costhetasum_cmf_two_taggable_tops = new TH1F( "h_toptop_costhetasum_cmf_two_taggable_tops", "top,top cos theta sum, CMF", 110, -1.1, 1.1 ) ;
   TH2F* h_toptop_cosdphi_vs_costhetasum_cmf = new TH2F( "h_toptop_cosdphi_vs_costhetasum_cmf", "top,top cos(dphi) vs cos(thetasum) cmf", 100, -1.1, 1.1,  100, -1.1, 1.1 ) ;
   TH2F* h_toptop_cosdphi_vs_costhetasum_cmf_two_taggable_tops = new TH2F( "h_toptop_cosdphi_vs_costhetasum_cmf_two_taggable_tops", "top,top cos(dphi) vs cos(thetasum) cmf, two taggable tops", 100, -1.1, 1.1,  100, -1.1, 1.1 ) ;

   TH1F* h_stop_sum_pz = new TH1F( "h_stop_sum_pz", "stop stop sum pz", 100, -6000, 6000 ) ;

   TH1F* h_n_taggable_tops = new TH1F( "h_n_taggable_tops", "Number of taggable tops in the event", 3, -0.5, 2.5 ) ;

   TH1F* h_top_pt_taggable = new TH1F( "h_top_pt_taggable", "Top pt, taggable", 100, 0., 1000. ) ;

   TH1F* h_extra_jet_pt = new TH1F( "h_extra_jet_pt", "Extra jet, pt", 100, 0., 1000. ) ;
   TH1F* h_extra_jet_eta = new TH1F( "h_extra_jet_eta", "Extra jet, eta", 100, -6., 6. ) ;
   TH1F* h_extra_njets = new TH1F( "h_extra_njets", "Extra number of jets", 11, -0.5, 10.5 ) ;

   TH1F* h_extra_jet_rank = new TH1F( "h_extra_jet_rank", "Pt rank of extra jet(s)", 15, 0.5, 15.5 ) ;
   TH2F* h_extra_jet_rank_vs_njets = new TH2F( "h_extra_jet_rank_vs_njets", "Pt rank of extra jet(s) vs Njets", 15, 0.5, 15.5, 15, 0.5, 15.5 ) ;
   TH1F* h_extra_lead_jet_rank = new TH1F( "h_extra_lead_jet_rank", "Pt rank of leading extra jet", 15, 0.5, 15.5 ) ;



  //+++++++++ gen only above here, reco below here

   TH1F* h_rec_top_mass = new TH1F( "h_rec_top_mass", "Reconstructed top mass", 100, 0., 350. ) ;
   TH1F* h_rec_w_mass = new TH1F( "h_rec_w_mass", "Reconstructed W mass", 100, 0., 350. ) ;

   TH1F* h_rec_chi_mass_3dau = new TH1F( "h_rec_chi_mass_3dau", "Reconstructed chi mass, 3 daughters", 100, 0., 200. ) ;
   TH1F* h_rec_chi_mass_2dau = new TH1F( "h_rec_chi_mass_2dau", "Reconstructed chi mass, 2 daughters", 100, 0., 200. ) ;

   TH1F* h_rec_stop_mass_complete = new TH1F( "h_rec_stop_mass_complete", "Reconstructed stop mass, complete", 100, 0., 700. ) ;
   TH1F* h_rec_stop_mass_partial  = new TH1F( "h_rec_stop_mass_partial", "Reconstructed stop mass, partial", 100, 0., 700. ) ;
   TH1F* h_rec_stop_mass_top_complete_chi_partial  = new TH1F( "h_rec_stop_mass_top_complete_chi_partial", "Reconstructed stop mass, top_complete_chi_partial", 100, 0., 700. ) ;

   TH2F* h_rec_ndau_top1_vs_ndau_top2 = new TH2F( "h_rec_ndau_top1_vs_ndau_top2", "Reconstructed Ndaughters top1 vs top2", 4, -0.5, 3.5,  4, -0.5, 3.5 ) ;
   TH2F* h_rec_ndau_chi1_vs_ndau_chi2 = new TH2F( "h_rec_ndau_chi1_vs_ndau_chi2", "Reconstructed Ndaughters chi1 vs chi2", 4, -0.5, 3.5,  4, -0.5, 3.5 ) ;
   TH2F* h_rec_ndau_stop1_vs_ndau_stop2 = new TH2F( "h_rec_ndau_stop1_vs_ndau_stop2", "Reconstructed Ndaughters stop1 vs stop2", 7, -0.5, 6.5,  7, -0.5, 6.5 ) ;

   TH1F* h_all_jets_sum_pz = new TH1F( "h_all_jets_sum_pz", "Sum of all jets, pz", 100, -6000, 6000 ) ;

   TH2F* h_sumpz_alljets_vs_genstops = new TH2F( "h_sumpz_alljets_vs_genstops", "Sum pz, all rec jets vs two gen stops", 100, -6000., 6000., 100, -6000., 6000. ) ;

   TH1F* h_rec_w_dau_pt = new TH1F( "h_rec_w_dau_pt", "W daughter pt", 100., 0., 600. ) ;
   TH1F* h_rec_chi_dau_pt = new TH1F( "h_rec_chi_dau_pt", "Chi daughter pt", 100., 0., 600. ) ;
   TH1F* h_rec_chi_dau_min_pt = new TH1F( "h_rec_chi_dau_min_pt", "Chi daughter min pt", 100., 0., 600. ) ;
   TH1F* h_rec_chi_dau_max_pt = new TH1F( "h_rec_chi_dau_max_pt", "Chi daughter max pt", 100., 0., 600. ) ;

   TH2F* h_ltbeta_recjets_vs_genstops = new TH2F( "h_ltbeta_recjets_vs_genstops", "LT beta: rec jets vs gen stops", 120, -1.1, 1.1,  120, -1.1, 1.1 ) ;
   TH2F* h_ltbeta_recjets_vs_genstops_dsumpzcutcheck = new TH2F( "h_ltbeta_recjets_vs_genstops_dsumpzcutcheck", "LT beta: rec jets vs gen stops, delta sumpz < 100", 120, -1.1, 1.1,  120, -1.1, 1.1 ) ;
   TH2F* h_ltbeta_recjets_vs_genstops_two_taggable_tops = new TH2F( "h_ltbeta_recjets_vs_genstops_two_taggable_tops", "LT beta: rec jets vs gen stops, two taggable tops", 120, -1.1, 1.1,  120, -1.1, 1.1 ) ;

   TH1F* h_rec_unmatched_jet_pt = new TH1F( "h_rec_unmatched_jet_pt", "Unmatched jet, pt", 100, 0., 1000. ) ;
   TH1F* h_rec_matched_jet_pt = new TH1F( "h_rec_matched_jet_pt", "matched jet, pt", 100, 0., 1000. ) ;
   TH1F* h_rec_unmatched_jet_eta = new TH1F( "h_rec_unmatched_jet_eta", "Unmatched jet, eta", 100, -5., 5. ) ;
   TH1F* h_rec_matched_jet_eta = new TH1F( "h_rec_matched_jet_eta", "matched jet, eta", 100, -5., 5. ) ;

   TH2F* h_rec_njet_unmatched_vs_matched = new TH2F( "h_rec_njet_unmatched_vs_matched", "Reconstructed jets, N unmatched vs N matched", 16, -0.5, 15.5,  16, -0.5, 15.5 ) ;
   TH2F* h_rec_njet_unmatched_vs_matched_pt40 = new TH2F( "h_rec_njet_unmatched_vs_matched_pt40", "Reconstructed jets, N unmatched vs N matched, pt>40", 16, -0.5, 15.5,  16, -0.5, 15.5 ) ;
   TH2F* h_rec_njet_unmatched_vs_matched_pt60 = new TH2F( "h_rec_njet_unmatched_vs_matched_pt60", "Reconstructed jets, N unmatched vs N matched, pt>60", 16, -0.5, 15.5,  16, -0.5, 15.5 ) ;
   TH2F* h_rec_njet_unmatched_vs_matched_eta24 = new TH2F( "h_rec_njet_unmatched_vs_matched_eta24", "Reconstructed jets, N unmatched vs N matched, |eta|<2.4", 16, -0.5, 15.5,  16, -0.5, 15.5 ) ;
   TH2F* h_rec_njet_unmatched_vs_matched_eta24_pt40 = new TH2F( "h_rec_njet_unmatched_vs_matched_eta24_pt40", "Reconstructed jets, N unmatched vs N matched, |eta|<2.4, pt>40", 16, -0.5, 15.5,  16, -0.5, 15.5 ) ;

   TH1F* h_rec_njet20 = new TH1F( "h_rec_njet20", "Njets, pt>20", 21, -0.5, 20.5 ) ;
   TH1F* h_rec_njet30 = new TH1F( "h_rec_njet30", "Njets, pt>30", 21, -0.5, 20.5 ) ;
   TH1F* h_rec_njet40 = new TH1F( "h_rec_njet40", "Njets, pt>40", 21, -0.5, 20.5 ) ;
   TH1F* h_rec_njet32 = new TH1F( "h_rec_njet32", "Njets, pt>32", 21, -0.5, 20.5 ) ;

   TH1F* h_rec_njet40_ht450 = new TH1F( "h_rec_njet40_ht450", "Njets, pt>40, HT>450", 21, -0.5, 20.5 ) ;
   TH1F* h_rec_njet32_ht380 = new TH1F( "h_rec_njet32_ht380", "Njets, pt>32, HT>380", 21, -0.5, 20.5 ) ;

   TH1F* h_rec_njet40_ht1100 = new TH1F( "h_rec_njet40_ht1100", "Njets, pt>40, HT>1100", 21, -0.5, 20.5 ) ;
   TH1F* h_rec_njet32_ht1100 = new TH1F( "h_rec_njet32_ht1100", "Njets, pt>32, HT>1100", 21, -0.5, 20.5 ) ;

   TH1F* h_rec_njet20_nottopdau = new TH1F( "h_rec_njet20_nottopdau", "Njets, pt>20, top daughters not included", 21, -0.5, 20.5 ) ;
   TH1F* h_rec_njet30_nottopdau = new TH1F( "h_rec_njet30_nottopdau", "Njets, pt>30, top daughters not included", 21, -0.5, 20.5 ) ;
   TH1F* h_rec_njet40_nottopdau = new TH1F( "h_rec_njet40_nottopdau", "Njets, pt>40, top daughters not included", 21, -0.5, 20.5 ) ;

   TH1F* h_rec_ht = new TH1F( "h_rec_ht", "HT", 100, 0., 6000. ) ;
   TH1F* h_rec_ht_njet40ge6 = new TH1F( "h_rec_ht_njet40ge6", "HT, Njet40>=6", 100, 0., 6000. ) ;
   TH1F* h_rec_ht_njet32ge6 = new TH1F( "h_rec_ht_njet32ge6", "HT, Njet32>=6", 100, 0., 6000. ) ;


   TH2F* h_rec_chi_dalitz = new TH2F( "h_rec_chi_dalitz", "Reconstructed jets, chi dalitz plot", 100, 0., 150.*150., 100, 0., 150.*150. ) ;

   TH2F* h_jet_rec_pt_vs_gen_pt = new TH2F( "h_jet_rec_pt_vs_gen_pt", "Jets rec pt vs gen pt", 100, 0., 500.,  100., 0., 500 ) ;

   TH1F* h_rec_toptop_cosdphi = new TH1F( "h_rec_toptop_cosdphi", "Reconstructed top top dphi", 110, -1.1, 1.1 ) ;
   TH1F* h_rec_toptop_cosdphi_two_taggable_tops = new TH1F( "h_rec_toptop_cosdphi_two_taggable_tops", "Reconstructed top top dphi, two taggable tops", 110, -1.1, 1.1 ) ;

   TH1F* h_rec_toptop_costhetasum_rec_cmf = new TH1F( "h_rec_toptop_costhetasum_rec_cmf", "Reconstructed top top cos(thetasum), rec CMF", 110, -1.1, 1.1 ) ;
   TH1F* h_rec_toptop_costhetasum_rec_cmf_two_taggable_tops = new TH1F( "h_rec_toptop_costhetasum_rec_cmf_two_taggable_tops", "Reconstructed top top cos(thetasum), rec CMF, two taggable tops", 110, -1.1, 1.1 ) ;
   TH2F* h_rec_toptop_cosdphi_vs_costhetasum_rec_cmf = new TH2F( "h_rec_toptop_cosdphi_vs_costhetasum_rec_cmf", "Reconstructed top,top cos(dphi) vs cos(thetasum) rec cmf", 100, -1.1, 1.1,  100, -1.1, 1.1 ) ;
   TH2F* h_rec_toptop_cosdphi_vs_costhetasum_rec_cmf_two_taggable_tops = new TH2F( "h_rec_toptop_cosdphi_vs_costhetasum_rec_cmf_two_taggable_tops", "Reconstructed top,top cos(dphi) vs cos(thetasum) rec cmf, two taggable tops", 100, -1.1, 1.1,  100, -1.1, 1.1 ) ;




   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   if ( max_events > 0 ) nentries = max_events ;

   Long64_t nbytes = 0, nb = 0;

   int n_missing_b(0) ;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


      if ( jentry % (nentries/10) == 0 ) printf("  Event %9llu / %9llu  (%2.0f%%)\n", jentry, nentries, 100*(jentry*1.)/(nentries*1.) ) ;










      if ( verb ) printf("\n\n =========== number %9llu : Run %9u , Lumi %9u , Event %9llu\n", jentry, RunNum, LumiBlockNum, EvtNum ) ;

      if ( verb ) printf("\n\n GenParticles: %lu\n", GenParticles->size() ) ;

      int stop1_gpi(-1) ;
      int stop2_gpi(-1) ;
      int chi1_gpi(-1) ;
      int chi2_gpi(-1) ;
      int top1_gpi(-1) ;
      int top2_gpi(-1) ;
      int w1_gpi(-1) ;
      int w2_gpi(-1) ;
      int b1_gpi(-1) ;
      int b2_gpi(-1) ;
      int w1_d1_gpi(-1) ;
      int w1_d2_gpi(-1) ;
      int w2_d1_gpi(-1) ;
      int w2_d2_gpi(-1) ;
      int chi1_d1_gpi(-1) ;
      int chi1_d2_gpi(-1) ;
      int chi1_d3_gpi(-1) ;
      int chi2_d1_gpi(-1) ;
      int chi2_d2_gpi(-1) ;
      int chi2_d3_gpi(-1) ;

      vector<int> extras_gpi ;
      vector<int> extras_rji ;

      bool skip_event(false) ;
      for ( unsigned int gpi=0; gpi < GenParticles->size() ; gpi++ ) {

         int spdgid = GenParticles_PdgId->at(gpi) ;
         int pdgid = abs( GenParticles_PdgId->at(gpi) ) ;
         int smomid = GenParticles_ParentId->at(gpi) ;
         int momid = abs( GenParticles_ParentId->at(gpi) ) ;
         int momidx = GenParticles_ParentIdx->at(gpi) ;

         if ( stop1_gpi < 0 && (pdgid < 6 || pdgid == 21 ) && GenParticles->at(gpi).Pt() > 1. ) { extras_gpi.emplace_back( gpi ) ; extras_rji.emplace_back(-1) ; }

         if ( spdgid ==  1000006 && stop1_gpi < 0 ) stop1_gpi = gpi ;
         if ( spdgid == -1000006 && stop2_gpi < 0 ) stop2_gpi = gpi ;

         if ( spdgid == 1000022 && momidx >= 0 && stop1_gpi >= 0 && momidx == stop1_gpi && chi1_gpi < 0 ) chi1_gpi = gpi ;
         if ( spdgid == 1000022 && momidx >= 0 && stop2_gpi >= 0 && momidx == stop2_gpi && chi2_gpi < 0 ) chi2_gpi = gpi ;

         if ( spdgid ==  6 && momidx >= 0 && stop1_gpi >= 0 && top1_gpi < 0 ) top1_gpi = gpi ;
         if ( spdgid == -6 && momidx >= 0 && stop2_gpi >= 0 && top2_gpi < 0 ) top2_gpi = gpi ;

         if ( pdgid < 6 && chi1_gpi >= 0 && momidx == chi1_gpi ) {
            if ( chi1_d1_gpi < 0 ) {
               chi1_d1_gpi = gpi ;
            } else if ( chi1_d2_gpi < 0 ) {
               chi1_d2_gpi = gpi ;
            } else if ( chi1_d3_gpi < 0 ) {
               chi1_d3_gpi = gpi ;
            }
         }

         if ( pdgid < 6 && chi2_gpi >= 0 && momidx >= 0 && momidx == chi2_gpi ) {
            if ( chi2_d1_gpi < 0 ) {
               chi2_d1_gpi = gpi ;
            } else if ( chi2_d2_gpi < 0 ) {
               chi2_d2_gpi = gpi ;
            } else if ( chi2_d3_gpi < 0 ) {
               chi2_d3_gpi = gpi ;
            }
         }


         if ( spdgid ==  24 && smomid == 6 && w1_gpi < 0 ) { w1_gpi = gpi ; }
         if ( spdgid == -24 && smomid ==-6 && w2_gpi < 0 ) { w2_gpi = gpi ; }


         if ( spdgid ==  5 && b1_gpi < 0 && GenParticles->at(gpi).Pt() > 1. ) { b1_gpi = gpi ; }
         if ( spdgid == -5 && b2_gpi < 0 && GenParticles->at(gpi).Pt() > 1. ) { b2_gpi = gpi ; }


         if ( pdgid < 6 && smomid == 24 ) {
            if ( w1_d1_gpi < 0 ) {
               w1_d1_gpi = gpi ;
            } else if ( w1_d2_gpi < 0 ) {
               w1_d2_gpi = gpi ;
            }
         }

         if ( pdgid < 6 && smomid == -24 ) {
            if ( w2_d1_gpi < 0 ) {
               w2_d1_gpi = gpi ;
            } else if ( w2_d2_gpi < 0 ) {
               w2_d2_gpi = gpi ;
            }
         }

         if ( skip_sl_top_events && ( ( pdgid==11 || pdgid==13 || pdgid==15 ) && momid == 24 ) ) { skip_event = true ; }

         if ( verb ) {
            char pname[100] ;
            char mname[100] ;
            sprintf( pname, "%s", mcname( GenParticles_PdgId->at(gpi) ) ) ;
            sprintf( mname, "%s", mcname( GenParticles_ParentId->at(gpi) ) ) ;
            double eta = 99. ;
            if ( GenParticles->at(gpi).Pt() > 0 ) eta = GenParticles->at(gpi).Eta() ;
            double phi = GenParticles->at(gpi).Phi() ;
            double pt = GenParticles->at(gpi).Pt() ;
            printf("  %3u :  ID=%9d %10s : MomID=%9d %10s MomIdx=%3d status=%2d :  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f,  px,py,pz,E = %6.1f, %6.1f, %6.1f,   %6.1f\n",
                gpi,
                GenParticles_PdgId->at(gpi), pname,
                GenParticles_ParentId->at(gpi), mname, GenParticles_ParentIdx->at(gpi),
                GenParticles_Status->at(gpi),
                GenParticles->at(gpi).Pt(),
                eta,
                phi,
                GenParticles->at(gpi).Px(),
                GenParticles->at(gpi).Py(),
                GenParticles->at(gpi).Pz(),
                GenParticles->at(gpi).E()
                ) ;
         } // verbose?

      } // gpi



      if ( verb ) {
         printf( "\n" ) ;
         printf( " stop1 = %3d ,  top1 = %3d ,  w1 = %3d , chi1 = %d, b1 = %d, w1,d1 = %d w1,d2 = %d\n", stop1_gpi, top1_gpi, w1_gpi, chi1_gpi, b1_gpi, w1_d1_gpi, w1_d2_gpi ) ;
         printf( " stop2 = %3d ,  top2 = %3d ,  w2 = %3d , chi2 = %d, b2 = %d, w2,d1 = %d w2,d2 = %d\n", stop2_gpi, top2_gpi, w2_gpi, chi2_gpi, b2_gpi, w2_d1_gpi, w2_d2_gpi ) ;
         printf( " chi1 daughters = %3d, %3d, %3d\n", chi1_d1_gpi, chi1_d2_gpi, chi1_d3_gpi ) ;
         printf( " chi2 daughters = %3d, %3d, %3d\n", chi2_d1_gpi, chi2_d2_gpi, chi2_d3_gpi ) ;
      }
      if ( skip_event ) {
         if ( verb ) printf("\n\n *** skipping SL W decay event.\n\n") ;
         continue ;
      }

      if ( stop1_gpi < 0 || stop2_gpi < 0 ) {
         printf("\n\n *** Couldn't find one of the stops.\n") ;
         printf( " stop1 = %3d ,  top1 = %3d ,  w1 = %3d , chi1 = %d\n", stop1_gpi, top1_gpi, w1_gpi, chi1_gpi ) ;
         printf( " stop2 = %3d ,  top2 = %3d ,  w2 = %3d , chi2 = %d\n", stop2_gpi, top2_gpi, w2_gpi, chi2_gpi ) ;
         continue ;
      }
      if ( top1_gpi < 0 || top2_gpi < 0 ) {
         printf("\n\n *** Couldn't find one of the tops.\n") ;
         printf( " stop1 = %3d ,  top1 = %3d ,  w1 = %3d , chi1 = %d\n", stop1_gpi, top1_gpi, w1_gpi, chi1_gpi ) ;
         printf( " stop2 = %3d ,  top2 = %3d ,  w2 = %3d , chi2 = %d\n", stop2_gpi, top2_gpi, w2_gpi, chi2_gpi ) ;
         continue ;
      }
      if ( w1_gpi < 0 || w2_gpi < 0 ) {
         printf("\n\n *** Couldn't find one of the Ws.\n") ;
         printf( " stop1 = %3d ,  top1 = %3d ,  w1 = %3d , chi1 = %d\n", stop1_gpi, top1_gpi, w1_gpi, chi1_gpi ) ;
         printf( " stop2 = %3d ,  top2 = %3d ,  w2 = %3d , chi2 = %d\n", stop2_gpi, top2_gpi, w2_gpi, chi2_gpi ) ;
         for ( unsigned int gpi=0; gpi < GenParticles->size() ; gpi++ ) {
            char pname[100] ;
            char mname[100] ;
            sprintf( pname, "%s", mcname( GenParticles_PdgId->at(gpi) ) ) ;
            sprintf( mname, "%s", mcname( GenParticles_ParentId->at(gpi) ) ) ;
            double eta = 99. ;
            if ( GenParticles->at(gpi).Pt() > 0 ) eta = GenParticles->at(gpi).Eta() ;
            double phi = GenParticles->at(gpi).Phi() ;
            double pt = GenParticles->at(gpi).Pt() ;
            printf("  %3u :  ID=%9d %10s : MomID=%9d %10s MomIdx=%3d status=%2d :  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n",
                gpi,
                GenParticles_PdgId->at(gpi), pname,
                GenParticles_ParentId->at(gpi), mname, GenParticles_ParentIdx->at(gpi),
                GenParticles_Status->at(gpi),
                GenParticles->at(gpi).Pt(),
                eta,
                phi
                ) ;
         } // gpi
         continue ;
      }
      if ( chi1_gpi < 0 || chi2_gpi < 0 ) {
         printf("\n\n *** Couldn't find one of the chis.\n") ;
         printf( " stop1 = %3d ,  top1 = %3d ,  w1 = %3d , chi1 = %d\n", stop1_gpi, top1_gpi, w1_gpi, chi1_gpi ) ;
         printf( " stop2 = %3d ,  top2 = %3d ,  w2 = %3d , chi2 = %d\n", stop2_gpi, top2_gpi, w2_gpi, chi2_gpi ) ;
         continue ;
      }
      if ( b1_gpi < 0 || b2_gpi < 0 ) {
     /// printf("\n\n *** Couldn't find one of the bs.\n") ;
     /// printf( " stop1 = %3d ,  top1 = %3d ,  w1 = %3d , chi1 = %d, b1 = %d\n", stop1_gpi, top1_gpi, w1_gpi, chi1_gpi, b1_gpi ) ;
     /// printf( " stop2 = %3d ,  top2 = %3d ,  w2 = %3d , chi2 = %d, b2 = %d\n", stop2_gpi, top2_gpi, w2_gpi, chi2_gpi, b2_gpi ) ;
     /// for ( unsigned int gpi=0; gpi < GenParticles->size() ; gpi++ ) {
     ///    char pname[100] ;
     ///    char mname[100] ;
     ///    sprintf( pname, "%s", mcname( GenParticles_PdgId->at(gpi) ) ) ;
     ///    sprintf( mname, "%s", mcname( GenParticles_ParentId->at(gpi) ) ) ;
     ///    double eta = 99. ;
     ///    if ( GenParticles->at(gpi).Pt() > 0 ) eta = GenParticles->at(gpi).Eta() ;
     ///    double phi = GenParticles->at(gpi).Phi() ;
     ///    double pt = GenParticles->at(gpi).Pt() ;
     ///    printf("  %3u :  ID=%9d %10s : MomID=%9d %10s MomIdx=%3d status=%2d :  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n",
     ///        gpi,
     ///        GenParticles_PdgId->at(gpi), pname,
     ///        GenParticles_ParentId->at(gpi), mname, GenParticles_ParentIdx->at(gpi),
     ///        GenParticles_Status->at(gpi),
     ///        GenParticles->at(gpi).Pt(),
     ///        eta,
     ///        phi
     ///        ) ;
     /// } // gpi
         n_missing_b++ ;
         continue ;
      }
      if ( w1_d1_gpi < 0 || w1_d2_gpi < 0 ) {
         printf("\n\n *** can't find one of the W1 daughters.\n\n") ;
         printf( " stop1 = %3d ,  top1 = %3d ,  w1 = %3d , chi1 = %d, b1 = %d, w1,d1 = %d w1,d2 = %d\n", stop1_gpi, top1_gpi, w1_gpi, chi1_gpi, b1_gpi, w1_d1_gpi, w1_d2_gpi ) ;
         printf( " stop2 = %3d ,  top2 = %3d ,  w2 = %3d , chi2 = %d, b2 = %d, w2,d1 = %d w2,d2 = %d\n", stop2_gpi, top2_gpi, w2_gpi, chi2_gpi, b2_gpi, w2_d1_gpi, w2_d2_gpi ) ;
         continue ;
      }

      if ( w2_d1_gpi < 0 || w2_d2_gpi < 0 ) {
         printf("\n\n *** can't find one of the W2 daughters.\n\n") ;
         printf( " stop1 = %3d ,  top1 = %3d ,  w1 = %3d , chi1 = %d, b1 = %d, w1,d1 = %d w1,d2 = %d\n", stop1_gpi, top1_gpi, w1_gpi, chi1_gpi, b1_gpi, w1_d1_gpi, w1_d2_gpi ) ;
         printf( " stop2 = %3d ,  top2 = %3d ,  w2 = %3d , chi2 = %d, b2 = %d, w2,d1 = %d w2,d2 = %d\n", stop2_gpi, top2_gpi, w2_gpi, chi2_gpi, b2_gpi, w2_d1_gpi, w2_d2_gpi ) ;
         continue ;
      }

      if ( chi1_d1_gpi < 0 || chi1_d2_gpi < 0 || chi1_d3_gpi < 0 ) {
         printf("\n\n *** can't find one of the expected 3 daughters of chi1.\n\n") ;
         printf( " chi1 daughters = %3d, %3d, %3d\n", chi1_d1_gpi, chi1_d2_gpi, chi1_d3_gpi ) ;
         continue ;
      }
      if ( chi2_d1_gpi < 0 || chi2_d2_gpi < 0 || chi2_d3_gpi < 0 ) {
         printf("\n\n *** can't find one of the expected 3 daughters of chi2.\n\n") ;
         printf( " chi2 daughters = %3d, %3d, %3d\n", chi2_d1_gpi, chi2_d2_gpi, chi2_d3_gpi ) ;
         continue ;
      }





    //--- Approximate a hadronic trigger

      int rec_njet_pt20(0) ;
      int rec_njet_pt30(0) ;
      int rec_njet_pt40(0) ;
      int rec_njet_pt32(0) ;

      for ( unsigned int rji=0; rji < Jets->size() ; rji++ ) {

            TLorentzVector jlv( Jets->at(rji) ) ;

            if ( jlv.Pt() > 20 ) rec_njet_pt20++ ;
            if ( jlv.Pt() > 30 ) rec_njet_pt30++ ;
            if ( jlv.Pt() > 40 ) rec_njet_pt40++ ;
            if ( jlv.Pt() > 32 ) rec_njet_pt32++ ;

      } // rji

      h_rec_njet20 -> Fill( rec_njet_pt20 ) ;
      h_rec_njet30 -> Fill( rec_njet_pt30 ) ;

      h_rec_njet40 -> Fill( rec_njet_pt40 ) ;
      if ( HT > 450 ) h_rec_njet40_ht450 -> Fill( rec_njet_pt40 ) ;
      if ( HT > 1100 ) h_rec_njet40_ht1100 -> Fill( rec_njet_pt40 ) ;

      h_rec_njet32 -> Fill( rec_njet_pt32 ) ;
      if ( HT > 380 ) h_rec_njet32_ht380 -> Fill( rec_njet_pt32 ) ;
      if ( HT > 1100 ) h_rec_njet32_ht1100 -> Fill( rec_njet_pt32 ) ;

      h_rec_ht -> Fill( HT ) ;
      if ( rec_njet_pt40 >= 6 ) h_rec_ht_njet40ge6 -> Fill( HT ) ;
      if ( rec_njet_pt32 >= 6 ) h_rec_ht_njet32ge6 -> Fill( HT ) ;

      if ( !( HT>450 && rec_njet_pt40>=6 ) ) continue ;

      //if ( HT < 1100 ) continue ; //************ HT cut *****************








      TLorentzVector gplv_stop1( GenParticles->at(stop1_gpi) ) ;
      TLorentzVector gplv_stop2( GenParticles->at(stop2_gpi) ) ;

      TLorentzVector gplv_top1( GenParticles->at(top1_gpi) ) ;
      TLorentzVector gplv_top2( GenParticles->at(top2_gpi) ) ;

      TLorentzVector gplv_chi1( GenParticles->at(chi1_gpi) ) ;
      TLorentzVector gplv_chi2( GenParticles->at(chi2_gpi) ) ;

      TLorentzVector gplv_b1( GenParticles->at(b1_gpi) ) ;
      TLorentzVector gplv_b2( GenParticles->at(b2_gpi) ) ;

      TLorentzVector gplv_w1( GenParticles->at(w1_gpi) ) ;
      TLorentzVector gplv_w2( GenParticles->at(w2_gpi) ) ;

      TLorentzVector gplv_w1_d1( GenParticles->at(w1_d1_gpi) ) ;
      TLorentzVector gplv_w1_d2( GenParticles->at(w1_d2_gpi) ) ;

      TLorentzVector gplv_w2_d1( GenParticles->at(w2_d1_gpi) ) ;
      TLorentzVector gplv_w2_d2( GenParticles->at(w2_d2_gpi) ) ;

      TLorentzVector gplv_two_stop_system = gplv_stop1 + gplv_stop2 ;

      TLorentzVector gplv_w1_dsum = gplv_w1_d1 + gplv_w1_d2 ; // for sanity check
      TLorentzVector gplv_w2_dsum = gplv_w2_d1 + gplv_w2_d2 ; // for sanity check

      TLorentzVector gplv_chi1_d1( GenParticles->at(chi1_d1_gpi) ) ;
      TLorentzVector gplv_chi1_d2( GenParticles->at(chi1_d2_gpi) ) ;
      TLorentzVector gplv_chi1_d3( GenParticles->at(chi1_d3_gpi) ) ;

      TLorentzVector gplv_chi2_d1( GenParticles->at(chi2_d1_gpi) ) ;
      TLorentzVector gplv_chi2_d2( GenParticles->at(chi2_d2_gpi) ) ;
      TLorentzVector gplv_chi2_d3( GenParticles->at(chi2_d3_gpi) ) ;

      TLorentzVector gplv_chi1_dsum = gplv_chi1_d1 + gplv_chi1_d2 + gplv_chi1_d3 ;
      TLorentzVector gplv_chi2_dsum = gplv_chi2_d1 + gplv_chi2_d2 + gplv_chi2_d3 ;

      TLorentzVector gplv_top1_dsum = gplv_b1 + gplv_w1 ;
      TLorentzVector gplv_top2_dsum = gplv_b2 + gplv_w2 ;


      TLorentzVector gplv_w1_dh ;
      TLorentzVector gplv_w1_ds ;
      if ( gplv_w1_d1.Pt() > gplv_w1_d2.Pt() ) {
         gplv_w1_dh = gplv_w1_d1 ;
         gplv_w1_ds = gplv_w1_d2 ;
      } else {
         gplv_w1_dh = gplv_w1_d2 ;
         gplv_w1_ds = gplv_w1_d1 ;
      }

      TLorentzVector gplv_w2_dh ;
      TLorentzVector gplv_w2_ds ;
      if ( gplv_w2_d1.Pt() > gplv_w2_d2.Pt() ) {
         gplv_w2_dh = gplv_w2_d1 ;
         gplv_w2_ds = gplv_w2_d2 ;
      } else {
         gplv_w2_dh = gplv_w2_d2 ;
         gplv_w2_ds = gplv_w2_d1 ;
      }




      if (verb) {
         if ( gplv_top1.Pt() > 200 && gplv_top2.Pt() > 200 ) { printf(" *** two hard tops.\n") ; }
      }

      if (verb) {
       //-- sanity checks
         printf("  stop1 mass = %.1f , chi1 mass = %.1f\n", gplv_stop1.M(), gplv_chi1.M() ) ;
         printf("  stop2 mass = %.1f , chi2 mass = %.1f\n", gplv_stop2.M(), gplv_chi2.M() ) ;
         printf("    Two stop system:  pt = %.1f ,  pz = %.1f\n", gplv_two_stop_system.Pt(), gplv_two_stop_system.Pz() ) ;
         printf("  W1 mass = %.1f, W2 mass = %.1f\n", gplv_w1.M(), gplv_w2.M() ) ;
         printf("  W1 daughter sum, M = %.1f\n", gplv_w1_dsum.M() ) ;
         printf("  W2 daughter sum, M = %.1f\n", gplv_w2_dsum.M() ) ;
         printf("  W1 genparticle  px,py,pz,E: %6.1f, %6.1f, %6.1f, %6.1f\n", gplv_w1.Px(), gplv_w1.Py(), gplv_w1.Pz(), gplv_w1.E() ) ;
         printf("  W1 daughter sum px,py,pz,E: %6.1f, %6.1f, %6.1f, %6.1f\n", gplv_w1_dsum.Px(), gplv_w1_dsum.Py(), gplv_w1_dsum.Pz(), gplv_w1_dsum.E() ) ;
         printf("  W2 genparticle  px,py,pz,E: %6.1f, %6.1f, %6.1f, %6.1f\n", gplv_w2.Px(), gplv_w2.Py(), gplv_w2.Pz(), gplv_w2.E() ) ;
         printf("  W2 daughter sum px,py,pz,E: %6.1f, %6.1f, %6.1f, %6.1f\n", gplv_w2_dsum.Px(), gplv_w2_dsum.Py(), gplv_w2_dsum.Pz(), gplv_w2_dsum.E() ) ;
         printf("  chi genparticle daughter sum masses : %.1f , %.1f\n", gplv_chi1_dsum.M(), gplv_chi2_dsum.M() ) ;
      }


      h_stop_sum_pz -> Fill( gplv_stop1.Pz() + gplv_stop2.Pz() ) ;

      h_stop_pt -> Fill( gplv_stop1.Pt() ) ;
      h_stop_pt -> Fill( gplv_stop2.Pt() ) ;

      h_top_pt -> Fill( gplv_top1.Pt() ) ;
      h_top_pt -> Fill( gplv_top2.Pt() ) ;
      h_top1_vs_top2_pt -> Fill( gplv_top2.Pt(), gplv_top1.Pt() ) ;

      h_chi_pt -> Fill( gplv_chi1.Pt() ) ;
      h_chi_pt -> Fill( gplv_chi2.Pt() ) ;

      h_b_pt -> Fill( gplv_b1.Pt() ) ;
      h_b_pt -> Fill( gplv_b2.Pt() ) ;

      h_w_pt -> Fill( gplv_w1.Pt() ) ;
      h_w_pt -> Fill( gplv_w2.Pt() ) ;

      h_w_dh_pt -> Fill( gplv_w1_dh.Pt() ) ;
      h_w_dh_pt -> Fill( gplv_w2_dh.Pt() ) ;

      h_w_ds_pt -> Fill( gplv_w1_ds.Pt() ) ;
      h_w_ds_pt -> Fill( gplv_w2_ds.Pt() ) ;


      TLorentzVector gplv_extras ;
      h_extra_njets -> Fill( extras_gpi.size() ) ;
      for ( int ei=0; ei<extras_gpi.size(); ei++ ) {
         h_extra_jet_pt -> Fill( GenParticles->at( extras_gpi[ei] ).Pt() ) ;
         h_extra_jet_eta -> Fill( GenParticles->at( extras_gpi[ei] ).Eta() ) ;
         gplv_extras += GenParticles->at( extras_gpi[ei] ) ;
      } // ei

      if ( verb ) {
         printf("  extra jets (%lu) 4-vector :  Pt = %6.1f , Phi = %7.3f\n", extras_gpi.size(), gplv_extras.Pt(), gplv_extras.Phi() ) ;
         printf("  sum of stops            :  Pt = %6.1f , Phi = %7.3f\n", gplv_two_stop_system.Pt(), gplv_two_stop_system.Phi() ) ;
         if ( extras_gpi.size() > 0 ) {
            double dphi = gplv_extras.Phi() - gplv_two_stop_system.Phi() ;
            if ( dphi >  3.14159265 ) dphi -= 2*3.14159265 ;
            if ( dphi < -3.14159265 ) dphi += 2*3.14159265 ;
            printf("  delta phi : %7.3f, cos(dphi) = %6.3f\n", dphi, cos(dphi) ) ;
         }
      }



      h_stop_eta -> Fill( gplv_stop1.Eta() ) ;
      h_stop_eta -> Fill( gplv_stop2.Eta() ) ;

      h_top_eta -> Fill( gplv_top1.Eta() ) ;
      h_top_eta -> Fill( gplv_top2.Eta() ) ;

      h_chi_eta -> Fill( gplv_chi1.Eta() ) ;
      h_chi_eta -> Fill( gplv_chi2.Eta() ) ;

      h_b_eta -> Fill( gplv_b1.Eta() ) ;
      h_b_eta -> Fill( gplv_b2.Eta() ) ;

      h_w_eta -> Fill( gplv_w1.Eta() ) ;
      h_w_eta -> Fill( gplv_w2.Eta() ) ;




      h_twostop_pt -> Fill( gplv_two_stop_system.Pt() ) ;
      h_twostop_pz -> Fill( gplv_two_stop_system.Pz() ) ;


      double stopstop_dphi = gplv_stop2.Phi() - gplv_stop1.Phi() ;
      if ( stopstop_dphi > 3.1415926 ) stopstop_dphi -= 2*3.1415926 ;
      if ( stopstop_dphi <-3.1415926 ) stopstop_dphi += 2*3.1415926 ;
      double stopstop_cosdphi = cos( stopstop_dphi ) ;
      h_stopstop_cosdphi -> Fill( stopstop_cosdphi ) ;
      h_stopstop_pt_vs_cosdphi -> Fill( stopstop_cosdphi, gplv_two_stop_system.Pt() ) ;

      if ( verb ) {
         printf(" stop,stop : pt = %.1f , phi = %.3f, cosdphi = %.3f\n", gplv_two_stop_system.Pt(), gplv_two_stop_system.Phi(), stopstop_cosdphi ) ;
      }




      double toptop_dphi = gplv_top2.Phi() - gplv_top1.Phi() ;
      if ( toptop_dphi > 3.1415926 ) toptop_dphi -= 2*3.1415926 ;
      if ( toptop_dphi <-3.1415926 ) toptop_dphi += 2*3.1415926 ;
      double toptop_cosdphi = cos( toptop_dphi ) ;

      h_toptop_cosdphi -> Fill( toptop_cosdphi ) ;




      h_w_gen_mass_check -> Fill( gplv_w1.M() ) ;
      h_w_gen_mass_check -> Fill( gplv_w2.M() ) ;

      h_w_gen_dsum_mass_check -> Fill( gplv_w1_dsum.M() ) ;
      h_w_gen_dsum_mass_check -> Fill( gplv_w2_dsum.M() ) ;


      h_chi_gen_mass_check -> Fill( gplv_chi1.M() ) ;
      h_chi_gen_mass_check -> Fill( gplv_chi2.M() ) ;

      h_chi_gen_dsum_mass_check -> Fill( gplv_chi1_dsum.M() ) ;
      h_chi_gen_dsum_mass_check -> Fill( gplv_chi2_dsum.M() ) ;


      h_top_gen_mass_check -> Fill( gplv_top1.M() ) ;
      h_top_gen_mass_check -> Fill( gplv_top2.M() ) ;

      h_top_gen_dsum_mass_check -> Fill( gplv_top1_dsum.M() ) ;
      h_top_gen_dsum_mass_check -> Fill( gplv_top2_dsum.M() ) ;



      double top1_drmin(99.) ;
      double top1_drmax(0.) ;
      double top1_drmax_wrt_top_axis(0.) ;
      {
         double dr ;
         dr = gplv_w1_d1.DeltaR( gplv_w1_d2 ) ;
         if ( dr > top1_drmax ) top1_drmax = dr ;
         if ( dr < top1_drmin ) top1_drmin = dr ;
         dr = gplv_w1_d1.DeltaR( gplv_b1 ) ;
         if ( dr > top1_drmax ) top1_drmax = dr ;
         if ( dr < top1_drmin ) top1_drmin = dr ;
         dr = gplv_w1_d2.DeltaR( gplv_b1 ) ;
         if ( dr > top1_drmax ) top1_drmax = dr ;
         if ( dr < top1_drmin ) top1_drmin = dr ;

         dr = gplv_top1.DeltaR( gplv_b1 ) ;
         if ( verb ) printf("  top1 pt=%6.1f : b quark pt=%6.1f dR w.r.t. top axis: %.3f\n", gplv_top1.Pt(), gplv_b1.Pt(), dr ) ;
         if ( dr > top1_drmax_wrt_top_axis ) top1_drmax_wrt_top_axis = dr ;
         dr = gplv_top1.DeltaR( gplv_w1_d1 ) ;
         if ( verb ) printf("  top1 pt=%6.1f : Wd1     pt=%6.1f dR w.r.t. top axis: %.3f\n", gplv_top1.Pt(), gplv_w1_d1.Pt(), dr ) ;
         if ( dr > top1_drmax_wrt_top_axis ) top1_drmax_wrt_top_axis = dr ;
         dr = gplv_top1.DeltaR( gplv_w1_d2 ) ;
         if ( verb ) printf("  top1 pt=%6.1f : Wd2     pt=%6.1f dR w.r.t. top axis: %.3f\n", gplv_top1.Pt(), gplv_w1_d2.Pt(), dr ) ;
         if ( dr > top1_drmax_wrt_top_axis ) top1_drmax_wrt_top_axis = dr ;
      }


      h_drmin_top_dau -> Fill( top1_drmin ) ;
      h_drmax_top_dau -> Fill( top1_drmax ) ;
      h_top_dau_drmax_vs_top_pt -> Fill( gplv_top1.Pt(), top1_drmax ) ;
      h_drmin_vs_drmax_top_dau -> Fill( top1_drmax, top1_drmin ) ;
      h_drmax_wrt_top_axis -> Fill( top1_drmax_wrt_top_axis ) ;
      h_drmax_wrt_top_axis_vs_top_pt -> Fill(  gplv_top1.Pt(), top1_drmax_wrt_top_axis ) ;
      if ( gplv_b1.Pt() > 20 && gplv_w1_d1.Pt() > 20 && gplv_w1_d2.Pt() > 20 ) {
         h_drmax_wrt_top_axis_vs_top_pt_dpt20 -> Fill(  gplv_top1.Pt(), top1_drmax_wrt_top_axis ) ;
      }
      if ( gplv_b1.Pt() > 30 && gplv_w1_d1.Pt() > 30 && gplv_w1_d2.Pt() > 30 ) {
         h_drmax_wrt_top_axis_vs_top_pt_dpt30 -> Fill(  gplv_top1.Pt(), top1_drmax_wrt_top_axis ) ;
      }


      double top2_drmin(99.) ;
      double top2_drmax(0.) ;
      double top2_drmax_wrt_top_axis(0.) ;
      {
         double dr ;
         dr = gplv_w2_d1.DeltaR( gplv_w2_d2 ) ;
         if ( dr > top2_drmax ) top2_drmax = dr ;
         if ( dr < top2_drmin ) top2_drmin = dr ;
         dr = gplv_w2_d1.DeltaR( gplv_b2 ) ;
         if ( dr > top2_drmax ) top2_drmax = dr ;
         if ( dr < top2_drmin ) top2_drmin = dr ;
         dr = gplv_w2_d2.DeltaR( gplv_b2 ) ;
         if ( dr > top2_drmax ) top2_drmax = dr ;
         if ( dr < top2_drmin ) top2_drmin = dr ;

         dr = gplv_top2.DeltaR( gplv_b2 ) ;
         if ( verb ) printf("  top2 pt=%6.1f : b quark pt=%6.1f dR w.r.t. top axis: %.3f\n", gplv_top2.Pt(), gplv_b2.Pt(), dr ) ;
         if ( dr > top2_drmax_wrt_top_axis ) top2_drmax_wrt_top_axis = dr ;
         dr = gplv_top2.DeltaR( gplv_w2_d1 ) ;
         if ( verb ) printf("  top2 pt=%6.1f : Wd1     pt=%6.1f dR w.r.t. top axis: %.3f\n", gplv_top2.Pt(), gplv_w2_d1.Pt(), dr ) ;
         if ( dr > top2_drmax_wrt_top_axis ) top2_drmax_wrt_top_axis = dr ;
         dr = gplv_top2.DeltaR( gplv_w2_d2 ) ;
         if ( verb ) printf("  top2 pt=%6.1f : Wd2     pt=%6.1f dR w.r.t. top axis: %.3f\n", gplv_top2.Pt(), gplv_w2_d2.Pt(), dr ) ;
         if ( dr > top2_drmax_wrt_top_axis ) top2_drmax_wrt_top_axis = dr ;
      }

      h_drmin_top_dau -> Fill( top2_drmin ) ;
      h_drmax_top_dau -> Fill( top2_drmax ) ;
      h_top_dau_drmax_vs_top_pt -> Fill( gplv_top2.Pt(), top2_drmax ) ;
      h_drmin_vs_drmax_top_dau -> Fill( top2_drmax, top2_drmin ) ;
      h_drmax_wrt_top_axis -> Fill( top2_drmax_wrt_top_axis ) ;
      h_drmax_wrt_top_axis_vs_top_pt -> Fill(  gplv_top2.Pt(), top2_drmax_wrt_top_axis ) ;
      if ( gplv_b2.Pt() > 20 && gplv_w2_d1.Pt() > 20 && gplv_w2_d2.Pt() > 20 ) {
         h_drmax_wrt_top_axis_vs_top_pt_dpt20 -> Fill(  gplv_top2.Pt(), top2_drmax_wrt_top_axis ) ;
      }
      if ( gplv_b2.Pt() > 30 && gplv_w2_d1.Pt() > 30 && gplv_w2_d2.Pt() > 30 ) {
         h_drmax_wrt_top_axis_vs_top_pt_dpt30 -> Fill(  gplv_top2.Pt(), top2_drmax_wrt_top_axis ) ;
      }


      double true_stopstop_beta = ( gplv_stop1.Pz() + gplv_stop2.Pz() ) / ( gplv_stop1.E() + gplv_stop2.E() ) ;
      TVector3 boost_beta_vec( 0., 0., -1.*true_stopstop_beta ) ;
      TLorentzVector gplv_stop1_cmf( gplv_stop1 ) ;
      TLorentzVector gplv_stop2_cmf( gplv_stop2 ) ;
      gplv_stop1_cmf.Boost( boost_beta_vec ) ;
      gplv_stop2_cmf.Boost( boost_beta_vec ) ;
      double stopstop_thetasum_cmf = gplv_stop1_cmf.Theta() + gplv_stop2_cmf.Theta() ;
      h_stopstop_costhetasum_cmf -> Fill( cos( stopstop_thetasum_cmf ) ) ;
      if ( verb ) {
         printf("  Stop pz1, pz2 : %6.1f , %6.1f\n", gplv_stop1.Pz(), gplv_stop2.Pz() ) ;
         printf("  beta z : %.3f\n", true_stopstop_beta ) ;
         printf("  Lab frame, stop1, px, py, pz:  %6.1f, %6.1f, %6.1f\n", gplv_stop1.Px(), gplv_stop1.Py(), gplv_stop1.Pz() ) ;
         printf("  CM  frame, stop1, px, py, pz:  %6.1f, %6.1f, %6.1f\n", gplv_stop1_cmf.Px(), gplv_stop1_cmf.Py(), gplv_stop1_cmf.Pz() ) ;
         printf("  Lab frame, stop2, px, py, pz:  %6.1f, %6.1f, %6.1f\n", gplv_stop2.Px(), gplv_stop2.Py(), gplv_stop2.Pz() ) ;
         printf("  CM  frame, stop2, px, py, pz:  %6.1f, %6.1f, %6.1f\n", gplv_stop2_cmf.Px(), gplv_stop2_cmf.Py(), gplv_stop2_cmf.Pz() ) ;
         printf("  stop,stop CMF theta1, theta2, theta1+theta2, cos(thetasum) : %7.3f, %7.3f, %7.3f, %6.3f\n",
             gplv_stop1_cmf.Theta(), gplv_stop2_cmf.Theta(), gplv_stop1_cmf.Theta()+gplv_stop2_cmf.Theta(), cos( stopstop_thetasum_cmf ) ) ;
      }


      TLorentzVector gplv_top1_cmf( gplv_top1 ) ;
      TLorentzVector gplv_top2_cmf( gplv_top2 ) ;
      gplv_top1_cmf.Boost( boost_beta_vec ) ;
      gplv_top2_cmf.Boost( boost_beta_vec ) ;
      double toptop_thetasum_cmf = gplv_top1_cmf.Theta() + gplv_top2_cmf.Theta() ;
      h_toptop_costhetasum_cmf -> Fill( cos(toptop_thetasum_cmf) ) ;
      h_toptop_cosdphi_vs_costhetasum_cmf -> Fill( cos( toptop_thetasum_cmf ), cos( toptop_dphi ) ) ;
      if ( verb ) {
         printf("  top,top CMS theta1, theta2, theta1+theta2, cos(thetasum) : %7.3f, %7.3f, %7.3f, %6.3f\n",
             gplv_top1_cmf.Theta(), gplv_top2_cmf.Theta(), gplv_top1_cmf.Theta()+gplv_top2_cmf.Theta(), cos( toptop_thetasum_cmf ) ) ;
      }





      bool top1_taggable(true) ;
      if ( gplv_b1.Pt() < 20 || fabs(gplv_b1.Eta()) > 2.4 ) top1_taggable = false ;
      if ( gplv_w1_d1.Pt() < 20 ) top1_taggable = false ;
      if ( gplv_w1_d2.Pt() < 20 ) top1_taggable = false ;
      if ( gplv_b1.DeltaR( gplv_top1 ) > 1.5 ) top1_taggable = false ;
      if ( gplv_w1_d1.DeltaR( gplv_top1 ) > 1.5 ) top1_taggable = false ;
      if ( gplv_w1_d2.DeltaR( gplv_top1 ) > 1.5 ) top1_taggable = false ;

      bool top2_taggable(true) ;
      if ( gplv_b2.Pt() < 20 || fabs(gplv_b2.Eta()) > 2.4 ) top2_taggable = false ;
      if ( gplv_w2_d1.Pt() < 20 ) top2_taggable = false ;
      if ( gplv_w2_d2.Pt() < 20 ) top2_taggable = false ;
      if ( gplv_b2.DeltaR( gplv_top2 ) > 1.5 ) top2_taggable = false ;
      if ( gplv_w2_d1.DeltaR( gplv_top2 ) > 1.5 ) top2_taggable = false ;
      if ( gplv_w2_d2.DeltaR( gplv_top2 ) > 1.5 ) top2_taggable = false ;

      int n_taggable_tops(0) ;
      if ( top1_taggable ) n_taggable_tops++ ;
      if ( top2_taggable ) n_taggable_tops++ ;

      h_n_taggable_tops -> Fill( n_taggable_tops ) ;

      if ( verb ) printf("\n Number %9llu : Number of taggable tops : %d\n\n", jentry, n_taggable_tops ) ;


      if ( top1_taggable ) h_top_pt_taggable -> Fill( gplv_top1.Pt() ) ;
      if ( top2_taggable ) h_top_pt_taggable -> Fill( gplv_top2.Pt() ) ;

      if ( n_taggable_tops == 2 ) {
         h_toptop_cosdphi_two_taggable_tops -> Fill( toptop_cosdphi ) ;
         h_toptop_costhetasum_cmf_two_taggable_tops -> Fill( cos(toptop_thetasum_cmf) ) ;
         h_toptop_cosdphi_vs_costhetasum_cmf_two_taggable_tops -> Fill( cos( toptop_thetasum_cmf ), cos( toptop_dphi ) ) ;
         h_stopstop_cosdphi_two_taggable_tops -> Fill( stopstop_cosdphi ) ;
         h_stopstop_costhetasum_cmf_two_taggable_tops -> Fill( cos( stopstop_thetasum_cmf ) ) ;
      }


    //+++++++++ Above here is gen-level only.  Below adds reco jets.  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++








      int b1_rji(-1) ;
      int b2_rji(-1) ;
      int w1_d1_rji(-1) ;
      int w1_d2_rji(-1) ;
      int w2_d1_rji(-1) ;
      int w2_d2_rji(-1) ;

      int chi1_d1_rji(-1) ;
      int chi1_d2_rji(-1) ;
      int chi1_d3_rji(-1) ;

      int chi2_d1_rji(-1) ;
      int chi2_d2_rji(-1) ;
      int chi2_d3_rji(-1) ;



      for ( unsigned int rji=0; rji < Jets->size() ; rji++ ) {

         double dr ;

         TLorentzVector jlv( Jets->at(rji) ) ;

         dr = gplv_b1.DeltaR( jlv ) ;
         if ( dr < 0.1 ) b1_rji = rji ;

         dr = gplv_b2.DeltaR( jlv ) ;
         if ( dr < 0.1 ) b2_rji = rji ;

         dr = gplv_w1_d1.DeltaR( jlv ) ;
         if ( dr < 0.1 ) w1_d1_rji = rji ;

         dr = gplv_w1_d2.DeltaR( jlv ) ;
         if ( dr < 0.1 ) w1_d2_rji = rji ;

         dr = gplv_w2_d1.DeltaR( jlv ) ;
         if ( dr < 0.1 ) w2_d1_rji = rji ;

         dr = gplv_w2_d2.DeltaR( jlv ) ;
         if ( dr < 0.1 ) w2_d2_rji = rji ;

         dr = gplv_chi1_d1.DeltaR( jlv ) ;
         if ( dr < 0.1 ) chi1_d1_rji = rji ;

         dr = gplv_chi1_d2.DeltaR( jlv ) ;
         if ( dr < 0.1 ) chi1_d2_rji = rji ;

         dr = gplv_chi1_d3.DeltaR( jlv ) ;
         if ( dr < 0.1 ) chi1_d3_rji = rji ;

         dr = gplv_chi2_d1.DeltaR( jlv ) ;
         if ( dr < 0.1 ) chi2_d1_rji = rji ;

         dr = gplv_chi2_d2.DeltaR( jlv ) ;
         if ( dr < 0.1 ) chi2_d2_rji = rji ;

         dr = gplv_chi2_d3.DeltaR( jlv ) ;
         if ( dr < 0.1 ) chi2_d3_rji = rji ;

         for ( int ei=0; ei < extras_gpi.size() ; ei++ ) {
            dr = GenParticles->at( extras_gpi[ei] ).DeltaR( jlv ) ;
            if ( dr < 0.1 ) extras_rji[ei] = rji ;
         } // ei


      } // rji.


      int rec_njets_matched(0) ;
      int rec_njets_unmatched(0) ;

      int rec_njets_matched_pt40(0) ;
      int rec_njets_unmatched_pt40(0) ;

      int rec_njets_matched_pt60(0) ;
      int rec_njets_unmatched_pt60(0) ;

      int rec_njets_matched_eta24(0) ;
      int rec_njets_unmatched_eta24(0) ;

      int rec_njets_matched_eta24_pt40(0) ;
      int rec_njets_unmatched_eta24_pt40(0) ;

      int njets_pt20_nottopdau(0) ;
      int njets_pt30_nottopdau(0) ;
      int njets_pt40_nottopdau(0) ;

      for ( unsigned int rji=0; rji < Jets->size() ; rji++ ) {

         TLorentzVector jlv( Jets->at(rji) ) ;

         bool unmatched(true) ;
         bool istopdau(false) ;
         if ( rji == b1_rji      ) { unmatched = false ; istopdau = true ; }
         if ( rji == b2_rji      ) { unmatched = false ; istopdau = true ; }
         if ( rji == w1_d1_rji   ) { unmatched = false ; istopdau = true ; }
         if ( rji == w1_d2_rji   ) { unmatched = false ; istopdau = true ; }
         if ( rji == w2_d1_rji   ) { unmatched = false ; istopdau = true ; }
         if ( rji == w2_d2_rji   ) { unmatched = false ; istopdau = true ; }
         if ( rji == chi1_d1_rji ) unmatched = false ;
         if ( rji == chi1_d2_rji ) unmatched = false ;
         if ( rji == chi1_d3_rji ) unmatched = false ;
         if ( rji == chi2_d1_rji ) unmatched = false ;
         if ( rji == chi2_d2_rji ) unmatched = false ;
         if ( rji == chi2_d3_rji ) unmatched = false ;

         if ( !istopdau ) {
            if ( jlv.Pt() > 20 ) njets_pt20_nottopdau++ ;
            if ( jlv.Pt() > 30 ) njets_pt30_nottopdau++ ;
            if ( jlv.Pt() > 40 ) njets_pt40_nottopdau++ ;
         }

         if ( unmatched ) {
            rec_njets_unmatched ++ ;
            h_rec_unmatched_jet_pt -> Fill( jlv.Pt() ) ;
            h_rec_unmatched_jet_eta -> Fill( jlv.Eta() ) ;
            if ( jlv.Pt() > 40 ) rec_njets_unmatched_pt40++ ;
            if ( jlv.Pt() > 60 ) rec_njets_unmatched_pt60++ ;
            if ( fabs(jlv.Eta()) < 2.4 ) rec_njets_unmatched_eta24++ ;
            if ( fabs(jlv.Eta()) < 2.4 && jlv.Pt() > 40 ) rec_njets_unmatched_eta24_pt40++ ;
         } else {
            rec_njets_matched ++ ;
            h_rec_matched_jet_pt -> Fill( jlv.Pt() ) ;
            h_rec_matched_jet_eta -> Fill( jlv.Eta() ) ;
            if ( jlv.Pt() > 40 ) rec_njets_matched_pt40++ ;
            if ( jlv.Pt() > 60 ) rec_njets_matched_pt60++ ;
            if ( fabs(jlv.Eta()) < 2.4 ) rec_njets_matched_eta24++ ;
            if ( fabs(jlv.Eta()) < 2.4 && jlv.Pt() > 40 ) rec_njets_matched_eta24_pt40++ ;
         }

      } // rji

      h_rec_njet20_nottopdau -> Fill( njets_pt20_nottopdau ) ;
      h_rec_njet30_nottopdau -> Fill( njets_pt30_nottopdau ) ;
      h_rec_njet40_nottopdau -> Fill( njets_pt40_nottopdau ) ;

      h_rec_njet_unmatched_vs_matched -> Fill( rec_njets_matched, rec_njets_unmatched ) ;
      h_rec_njet_unmatched_vs_matched_pt40 -> Fill( rec_njets_matched_pt40, rec_njets_unmatched_pt40 ) ;
      h_rec_njet_unmatched_vs_matched_pt60 -> Fill( rec_njets_matched_pt60, rec_njets_unmatched_pt60 ) ;
      h_rec_njet_unmatched_vs_matched_eta24 -> Fill( rec_njets_matched_eta24, rec_njets_unmatched_eta24 ) ;
      h_rec_njet_unmatched_vs_matched_eta24_pt40 -> Fill( rec_njets_matched_eta24_pt40, rec_njets_unmatched_eta24_pt40 ) ;





      for ( unsigned int rji=0; rji < Jets->size() ; rji++ ) {

         TLorentzVector jlv( Jets->at(rji) ) ;

         if ( rji == b1_rji ) h_jet_rec_pt_vs_gen_pt -> Fill( gplv_b1.Pt(), jlv.Pt() ) ;
         if ( rji == b2_rji ) h_jet_rec_pt_vs_gen_pt -> Fill( gplv_b2.Pt(), jlv.Pt() ) ;
         if ( rji == w1_d1_rji ) h_jet_rec_pt_vs_gen_pt -> Fill( gplv_w1_d1.Pt(), jlv.Pt() ) ;
         if ( rji == w1_d2_rji ) h_jet_rec_pt_vs_gen_pt -> Fill( gplv_w1_d2.Pt(), jlv.Pt() ) ;
         if ( rji == w2_d1_rji ) h_jet_rec_pt_vs_gen_pt -> Fill( gplv_w2_d1.Pt(), jlv.Pt() ) ;
         if ( rji == w2_d2_rji ) h_jet_rec_pt_vs_gen_pt -> Fill( gplv_w2_d2.Pt(), jlv.Pt() ) ;
         if ( rji == chi1_d1_rji ) h_jet_rec_pt_vs_gen_pt -> Fill( gplv_chi1_d1.Pt(), jlv.Pt() ) ;
         if ( rji == chi1_d2_rji ) h_jet_rec_pt_vs_gen_pt -> Fill( gplv_chi1_d2.Pt(), jlv.Pt() ) ;
         if ( rji == chi1_d3_rji ) h_jet_rec_pt_vs_gen_pt -> Fill( gplv_chi1_d3.Pt(), jlv.Pt() ) ;
         if ( rji == chi2_d1_rji ) h_jet_rec_pt_vs_gen_pt -> Fill( gplv_chi2_d1.Pt(), jlv.Pt() ) ;
         if ( rji == chi2_d2_rji ) h_jet_rec_pt_vs_gen_pt -> Fill( gplv_chi2_d2.Pt(), jlv.Pt() ) ;
         if ( rji == chi2_d3_rji ) h_jet_rec_pt_vs_gen_pt -> Fill( gplv_chi2_d3.Pt(), jlv.Pt() ) ;

      } // rji



      for ( unsigned int rji=0; rji < Jets->size() ; rji++ ) {
         for ( int ei=0; ei<extras_rji.size() ; ei++ ) {
            if ( rji == extras_rji[ei] ) {
               if ( ei==0 ) {
                  h_extra_lead_jet_rank -> Fill( rji+1 ) ;
               }
               h_extra_jet_rank -> Fill( rji+1 ) ;
               h_extra_jet_rank_vs_njets -> Fill( Jets->size(), rji+1 ) ;
            }
         } // ei
      } // rji



      int w1_fji(-1) ;
      int w2_fji(-1) ;
      int top1_fji(-1) ;
      int top2_fji(-1) ;
      int chi1_fji(-1) ;
      int chi2_fji(-1) ;

      for ( int fji=0; fji < JetsAK8->size(); fji++ ) {

         double dr ;

         TLorentzVector fjlv( JetsAK8->at(fji) ) ;

         dr = gplv_w1.DeltaR( fjlv ) ;
         if ( dr < 0.2 ) w1_fji = fji ;

         dr = gplv_w2.DeltaR( fjlv ) ;
         if ( dr < 0.2 ) w2_fji = fji ;

         dr = gplv_top1.DeltaR( fjlv ) ;
         if ( dr < 0.2 ) top1_fji = fji ;

         dr = gplv_top2.DeltaR( fjlv ) ;
         if ( dr < 0.2 ) top2_fji = fji ;

         dr = gplv_chi1.DeltaR( fjlv ) ;
         if ( dr < 0.2 ) chi1_fji = fji ;

         dr = gplv_chi2.DeltaR( fjlv ) ;
         if ( dr < 0.2 ) chi2_fji = fji ;

      } // fji








      if ( verb ) {

         printf( "\n top tags:\n") ;
         for ( unsigned int ti=0; ti < toptag_tlv->size() ; ti ++ ) {
            TLorentzVector tlv( toptag_tlv->at(ti) ) ;
            printf("  tt %3d :  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f, Mass = %6.1f,  Nconstituents = %d\n",
              ti,
                tlv.Pt(),
                tlv.Eta(),
                tlv.Phi(),
                tlv.M(),
                toptag_nconstituents->at(ti)
            ) ;
         } // ti

         printf("\n\n AK4 jets:\n") ;
         for ( unsigned int rji=0; rji < Jets->size() ; rji++ ) {
            TLorentzVector jlv( Jets->at(rji) ) ;
            printf("  %3d :  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f",
              rji,
                Jets->at(rji).Pt(),
                Jets->at(rji).Eta(),
                Jets->at(rji).Phi()
            ) ;
            if ( Jets_toptag_index->at(rji) >= 0 ) {
               printf(" : toptag %d ", Jets_toptag_index->at(rji) ) ;
            } else {
               printf("            ") ;
            }
            if ( rji == b1_rji      ) printf(" : b1     %.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_b1.DeltaR( jlv ), gplv_b1.Pt(), gplv_b1.Eta(), gplv_b1.Phi() ) ;
            if ( rji == b2_rji      ) printf(" : b2     %.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_b2.DeltaR( jlv ), gplv_b2.Pt(), gplv_b2.Eta(), gplv_b2.Phi() ) ;
            if ( rji == w1_d1_rji   ) printf(" : W1d1   %.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_w1_d1.DeltaR( jlv ), gplv_w1_d1.Pt(), gplv_w1_d1.Eta(), gplv_w1_d1.Phi() ) ;
            if ( rji == w1_d2_rji   ) printf(" : W1d2   %.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_w1_d2.DeltaR( jlv ), gplv_w1_d2.Pt(), gplv_w1_d2.Eta(), gplv_w1_d2.Phi() ) ;
            if ( rji == w2_d1_rji   ) printf(" : W2d1   %.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_w2_d1.DeltaR( jlv ), gplv_w2_d1.Pt(), gplv_w2_d1.Eta(), gplv_w2_d1.Phi() ) ;
            if ( rji == w2_d2_rji   ) printf(" : W2d2   %.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_w2_d2.DeltaR( jlv ), gplv_w2_d2.Pt(), gplv_w2_d2.Eta(), gplv_w2_d2.Phi() ) ;
            if ( rji == chi1_d1_rji ) printf(" : Chi1d1 %.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_chi1_d1.DeltaR( jlv ), gplv_chi1_d1.Pt(), gplv_chi1_d1.Eta(), gplv_chi1_d1.Phi() ) ;
            if ( rji == chi1_d2_rji ) printf(" : Chi1d2 %.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_chi1_d2.DeltaR( jlv ), gplv_chi1_d2.Pt(), gplv_chi1_d2.Eta(), gplv_chi1_d2.Phi() ) ;
            if ( rji == chi1_d3_rji ) printf(" : Chi1d3 %.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_chi1_d3.DeltaR( jlv ), gplv_chi1_d3.Pt(), gplv_chi1_d3.Eta(), gplv_chi1_d3.Phi() ) ;
            if ( rji == chi2_d1_rji ) printf(" : Chi2d1 %.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_chi2_d1.DeltaR( jlv ), gplv_chi2_d1.Pt(), gplv_chi2_d1.Eta(), gplv_chi2_d1.Phi() ) ;
            if ( rji == chi2_d2_rji ) printf(" : Chi2d2 %.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_chi2_d2.DeltaR( jlv ), gplv_chi2_d2.Pt(), gplv_chi2_d2.Eta(), gplv_chi2_d2.Phi() ) ;
            if ( rji == chi2_d3_rji ) printf(" : Chi2d3 %.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_chi2_d3.DeltaR( jlv ), gplv_chi2_d3.Pt(), gplv_chi2_d3.Eta(), gplv_chi2_d3.Phi() ) ;
            for ( int ei=0; ei<extras_rji.size() ; ei++ ) {
               if ( rji == extras_rji[ei] ) printf(" : extra  %.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f",
                   GenParticles->at( extras_gpi[ei] ).DeltaR( jlv ),
                   GenParticles->at( extras_gpi[ei] ).Pt(),
                   GenParticles->at( extras_gpi[ei] ).Eta(),
                   GenParticles->at( extras_gpi[ei] ).Phi() ) ;
            }
            printf("\n") ;
         } // rji

         printf("\n\n AK8 jets:\n") ;
         for ( int fji=0; fji < JetsAK8->size(); fji++ ) {

            TLorentzVector fjlv( JetsAK8->at(fji) ) ;

            printf("  %3d :  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f, pruned mass = %6.1f, softdrop mass = %6.1f",
              fji, fjlv.Pt(), fjlv.Eta(), fjlv.Phi(), JetsAK8_prunedMass->at(fji), JetsAK8_softDropMass->at(fji)
            ) ;

            if ( JetsAK8_toptag_index->at(fji) >= 0 ) {
               printf(" : toptag %d ", JetsAK8_toptag_index->at(fji) ) ;
            }
            if ( fji == w1_fji   ) printf(" : W1     dR=%.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_w1.DeltaR( fjlv ), gplv_w1.Pt(), gplv_w1.Eta(), gplv_w1.Phi() ) ;
            if ( fji == w2_fji   ) printf(" : W2     dR=%.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_w2.DeltaR( fjlv ), gplv_w2.Pt(), gplv_w2.Eta(), gplv_w2.Phi() ) ;
            if ( fji == top1_fji ) printf(" : top1   dR=%.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_top1.DeltaR( fjlv ), gplv_top1.Pt(), gplv_top1.Eta(), gplv_top1.Phi() ) ;
            if ( fji == top2_fji ) printf(" : top2   dR=%.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_top2.DeltaR( fjlv ), gplv_top2.Pt(), gplv_top2.Eta(), gplv_top2.Phi() ) ;
            if ( fji == chi1_fji ) printf(" : chi1   dR=%.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_chi1.DeltaR( fjlv ), gplv_chi1.Pt(), gplv_chi1.Eta(), gplv_chi1.Phi() ) ;
            if ( fji == chi2_fji ) printf(" : chi2   dR=%.3f  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f", gplv_chi2.DeltaR( fjlv ), gplv_chi2.Pt(), gplv_chi2.Eta(), gplv_chi2.Phi() ) ;

            printf("\n") ;

         } // fji

         printf("\n  Unmatched objects:\n") ;
         if ( b1_rji < 0 && w1_d1_rji < 0 &&  w1_d2_rji < 0 ) printf("  t1       Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_top1.Pt(), gplv_top1.Eta(), gplv_top1.Phi() ) ;
         if ( b2_rji < 0 && w2_d1_rji < 0 &&  w2_d2_rji < 0 ) printf("  t2       Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_top2.Pt(), gplv_top2.Eta(), gplv_top2.Phi() ) ;
         if ( w1_d1_rji < 0 &&  w1_d2_rji < 0 ) printf("  W1       Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_w1.Pt(), gplv_w1.Eta(), gplv_w1.Phi() ) ;
         if ( w2_d1_rji < 0 &&  w2_d2_rji < 0 ) printf("  W2       Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_w2.Pt(), gplv_w2.Eta(), gplv_w2.Phi() ) ;
         if ( chi1_d1_rji < 0 && chi1_d2_rji < 0 && chi1_d3_rji < 0 ) printf("  Chi1     Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_chi1.Pt(), gplv_chi1.Eta(), gplv_chi1.Phi() ) ;
         if ( chi2_d1_rji < 0 && chi2_d2_rji < 0 && chi2_d3_rji < 0 ) printf("  Chi2     Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_chi2.Pt(), gplv_chi2.Eta(), gplv_chi2.Phi() ) ;
         if ( b1_rji      < 0 ) printf("  b1       Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_b1.Pt(), gplv_b1.Eta(), gplv_b1.Phi() ) ;
         if ( b2_rji      < 0 ) printf("  b2       Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_b2.Pt(), gplv_b2.Eta(), gplv_b2.Phi() ) ;
         if ( w1_d1_rji   < 0 ) printf("  W1d1     Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_w1_d1.Pt(), gplv_w1_d1.Eta(), gplv_w1_d1.Phi() ) ;
         if ( w1_d2_rji   < 0 ) printf("  W1d2     Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_w1_d2.Pt(), gplv_w1_d2.Eta(), gplv_w1_d2.Phi() ) ;
         if ( w2_d1_rji   < 0 ) printf("  W2d1     Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_w2_d1.Pt(), gplv_w2_d1.Eta(), gplv_w2_d1.Phi() ) ;
         if ( w2_d2_rji   < 0 ) printf("  W2d2     Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_w2_d2.Pt(), gplv_w2_d2.Eta(), gplv_w2_d2.Phi() ) ;
         if ( chi1_d1_rji < 0 ) printf("  Chi1d1   Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_chi1_d1.Pt(), gplv_chi1_d1.Eta(), gplv_chi1_d1.Phi() ) ;
         if ( chi1_d2_rji < 0 ) printf("  Chi1d2   Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_chi1_d2.Pt(), gplv_chi1_d2.Eta(), gplv_chi1_d2.Phi() ) ;
         if ( chi1_d3_rji < 0 ) printf("  Chi1d3   Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_chi1_d3.Pt(), gplv_chi1_d3.Eta(), gplv_chi1_d3.Phi() ) ;
         if ( chi2_d1_rji < 0 ) printf("  chi2d1   Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_chi2_d1.Pt(), gplv_chi2_d1.Eta(), gplv_chi2_d1.Phi() ) ;
         if ( chi2_d2_rji < 0 ) printf("  chi2d2   Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_chi2_d2.Pt(), gplv_chi2_d2.Eta(), gplv_chi2_d2.Phi() ) ;
         if ( chi2_d3_rji < 0 ) printf("  chi2d3   Pt = %7.1f , Eta = %6.3f, Phi = %6.3f\n", gplv_chi2_d3.Pt(), gplv_chi2_d3.Eta(), gplv_chi2_d3.Phi() ) ;
         printf("\n") ;



      } // verbose?




      TLorentzVector rlv_top1 ;
      TLorentzVector rlv_top2 ;
      TLorentzVector rlv_w1 ;
      TLorentzVector rlv_w2 ;
      TLorentzVector rlv_chi1 ;
      TLorentzVector rlv_chi2 ;

      TLorentzVector rlv_all_jets ;

      for ( unsigned int rji=0; rji < Jets->size() ; rji++ ) {

         TLorentzVector jlv( Jets->at(rji) ) ;

         rlv_all_jets += jlv ;

         if ( rji == b1_rji    ) rlv_top1 += jlv ;
         if ( rji == w1_d1_rji ) rlv_top1 += jlv ;
         if ( rji == w1_d2_rji ) rlv_top1 += jlv ;

         if ( rji == w1_d1_rji ) rlv_w1 += jlv ;
         if ( rji == w1_d2_rji ) rlv_w1 += jlv ;

         if ( rji == b2_rji    ) rlv_top2 += jlv ;
         if ( rji == w2_d1_rji ) rlv_top2 += jlv ;
         if ( rji == w2_d2_rji ) rlv_top2 += jlv ;

         if ( rji == w2_d1_rji ) rlv_w2 += jlv ;
         if ( rji == w2_d2_rji ) rlv_w2 += jlv ;

         if ( rji == chi1_d1_rji ) rlv_chi1 += jlv ;
         if ( rji == chi1_d2_rji ) rlv_chi1 += jlv ;
         if ( rji == chi1_d3_rji ) rlv_chi1 += jlv ;

         if ( rji == chi2_d1_rji ) rlv_chi2 += jlv ;
         if ( rji == chi2_d2_rji ) rlv_chi2 += jlv ;
         if ( rji == chi2_d3_rji ) rlv_chi2 += jlv ;

      } // rji

      TLorentzVector rlv_stop1 = rlv_top1 + rlv_chi1 ;
      TLorentzVector rlv_stop2 = rlv_top2 + rlv_chi2 ;

      double reco_jets_beta = rlv_all_jets.Pz() / rlv_all_jets.E() ;

      TVector3 rec_boost_beta_vec( 0., 0., -1.*reco_jets_beta ) ;
      TLorentzVector rlv_top1_rec_cmf( rlv_top1 ) ;
      TLorentzVector rlv_top2_rec_cmf( rlv_top2 ) ;
      rlv_top1_rec_cmf.Boost( rec_boost_beta_vec ) ;
      rlv_top2_rec_cmf.Boost( rec_boost_beta_vec ) ;


      if ( verb ) {
         printf("    Reconstructed top1 mass  : %6.1f ,  W1 mass %6.1f\n", rlv_top1.M(), rlv_w1.M() ) ;
         printf("    Reconstructed top2 mass  : %6.1f ,  W2 mass %6.1f\n", rlv_top2.M(), rlv_w2.M() ) ;
         printf("    Reconstructed chi1 mass  : %6.1f\n", rlv_chi1.M() ) ;
         printf("    Reconstructed chi2 mass  : %6.1f\n", rlv_chi2.M() ) ;
         printf("    Reconstructed stop1 mass : %6.1f\n", rlv_stop1.M() ) ;
         printf("    Reconstructed stop2 mass : %6.1f\n", rlv_stop2.M() ) ;
         printf("      Lorentz transformation beta to CM from jets:  %.3f ,  from gen stops %.3f\n", reco_jets_beta, true_stopstop_beta ) ;
      }

      h_ltbeta_recjets_vs_genstops -> Fill( true_stopstop_beta, reco_jets_beta ) ;
      if ( fabs( (gplv_stop1.Pz() + gplv_stop2.Pz()) - rlv_all_jets.Pz() ) < 100. ) h_ltbeta_recjets_vs_genstops_dsumpzcutcheck -> Fill( true_stopstop_beta, reco_jets_beta ) ;
      if ( n_taggable_tops == 2 ) h_ltbeta_recjets_vs_genstops_two_taggable_tops -> Fill( true_stopstop_beta, reco_jets_beta ) ;

      if ( w1_d1_rji >= 0 && w1_d2_rji >= 0 && b1_rji >= 0 ) h_rec_top_mass -> Fill( rlv_top1.M() ) ;
      if ( w1_d1_rji >= 0 && w1_d2_rji >= 0                ) h_rec_w_mass   -> Fill( rlv_w1.M() ) ;
      if ( w2_d1_rji >= 0 && w2_d2_rji >= 0 && b2_rji >= 0 ) h_rec_top_mass -> Fill( rlv_top2.M() ) ;
      if ( w2_d1_rji >= 0 && w2_d2_rji >= 0                ) h_rec_w_mass   -> Fill( rlv_w2.M() ) ;

      int nrd_chi1(0) ;
      if ( chi1_d1_rji >= 0 ) nrd_chi1++ ;
      if ( chi1_d2_rji >= 0 ) nrd_chi1++ ;
      if ( chi1_d3_rji >= 0 ) nrd_chi1++ ;

      int nrd_chi2(0) ;
      if ( chi2_d1_rji >= 0 ) nrd_chi2++ ;
      if ( chi2_d2_rji >= 0 ) nrd_chi2++ ;
      if ( chi2_d3_rji >= 0 ) nrd_chi2++ ;

      int nrd_top1(0) ;
      if ( w1_d1_rji >= 0 ) nrd_top1++ ;
      if ( w1_d2_rji >= 0 ) nrd_top1++ ;
      if ( b1_rji    >= 0 ) nrd_top1++ ;

      int nrd_top2(0) ;
      if ( w2_d1_rji >= 0 ) nrd_top2++ ;
      if ( w2_d2_rji >= 0 ) nrd_top2++ ;
      if ( b2_rji    >= 0 ) nrd_top2++ ;

      if ( nrd_chi1==3 ) h_rec_chi_mass_3dau -> Fill( rlv_chi1.M() ) ;
      if ( nrd_chi2==3 ) h_rec_chi_mass_3dau -> Fill( rlv_chi2.M() ) ;

      if ( nrd_chi1==2 ) h_rec_chi_mass_2dau -> Fill( rlv_chi1.M() ) ;
      if ( nrd_chi2==2 ) h_rec_chi_mass_2dau -> Fill( rlv_chi2.M() ) ;

      if ( nrd_top1==3 && nrd_chi1==3 ) h_rec_stop_mass_complete -> Fill( rlv_stop1.M() ) ;
      if ( nrd_top2==3 && nrd_chi2==3 ) h_rec_stop_mass_complete -> Fill( rlv_stop2.M() ) ;

      if ( !(nrd_top1==3 && nrd_chi1==3) ) h_rec_stop_mass_partial -> Fill( rlv_stop1.M() ) ;
      if ( !(nrd_top2==3 && nrd_chi2==3) ) h_rec_stop_mass_partial -> Fill( rlv_stop2.M() ) ;

      if ( nrd_top1==3 && nrd_chi1<3 ) h_rec_stop_mass_top_complete_chi_partial -> Fill( rlv_stop1.M() ) ;
      if ( nrd_top2==3 && nrd_chi2<3 ) h_rec_stop_mass_top_complete_chi_partial -> Fill( rlv_stop2.M() ) ;

      h_rec_ndau_top1_vs_ndau_top2 -> Fill( nrd_top2, nrd_top1 ) ;
      h_rec_ndau_chi1_vs_ndau_chi2 -> Fill( nrd_chi2, nrd_chi1 ) ;
      h_rec_ndau_stop1_vs_ndau_stop2 -> Fill( nrd_top2+nrd_chi2, nrd_top1+nrd_chi1 ) ;

      h_all_jets_sum_pz -> Fill( rlv_all_jets.Pz() ) ;

      h_sumpz_alljets_vs_genstops -> Fill( gplv_stop1.Pz() + gplv_stop2.Pz(), rlv_all_jets.Pz() ) ;

      if ( w1_d1_rji >= 0 ) h_rec_w_dau_pt -> Fill( Jets->at( w1_d1_rji ).Pt() ) ;
      if ( w1_d2_rji >= 0 ) h_rec_w_dau_pt -> Fill( Jets->at( w1_d2_rji ).Pt() ) ;
      if ( w2_d1_rji >= 0 ) h_rec_w_dau_pt -> Fill( Jets->at( w2_d1_rji ).Pt() ) ;
      if ( w2_d2_rji >= 0 ) h_rec_w_dau_pt -> Fill( Jets->at( w2_d2_rji ).Pt() ) ;


      double chi1_d_rec_minpt(1000.) ;
      double chi1_d_rec_maxpt(0.) ;
      int chi1_d_rec_minpt_rji(-1) ;
      int chi1_d_rec_maxpt_rji(-1) ;

      if ( chi1_d1_rji >= 0 ) {
         double pt = Jets->at( chi1_d1_rji ).Pt() ;
         h_rec_chi_dau_pt -> Fill( pt ) ;
         if ( pt > chi1_d_rec_maxpt ) {
            chi1_d_rec_maxpt = pt ;
            chi1_d_rec_maxpt_rji = chi1_d1_rji ;
         }
         if ( pt < chi1_d_rec_minpt ) {
            chi1_d_rec_minpt = pt ;
            chi1_d_rec_minpt_rji = chi1_d1_rji ;
         }
      }
      if ( chi1_d2_rji >= 0 ) {
         double pt = Jets->at( chi1_d2_rji ).Pt() ;
         h_rec_chi_dau_pt -> Fill( pt ) ;
         if ( pt > chi1_d_rec_maxpt ) {
            chi1_d_rec_maxpt = pt ;
            chi1_d_rec_maxpt_rji = chi1_d2_rji ;
         }
         if ( pt < chi1_d_rec_minpt ) {
            chi1_d_rec_minpt = pt ;
            chi1_d_rec_minpt_rji = chi1_d2_rji ;
         }
      }
      if ( chi1_d3_rji >= 0 ) {
         double pt = Jets->at( chi1_d3_rji ).Pt() ;
         h_rec_chi_dau_pt -> Fill( pt ) ;
         if ( pt > chi1_d_rec_maxpt ) {
            chi1_d_rec_maxpt = pt ;
            chi1_d_rec_maxpt_rji = chi1_d3_rji ;
         }
         if ( pt < chi1_d_rec_minpt ) {
            chi1_d_rec_minpt = pt ;
            chi1_d_rec_minpt_rji = chi1_d3_rji ;
         }
      }
      if ( nrd_chi1 > 0 ) {
         h_rec_chi_dau_min_pt -> Fill( chi1_d_rec_minpt ) ;
         h_rec_chi_dau_max_pt -> Fill( chi1_d_rec_maxpt ) ;
      }



      double chi2_d_rec_minpt(1000.) ;
      double chi2_d_rec_maxpt(0.) ;
      int chi2_d_rec_minpt_rji(-1) ;
      int chi2_d_rec_maxpt_rji(-1) ;

      if ( chi2_d1_rji >= 0 ) {
         double pt = Jets->at( chi2_d1_rji ).Pt() ;
         h_rec_chi_dau_pt -> Fill( pt ) ;
         if ( pt > chi2_d_rec_maxpt ) {
            chi2_d_rec_maxpt = pt ;
            chi2_d_rec_maxpt_rji = chi2_d1_rji ;
         }
         if ( pt < chi2_d_rec_minpt ) {
            chi2_d_rec_minpt = pt ;
            chi2_d_rec_minpt_rji = chi2_d1_rji ;
         }
      }
      if ( chi2_d2_rji >= 0 ) {
         double pt = Jets->at( chi2_d2_rji ).Pt() ;
         h_rec_chi_dau_pt -> Fill( pt ) ;
         if ( pt > chi2_d_rec_maxpt ) {
            chi2_d_rec_maxpt = pt ;
            chi2_d_rec_maxpt_rji = chi2_d2_rji ;
         }
         if ( pt < chi2_d_rec_minpt ) {
            chi2_d_rec_minpt = pt ;
            chi2_d_rec_minpt_rji = chi2_d2_rji ;
         }
      }
      if ( chi2_d3_rji >= 0 ) {
         double pt = Jets->at( chi2_d3_rji ).Pt() ;
         h_rec_chi_dau_pt -> Fill( pt ) ;
         if ( pt > chi2_d_rec_maxpt ) {
            chi2_d_rec_maxpt = pt ;
            chi2_d_rec_maxpt_rji = chi2_d3_rji ;
         }
         if ( pt < chi2_d_rec_minpt ) {
            chi2_d_rec_minpt = pt ;
            chi2_d_rec_minpt_rji = chi2_d3_rji ;
         }
      }
      if ( nrd_chi2 > 0 ) {
         h_rec_chi_dau_min_pt -> Fill( chi2_d_rec_minpt ) ;
         h_rec_chi_dau_max_pt -> Fill( chi2_d_rec_maxpt ) ;
      }





      if ( nrd_chi1==3 && chi1_d1_rji != chi1_d2_rji && chi1_d1_rji != chi1_d3_rji && chi1_d2_rji != chi1_d3_rji) {
         TLorentzVector rlv_chi1_d1 = Jets->at( chi1_d1_rji ) ;
         TLorentzVector rlv_chi1_d2 = Jets->at( chi1_d2_rji ) ;
         TLorentzVector rlv_chi1_d3 = Jets->at( chi1_d3_rji ) ;
         TLorentzVector rlv_chi1_d1d2 = rlv_chi1_d1 + rlv_chi1_d2 ;
         TLorentzVector rlv_chi1_d1d3 = rlv_chi1_d1 + rlv_chi1_d3 ;
         h_rec_chi_dalitz -> Fill( rlv_chi1_d1d2.M2(), rlv_chi1_d1d3.M2() ) ;
      }
      if ( nrd_chi2==3 && chi2_d1_rji != chi2_d2_rji && chi2_d1_rji != chi2_d3_rji && chi2_d2_rji != chi2_d3_rji) {
         TLorentzVector rlv_chi2_d1 = Jets->at( chi2_d1_rji ) ;
         TLorentzVector rlv_chi2_d2 = Jets->at( chi2_d2_rji ) ;
         TLorentzVector rlv_chi2_d3 = Jets->at( chi2_d3_rji ) ;
         TLorentzVector rlv_chi2_d1d2 = rlv_chi2_d1 + rlv_chi2_d2 ;
         TLorentzVector rlv_chi2_d1d3 = rlv_chi2_d1 + rlv_chi2_d3 ;
         h_rec_chi_dalitz -> Fill( rlv_chi2_d1d2.M2(), rlv_chi2_d1d3.M2() ) ;
      }



      if ( nrd_top1 > 0 && nrd_top2 > 0 ) {
         double rec_toptop_dphi = rlv_top1.Phi() - rlv_top2.Phi() ;
         if ( rec_toptop_dphi >  3.1415926 ) rec_toptop_dphi -= 2*3.1415926 ;
         if ( rec_toptop_dphi < -3.1415926 ) rec_toptop_dphi += 2*3.1415926 ;
         h_rec_toptop_cosdphi -> Fill( cos( rec_toptop_dphi ) ) ;
         double rec_toptop_thetasum = rlv_top1_rec_cmf.Theta() + rlv_top2_rec_cmf.Theta() ;
         h_rec_toptop_costhetasum_rec_cmf -> Fill( cos(rec_toptop_thetasum) ) ;
         h_rec_toptop_cosdphi_vs_costhetasum_rec_cmf -> Fill( cos(rec_toptop_thetasum), cos( rec_toptop_dphi ) ) ;
         if ( n_taggable_tops == 2 ) {
            h_rec_toptop_cosdphi_two_taggable_tops -> Fill( cos( rec_toptop_dphi ) ) ;
            h_rec_toptop_costhetasum_rec_cmf_two_taggable_tops -> Fill( cos(rec_toptop_thetasum) ) ;
            h_rec_toptop_cosdphi_vs_costhetasum_rec_cmf_two_taggable_tops -> Fill( cos(rec_toptop_thetasum), cos( rec_toptop_dphi ) ) ;
         }
      }











   } // jentry

   printf("\n\n") ;
   gDirectory -> ls() ;
   printf("\n\n") ;

   saveHist( hist_file, "h*" ) ;

   printf("\n\n") ;
   printf("  Number of events with missing b in gen particles: %d\n", n_missing_b ) ;
   printf("\n\n") ;

} // Loop

//=================================================================================
