#define dump1_cxx
#include "dump1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

char pname[100] ;

const char* mcname( int pdgid ) ;

double calcDr( double eta1, double eta2, double phi1, double phi2 ) ;

void dump1::Loop( int display_event, float rhophi_scale )
{

   bool skip_sl_top_events(true) ;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   TCanvas* can1 = (TCanvas*) gDirectory -> FindObject( "can1" ) ;
   if ( can1 == 0x0 ) {
      can1 = new TCanvas( "can1", "phi vs eta plane", 1200, 800 ) ;
   }

   TCanvas* can2 = (TCanvas*) gDirectory -> FindObject( "can2" ) ;
   if ( can2 == 0x0 ) {
      can2 = new TCanvas( "can2", "GenParticle : rho phi plane", 10, 10, 700, 700 ) ;
   }

   TCanvas* can3 = (TCanvas*) gDirectory -> FindObject( "can3" ) ;
   if ( can3 == 0x0 ) {
      can3 = new TCanvas( "can3", "Reco jets : rho phi plane", 710, 10, 700, 700 ) ;
   }

   TH2F* h_dummy = new TH2F( "h_dummy", "Phi vs Eta", 200, -5., 5.,  200, -3.1415926, 3.1415926 ) ;
   h_dummy -> SetXTitle( "Eta" ) ;
   h_dummy -> SetYTitle( "Phi" ) ;

   TH2F* h_dummy_rhophi = new TH2F( "h_dummy_rhophi", "", 200, -1*rhophi_scale, rhophi_scale,  200, -1*rhophi_scale, rhophi_scale ) ;

   char label[100] ;
   TText* text = new TText() ;
   text -> SetTextSize(0.03) ;
   text -> SetTextFont( 82 ) ;

   TMarker* marker_gp = new TMarker() ;
   marker_gp -> SetMarkerStyle( 22 ) ;
   marker_gp -> SetMarkerSize( 2.5 ) ;
   marker_gp -> SetMarkerColor( 4 ) ;

   TMarker* marker_rj = new TMarker() ;
   marker_rj -> SetMarkerStyle(24) ;
   marker_rj -> SetMarkerSize( 2.5 ) ;
   marker_rj -> SetMarkerColor(2) ;

   TMarker* marker_fj = new TMarker() ;
   marker_fj -> SetMarkerStyle(24) ;
   marker_fj -> SetMarkerSize( 4.5 ) ;
   marker_fj -> SetMarkerColor(3) ;

   TEllipse* circle_rj = new TEllipse() ;
   circle_rj -> SetLineColor(2) ;
   circle_rj -> SetLineWidth(2) ;
   circle_rj -> SetFillColor(0) ;
   circle_rj -> SetFillStyle(0) ;
   float circle_rj_radius = 0.4 ;

   TEllipse* circle_fj = new TEllipse() ;
   circle_fj -> SetLineColor(3) ;
   circle_fj -> SetLineWidth(2) ;
   circle_fj -> SetFillColor(0) ;
   circle_fj -> SetFillStyle(0) ;
   float circle_fj_radius = 0.8 ;

   TLine* jet_pt_bar = new TLine() ;
   jet_pt_bar -> SetLineWidth(3) ;
   jet_pt_bar -> SetLineColor(2) ;
   double jet_pt_bar_scale = 0.2/20. ; // deta/pt
   double jet_pt_bar_deta = 0.1 ;
   double jet_pt_bar_dphi = -0.20 ;

   TLine* fatjet_pt_bar = new TLine() ;
   fatjet_pt_bar -> SetLineWidth(3) ;
   fatjet_pt_bar -> SetLineColor(3) ;
   double fatjet_pt_bar_scale = 0.2/20. ; // deta/pt
   double fatjet_pt_bar_deta = 0.1 ;
   double fatjet_pt_bar_dphi = -0.15 ;

   TArrow* arrow = new TArrow() ;
   arrow -> SetLineWidth(3) ;
   arrow -> SetFillStyle(0) ;

   Long64_t nbytes = 0, nb = 0;

   int first_event ;
   int last_event ;
   if ( display_event > 0 ) {
      first_event = display_event ;
      last_event = display_event ;
   } else {
      first_event = 1 ;
      last_event = nentries-1 ;
      printf("\n\n Running in loop mode.  %d to %d\n\n", first_event, last_event) ;
   }

   for ( int jentry=first_event; jentry<=last_event; jentry++ ) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      can1 -> cd() ;
      can1 -> Clear() ;
      h_dummy -> Draw() ;
      gPad -> SetGridx(1) ;
      gPad -> SetGridy(1) ;
      can1 -> Update() ; can1 -> Draw() ;

      can2 -> cd() ;
      can2 -> Clear() ;
      h_dummy_rhophi -> Draw() ;
      can2 -> Update() ; can2 -> Draw() ;

      can3 -> cd() ;
      can3 -> Clear() ;
      h_dummy_rhophi -> Draw() ;
      can3 -> Update() ; can3 -> Draw() ;



      printf("\n\n =========== Run %9u , Lumi %9u , Event %9llu\n", RunNum, LumiBlockNum, EvtNum ) ;

      printf("\n") ;
      printf("  GenParticles: %lu\n", GenParticles->size() ) ;
      bool skip_event(false) ;
      int top1_gpi(-1) ;
      int top2_gpi(-1) ;
      int w1_gpi(-1) ;
      int w2_gpi(-1) ;
      int w1_d1_gpi(-1) ;
      int w1_d2_gpi(-1) ;
      int w2_d1_gpi(-1) ;
      int w2_d2_gpi(-1) ;
      int stop1_gpi(-1) ;
      int stop2_gpi(-1) ;
      int chi1_gpi(-1) ;
      int chi2_gpi(-1) ;
      for ( unsigned int gpi=0; gpi < GenParticles->size() ; gpi++ ) {
         int spdgid = GenParticles_PdgId->at(gpi) ;
         int pdgid = abs( GenParticles_PdgId->at(gpi) ) ;
         int smomid = GenParticles_ParentId->at(gpi) ;
         int momid = abs( GenParticles_ParentId->at(gpi) ) ;
         int momidx = GenParticles_ParentIdx->at(gpi) ;

         TLorentzVector gplv( GenParticles->at(gpi) ) ;

         if ( spdgid ==  1000006 && stop1_gpi < 0 ) stop1_gpi = gpi ;
         if ( spdgid == -1000006 && stop2_gpi < 0 ) stop2_gpi = gpi ;

         if ( spdgid == 1000022 && momidx >= 0 && stop1_gpi >= 0 && momidx == stop1_gpi && chi1_gpi < 0 ) chi1_gpi = gpi ;
         if ( spdgid == 1000022 && momidx >= 0 && stop2_gpi >= 0 && momidx == stop2_gpi && chi2_gpi < 0 ) chi2_gpi = gpi ;

         if ( spdgid ==  6 && momidx >= 0 && stop1_gpi >= 0 && top1_gpi < 0 ) top1_gpi = gpi ;
         if ( spdgid == -6 && momidx >= 0 && stop2_gpi >= 0 && top2_gpi < 0 ) top2_gpi = gpi ;

         if ( spdgid ==  24 && smomid == 6 && w1_gpi < 0 ) { w1_gpi = gpi ; }
         if ( spdgid == -24 && smomid ==-6 && w2_gpi < 0 ) { w2_gpi = gpi ; }

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
         arrow -> SetLineColor(1) ;
         arrow -> SetLineWidth(1) ;
         arrow -> SetLineStyle(2) ;
         if ( ! (GenParticles_Status->at(gpi) == 1 && GenParticles_ParentIdx->at(gpi) == -1)
              && (GenParticles_Status->at(gpi) < 30 || GenParticles_Status->at(gpi) == 62 ) ) {
            if ( pdgid==1 || pdgid==2 || pdgid==3 || pdgid==4 || pdgid==5 || pdgid==21 ) {
               marker_gp -> SetMarkerColor(4) ;
            } else {
               marker_gp -> SetMarkerColor(6) ;
            }
            if ( momidx >= 0 && momidx == chi1_gpi ) {
               arrow -> SetLineColor(2) ;
               arrow -> SetLineWidth(3) ;
               arrow -> SetLineStyle(1) ;
               if ( pt > 20 ) {
                  marker_gp -> SetMarkerStyle( 20 ) ;
               } else {
                  marker_gp -> SetMarkerStyle( 24 ) ;
               }
            } else if ( momidx >= 0 && momidx == chi2_gpi ) {
               arrow -> SetLineColor(4) ;
               arrow -> SetLineWidth(3) ;
               arrow -> SetLineStyle(1) ;
               if ( pt > 20 ) {
                  marker_gp -> SetMarkerStyle( 21 ) ;
               } else {
                  marker_gp -> SetMarkerStyle( 25 ) ;
               }
            } else if ( (spdgid == 5 ) || gpi == w1_d1_gpi || gpi == w1_d2_gpi  ) {
               arrow -> SetLineColor(3) ;
               arrow -> SetLineWidth(3) ;
               arrow -> SetLineStyle(1) ;
               if ( pt > 20 ) {
                  marker_gp -> SetMarkerStyle( 22 ) ;
               } else {
                  marker_gp -> SetMarkerStyle( 26 ) ;
               }
            } else if ( (spdgid ==-5 ) || gpi == w2_d1_gpi || gpi == w2_d2_gpi ) {
               arrow -> SetLineColor(6) ;
               arrow -> SetLineWidth(3) ;
               arrow -> SetLineStyle(1) ;
               if ( pt > 20 ) {
                  marker_gp -> SetMarkerStyle( 23 ) ;
               } else {
                  marker_gp -> SetMarkerStyle( 32 ) ;
               }
            } else {
               marker_gp -> SetMarkerStyle( 27 ) ;
            }
            if ( pdgid == 1000006 ) {
               marker_gp -> SetMarkerStyle( 27 ) ;
            }
            if ( pdgid==6 || pdgid == 1000006 ) {
               text -> SetTextSize(0.035) ;
               text -> SetTextFont( 102 ) ;
            } else {
               text -> SetTextSize(0.025) ;
               text -> SetTextFont( 82 ) ;
            }
            can1->cd() ;
            marker_gp -> DrawMarker( eta, phi ) ;
            text -> DrawText( eta+0.1, phi+0.1, pname ) ;
            can2->cd() ;
            arrow -> DrawArrow( 0., 0., gplv.Px(), gplv.Py(), 0.02 ) ;
            text -> DrawText( gplv.Px(), gplv.Py(), pname ) ;
         }
      } // gpi
      printf( "\n" ) ;
      printf( " stop1 = %3d ,  top1 = %3d ,  w1 = %3d , chi1 = %d\n", stop1_gpi, top1_gpi, w1_gpi, chi1_gpi ) ;
      printf( " stop2 = %3d ,  top2 = %3d ,  w2 = %3d , chi2 = %d\n", stop2_gpi, top2_gpi, w2_gpi, chi2_gpi ) ;
      if ( skip_event ) {
         printf("\n\n *** skipping SL W decay event.\n\n") ;
         continue ;
      }

  //*********************************************************************************************************
  //*********************************************************************************************************
  //       Only look at events with two tops with pt>200
  //*********************************************************************************************************
  //*********************************************************************************************************
  ///////////////    if ( GenParticles->at(top1_gpi).Pt()<200 || GenParticles->at(top2_gpi).Pt()<200 ) continue ;
  //*********************************************************************************************************
  //*********************************************************************************************************
  //*********************************************************************************************************

      can1 -> cd() ;

      printf("\n") ;
      printf("  GenJets: %lu\n", GenJets->size() ) ;
      for ( unsigned int gji=0; gji < GenJets->size() ; gji++ ) {

         int best_match_gpi(-1) ;
         double best_match_dr(99.) ;

         double gj_eta = GenJets->at(gji).Eta() ;
         double gj_phi = GenJets->at(gji).Phi() ;

         for ( unsigned int gpi=0; gpi < GenParticles->size() ; gpi++ ) {
            int pdgid = abs( GenParticles_PdgId->at(gpi) ) ;
            if ( !( pdgid == 1 || pdgid == 2 || pdgid == 3 || pdgid == 4 || pdgid == 5
                 || pdgid == 11 || pdgid == 13 || pdgid == 21 || pdgid == 22 ) ) continue ;
            double gp_eta = 99. ;
            if ( GenParticles->at(gpi).Pt() > 0 ) gp_eta = GenParticles->at(gpi).Eta() ;
            double ptr = (GenParticles->at(gpi).Pt())/(GenJets->at(gji).Pt()) ;
            double dr = calcDr( gj_eta, gp_eta, gj_phi, GenParticles->at(gpi).Phi() ) ;
            //printf("  gji=%3d gpi=%3d : dr=%6.3f\n", gji, gpi, dr ) ;
            if ( dr < best_match_dr && ptr > 0.5 ) {
               best_match_dr = dr ;
               best_match_gpi = gpi ;
            }
         } // gpi

         printf("  %3d :  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f",
           gji,
             GenJets->at(gji).Pt(),
             GenJets->at(gji).Eta(),
             GenJets->at(gji).Phi()
         ) ;
         if ( best_match_dr < 0.1 && best_match_gpi >= 0 ) {
            char pname[100] ;
            sprintf( pname, "%s", mcname( GenParticles_PdgId->at(best_match_gpi) ) ) ;
            double ptr = (GenParticles->at(best_match_gpi).Pt())/(GenJets->at(gji).Pt()) ;
            printf(" : %3d dr=%5.3f ptr=%5.3f %10s", best_match_gpi, best_match_dr, ptr, pname ) ;
         }
         printf("\n") ;
      } // gji

      printf("\n") ;
      printf("  Jets: %lu\n", Jets->size() ) ;
      for ( unsigned int ji=0; ji < Jets->size() ; ji++ ) {
         printf("  %3d :  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f",
           ji,
             Jets->at(ji).Pt(),
             Jets->at(ji).Eta(),
             Jets->at(ji).Phi()
         ) ;
         printf("\n") ;
         ///marker_rj -> DrawMarker( Jets->at(ji).Eta(), Jets->at(ji).Phi() ) ;
         if ( Jets->at(ji).Pt() > 30 ) {
            circle_rj -> SetLineStyle(1) ;
         } else {
            circle_rj -> SetLineStyle(3) ;
         }
         double eta = Jets->at(ji).Eta() ;
         double phi = Jets->at(ji).Phi() ;
         double pt = Jets->at(ji).Pt() ;

         can1->cd() ;
         circle_rj -> DrawEllipse( eta, phi, circle_rj_radius, circle_rj_radius, 0., 360., 0. ) ;
         jet_pt_bar -> DrawLine( eta+jet_pt_bar_deta, phi+jet_pt_bar_dphi,
                         eta+jet_pt_bar_deta + jet_pt_bar_scale*pt, phi+jet_pt_bar_dphi ) ;

         TLorentzVector rjlv( Jets->at(ji) ) ;
         can3->cd() ;
         arrow->SetLineColor(4) ;
         arrow->SetLineStyle(1) ;
         arrow->SetLineWidth(3) ;
         arrow -> DrawArrow( 0., 0., rjlv.Px(), rjlv.Py(), 0.02 ) ;

      } // ji

      printf("\n") ;
      printf("  Jets AK8: %lu\n", JetsAK8->size() ) ;
      for ( unsigned int ji=0; ji < JetsAK8->size() ; ji++ ) {
         double tau21 = -1;
         if ( JetsAK8_NsubjettinessTau2->at(ji) > 0 && JetsAK8_NsubjettinessTau1->at(ji) > 0 ) {
            tau21 = (JetsAK8_NsubjettinessTau2->at(ji))/(JetsAK8_NsubjettinessTau1->at(ji)) ;
         }
         double tau31 = -1;
         if ( JetsAK8_NsubjettinessTau3->at(ji) > 0 && JetsAK8_NsubjettinessTau1->at(ji) > 0 ) {
            tau31 = (JetsAK8_NsubjettinessTau3->at(ji))/(JetsAK8_NsubjettinessTau1->at(ji)) ;
         }
         printf("  %3d :  Pt = %7.1f , Eta = %6.3f, Phi = %6.3f , Mass = %6.1f, tau21 = %5.3f, tau31 = %5.3f",
           ji,
             JetsAK8->at(ji).Pt(),
             JetsAK8->at(ji).Eta(),
             JetsAK8->at(ji).Phi(),
             JetsAK8_prunedMass->at(ji),
             tau21,
             tau31
         ) ;
         printf("\n") ;
         double eta = JetsAK8->at(ji).Eta() ;
         double phi = JetsAK8->at(ji).Phi() ;
         double pt = JetsAK8->at(ji).Pt() ;
         can1->cd() ;
         circle_fj -> DrawEllipse( eta, phi, circle_fj_radius, circle_fj_radius, 0., 360., 0. ) ;
         fatjet_pt_bar -> DrawLine( eta+fatjet_pt_bar_deta, phi+fatjet_pt_bar_dphi,
                         eta+fatjet_pt_bar_deta + fatjet_pt_bar_scale*pt, phi+fatjet_pt_bar_dphi ) ;

         TLorentzVector rjlv( JetsAK8->at(ji) ) ;
         can3->cd() ;
         arrow->SetLineColor(2) ;
         arrow->SetLineStyle(3) ;
         arrow->SetLineWidth(5) ;
         arrow -> DrawArrow( 0., 0., rjlv.Px(), rjlv.Py(), 0.02 ) ;

      } // ji

      can1 -> Update() ; can1 -> Draw() ;
      can2 -> Update() ; can2 -> Draw() ;
      can3 -> Update() ; can3 -> Draw() ;

      char answer[100] ;
      printf("\n\n Type n for next, q to quit: ") ;
      scanf( "%s", answer ) ;
      if ( strcmp( answer, "q") == 0 ) return ;

   } // jentry

} // Loop

//=================================================================================


const char* mcname( int pdgid ) {

   sprintf( pname, "" ) ;

   if ( pdgid == 1 ) sprintf( pname, "d" ) ;
   if ( pdgid == 2 ) sprintf( pname, "u" ) ;
   if ( pdgid == 3 ) sprintf( pname, "s" ) ;
   if ( pdgid == 4 ) sprintf( pname, "c" ) ;
   if ( pdgid == 5 ) sprintf( pname, "b" ) ;
   if ( pdgid == 6 ) sprintf( pname, "t" ) ;

   if ( pdgid == -1 ) sprintf( pname, "d-bar" ) ;
   if ( pdgid == -2 ) sprintf( pname, "u-bar" ) ;
   if ( pdgid == -3 ) sprintf( pname, "s-bar" ) ;
   if ( pdgid == -4 ) sprintf( pname, "c-bar" ) ;
   if ( pdgid == -5 ) sprintf( pname, "b-bar" ) ;
   if ( pdgid == -6 ) sprintf( pname, "t-bar" ) ;

   if ( pdgid == 11 ) sprintf( pname, "e-" ) ;
   if ( pdgid == 12 ) sprintf( pname, "nu_e" ) ;
   if ( pdgid == 13 ) sprintf( pname, "mu-" ) ;
   if ( pdgid == 14 ) sprintf( pname, "nu_mu" ) ;
   if ( pdgid == 15 ) sprintf( pname, "tau-" ) ;
   if ( pdgid == 16 ) sprintf( pname, "nu_tau" ) ;

   if ( pdgid == -11 ) sprintf( pname, "e+" ) ;
   if ( pdgid == -12 ) sprintf( pname, "nu_e-bar" ) ;
   if ( pdgid == -13 ) sprintf( pname, "mu+" ) ;
   if ( pdgid == -14 ) sprintf( pname, "nu_mu-bar" ) ;
   if ( pdgid == -15 ) sprintf( pname, "tau+" ) ;
   if ( pdgid == -16 ) sprintf( pname, "nu_tau-bar" ) ;

   if ( pdgid == 21 ) sprintf( pname, "gluon" ) ;
   if ( pdgid == 22 ) sprintf( pname, "photon" ) ;
   if ( pdgid == 23 ) sprintf( pname, "Z0" ) ;
   if ( pdgid == 24 ) sprintf( pname, "W+" ) ;
   if ( pdgid ==-24 ) sprintf( pname, "W-" ) ;
   if ( pdgid == 25 ) sprintf( pname, "h" ) ;
   if ( pdgid == 35 ) sprintf( pname, "H" ) ;
   if ( pdgid == 36 ) sprintf( pname, "a" ) ;


   if ( pdgid == 1000001 ) sprintf( pname, "~dL" ) ;
   if ( pdgid == 1000002 ) sprintf( pname, "~uL" ) ;
   if ( pdgid == 1000003 ) sprintf( pname, "~sL" ) ;
   if ( pdgid == 1000004 ) sprintf( pname, "~cL" ) ;
   if ( pdgid == 1000005 ) sprintf( pname, "~b1" ) ;
   if ( pdgid == 1000006 ) sprintf( pname, "~t1" ) ;

   if ( pdgid == -1000001 ) sprintf( pname, "~dL*" ) ;
   if ( pdgid == -1000002 ) sprintf( pname, "~uL*" ) ;
   if ( pdgid == -1000003 ) sprintf( pname, "~sL*" ) ;
   if ( pdgid == -1000004 ) sprintf( pname, "~cL*" ) ;
   if ( pdgid == -1000005 ) sprintf( pname, "~b1*" ) ;
   if ( pdgid == -1000006 ) sprintf( pname, "~t1*" ) ;

   if ( pdgid == 1000022 ) sprintf( pname, "~chi01" ) ;

   return pname ;


} // mcname

//=====================================================================

double calcDr( double eta1, double eta2, double phi1, double phi2 ) {

   double deta = fabs( eta1 - eta2 ) ;

   double dphi = phi1 - phi2 ;
   if ( dphi > 3.1415926 ) dphi -= 2*3.1415926 ;
   if ( dphi <-3.1415926 ) dphi += 2*3.1415926 ;

   return sqrt( dphi*dphi + deta*deta ) ;

} // calcDr





