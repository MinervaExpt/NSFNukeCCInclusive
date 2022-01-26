#ifndef MasterAnaDevvars_h
#define MasterAnaDevvars_h
//#ifndef NukeCCvars_h
//#define NukeCCvars_h
TTree          *fChain;   //!pointer to the analyzed TTree or TChain
Int_t           fCurrent; //!current Tree number in a TChain

static const size_t MAX_USACT_EXTENT = 250;
static const size_t MAX_MC_FR_NNUANCESTORIDS = 4*10;
static const size_t MAX_MC_ER_NPART = 50;
static const size_t MAX_MC_NFSPART = 25;
static const size_t MAX_N_PRONGS = 15;
static const size_t MAX_N_SLICES = 3;
static const size_t MAX_GENIE_WGT_N_SHIFTS = 7;
static const size_t MAX_MC_WGT_GENIE_SZ = 100;
static const size_t MAX_MC_WGT_FLUX_TERTIARY_SZ = 100;
static const size_t MAX_MC_WGT_FLUX_NA49_SZ = 100;
static const size_t MAX_MC_WGT_FLUX_BEAMFOCUS_SZ = 100;


// Declaration of leaf types
Double_t        eventID;
Int_t           physEvtNum;
Int_t           n_hyps;
Int_t           processType;
Int_t           primaryPart;
Int_t           n_slices;
Int_t           slice_numbers[MAX_N_SLICES];   //[n_slices]
Int_t           shared_slice;
Double_t        vtx[4];
Double_t        vtxErr[4];
Double_t        E[4];
Bool_t          found_truth;
Bool_t          phys_front_activity;
Bool_t          phys_energy_in_road_upstream_is_rockmuon_consistent;
Bool_t          rock_muons_removed;
Bool_t          minos_track_match;
Bool_t          minos_stub_match;
Bool_t          unknown_helicity;
Bool_t          minos_track_inside_partial_plane;
Bool_t          prim_vtx_has_misassigned_track_direction;
Bool_t          prim_vtx_has_broken_track;
Bool_t          pass_MasterAnaDev;
Bool_t          short_track_vtx_used;
Bool_t          muon_sp_moved;
Bool_t          vtx_fit_converged;
Bool_t          muon_is_correct;
Bool_t          has_int_vtx;
Bool_t          has_bad_object;
Bool_t          has_muon;
Bool_t          muon_has_charge;
Bool_t          has_good_vtx;
Bool_t          is_rock_muon;
Int_t           MuonTaggedAsVetoButNotMatched;
Int_t           NonVetoMuonExtrpToVeto;
Int_t           NonVetoMuonWallOnePaddle;
Int_t           NonVetoMuonWallOnePaddleOverlap;
Int_t           NonVetoMuonWallOneSector;
Int_t           NonVetoMuonWallTwoPaddle;
Int_t           NonVetoMuonWallTwoPaddleOverlap;
Int_t           NonVetoMuonWallTwoSector;
Int_t           VetoMuonWallOnePaddle;
Int_t           VetoMuonWallOnePaddleOverlap;
Int_t           VetoMuonWallOneSector;
Int_t           VetoMuonWallOneTypeOfMatchNonOverlap;
Int_t           VetoMuonWallOneTypeOfMatchOverlap;
Int_t           VetoMuonWallTwoPaddle;
Int_t           VetoMuonWallTwoPaddleOverlap;
Int_t           VetoMuonWallTwoSector;
Int_t           VetoMuonWallTwoTypeOfMatchNonOverlap;
Int_t           VetoMuonWallTwoTypeOfMatchOverlap;
Int_t           blob_disp_nBlobs;
Int_t           blob_disp_nClus;
Int_t           blob_disp_nClus_ecal;
Int_t           blob_disp_nClus_hcal;
Int_t           blob_disp_nClus_nucl;
Int_t           blob_disp_nClus_od;
Int_t           blob_disp_nClus_tracker;
Int_t           blob_iso_nBlobs;
Int_t           blob_iso_nClus;
Int_t           blob_iso_nClus_ecal;
Int_t           blob_iso_nClus_hcal;
Int_t           blob_iso_nClus_nucl;
Int_t           blob_iso_nClus_od;
Int_t           blob_iso_nClus_tracker;
Int_t           blob_mufuzz_nBlobs;
Int_t           blob_mufuzz_nClus;
Int_t           blob_mufuzz_nClus_ecal;
Int_t           blob_mufuzz_nClus_hcal;
Int_t           blob_mufuzz_nClus_nucl;
Int_t           blob_mufuzz_nClus_od;
Int_t           blob_mufuzz_nClus_tracker;
Int_t           blob_recoil_nBlobs;
Int_t           blob_recoil_nClus;
Int_t           blob_recoil_nClus_ecal;
Int_t           blob_recoil_nClus_hcal;
Int_t           blob_recoil_nClus_nucl;
Int_t           blob_recoil_nClus_od;
Int_t           blob_recoil_nClus_tracker;
Int_t           blob_vtx_nBlobs;
Int_t           blob_vtx_nClus;
Int_t           blob_vtx_nClus_ecal;
Int_t           blob_vtx_nClus_hcal;
Int_t           blob_vtx_nClus_nucl;
Int_t           blob_vtx_nClus_od;
Int_t           blob_vtx_nClus_tracker;
Int_t           broken_track_most_us_plane;
Int_t           muon_n_USclusters;
Int_t           muon_truthMatch_track_id;
Int_t           n_prim_long_tracks;
Int_t           n_prim_short_tracks;
Int_t           n_start_vertices;
Int_t           n_tracks;
Int_t           n_tracks_non_prim;
Int_t           n_tracks_prim;
Int_t           n_tracks_prim_forked;
Int_t           n_tracks_prim_kinked;
Int_t           n_vertices_startpoint;
Int_t           passVetoMuonCut;
Int_t           phys_energy_in_road_downstream_nplanes;
Int_t           phys_energy_in_road_upstream_nplanes;
Int_t           phys_n_dead_discr_pair;
Int_t           phys_n_dead_discr_pair_in_prim_track_region;
Int_t           phys_n_dead_discr_pair_two_mod_downstream_prim_track;
Int_t           phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;
Int_t           phys_n_dead_discr_pair_upstream_prim_track_proj;
Int_t           phys_vertex_is_fiducial;
Int_t           rock_muon_code;
Int_t           truth_has_michel_electron;
Int_t           usact_extent_huge;
Int_t           usact_extent_inf;
Int_t           usact_extent_large;
Int_t           usact_extent_normal;
Int_t           usact_extent_small;
Int_t           usact_extent_texas_sized;
Int_t           usact_extent_tiny;
Int_t           usact_n_consecutive_huge;
Int_t           usact_n_consecutive_inf;
Int_t           usact_n_consecutive_large;
Int_t           usact_n_consecutive_normal;
Int_t           usact_n_consecutive_small;
Int_t           usact_n_consecutive_texas_sized;
Int_t           usact_n_consecutive_tiny;
Int_t           usact_n_planes_huge;
Int_t           usact_n_planes_inf;
Int_t           usact_n_planes_large;
Int_t           usact_n_planes_normal;
Int_t           usact_n_planes_small;
Int_t           usact_n_planes_texas_sized;
Int_t           usact_n_planes_tiny;
Double_t        NonVetoMuonWallOneBadPosX;
Double_t        NonVetoMuonWallOneBadPosY;
Double_t        NonVetoMuonWallOnePosX;
Double_t        NonVetoMuonWallOnePosY;
Double_t        NonVetoMuonWallOne_ANDEfficiency_Central;
Double_t        NonVetoMuonWallOne_ANDEfficiency_Overlap;
Double_t        NonVetoMuonWallOne_ANDError_Central;
Double_t        NonVetoMuonWallOne_ANDError_Overlap;
Double_t        NonVetoMuonWallOne_AccRatesError_PaddleAbove;
Double_t        NonVetoMuonWallOne_AccRatesError_PaddleBelow;
Double_t        NonVetoMuonWallOne_AccRates_PaddleAbove;
Double_t        NonVetoMuonWallOne_AccRates_PaddleBelow;
Double_t        NonVetoMuonWallOne_OREfficiency_Central;
Double_t        NonVetoMuonWallOne_OREfficiency_Overlap;
Double_t        NonVetoMuonWallOne_ORError_Central;
Double_t        NonVetoMuonWallOne_ORError_Overlap;
Double_t        NonVetoMuonWallTwoBadPosX;
Double_t        NonVetoMuonWallTwoBadPosY;
Double_t        NonVetoMuonWallTwoPosX;
Double_t        NonVetoMuonWallTwoPosY;
Double_t        NonVetoMuonWallTwo_ANDEfficiency_Central;
Double_t        NonVetoMuonWallTwo_ANDEfficiency_Overlap;
Double_t        NonVetoMuonWallTwo_ANDError_Central;
Double_t        NonVetoMuonWallTwo_ANDError_Overlap;
Double_t        NonVetoMuonWallTwo_AccRatesError_PaddleAbove;
Double_t        NonVetoMuonWallTwo_AccRatesError_PaddleBelow;
Double_t        NonVetoMuonWallTwo_AccRates_PaddleAbove;
Double_t        NonVetoMuonWallTwo_AccRates_PaddleBelow;
Double_t        NonVetoMuonWallTwo_OREfficiency_Central;
Double_t        NonVetoMuonWallTwo_OREfficiency_Overlap;
Double_t        NonVetoMuonWallTwo_ORError_Central;
Double_t        NonVetoMuonWallTwo_ORError_Overlap;
Double_t        VetoMuonWallOneDeltaTime;
Double_t        VetoMuonWallOnePosX;
Double_t        VetoMuonWallOnePosY;
Double_t        VetoMuonWallOne_ANDEfficiency_Central;
Double_t        VetoMuonWallOne_ANDEfficiency_Overlap;
Double_t        VetoMuonWallOne_ANDError_Central;
Double_t        VetoMuonWallOne_ANDError_Overlap;
Double_t        VetoMuonWallOne_AccRatesError_PaddleAbove;
Double_t        VetoMuonWallOne_AccRatesError_PaddleBelow;
Double_t        VetoMuonWallOne_AccRates_PaddleAbove;
Double_t        VetoMuonWallOne_AccRates_PaddleBelow;
Double_t        VetoMuonWallOne_OREfficiency_Central;
Double_t        VetoMuonWallOne_OREfficiency_Overlap;
Double_t        VetoMuonWallOne_ORError_Central;
Double_t        VetoMuonWallOne_ORError_Overlap;
Double_t        VetoMuonWallTwoDeltaTime;
Double_t        VetoMuonWallTwoPosX;
Double_t        VetoMuonWallTwoPosY;
Double_t        VetoMuonWallTwo_ANDEfficiency_Central;
Double_t        VetoMuonWallTwo_ANDEfficiency_Overlap;
Double_t        VetoMuonWallTwo_ANDError_Central;
Double_t        VetoMuonWallTwo_ANDError_Overlap;
Double_t        VetoMuonWallTwo_AccRatesError_PaddleAbove;
Double_t        VetoMuonWallTwo_AccRatesError_PaddleBelow;
Double_t        VetoMuonWallTwo_AccRates_PaddleAbove;
Double_t        VetoMuonWallTwo_AccRates_PaddleBelow;
Double_t        VetoMuonWallTwo_OREfficiency_Central;
Double_t        VetoMuonWallTwo_OREfficiency_Overlap;
Double_t        VetoMuonWallTwo_ORError_Central;
Double_t        VetoMuonWallTwo_ORError_Overlap;
Double_t        blob_ccqe_recoil_E;
Double_t        blob_disp_E;
Double_t        blob_disp_E_ecal;
Double_t        blob_disp_E_ecal_em;
Double_t        blob_disp_E_ecal_highn;
Double_t        blob_disp_E_ecal_lown;
Double_t        blob_disp_E_ecal_meson;
Double_t        blob_disp_E_ecal_midn;
Double_t        blob_disp_E_ecal_mu;
Double_t        blob_disp_E_ecal_other;
Double_t        blob_disp_E_ecal_p;
Double_t        blob_disp_E_ecal_xtalk;
Double_t        blob_disp_E_em;
Double_t        blob_disp_E_hcal;
Double_t        blob_disp_E_hcal_em;
Double_t        blob_disp_E_hcal_highn;
Double_t        blob_disp_E_hcal_lown;
Double_t        blob_disp_E_hcal_meson;
Double_t        blob_disp_E_hcal_midn;
Double_t        blob_disp_E_hcal_mu;
Double_t        blob_disp_E_hcal_other;
Double_t        blob_disp_E_hcal_p;
Double_t        blob_disp_E_hcal_xtalk;
Double_t        blob_disp_E_highn;
Double_t        blob_disp_E_lown;
Double_t        blob_disp_E_meson;
Double_t        blob_disp_E_midn;
Double_t        blob_disp_E_mu;
Double_t        blob_disp_E_nucl;
Double_t        blob_disp_E_nucl_em;
Double_t        blob_disp_E_nucl_highn;
Double_t        blob_disp_E_nucl_lown;
Double_t        blob_disp_E_nucl_meson;
Double_t        blob_disp_E_nucl_midn;
Double_t        blob_disp_E_nucl_mu;
Double_t        blob_disp_E_nucl_other;
Double_t        blob_disp_E_nucl_p;
Double_t        blob_disp_E_nucl_xtalk;
Double_t        blob_disp_E_od;
Double_t        blob_disp_E_od_em;
Double_t        blob_disp_E_od_highn;
Double_t        blob_disp_E_od_lown;
Double_t        blob_disp_E_od_meson;
Double_t        blob_disp_E_od_midn;
Double_t        blob_disp_E_od_mu;
Double_t        blob_disp_E_od_other;
Double_t        blob_disp_E_od_p;
Double_t        blob_disp_E_od_xtalk;
Double_t        blob_disp_E_other;
Double_t        blob_disp_E_p;
Double_t        blob_disp_E_tracker;
Double_t        blob_disp_E_tracker_em;
Double_t        blob_disp_E_tracker_highn;
Double_t        blob_disp_E_tracker_lown;
Double_t        blob_disp_E_tracker_meson;
Double_t        blob_disp_E_tracker_midn;
Double_t        blob_disp_E_tracker_mu;
Double_t        blob_disp_E_tracker_other;
Double_t        blob_disp_E_tracker_p;
Double_t        blob_disp_E_tracker_xtalk;
Double_t        blob_disp_E_xtalk;
Double_t        blob_iso_E;
Double_t        blob_iso_E_ecal;
Double_t        blob_iso_E_ecal_em;
Double_t        blob_iso_E_ecal_highn;
Double_t        blob_iso_E_ecal_lown;
Double_t        blob_iso_E_ecal_meson;
Double_t        blob_iso_E_ecal_midn;
Double_t        blob_iso_E_ecal_mu;
Double_t        blob_iso_E_ecal_other;
Double_t        blob_iso_E_ecal_p;
Double_t        blob_iso_E_ecal_xtalk;
Double_t        blob_iso_E_em;
Double_t        blob_iso_E_hcal;
Double_t        blob_iso_E_hcal_em;
Double_t        blob_iso_E_hcal_highn;
Double_t        blob_iso_E_hcal_lown;
Double_t        blob_iso_E_hcal_meson;
Double_t        blob_iso_E_hcal_midn;
Double_t        blob_iso_E_hcal_mu;
Double_t        blob_iso_E_hcal_other;
Double_t        blob_iso_E_hcal_p;
Double_t        blob_iso_E_hcal_xtalk;
Double_t        blob_iso_E_highn;
Double_t        blob_iso_E_lown;
Double_t        blob_iso_E_meson;
Double_t        blob_iso_E_midn;
Double_t        blob_iso_E_mu;
Double_t        blob_iso_E_nucl;
Double_t        blob_iso_E_nucl_em;
Double_t        blob_iso_E_nucl_highn;
Double_t        blob_iso_E_nucl_lown;
Double_t        blob_iso_E_nucl_meson;
Double_t        blob_iso_E_nucl_midn;
Double_t        blob_iso_E_nucl_mu;
Double_t        blob_iso_E_nucl_other;
Double_t        blob_iso_E_nucl_p;
Double_t        blob_iso_E_nucl_xtalk;
Double_t        blob_iso_E_od;
Double_t        blob_iso_E_od_em;
Double_t        blob_iso_E_od_highn;
Double_t        blob_iso_E_od_lown;
Double_t        blob_iso_E_od_meson;
Double_t        blob_iso_E_od_midn;
Double_t        blob_iso_E_od_mu;
Double_t        blob_iso_E_od_other;
Double_t        blob_iso_E_od_p;
Double_t        blob_iso_E_od_xtalk;
Double_t        blob_iso_E_other;
Double_t        blob_iso_E_p;
Double_t        blob_iso_E_tracker;
Double_t        blob_iso_E_tracker_em;
Double_t        blob_iso_E_tracker_highn;
Double_t        blob_iso_E_tracker_lown;
Double_t        blob_iso_E_tracker_meson;
Double_t        blob_iso_E_tracker_midn;
Double_t        blob_iso_E_tracker_mu;
Double_t        blob_iso_E_tracker_other;
Double_t        blob_iso_E_tracker_p;
Double_t        blob_iso_E_tracker_xtalk;
Double_t        blob_iso_E_xtalk;
Double_t        blob_mufuzz_E;
Double_t        blob_mufuzz_E_ecal;
Double_t        blob_mufuzz_E_ecal_em;
Double_t        blob_mufuzz_E_ecal_highn;
Double_t        blob_mufuzz_E_ecal_lown;
Double_t        blob_mufuzz_E_ecal_meson;
Double_t        blob_mufuzz_E_ecal_midn;
Double_t        blob_mufuzz_E_ecal_mu;
Double_t        blob_mufuzz_E_ecal_other;
Double_t        blob_mufuzz_E_ecal_p;
Double_t        blob_mufuzz_E_ecal_xtalk;
Double_t        blob_mufuzz_E_em;
Double_t        blob_mufuzz_E_hcal;
Double_t        blob_mufuzz_E_hcal_em;
Double_t        blob_mufuzz_E_hcal_highn;
Double_t        blob_mufuzz_E_hcal_lown;
Double_t        blob_mufuzz_E_hcal_meson;
Double_t        blob_mufuzz_E_hcal_midn;
Double_t        blob_mufuzz_E_hcal_mu;
Double_t        blob_mufuzz_E_hcal_other;
Double_t        blob_mufuzz_E_hcal_p;
Double_t        blob_mufuzz_E_hcal_xtalk;
Double_t        blob_mufuzz_E_highn;
Double_t        blob_mufuzz_E_lown;
Double_t        blob_mufuzz_E_meson;
Double_t        blob_mufuzz_E_midn;
Double_t        blob_mufuzz_E_mu;
Double_t        blob_mufuzz_E_nucl;
Double_t        blob_mufuzz_E_nucl_em;
Double_t        blob_mufuzz_E_nucl_highn;
Double_t        blob_mufuzz_E_nucl_lown;
Double_t        blob_mufuzz_E_nucl_meson;
Double_t        blob_mufuzz_E_nucl_midn;
Double_t        blob_mufuzz_E_nucl_mu;
Double_t        blob_mufuzz_E_nucl_other;
Double_t        blob_mufuzz_E_nucl_p;
Double_t        blob_mufuzz_E_nucl_xtalk;
Double_t        blob_mufuzz_E_od;
Double_t        blob_mufuzz_E_od_em;
Double_t        blob_mufuzz_E_od_highn;
Double_t        blob_mufuzz_E_od_lown;
Double_t        blob_mufuzz_E_od_meson;
Double_t        blob_mufuzz_E_od_midn;
Double_t        blob_mufuzz_E_od_mu;
Double_t        blob_mufuzz_E_od_other;
Double_t        blob_mufuzz_E_od_p;
Double_t        blob_mufuzz_E_od_xtalk;
Double_t        blob_mufuzz_E_other;
Double_t        blob_mufuzz_E_p;
Double_t        blob_mufuzz_E_tracker;
Double_t        blob_mufuzz_E_tracker_em;
Double_t        blob_mufuzz_E_tracker_highn;
Double_t        blob_mufuzz_E_tracker_lown;
Double_t        blob_mufuzz_E_tracker_meson;
Double_t        blob_mufuzz_E_tracker_midn;
Double_t        blob_mufuzz_E_tracker_mu;
Double_t        blob_mufuzz_E_tracker_other;
Double_t        blob_mufuzz_E_tracker_p;
Double_t        blob_mufuzz_E_tracker_xtalk;
Double_t        blob_mufuzz_E_xtalk;
Double_t        blob_recoil_E;
Double_t        blob_recoil_E_ecal;
Double_t        blob_recoil_E_ecal_em;
Double_t        blob_recoil_E_ecal_highn;
Double_t        blob_recoil_E_ecal_lown;
Double_t        blob_recoil_E_ecal_meson;
Double_t        blob_recoil_E_ecal_midn;
Double_t        blob_recoil_E_ecal_mu;
Double_t        blob_recoil_E_ecal_other;
Double_t        blob_recoil_E_ecal_p;
Double_t        blob_recoil_E_ecal_xtalk;
Double_t        blob_recoil_E_em;
Double_t        blob_recoil_E_hcal;
Double_t        blob_recoil_E_hcal_em;
Double_t        blob_recoil_E_hcal_highn;
Double_t        blob_recoil_E_hcal_lown;
Double_t        blob_recoil_E_hcal_meson;
Double_t        blob_recoil_E_hcal_midn;
Double_t        blob_recoil_E_hcal_mu;
Double_t        blob_recoil_E_hcal_other;
Double_t        blob_recoil_E_hcal_p;
Double_t        blob_recoil_E_hcal_xtalk;
Double_t        blob_recoil_E_highn;
Double_t        blob_recoil_E_lown;
Double_t        blob_recoil_E_meson;
Double_t        blob_recoil_E_midn;
Double_t        blob_recoil_E_mu;
Double_t        blob_recoil_E_nucl;
Double_t        blob_recoil_E_nucl_em;
Double_t        blob_recoil_E_nucl_highn;
Double_t        blob_recoil_E_nucl_lown;
Double_t        blob_recoil_E_nucl_meson;
Double_t        blob_recoil_E_nucl_midn;
Double_t        blob_recoil_E_nucl_mu;
Double_t        blob_recoil_E_nucl_other;
Double_t        blob_recoil_E_nucl_p;
Double_t        blob_recoil_E_nucl_xtalk;
Double_t        blob_recoil_E_od;
Double_t        blob_recoil_E_od_em;
Double_t        blob_recoil_E_od_highn;
Double_t        blob_recoil_E_od_lown;
Double_t        blob_recoil_E_od_meson;
Double_t        blob_recoil_E_od_midn;
Double_t        blob_recoil_E_od_mu;
Double_t        blob_recoil_E_od_other;
Double_t        blob_recoil_E_od_p;
Double_t        blob_recoil_E_od_xtalk;
Double_t        blob_recoil_E_other;
Double_t        blob_recoil_E_p;
Double_t        blob_recoil_E_tracker;
Double_t        blob_recoil_E_tracker_em;
Double_t        blob_recoil_E_tracker_highn;
Double_t        blob_recoil_E_tracker_lown;
Double_t        blob_recoil_E_tracker_meson;
Double_t        blob_recoil_E_tracker_midn;
Double_t        blob_recoil_E_tracker_mu;
Double_t        blob_recoil_E_tracker_other;
Double_t        blob_recoil_E_tracker_p;
Double_t        blob_recoil_E_tracker_xtalk;
Double_t        blob_recoil_E_xtalk;
Double_t        blob_vtx_E;
Double_t        blob_vtx_E_ecal;
Double_t        blob_vtx_E_ecal_em;
Double_t        blob_vtx_E_ecal_highn;
Double_t        blob_vtx_E_ecal_lown;
Double_t        blob_vtx_E_ecal_meson;
Double_t        blob_vtx_E_ecal_midn;
Double_t        blob_vtx_E_ecal_mu;
Double_t        blob_vtx_E_ecal_other;
Double_t        blob_vtx_E_ecal_p;
Double_t        blob_vtx_E_ecal_xtalk;
Double_t        blob_vtx_E_em;
Double_t        blob_vtx_E_hcal;
Double_t        blob_vtx_E_hcal_em;
Double_t        blob_vtx_E_hcal_highn;
Double_t        blob_vtx_E_hcal_lown;
Double_t        blob_vtx_E_hcal_meson;
Double_t        blob_vtx_E_hcal_midn;
Double_t        blob_vtx_E_hcal_mu;
Double_t        blob_vtx_E_hcal_other;
Double_t        blob_vtx_E_hcal_p;
Double_t        blob_vtx_E_hcal_xtalk;
Double_t        blob_vtx_E_highn;
Double_t        blob_vtx_E_lown;
Double_t        blob_vtx_E_meson;
Double_t        blob_vtx_E_midn;
Double_t        blob_vtx_E_mu;
Double_t        blob_vtx_E_nucl;
Double_t        blob_vtx_E_nucl_em;
Double_t        blob_vtx_E_nucl_highn;
Double_t        blob_vtx_E_nucl_lown;
Double_t        blob_vtx_E_nucl_meson;
Double_t        blob_vtx_E_nucl_midn;
Double_t        blob_vtx_E_nucl_mu;
Double_t        blob_vtx_E_nucl_other;
Double_t        blob_vtx_E_nucl_p;
Double_t        blob_vtx_E_nucl_xtalk;
Double_t        blob_vtx_E_od;
Double_t        blob_vtx_E_od_em;
Double_t        blob_vtx_E_od_highn;
Double_t        blob_vtx_E_od_lown;
Double_t        blob_vtx_E_od_meson;
Double_t        blob_vtx_E_od_midn;
Double_t        blob_vtx_E_od_mu;
Double_t        blob_vtx_E_od_other;
Double_t        blob_vtx_E_od_p;
Double_t        blob_vtx_E_od_xtalk;
Double_t        blob_vtx_E_other;
Double_t        blob_vtx_E_p;
Double_t        blob_vtx_E_tracker;
Double_t        blob_vtx_E_tracker_em;
Double_t        blob_vtx_E_tracker_highn;
Double_t        blob_vtx_E_tracker_lown;
Double_t        blob_vtx_E_tracker_meson;
Double_t        blob_vtx_E_tracker_midn;
Double_t        blob_vtx_E_tracker_mu;
Double_t        blob_vtx_E_tracker_other;
Double_t        blob_vtx_E_tracker_p;
Double_t        blob_vtx_E_tracker_xtalk;
Double_t        blob_vtx_E_xtalk;
Double_t        energy_from_mc;
Double_t        energy_from_mc_fraction;
Double_t        energy_from_mc_fraction_of_highest;
Double_t        muon_phi;
Double_t        muon_theta;
Double_t        muon_thetaX;
Double_t        muon_thetaY;
Double_t        muon_truthMatch_fraction;
Double_t        numi_horn_curr;
Double_t        numi_pot;
Double_t        numi_x;
Double_t        numi_x_width;
Double_t        numi_y;
Double_t        numi_y_width;
Double_t        phys_energy_dispersed;
Double_t        phys_energy_in_road_downstream;
Double_t        phys_energy_in_road_upstream;
Double_t        phys_energy_unattached;
Double_t        prim_vtx_smallest_opening_angle;
Double_t        primary_track_minerva_energy;
Double_t        primary_track_minerva_phi;
Double_t        primary_track_minerva_theta;
Double_t        usact_avg_E_consecutive_huge;
Double_t        usact_avg_E_consecutive_inf;
Double_t        usact_avg_E_consecutive_large;
Double_t        usact_avg_E_consecutive_normal;
Double_t        usact_avg_E_consecutive_small;
Double_t        usact_avg_E_consecutive_texas_sized;
Double_t        usact_avg_E_consecutive_tiny;
Double_t        usact_avg_E_huge;
Double_t        usact_avg_E_inf;
Double_t        usact_avg_E_large;
Double_t        usact_avg_E_normal;
Double_t        usact_avg_E_small;
Double_t        usact_avg_E_texas_sized;
Double_t        usact_avg_E_tiny;
Double_t        usact_frac_withE_huge;
Double_t        usact_frac_withE_inf;
Double_t        usact_frac_withE_large;
Double_t        usact_frac_withE_normal;
Double_t        usact_frac_withE_small;
Double_t        usact_frac_withE_texas_sized;
Double_t        usact_frac_withE_tiny;
Double_t        vetoMuonTime;
Double_t        vtx_fit_chi2;
Int_t           n_vetoDigits;
Int_t           discrFired[1];   //[n_vetoDigits]
Int_t           has_michel_at_vertex_sz;
Int_t           has_michel_at_vertex[2];   //[has_michel_at_vertex_sz]
Int_t           has_michel_beginModule_sz;
Int_t           has_michel_beginModule[2];   //[has_michel_beginModule_sz]
Int_t           has_michel_category_sz;
Int_t           has_michel_category[2];   //[has_michel_category_sz]
Int_t           has_michel_endModule_sz;
Int_t           has_michel_endModule[2];   //[has_michel_endModule_sz]
Int_t           has_michel_numDigits_sz;
Int_t           has_michel_numDigits[2];   //[has_michel_numDigits_sz]
Int_t           has_michel_numModules_sz;
Int_t           has_michel_numModules[2];   //[has_michel_numModules_sz]
Int_t           has_michel_numPlanes_sz;
Int_t           has_michel_numPlanes[2];   //[has_michel_numPlanes_sz]
Int_t           has_michel_numTracks_sz;
Int_t           has_michel_numTracks[2];   //[has_michel_numTracks_sz]
Int_t           has_michel_slice_number_sz;
Int_t           has_michel_slice_number[1];   //[has_michel_slice_number_sz]
Int_t           has_michel_vertex_type_sz;
Int_t           has_michel_vertex_type[2];   //[has_michel_vertex_type_sz]
Int_t           has_michel_view_sum_sz;
Int_t           has_michel_view_sum[2];   //[has_michel_view_sum_sz]
Int_t           paddle[1];   //[n_vetoDigits]
Int_t           pmt[1];   //[n_vetoDigits]
Int_t           pmt_occupancy[1];   //[n_vetoDigits]
Int_t           wall[1];   //[n_vetoDigits]
Int_t           has_michel_distance_sz;
Double_t        has_michel_distance[2];   //[has_michel_distance_sz]
Int_t           has_michel_energy_sz;
Double_t        has_michel_energy[2];   //[has_michel_energy_sz]
Int_t           has_michel_slice_energy_sz;
Double_t        has_michel_slice_energy[2];   //[has_michel_slice_energy_sz]
Int_t           has_michel_time_diff_sz;
Double_t        has_michel_time_diff[2];   //[has_michel_time_diff_sz]
Double_t        muon_corrected_p[4];
Double_t        muon_sp[3];
Double_t        muon_sp_orig[3];
Double_t        orig_short_vtx[4];
Double_t        orig_vtx[4];
Double_t        primary_track_minerva_end_position[3];
Double_t        primary_track_minerva_start_position[3];
Double_t        time_distribution[1];   //[n_vetoDigits]
Int_t           truth_has_michel_from_pion_minus_momentum_sz;
Double_t        truth_has_michel_from_pion_minus_momentum[1];   //[truth_has_michel_from_pion_minus_momentum_sz]
Int_t           truth_has_michel_from_pion_plus_momentum_sz;
Double_t        truth_has_michel_from_pion_plus_momentum[4];   //[truth_has_michel_from_pion_plus_momentum_sz]
Double_t        usact_E_per_plane_huge[151];   //[usact_extent_huge]
Double_t        usact_E_per_plane_inf[164];   //[usact_extent_inf]
Double_t        usact_E_per_plane_large[151];   //[usact_extent_large]
Double_t        usact_E_per_plane_normal[147];   //[usact_extent_normal]
Double_t        usact_E_per_plane_small[147];   //[usact_extent_small]
Double_t        usact_E_per_plane_texas_sized[163];   //[usact_extent_texas_sized]
Double_t        usact_E_per_plane_tiny[139];   //[usact_extent_tiny]
Double_t        vtx_blob_radius;
Bool_t          truth_has_physics_event;
Bool_t          truth_muon_is_plausible;
Bool_t          truth_reco_has_int_vtx;
Bool_t          truth_reco_has_bad_object;
Bool_t          truth_reco_has_muon;
Bool_t          truth_reco_muon_has_charge;
Bool_t          truth_reco_has_good_vtx;
Bool_t          truth_reco_is_rock_muon;
Bool_t          truth_pass_MasterAnaDev;
Int_t           truth_in_analyzable_area;
Int_t           truth_in_fiducial_area;
Int_t           truth_muon_leaving_code;
Int_t           truth_muon_track_id;
Int_t           truth_reco_has_michel_electron;
Int_t           truth_targetID;
Int_t           truth_targetZ;
Int_t           truth_target_code;
Int_t           truth_vtx_module;
Int_t           truth_vtx_plane;
Double_t        truth_fs_energy_em;
Double_t        truth_fs_energy_highn;
Double_t        truth_fs_energy_hmidn;
Double_t        truth_fs_energy_lown;
Double_t        truth_fs_energy_meson;
Double_t        truth_fs_energy_mu;
Double_t        truth_fs_energy_o;
Double_t        truth_fs_energy_p;
Double_t        truth_muon_phi;
Double_t        truth_muon_theta;
Double_t        truth_muon_thetaX;
Double_t        truth_muon_thetaY;
Double_t        truth_target_dist_to_division;
Double_t        truth_target_zDist;
Int_t           truth_ref_targZ[5];
Int_t           genie_wgt_n_shifts;
Double_t        truth_genie_wgt_AGKYxF1pi[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_AhtBY[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_BhtBY[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_CCQEPauliSupViaKF[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_CV1uBY[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_CV2uBY[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_EtaNCEL[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_FrAbs_N[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_FrAbs_pi[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_FrCEx_N[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_FrCEx_pi[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_FrElas_N[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_FrElas_pi[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_FrInel_N[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_FrInel_pi[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_FrPiProd_N[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_FrPiProd_pi[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_MFP_N[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_MFP_pi[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_MaCCQE[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_MaCCQEshape[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_MaNCEL[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_MaRES[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_MvRES[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_NormCCQE[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_NormCCRES[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_NormDISCC[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_NormNCRES[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_RDecBR1gamma[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_Rvn1pi[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_Rvn2pi[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_Rvp1pi[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_Rvp2pi[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_Theta_Delta2Npi[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_VecFFCCQEshape[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_genie_wgt_shifts[MAX_GENIE_WGT_N_SHIFTS];   //[genie_wgt_n_shifts]
Double_t        truth_muon_end_momentum[4];
Double_t        truth_muon_end_position[4];
Double_t        truth_ref_dist_to_division[5];
Double_t        truth_ref_dist_to_target[5];
Int_t           MasterAnaDev_nuFlavor;
Int_t           MasterAnaDev_nuHelicity;
Int_t           MasterAnaDev_intCurrent;
Int_t           MasterAnaDev_intType;
Double_t        MasterAnaDev_E;
Double_t        MasterAnaDev_Q2;
Double_t        MasterAnaDev_x;
Double_t        MasterAnaDev_y;
Double_t        MasterAnaDev_W;
Double_t        MasterAnaDev_score;
//Double_t        MasterAnaDev_leptonE[4];
//std::copy(GetVecDouble("MasterAnaDev_leptonE").begin(), GetVecDouble("MasterAnaDev_leptonE").end(),MasterAnaDev_leptonE);
//MasterAnaDev_leptonE[3] = GetVecElem("MasterAnaDev_leptonE",3);//.begin(), GetVecDouble("MasterAnaDev_leptonE").end(),MasterAnaDev_leptonE);
Double_t        MasterAnaDev_vtx[4];
Bool_t          MasterAnaDev_minos_trk_is_contained;
Bool_t          MasterAnaDev_minos_trk_is_ok;
Bool_t          MasterAnaDev_minos_used_range;
Bool_t          MasterAnaDev_minos_used_curvature;
Bool_t          MasterAnaDev_pass_canonical_cut;
Bool_t          MasterAnaDev_is_cc;
Int_t           MasterAnaDev_in_analyzable_area;
Int_t           MasterAnaDev_in_fiducial_area;
Int_t           MasterAnaDev_minos_trk_end_plane;
Int_t           MasterAnaDev_minos_trk_quality;
Int_t           MasterAnaDev_r_minos_trk_vtx_plane;
Int_t           MasterAnaDev_t_minos_trk_numFSMuons;
Int_t           MasterAnaDev_t_minos_trk_primFSLeptonPDG;
Int_t           MasterAnaDev_targetID;
Int_t           MasterAnaDev_targetZ;
Int_t           MasterAnaDev_target_code;
Int_t           MasterAnaDev_vtx_module;
Int_t           MasterAnaDev_vtx_plane;
Double_t        MasterAnaDev_E_ccqe;
Double_t        MasterAnaDev_E_wide_window;
Double_t        MasterAnaDev_Q2_ccqe;
Double_t        MasterAnaDev_Q2_wide_window;
Double_t        MasterAnaDev_W_wide_window;
Double_t        MasterAnaDev_minos_trk_bave;
Double_t        MasterAnaDev_minos_trk_chi2;
Double_t        MasterAnaDev_minos_trk_end_u;
Double_t        MasterAnaDev_minos_trk_end_v;
Double_t        MasterAnaDev_minos_trk_end_x;
Double_t        MasterAnaDev_minos_trk_end_y;
Double_t        MasterAnaDev_minos_trk_end_z;
Double_t        MasterAnaDev_minos_trk_eqp;
Double_t        MasterAnaDev_minos_trk_eqp_qp;
Double_t        MasterAnaDev_minos_trk_fit_pass;
Double_t        MasterAnaDev_minos_trk_ndf;
Double_t        MasterAnaDev_minos_trk_p;
Double_t        MasterAnaDev_minos_trk_p_curvature;
Double_t        MasterAnaDev_minos_trk_p_range;
Double_t        MasterAnaDev_minos_trk_qp;
Double_t        MasterAnaDev_minos_trk_vtx_x;
Double_t        MasterAnaDev_minos_trk_vtx_y;
Double_t        MasterAnaDev_minos_trk_vtx_z;
Double_t        MasterAnaDev_nu_energy_recoil;
Double_t        MasterAnaDev_r_minos_trk_bdL;
Double_t        MasterAnaDev_r_minos_trk_end_dcosx;
Double_t        MasterAnaDev_r_minos_trk_end_dcosy;
Double_t        MasterAnaDev_r_minos_trk_end_dcosz;
Double_t        MasterAnaDev_r_minos_trk_vtx_dcosx;
Double_t        MasterAnaDev_r_minos_trk_vtx_dcosy;
Double_t        MasterAnaDev_r_minos_trk_vtx_dcosz;
Double_t        MasterAnaDev_recoil_E;
Double_t        MasterAnaDev_recoil_E_wide_window;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMinosInitProjPx;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMinosInitProjPy;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMinosInitProjPz;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMinosInitProjX;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMinosInitProjY;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMinosInitProjZ;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMnvFinalPx;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMnvFinalPy;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMnvFinalPz;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMnvFinalX;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMnvFinalY;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMnvFinalZ;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMnvInitPx;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMnvInitPy;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMnvInitPz;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMnvInitX;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMnvInitY;
Double_t        MasterAnaDev_t_minos_trk_primFSLepMnvInitZ;
Double_t        MasterAnaDev_target_dist_to_division;
Double_t        MasterAnaDev_target_zDist;
Double_t        MasterAnaDev_x_wide_window;
Double_t        MasterAnaDev_y_wide_window;
Int_t           MasterAnaDev_ref_targZ[5];
Int_t           MasterAnaDev_smeared_is_fiducial[100];
Int_t           MasterAnaDev_smeared_pass_dist_to_div[100];
Int_t           MasterAnaDev_smeared_targetID[100];
Int_t           MasterAnaDev_smeared_targetZ[100];
Int_t           MasterAnaDev_smeared_vtx_mod[100];
Double_t        MasterAnaDev_muon_vtx[4];
Double_t        MasterAnaDev_ref_dist_to_division[5];
Double_t        MasterAnaDev_ref_dist_to_target[5];
Double_t        MasterAnaDev_sys_muon_curve_energy_shift[2];
Double_t        MasterAnaDev_sys_muon_energy_shift[2];
Double_t        MasterAnaDev_sys_muon_minerva_energy_shift[2];
Double_t        MasterAnaDev_sys_muon_qSquared_shift[2];
Double_t        MasterAnaDev_sys_muon_range_energy_shift[2];
Double_t        MasterAnaDev_sys_muon_wSquared_shift[2];
Double_t        MasterAnaDev_sys_muon_xbj_shift[2];
Double_t        MasterAnaDev_sys_muon_y_shift[2];
Double_t        MasterAnaDev_sys_nu_energy_shift[2];
Double_t        MasterAnaDev_sys_recoil_energy_shift[2];
Double_t        MasterAnaDev_sys_recoil_qSquared_shift[2];
Double_t        MasterAnaDev_sys_recoil_wSquared_shift[2];
Double_t        MasterAnaDev_sys_recoil_xbj_shift[2];
Double_t        MasterAnaDev_sys_recoil_y_shift[2];
Double_t        MasterAnaDev_sys_total_qSquared_shift[2];
Double_t        MasterAnaDev_sys_total_wSquared_shift[2];
Double_t        MasterAnaDev_sys_total_xbj_shift[2];
Double_t        MasterAnaDev_sys_total_y_shift[2];
Int_t           ev_run;
Int_t           ev_subrun;
Int_t           ev_detector;
Int_t           ev_triggerType;
Int_t           ev_gate;
Int_t           ev_global_gate;
Int_t           ev_gps_time_sec;
Int_t           ev_gps_time_usec;
Int_t           mc_run;
Int_t           mc_subrun;
Int_t           mc_nInteractions;
Int_t           mc_MIState;
Double_t        mc_pot;
Int_t           mc_beamConfig;
Int_t           mc_processType;
Int_t           mc_nthEvtInSpill;
Int_t           mc_nthEvtInFile;
Int_t           mc_intType;
Int_t           mc_current;
Int_t           mc_charm;
Double_t        mc_weight;
Double_t        mc_XSec;
Double_t        mc_diffXSec;
Int_t           mc_incoming;
Double_t        mc_fluxDriverProb;
Int_t           mc_targetNucleus;
Int_t           mc_targetZ;
Int_t           mc_targetA;
Int_t           mc_targetNucleon;
Int_t           mc_struckQuark;
Int_t           mc_seaQuark;
Int_t           mc_resID;
Int_t           mc_primaryLepton;
Double_t        mc_incomingE;
Double_t        mc_Bjorkenx;
Double_t        mc_Bjorkeny;
Double_t        mc_Q2;
Double_t        mc_nuT;
Double_t        mc_w;
Double_t        mc_vtx[4];
Double_t        mc_incomingPartVec[4];
Double_t        mc_initNucVec[4];
Double_t        mc_primFSLepton[4];
Int_t           mc_nFSPart;
Double_t        mc_FSPartPx[MAX_MC_NFSPART];   //[mc_nFSPart]
Double_t        mc_FSPartPy[MAX_MC_NFSPART];   //[mc_nFSPart]
Double_t        mc_FSPartPz[MAX_MC_NFSPART];   //[mc_nFSPart]
Double_t        mc_FSPartE[MAX_MC_NFSPART];   //[mc_nFSPart]
Int_t           mc_FSPartPDG[MAX_MC_NFSPART];   //[mc_nFSPart]
Int_t           mc_er_nPart;
Int_t           mc_er_ID[MAX_MC_ER_NPART];   //[mc_er_nPart]
Int_t           mc_er_status[MAX_MC_ER_NPART];   //[mc_er_nPart]
Double_t        mc_er_posInNucX[MAX_MC_ER_NPART];   //[mc_er_nPart]
Double_t        mc_er_posInNucY[MAX_MC_ER_NPART];   //[mc_er_nPart]
Double_t        mc_er_posInNucZ[MAX_MC_ER_NPART];   //[mc_er_nPart]
Double_t        mc_er_Px[MAX_MC_ER_NPART];   //[mc_er_nPart]
Double_t        mc_er_Py[MAX_MC_ER_NPART];   //[mc_er_nPart]
Double_t        mc_er_Pz[MAX_MC_ER_NPART];   //[mc_er_nPart]
Double_t        mc_er_E[MAX_MC_ER_NPART];   //[mc_er_nPart]
Int_t           mc_er_FD[MAX_MC_ER_NPART];   //[mc_er_nPart]
Int_t           mc_er_LD[MAX_MC_ER_NPART];   //[mc_er_nPart]
Int_t           mc_er_mother[MAX_MC_ER_NPART];   //[mc_er_nPart]
Int_t           mc_fr_nNuAncestorIDs;
Int_t           mc_fr_nuAncestorIDs[MAX_MC_FR_NNUANCESTORIDS];   //[mc_fr_nNuAncestorIDs]
Int_t           mc_fr_nuParentID;
Int_t           mc_fr_decMode;
Double_t        mc_fr_primProtonVtx[3];
Double_t        mc_fr_primProtonP[4];
Double_t        mc_fr_nuParentDecVtx[3];
Double_t        mc_fr_nuParentProdVtx[3];
Double_t        mc_fr_nuParentProdP[4];
Double_t        mc_cvweight_total;
Double_t        wgt;
Double_t        mc_cvweight_totalFlux;
Double_t        mc_cvweight_totalXsec;
Double_t        mc_ppfx1_cvweight;
Double_t        mc_hornCurrent_cvweight;
Double_t        mc_gen1_cvweight_total;
Double_t        gen1_wgt;
Double_t        mc_gen1_cvweight_totalFlux;
Double_t        mc_gen1_cvweight_NA49;
Int_t           mc_wgt_Flux_BeamFocus_sz;
Double_t        mc_wgt_Flux_BeamFocus[MAX_MC_WGT_FLUX_BEAMFOCUS_SZ];   //[mc_wgt_Flux_BeamFocus_sz]
Int_t           mc_wgt_gen1_Flux_Tertiary_sz;
Double_t        mc_wgt_gen1_Flux_Tertiary[100];   //[mc_wgt_gen1_Flux_Tertiary_sz]
Int_t           mc_wgt_gen1_Flux_NA49_sz;
Double_t        mc_wgt_gen1_Flux_NA49[100];   //[mc_wgt_gen1_Flux_NA49_sz]
Int_t           mc_wgt_Norm_sz;
Double_t        mc_wgt_Norm[1];   //[mc_wgt_Norm_sz]
Int_t           mc_wgt_ppfx1_Total_sz;
Double_t        mc_wgt_ppfx1_Total[100];   //[mc_wgt_ppfx1_Total_sz]
Int_t           n_prongs;
Int_t           prong_nParticles[MAX_N_PRONGS];   //[n_prongs]
Double_t        prong_part_score[MAX_N_PRONGS];   //[n_prongs]
Double_t        prong_part_mass[MAX_N_PRONGS];   //[n_prongs]
Int_t           prong_part_charge[MAX_N_PRONGS];   //[n_prongs]
Int_t           prong_part_pid[MAX_N_PRONGS];   //[n_prongs]
vector<vector<double> > *prong_part_E;
vector<vector<double> > *prong_part_pos;

// List of branches
TBranch        *b_eventID;   //!
TBranch        *b_physEvtNum;   //!
TBranch        *b_n_hyps;   //!
TBranch        *b_processType;   //!
TBranch        *b_primaryPart;   //!
TBranch        *b_n_slices;   //!
TBranch        *b_slice_numbers;   //!
TBranch        *b_shared_slice;   //!
TBranch        *b_vtx;   //!
TBranch        *b_vtxErr;   //!
TBranch        *b_E;   //!
TBranch        *b_found_truth;   //!
TBranch        *b_phys_front_activity;   //!
TBranch        *b_phys_energy_in_road_upstream_is_rockmuon_consistent;   //!
TBranch        *b_rock_muons_removed;   //!
TBranch        *b_minos_track_match;   //!
TBranch        *b_minos_stub_match;   //!
TBranch        *b_unknown_helicity;   //!
TBranch        *b_minos_track_inside_partial_plane;   //!
TBranch        *b_prim_vtx_has_misassigned_track_direction;   //!
TBranch        *b_prim_vtx_has_broken_track;   //!
TBranch        *b_pass_MasterAnaDev;   //!
TBranch        *b_short_track_vtx_used;   //!
TBranch        *b_muon_sp_moved;   //!
TBranch        *b_vtx_fit_converged;   //!
TBranch        *b_muon_is_correct;   //!
TBranch        *b_has_int_vtx;   //!
TBranch        *b_has_bad_object;   //!
TBranch        *b_has_muon;   //!
TBranch        *b_muon_has_charge;   //!
TBranch        *b_has_good_vtx;   //!
TBranch        *b_is_rock_muon;   //!
TBranch        *b_MuonTaggedAsVetoButNotMatched;   //!
TBranch        *b_NonVetoMuonExtrpToVeto;   //!
TBranch        *b_NonVetoMuonWallOnePaddle;   //!
TBranch        *b_NonVetoMuonWallOnePaddleOverlap;   //!
TBranch        *b_NonVetoMuonWallOneSector;   //!
TBranch        *b_NonVetoMuonWallTwoPaddle;   //!
TBranch        *b_NonVetoMuonWallTwoPaddleOverlap;   //!
TBranch        *b_NonVetoMuonWallTwoSector;   //!
TBranch        *b_VetoMuonWallOnePaddle;   //!
TBranch        *b_VetoMuonWallOnePaddleOverlap;   //!
TBranch        *b_VetoMuonWallOneSector;   //!
TBranch        *b_VetoMuonWallOneTypeOfMatchNonOverlap;   //!
TBranch        *b_VetoMuonWallOneTypeOfMatchOverlap;   //!
TBranch        *b_VetoMuonWallTwoPaddle;   //!
TBranch        *b_VetoMuonWallTwoPaddleOverlap;   //!
TBranch        *b_VetoMuonWallTwoSector;   //!
TBranch        *b_VetoMuonWallTwoTypeOfMatchNonOverlap;   //!
TBranch        *b_VetoMuonWallTwoTypeOfMatchOverlap;   //!
TBranch        *b_blob_disp_nBlobs;   //!
TBranch        *b_blob_disp_nClus;   //!
TBranch        *b_blob_disp_nClus_ecal;   //!
TBranch        *b_blob_disp_nClus_hcal;   //!
TBranch        *b_blob_disp_nClus_nucl;   //!
TBranch        *b_blob_disp_nClus_od;   //!
TBranch        *b_blob_disp_nClus_tracker;   //!
TBranch        *b_blob_iso_nBlobs;   //!
TBranch        *b_blob_iso_nClus;   //!
TBranch        *b_blob_iso_nClus_ecal;   //!
TBranch        *b_blob_iso_nClus_hcal;   //!
TBranch        *b_blob_iso_nClus_nucl;   //!
TBranch        *b_blob_iso_nClus_od;   //!
TBranch        *b_blob_iso_nClus_tracker;   //!
TBranch        *b_blob_mufuzz_nBlobs;   //!
TBranch        *b_blob_mufuzz_nClus;   //!
TBranch        *b_blob_mufuzz_nClus_ecal;   //!
TBranch        *b_blob_mufuzz_nClus_hcal;   //!
TBranch        *b_blob_mufuzz_nClus_nucl;   //!
TBranch        *b_blob_mufuzz_nClus_od;   //!
TBranch        *b_blob_mufuzz_nClus_tracker;   //!
TBranch        *b_blob_recoil_nBlobs;   //!
TBranch        *b_blob_recoil_nClus;   //!
TBranch        *b_blob_recoil_nClus_ecal;   //!
TBranch        *b_blob_recoil_nClus_hcal;   //!
TBranch        *b_blob_recoil_nClus_nucl;   //!
TBranch        *b_blob_recoil_nClus_od;   //!
TBranch        *b_blob_recoil_nClus_tracker;   //!
TBranch        *b_blob_vtx_nBlobs;   //!
TBranch        *b_blob_vtx_nClus;   //!
TBranch        *b_blob_vtx_nClus_ecal;   //!
TBranch        *b_blob_vtx_nClus_hcal;   //!
TBranch        *b_blob_vtx_nClus_nucl;   //!
TBranch        *b_blob_vtx_nClus_od;   //!
TBranch        *b_blob_vtx_nClus_tracker;   //!
TBranch        *b_broken_track_most_us_plane;   //!
TBranch        *b_muon_n_USclusters;   //!
TBranch        *b_muon_truthMatch_track_id;   //!
TBranch        *b_n_prim_long_tracks;   //!
TBranch        *b_n_prim_short_tracks;   //!
TBranch        *b_n_start_vertices;   //!
TBranch        *b_n_tracks;   //!
TBranch        *b_n_tracks_non_prim;   //!
TBranch        *b_n_tracks_prim;   //!
TBranch        *b_n_tracks_prim_forked;   //!
TBranch        *b_n_tracks_prim_kinked;   //!
TBranch        *b_n_vertices_startpoint;   //!
TBranch        *b_passVetoMuonCut;   //!
TBranch        *b_phys_energy_in_road_downstream_nplanes;   //!
TBranch        *b_phys_energy_in_road_upstream_nplanes;   //!
TBranch        *b_phys_n_dead_discr_pair;   //!
TBranch        *b_phys_n_dead_discr_pair_in_prim_track_region;   //!
TBranch        *b_phys_n_dead_discr_pair_two_mod_downstream_prim_track;   //!
TBranch        *b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;   //!
TBranch        *b_phys_n_dead_discr_pair_upstream_prim_track_proj;   //!
TBranch        *b_phys_vertex_is_fiducial;   //!
TBranch        *b_rock_muon_code;   //!
TBranch        *b_truth_has_michel_electron;   //!
TBranch        *b_usact_extent_huge;   //!
TBranch        *b_usact_extent_inf;   //!
TBranch        *b_usact_extent_large;   //!
TBranch        *b_usact_extent_normal;   //!
TBranch        *b_usact_extent_small;   //!
TBranch        *b_usact_extent_texas_sized;   //!
TBranch        *b_usact_extent_tiny;   //!
TBranch        *b_usact_n_consecutive_huge;   //!
TBranch        *b_usact_n_consecutive_inf;   //!
TBranch        *b_usact_n_consecutive_large;   //!
TBranch        *b_usact_n_consecutive_normal;   //!
TBranch        *b_usact_n_consecutive_small;   //!
TBranch        *b_usact_n_consecutive_texas_sized;   //!
TBranch        *b_usact_n_consecutive_tiny;   //!
TBranch        *b_usact_n_planes_huge;   //!
TBranch        *b_usact_n_planes_inf;   //!
TBranch        *b_usact_n_planes_large;   //!
TBranch        *b_usact_n_planes_normal;   //!
TBranch        *b_usact_n_planes_small;   //!
TBranch        *b_usact_n_planes_texas_sized;   //!
TBranch        *b_usact_n_planes_tiny;   //!
TBranch        *b_NonVetoMuonWallOneBadPosX;   //!
TBranch        *b_NonVetoMuonWallOneBadPosY;   //!
TBranch        *b_NonVetoMuonWallOnePosX;   //!
TBranch        *b_NonVetoMuonWallOnePosY;   //!
TBranch        *b_NonVetoMuonWallOne_ANDEfficiency_Central;   //!
TBranch        *b_NonVetoMuonWallOne_ANDEfficiency_Overlap;   //!
TBranch        *b_NonVetoMuonWallOne_ANDError_Central;   //!
TBranch        *b_NonVetoMuonWallOne_ANDError_Overlap;   //!
TBranch        *b_NonVetoMuonWallOne_AccRatesError_PaddleAbove;   //!
TBranch        *b_NonVetoMuonWallOne_AccRatesError_PaddleBelow;   //!
TBranch        *b_NonVetoMuonWallOne_AccRates_PaddleAbove;   //!
TBranch        *b_NonVetoMuonWallOne_AccRates_PaddleBelow;   //!
TBranch        *b_NonVetoMuonWallOne_OREfficiency_Central;   //!
TBranch        *b_NonVetoMuonWallOne_OREfficiency_Overlap;   //!
TBranch        *b_NonVetoMuonWallOne_ORError_Central;   //!
TBranch        *b_NonVetoMuonWallOne_ORError_Overlap;   //!
TBranch        *b_NonVetoMuonWallTwoBadPosX;   //!
TBranch        *b_NonVetoMuonWallTwoBadPosY;   //!
TBranch        *b_NonVetoMuonWallTwoPosX;   //!
TBranch        *b_NonVetoMuonWallTwoPosY;   //!
TBranch        *b_NonVetoMuonWallTwo_ANDEfficiency_Central;   //!
TBranch        *b_NonVetoMuonWallTwo_ANDEfficiency_Overlap;   //!
TBranch        *b_NonVetoMuonWallTwo_ANDError_Central;   //!
TBranch        *b_NonVetoMuonWallTwo_ANDError_Overlap;   //!
TBranch        *b_NonVetoMuonWallTwo_AccRatesError_PaddleAbove;   //!
TBranch        *b_NonVetoMuonWallTwo_AccRatesError_PaddleBelow;   //!
TBranch        *b_NonVetoMuonWallTwo_AccRates_PaddleAbove;   //!
TBranch        *b_NonVetoMuonWallTwo_AccRates_PaddleBelow;   //!
TBranch        *b_NonVetoMuonWallTwo_OREfficiency_Central;   //!
TBranch        *b_NonVetoMuonWallTwo_OREfficiency_Overlap;   //!
TBranch        *b_NonVetoMuonWallTwo_ORError_Central;   //!
TBranch        *b_NonVetoMuonWallTwo_ORError_Overlap;   //!
TBranch        *b_VetoMuonWallOneDeltaTime;   //!
TBranch        *b_VetoMuonWallOnePosX;   //!
TBranch        *b_VetoMuonWallOnePosY;   //!
TBranch        *b_VetoMuonWallOne_ANDEfficiency_Central;   //!
TBranch        *b_VetoMuonWallOne_ANDEfficiency_Overlap;   //!
TBranch        *b_VetoMuonWallOne_ANDError_Central;   //!
TBranch        *b_VetoMuonWallOne_ANDError_Overlap;   //!
TBranch        *b_VetoMuonWallOne_AccRatesError_PaddleAbove;   //!
TBranch        *b_VetoMuonWallOne_AccRatesError_PaddleBelow;   //!
TBranch        *b_VetoMuonWallOne_AccRates_PaddleAbove;   //!
TBranch        *b_VetoMuonWallOne_AccRates_PaddleBelow;   //!
TBranch        *b_VetoMuonWallOne_OREfficiency_Central;   //!
TBranch        *b_VetoMuonWallOne_OREfficiency_Overlap;   //!
TBranch        *b_VetoMuonWallOne_ORError_Central;   //!
TBranch        *b_VetoMuonWallOne_ORError_Overlap;   //!
TBranch        *b_VetoMuonWallTwoDeltaTime;   //!
TBranch        *b_VetoMuonWallTwoPosX;   //!
TBranch        *b_VetoMuonWallTwoPosY;   //!
TBranch        *b_VetoMuonWallTwo_ANDEfficiency_Central;   //!
TBranch        *b_VetoMuonWallTwo_ANDEfficiency_Overlap;   //!
TBranch        *b_VetoMuonWallTwo_ANDError_Central;   //!
TBranch        *b_VetoMuonWallTwo_ANDError_Overlap;   //!
TBranch        *b_VetoMuonWallTwo_AccRatesError_PaddleAbove;   //!
TBranch        *b_VetoMuonWallTwo_AccRatesError_PaddleBelow;   //!
TBranch        *b_VetoMuonWallTwo_AccRates_PaddleAbove;   //!
TBranch        *b_VetoMuonWallTwo_AccRates_PaddleBelow;   //!
TBranch        *b_VetoMuonWallTwo_OREfficiency_Central;   //!
TBranch        *b_VetoMuonWallTwo_OREfficiency_Overlap;   //!
TBranch        *b_VetoMuonWallTwo_ORError_Central;   //!
TBranch        *b_VetoMuonWallTwo_ORError_Overlap;   //!
TBranch        *b_blob_ccqe_recoil_E;   //!
TBranch        *b_blob_disp_E;   //!
TBranch        *b_blob_disp_E_ecal;   //!
TBranch        *b_blob_disp_E_ecal_em;   //!
TBranch        *b_blob_disp_E_ecal_highn;   //!
TBranch        *b_blob_disp_E_ecal_lown;   //!
TBranch        *b_blob_disp_E_ecal_meson;   //!
TBranch        *b_blob_disp_E_ecal_midn;   //!
TBranch        *b_blob_disp_E_ecal_mu;   //!
TBranch        *b_blob_disp_E_ecal_other;   //!
TBranch        *b_blob_disp_E_ecal_p;   //!
TBranch        *b_blob_disp_E_ecal_xtalk;   //!
TBranch        *b_blob_disp_E_em;   //!
TBranch        *b_blob_disp_E_hcal;   //!
TBranch        *b_blob_disp_E_hcal_em;   //!
TBranch        *b_blob_disp_E_hcal_highn;   //!
TBranch        *b_blob_disp_E_hcal_lown;   //!
TBranch        *b_blob_disp_E_hcal_meson;   //!
TBranch        *b_blob_disp_E_hcal_midn;   //!
TBranch        *b_blob_disp_E_hcal_mu;   //!
TBranch        *b_blob_disp_E_hcal_other;   //!
TBranch        *b_blob_disp_E_hcal_p;   //!
TBranch        *b_blob_disp_E_hcal_xtalk;   //!
TBranch        *b_blob_disp_E_highn;   //!
TBranch        *b_blob_disp_E_lown;   //!
TBranch        *b_blob_disp_E_meson;   //!
TBranch        *b_blob_disp_E_midn;   //!
TBranch        *b_blob_disp_E_mu;   //!
TBranch        *b_blob_disp_E_nucl;   //!
TBranch        *b_blob_disp_E_nucl_em;   //!
TBranch        *b_blob_disp_E_nucl_highn;   //!
TBranch        *b_blob_disp_E_nucl_lown;   //!
TBranch        *b_blob_disp_E_nucl_meson;   //!
TBranch        *b_blob_disp_E_nucl_midn;   //!
TBranch        *b_blob_disp_E_nucl_mu;   //!
TBranch        *b_blob_disp_E_nucl_other;   //!
TBranch        *b_blob_disp_E_nucl_p;   //!
TBranch        *b_blob_disp_E_nucl_xtalk;   //!
TBranch        *b_blob_disp_E_od;   //!
TBranch        *b_blob_disp_E_od_em;   //!
TBranch        *b_blob_disp_E_od_highn;   //!
TBranch        *b_blob_disp_E_od_lown;   //!
TBranch        *b_blob_disp_E_od_meson;   //!
TBranch        *b_blob_disp_E_od_midn;   //!
TBranch        *b_blob_disp_E_od_mu;   //!
TBranch        *b_blob_disp_E_od_other;   //!
TBranch        *b_blob_disp_E_od_p;   //!
TBranch        *b_blob_disp_E_od_xtalk;   //!
TBranch        *b_blob_disp_E_other;   //!
TBranch        *b_blob_disp_E_p;   //!
TBranch        *b_blob_disp_E_tracker;   //!
TBranch        *b_blob_disp_E_tracker_em;   //!
TBranch        *b_blob_disp_E_tracker_highn;   //!
TBranch        *b_blob_disp_E_tracker_lown;   //!
TBranch        *b_blob_disp_E_tracker_meson;   //!
TBranch        *b_blob_disp_E_tracker_midn;   //!
TBranch        *b_blob_disp_E_tracker_mu;   //!
TBranch        *b_blob_disp_E_tracker_other;   //!
TBranch        *b_blob_disp_E_tracker_p;   //!
TBranch        *b_blob_disp_E_tracker_xtalk;   //!
TBranch        *b_blob_disp_E_xtalk;   //!
TBranch        *b_blob_iso_E;   //!
TBranch        *b_blob_iso_E_ecal;   //!
TBranch        *b_blob_iso_E_ecal_em;   //!
TBranch        *b_blob_iso_E_ecal_highn;   //!
TBranch        *b_blob_iso_E_ecal_lown;   //!
TBranch        *b_blob_iso_E_ecal_meson;   //!
TBranch        *b_blob_iso_E_ecal_midn;   //!
TBranch        *b_blob_iso_E_ecal_mu;   //!
TBranch        *b_blob_iso_E_ecal_other;   //!
TBranch        *b_blob_iso_E_ecal_p;   //!
TBranch        *b_blob_iso_E_ecal_xtalk;   //!
TBranch        *b_blob_iso_E_em;   //!
TBranch        *b_blob_iso_E_hcal;   //!
TBranch        *b_blob_iso_E_hcal_em;   //!
TBranch        *b_blob_iso_E_hcal_highn;   //!
TBranch        *b_blob_iso_E_hcal_lown;   //!
TBranch        *b_blob_iso_E_hcal_meson;   //!
TBranch        *b_blob_iso_E_hcal_midn;   //!
TBranch        *b_blob_iso_E_hcal_mu;   //!
TBranch        *b_blob_iso_E_hcal_other;   //!
TBranch        *b_blob_iso_E_hcal_p;   //!
TBranch        *b_blob_iso_E_hcal_xtalk;   //!
TBranch        *b_blob_iso_E_highn;   //!
TBranch        *b_blob_iso_E_lown;   //!
TBranch        *b_blob_iso_E_meson;   //!
TBranch        *b_blob_iso_E_midn;   //!
TBranch        *b_blob_iso_E_mu;   //!
TBranch        *b_blob_iso_E_nucl;   //!
TBranch        *b_blob_iso_E_nucl_em;   //!
TBranch        *b_blob_iso_E_nucl_highn;   //!
TBranch        *b_blob_iso_E_nucl_lown;   //!
TBranch        *b_blob_iso_E_nucl_meson;   //!
TBranch        *b_blob_iso_E_nucl_midn;   //!
TBranch        *b_blob_iso_E_nucl_mu;   //!
TBranch        *b_blob_iso_E_nucl_other;   //!
TBranch        *b_blob_iso_E_nucl_p;   //!
TBranch        *b_blob_iso_E_nucl_xtalk;   //!
TBranch        *b_blob_iso_E_od;   //!
TBranch        *b_blob_iso_E_od_em;   //!
TBranch        *b_blob_iso_E_od_highn;   //!
TBranch        *b_blob_iso_E_od_lown;   //!
TBranch        *b_blob_iso_E_od_meson;   //!
TBranch        *b_blob_iso_E_od_midn;   //!
TBranch        *b_blob_iso_E_od_mu;   //!
TBranch        *b_blob_iso_E_od_other;   //!
TBranch        *b_blob_iso_E_od_p;   //!
TBranch        *b_blob_iso_E_od_xtalk;   //!
TBranch        *b_blob_iso_E_other;   //!
TBranch        *b_blob_iso_E_p;   //!
TBranch        *b_blob_iso_E_tracker;   //!
TBranch        *b_blob_iso_E_tracker_em;   //!
TBranch        *b_blob_iso_E_tracker_highn;   //!
TBranch        *b_blob_iso_E_tracker_lown;   //!
TBranch        *b_blob_iso_E_tracker_meson;   //!
TBranch        *b_blob_iso_E_tracker_midn;   //!
TBranch        *b_blob_iso_E_tracker_mu;   //!
TBranch        *b_blob_iso_E_tracker_other;   //!
TBranch        *b_blob_iso_E_tracker_p;   //!
TBranch        *b_blob_iso_E_tracker_xtalk;   //!
TBranch        *b_blob_iso_E_xtalk;   //!
TBranch        *b_blob_mufuzz_E;   //!
TBranch        *b_blob_mufuzz_E_ecal;   //!
TBranch        *b_blob_mufuzz_E_ecal_em;   //!
TBranch        *b_blob_mufuzz_E_ecal_highn;   //!
TBranch        *b_blob_mufuzz_E_ecal_lown;   //!
TBranch        *b_blob_mufuzz_E_ecal_meson;   //!
TBranch        *b_blob_mufuzz_E_ecal_midn;   //!
TBranch        *b_blob_mufuzz_E_ecal_mu;   //!
TBranch        *b_blob_mufuzz_E_ecal_other;   //!
TBranch        *b_blob_mufuzz_E_ecal_p;   //!
TBranch        *b_blob_mufuzz_E_ecal_xtalk;   //!
TBranch        *b_blob_mufuzz_E_em;   //!
TBranch        *b_blob_mufuzz_E_hcal;   //!
TBranch        *b_blob_mufuzz_E_hcal_em;   //!
TBranch        *b_blob_mufuzz_E_hcal_highn;   //!
TBranch        *b_blob_mufuzz_E_hcal_lown;   //!
TBranch        *b_blob_mufuzz_E_hcal_meson;   //!
TBranch        *b_blob_mufuzz_E_hcal_midn;   //!
TBranch        *b_blob_mufuzz_E_hcal_mu;   //!
TBranch        *b_blob_mufuzz_E_hcal_other;   //!
TBranch        *b_blob_mufuzz_E_hcal_p;   //!
TBranch        *b_blob_mufuzz_E_hcal_xtalk;   //!
TBranch        *b_blob_mufuzz_E_highn;   //!
TBranch        *b_blob_mufuzz_E_lown;   //!
TBranch        *b_blob_mufuzz_E_meson;   //!
TBranch        *b_blob_mufuzz_E_midn;   //!
TBranch        *b_blob_mufuzz_E_mu;   //!
TBranch        *b_blob_mufuzz_E_nucl;   //!
TBranch        *b_blob_mufuzz_E_nucl_em;   //!
TBranch        *b_blob_mufuzz_E_nucl_highn;   //!
TBranch        *b_blob_mufuzz_E_nucl_lown;   //!
TBranch        *b_blob_mufuzz_E_nucl_meson;   //!
TBranch        *b_blob_mufuzz_E_nucl_midn;   //!
TBranch        *b_blob_mufuzz_E_nucl_mu;   //!
TBranch        *b_blob_mufuzz_E_nucl_other;   //!
TBranch        *b_blob_mufuzz_E_nucl_p;   //!
TBranch        *b_blob_mufuzz_E_nucl_xtalk;   //!
TBranch        *b_blob_mufuzz_E_od;   //!
TBranch        *b_blob_mufuzz_E_od_em;   //!
TBranch        *b_blob_mufuzz_E_od_highn;   //!
TBranch        *b_blob_mufuzz_E_od_lown;   //!
TBranch        *b_blob_mufuzz_E_od_meson;   //!
TBranch        *b_blob_mufuzz_E_od_midn;   //!
TBranch        *b_blob_mufuzz_E_od_mu;   //!
TBranch        *b_blob_mufuzz_E_od_other;   //!
TBranch        *b_blob_mufuzz_E_od_p;   //!
TBranch        *b_blob_mufuzz_E_od_xtalk;   //!
TBranch        *b_blob_mufuzz_E_other;   //!
TBranch        *b_blob_mufuzz_E_p;   //!
TBranch        *b_blob_mufuzz_E_tracker;   //!
TBranch        *b_blob_mufuzz_E_tracker_em;   //!
TBranch        *b_blob_mufuzz_E_tracker_highn;   //!
TBranch        *b_blob_mufuzz_E_tracker_lown;   //!
TBranch        *b_blob_mufuzz_E_tracker_meson;   //!
TBranch        *b_blob_mufuzz_E_tracker_midn;   //!
TBranch        *b_blob_mufuzz_E_tracker_mu;   //!
TBranch        *b_blob_mufuzz_E_tracker_other;   //!
TBranch        *b_blob_mufuzz_E_tracker_p;   //!
TBranch        *b_blob_mufuzz_E_tracker_xtalk;   //!
TBranch        *b_blob_mufuzz_E_xtalk;   //!
TBranch        *b_blob_recoil_E;   //!
TBranch        *b_blob_recoil_E_ecal;   //!
TBranch        *b_blob_recoil_E_ecal_em;   //!
TBranch        *b_blob_recoil_E_ecal_highn;   //!
TBranch        *b_blob_recoil_E_ecal_lown;   //!
TBranch        *b_blob_recoil_E_ecal_meson;   //!
TBranch        *b_blob_recoil_E_ecal_midn;   //!
TBranch        *b_blob_recoil_E_ecal_mu;   //!
TBranch        *b_blob_recoil_E_ecal_other;   //!
TBranch        *b_blob_recoil_E_ecal_p;   //!
TBranch        *b_blob_recoil_E_ecal_xtalk;   //!
TBranch        *b_blob_recoil_E_em;   //!
TBranch        *b_blob_recoil_E_hcal;   //!
TBranch        *b_blob_recoil_E_hcal_em;   //!
TBranch        *b_blob_recoil_E_hcal_highn;   //!
TBranch        *b_blob_recoil_E_hcal_lown;   //!
TBranch        *b_blob_recoil_E_hcal_meson;   //!
TBranch        *b_blob_recoil_E_hcal_midn;   //!
TBranch        *b_blob_recoil_E_hcal_mu;   //!
TBranch        *b_blob_recoil_E_hcal_other;   //!
TBranch        *b_blob_recoil_E_hcal_p;   //!
TBranch        *b_blob_recoil_E_hcal_xtalk;   //!
TBranch        *b_blob_recoil_E_highn;   //!
TBranch        *b_blob_recoil_E_lown;   //!
TBranch        *b_blob_recoil_E_meson;   //!
TBranch        *b_blob_recoil_E_midn;   //!
TBranch        *b_blob_recoil_E_mu;   //!
TBranch        *b_blob_recoil_E_nucl;   //!
TBranch        *b_blob_recoil_E_nucl_em;   //!
TBranch        *b_blob_recoil_E_nucl_highn;   //!
TBranch        *b_blob_recoil_E_nucl_lown;   //!
TBranch        *b_blob_recoil_E_nucl_meson;   //!
TBranch        *b_blob_recoil_E_nucl_midn;   //!
TBranch        *b_blob_recoil_E_nucl_mu;   //!
TBranch        *b_blob_recoil_E_nucl_other;   //!
TBranch        *b_blob_recoil_E_nucl_p;   //!
TBranch        *b_blob_recoil_E_nucl_xtalk;   //!
TBranch        *b_blob_recoil_E_od;   //!
TBranch        *b_blob_recoil_E_od_em;   //!
TBranch        *b_blob_recoil_E_od_highn;   //!
TBranch        *b_blob_recoil_E_od_lown;   //!
TBranch        *b_blob_recoil_E_od_meson;   //!
TBranch        *b_blob_recoil_E_od_midn;   //!
TBranch        *b_blob_recoil_E_od_mu;   //!
TBranch        *b_blob_recoil_E_od_other;   //!
TBranch        *b_blob_recoil_E_od_p;   //!
TBranch        *b_blob_recoil_E_od_xtalk;   //!
TBranch        *b_blob_recoil_E_other;   //!
TBranch        *b_blob_recoil_E_p;   //!
TBranch        *b_blob_recoil_E_tracker;   //!
TBranch        *b_blob_recoil_E_tracker_em;   //!
TBranch        *b_blob_recoil_E_tracker_highn;   //!
TBranch        *b_blob_recoil_E_tracker_lown;   //!
TBranch        *b_blob_recoil_E_tracker_meson;   //!
TBranch        *b_blob_recoil_E_tracker_midn;   //!
TBranch        *b_blob_recoil_E_tracker_mu;   //!
TBranch        *b_blob_recoil_E_tracker_other;   //!
TBranch        *b_blob_recoil_E_tracker_p;   //!
TBranch        *b_blob_recoil_E_tracker_xtalk;   //!
TBranch        *b_blob_recoil_E_xtalk;   //!
TBranch        *b_blob_vtx_E;   //!
TBranch        *b_blob_vtx_E_ecal;   //!
TBranch        *b_blob_vtx_E_ecal_em;   //!
TBranch        *b_blob_vtx_E_ecal_highn;   //!
TBranch        *b_blob_vtx_E_ecal_lown;   //!
TBranch        *b_blob_vtx_E_ecal_meson;   //!
TBranch        *b_blob_vtx_E_ecal_midn;   //!
TBranch        *b_blob_vtx_E_ecal_mu;   //!
TBranch        *b_blob_vtx_E_ecal_other;   //!
TBranch        *b_blob_vtx_E_ecal_p;   //!
TBranch        *b_blob_vtx_E_ecal_xtalk;   //!
TBranch        *b_blob_vtx_E_em;   //!
TBranch        *b_blob_vtx_E_hcal;   //!
TBranch        *b_blob_vtx_E_hcal_em;   //!
TBranch        *b_blob_vtx_E_hcal_highn;   //!
TBranch        *b_blob_vtx_E_hcal_lown;   //!
TBranch        *b_blob_vtx_E_hcal_meson;   //!
TBranch        *b_blob_vtx_E_hcal_midn;   //!
TBranch        *b_blob_vtx_E_hcal_mu;   //!
TBranch        *b_blob_vtx_E_hcal_other;   //!
TBranch        *b_blob_vtx_E_hcal_p;   //!
TBranch        *b_blob_vtx_E_hcal_xtalk;   //!
TBranch        *b_blob_vtx_E_highn;   //!
TBranch        *b_blob_vtx_E_lown;   //!
TBranch        *b_blob_vtx_E_meson;   //!
TBranch        *b_blob_vtx_E_midn;   //!
TBranch        *b_blob_vtx_E_mu;   //!
TBranch        *b_blob_vtx_E_nucl;   //!
TBranch        *b_blob_vtx_E_nucl_em;   //!
TBranch        *b_blob_vtx_E_nucl_highn;   //!
TBranch        *b_blob_vtx_E_nucl_lown;   //!
TBranch        *b_blob_vtx_E_nucl_meson;   //!
TBranch        *b_blob_vtx_E_nucl_midn;   //!
TBranch        *b_blob_vtx_E_nucl_mu;   //!
TBranch        *b_blob_vtx_E_nucl_other;   //!
TBranch        *b_blob_vtx_E_nucl_p;   //!
TBranch        *b_blob_vtx_E_nucl_xtalk;   //!
TBranch        *b_blob_vtx_E_od;   //!
TBranch        *b_blob_vtx_E_od_em;   //!
TBranch        *b_blob_vtx_E_od_highn;   //!
TBranch        *b_blob_vtx_E_od_lown;   //!
TBranch        *b_blob_vtx_E_od_meson;   //!
TBranch        *b_blob_vtx_E_od_midn;   //!
TBranch        *b_blob_vtx_E_od_mu;   //!
TBranch        *b_blob_vtx_E_od_other;   //!
TBranch        *b_blob_vtx_E_od_p;   //!
TBranch        *b_blob_vtx_E_od_xtalk;   //!
TBranch        *b_blob_vtx_E_other;   //!
TBranch        *b_blob_vtx_E_p;   //!
TBranch        *b_blob_vtx_E_tracker;   //!
TBranch        *b_blob_vtx_E_tracker_em;   //!
TBranch        *b_blob_vtx_E_tracker_highn;   //!
TBranch        *b_blob_vtx_E_tracker_lown;   //!
TBranch        *b_blob_vtx_E_tracker_meson;   //!
TBranch        *b_blob_vtx_E_tracker_midn;   //!
TBranch        *b_blob_vtx_E_tracker_mu;   //!
TBranch        *b_blob_vtx_E_tracker_other;   //!
TBranch        *b_blob_vtx_E_tracker_p;   //!
TBranch        *b_blob_vtx_E_tracker_xtalk;   //!
TBranch        *b_blob_vtx_E_xtalk;   //!
TBranch        *b_energy_from_mc;   //!
TBranch        *b_energy_from_mc_fraction;   //!
TBranch        *b_energy_from_mc_fraction_of_highest;   //!
TBranch        *b_muon_phi;   //!
TBranch        *b_muon_theta;   //!
TBranch        *b_muon_thetaX;   //!
TBranch        *b_muon_thetaY;   //!
TBranch        *b_muon_truthMatch_fraction;   //!
TBranch        *b_numi_horn_curr;   //!
TBranch        *b_numi_pot;   //!
TBranch        *b_numi_x;   //!
TBranch        *b_numi_x_width;   //!
TBranch        *b_numi_y;   //!
TBranch        *b_numi_y_width;   //!
TBranch        *b_phys_energy_dispersed;   //!
TBranch        *b_phys_energy_in_road_downstream;   //!
TBranch        *b_phys_energy_in_road_upstream;   //!
TBranch        *b_phys_energy_unattached;   //!
TBranch        *b_prim_vtx_smallest_opening_angle;   //!
TBranch        *b_primary_track_minerva_energy;   //!
TBranch        *b_primary_track_minerva_phi;   //!
TBranch        *b_primary_track_minerva_theta;   //!
TBranch        *b_usact_avg_E_consecutive_huge;   //!
TBranch        *b_usact_avg_E_consecutive_inf;   //!
TBranch        *b_usact_avg_E_consecutive_large;   //!
TBranch        *b_usact_avg_E_consecutive_normal;   //!
TBranch        *b_usact_avg_E_consecutive_small;   //!
TBranch        *b_usact_avg_E_consecutive_texas_sized;   //!
TBranch        *b_usact_avg_E_consecutive_tiny;   //!
TBranch        *b_usact_avg_E_huge;   //!
TBranch        *b_usact_avg_E_inf;   //!
TBranch        *b_usact_avg_E_large;   //!
TBranch        *b_usact_avg_E_normal;   //!
TBranch        *b_usact_avg_E_small;   //!
TBranch        *b_usact_avg_E_texas_sized;   //!
TBranch        *b_usact_avg_E_tiny;   //!
TBranch        *b_usact_frac_withE_huge;   //!
TBranch        *b_usact_frac_withE_inf;   //!
TBranch        *b_usact_frac_withE_large;   //!
TBranch        *b_usact_frac_withE_normal;   //!
TBranch        *b_usact_frac_withE_small;   //!
TBranch        *b_usact_frac_withE_texas_sized;   //!
TBranch        *b_usact_frac_withE_tiny;   //!
TBranch        *b_vetoMuonTime;   //!
TBranch        *b_vtx_fit_chi2;   //!
TBranch        *b_n_vetoDigits;   //!
TBranch        *b_discrFired;   //!
TBranch        *b_has_michel_at_vertex_sz;   //!
TBranch        *b_has_michel_at_vertex;   //!
TBranch        *b_has_michel_beginModule_sz;   //!
TBranch        *b_has_michel_beginModule;   //!
TBranch        *b_has_michel_category_sz;   //!
TBranch        *b_has_michel_category;   //!
TBranch        *b_has_michel_endModule_sz;   //!
TBranch        *b_has_michel_endModule;   //!
TBranch        *b_has_michel_numDigits_sz;   //!
TBranch        *b_has_michel_numDigits;   //!
TBranch        *b_has_michel_numModules_sz;   //!
TBranch        *b_has_michel_numModules;   //!
TBranch        *b_has_michel_numPlanes_sz;   //!
TBranch        *b_has_michel_numPlanes;   //!
TBranch        *b_has_michel_numTracks_sz;   //!
TBranch        *b_has_michel_numTracks;   //!
TBranch        *b_has_michel_slice_number_sz;   //!
TBranch        *b_has_michel_slice_number;   //!
TBranch        *b_has_michel_vertex_type_sz;   //!
TBranch        *b_has_michel_vertex_type;   //!
TBranch        *b_has_michel_view_sum_sz;   //!
TBranch        *b_has_michel_view_sum;   //!
TBranch        *b_paddle;   //!
TBranch        *b_pmt;   //!
TBranch        *b_pmt_occupancy;   //!
TBranch        *b_wall;   //!
TBranch        *b_has_michel_distance_sz;   //!
TBranch        *b_has_michel_distance;   //!
TBranch        *b_has_michel_energy_sz;   //!
TBranch        *b_has_michel_energy;   //!
TBranch        *b_has_michel_slice_energy_sz;   //!
TBranch        *b_has_michel_slice_energy;   //!
TBranch        *b_has_michel_time_diff_sz;   //!
TBranch        *b_has_michel_time_diff;   //!
TBranch        *b_muon_corrected_p;   //!
TBranch        *b_muon_sp;   //!
TBranch        *b_muon_sp_orig;   //!
TBranch        *b_orig_short_vtx;   //!
TBranch        *b_orig_vtx;   //!
TBranch        *b_primary_track_minerva_end_position;   //!
TBranch        *b_primary_track_minerva_start_position;   //!
TBranch        *b_time_distribution;   //!
TBranch        *b_truth_has_michel_from_pion_minus_momentum_sz;   //!
TBranch        *b_truth_has_michel_from_pion_minus_momentum;   //!
TBranch        *b_truth_has_michel_from_pion_plus_momentum_sz;   //!
TBranch        *b_truth_has_michel_from_pion_plus_momentum;   //!
TBranch        *b_usact_E_per_plane_huge;   //!
TBranch        *b_usact_E_per_plane_inf;   //!
TBranch        *b_usact_E_per_plane_large;   //!
TBranch        *b_usact_E_per_plane_normal;   //!
TBranch        *b_usact_E_per_plane_small;   //!
TBranch        *b_usact_E_per_plane_texas_sized;   //!
TBranch        *b_usact_E_per_plane_tiny;   //!
TBranch        *b_vtx_blob_radius;   //!
TBranch        *b_truth_has_physics_event;   //!
TBranch        *b_truth_muon_is_plausible;   //!
TBranch        *b_truth_reco_has_int_vtx;   //!
TBranch        *b_truth_reco_has_bad_object;   //!
TBranch        *b_truth_reco_has_muon;   //!
TBranch        *b_truth_reco_muon_has_charge;   //!
TBranch        *b_truth_reco_has_good_vtx;   //!
TBranch        *b_truth_reco_is_rock_muon;   //!
TBranch        *b_truth_pass_MasterAnaDev;   //!
TBranch        *b_truth_in_analyzable_area;   //!
TBranch        *b_truth_in_fiducial_area;   //!
TBranch        *b_truth_muon_leaving_code;   //!
TBranch        *b_truth_muon_track_id;   //!
TBranch        *b_truth_reco_has_michel_electron;   //!
TBranch        *b_truth_targetID;   //!
TBranch        *b_truth_targetZ;   //!
TBranch        *b_truth_target_code;   //!
TBranch        *b_truth_vtx_module;   //!
TBranch        *b_truth_vtx_plane;   //!
TBranch        *b_truth_fs_energy_em;   //!
TBranch        *b_truth_fs_energy_highn;   //!
TBranch        *b_truth_fs_energy_hmidn;   //!
TBranch        *b_truth_fs_energy_lown;   //!
TBranch        *b_truth_fs_energy_meson;   //!
TBranch        *b_truth_fs_energy_mu;   //!
TBranch        *b_truth_fs_energy_o;   //!
TBranch        *b_truth_fs_energy_p;   //!
TBranch        *b_truth_muon_phi;   //!
TBranch        *b_truth_muon_theta;   //!
TBranch        *b_truth_muon_thetaX;   //!
TBranch        *b_truth_muon_thetaY;   //!
TBranch        *b_truth_target_dist_to_division;   //!
TBranch        *b_truth_target_zDist;   //!
TBranch        *b_truth_ref_targZ;   //!
TBranch        *b_genie_wgt_n_shifts;   //!
TBranch        *b_truth_genie_wgt_AGKYxF1pi;   //!
TBranch        *b_truth_genie_wgt_AhtBY;   //!
TBranch        *b_truth_genie_wgt_BhtBY;   //!
TBranch        *b_truth_genie_wgt_CCQEPauliSupViaKF;   //!
TBranch        *b_truth_genie_wgt_CV1uBY;   //!
TBranch        *b_truth_genie_wgt_CV2uBY;   //!
TBranch        *b_truth_genie_wgt_EtaNCEL;   //!
TBranch        *b_truth_genie_wgt_FrAbs_N;   //!
TBranch        *b_truth_genie_wgt_FrAbs_pi;   //!
TBranch        *b_truth_genie_wgt_FrCEx_N;   //!
TBranch        *b_truth_genie_wgt_FrCEx_pi;   //!
TBranch        *b_truth_genie_wgt_FrElas_N;   //!
TBranch        *b_truth_genie_wgt_FrElas_pi;   //!
TBranch        *b_truth_genie_wgt_FrInel_N;   //!
TBranch        *b_truth_genie_wgt_FrInel_pi;   //!
TBranch        *b_truth_genie_wgt_FrPiProd_N;   //!
TBranch        *b_truth_genie_wgt_FrPiProd_pi;   //!
TBranch        *b_truth_genie_wgt_MFP_N;   //!
TBranch        *b_truth_genie_wgt_MFP_pi;   //!
TBranch        *b_truth_genie_wgt_MaCCQE;   //!
TBranch        *b_truth_genie_wgt_MaCCQEshape;   //!
TBranch        *b_truth_genie_wgt_MaNCEL;   //!
TBranch        *b_truth_genie_wgt_MaRES;   //!
TBranch        *b_truth_genie_wgt_MvRES;   //!
TBranch        *b_truth_genie_wgt_NormCCQE;   //!
TBranch        *b_truth_genie_wgt_NormCCRES;   //!
TBranch        *b_truth_genie_wgt_NormDISCC;   //!
TBranch        *b_truth_genie_wgt_NormNCRES;   //!
TBranch        *b_truth_genie_wgt_RDecBR1gamma;   //!
TBranch        *b_truth_genie_wgt_Rvn1pi;   //!
TBranch        *b_truth_genie_wgt_Rvn2pi;   //!
TBranch        *b_truth_genie_wgt_Rvp1pi;   //!
TBranch        *b_truth_genie_wgt_Rvp2pi;   //!
TBranch        *b_truth_genie_wgt_Theta_Delta2Npi;   //!
TBranch        *b_truth_genie_wgt_VecFFCCQEshape;   //!
TBranch        *b_truth_genie_wgt_shifts;   //!
TBranch        *b_truth_muon_end_momentum;   //!
TBranch        *b_truth_muon_end_position;   //!
TBranch        *b_truth_ref_dist_to_division;   //!
TBranch        *b_truth_ref_dist_to_target;   //!
TBranch        *b_MasterAnaDev_nuFlavor;   //!
TBranch        *b_MasterAnaDev_nuHelicity;   //!
TBranch        *b_MasterAnaDev_intCurrent;   //!
TBranch        *b_MasterAnaDev_intType;   //!
TBranch        *b_MasterAnaDev_E;   //!
TBranch        *b_MasterAnaDev_Q2;   //!
TBranch        *b_MasterAnaDev_x;   //!
TBranch        *b_MasterAnaDev_y;   //!
TBranch        *b_MasterAnaDev_W;   //!
TBranch        *b_MasterAnaDev_score;   //!
TBranch        *b_MasterAnaDev_leptonE;   //!
TBranch        *b_MasterAnaDev_vtx;   //!
TBranch        *b_MasterAnaDev_minos_trk_is_contained;   //!
TBranch        *b_MasterAnaDev_minos_trk_is_ok;   //!
TBranch        *b_MasterAnaDev_minos_used_range;   //!
TBranch        *b_MasterAnaDev_minos_used_curvature;   //!
TBranch        *b_MasterAnaDev_pass_canonical_cut;   //!
TBranch        *b_MasterAnaDev_is_cc;   //!
TBranch        *b_MasterAnaDev_in_analyzable_area;   //!
TBranch        *b_MasterAnaDev_in_fiducial_area;   //!
TBranch        *b_MasterAnaDev_minos_trk_end_plane;   //!
TBranch        *b_MasterAnaDev_minos_trk_quality;   //!
TBranch        *b_MasterAnaDev_r_minos_trk_vtx_plane;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_numFSMuons;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLeptonPDG;   //!
TBranch        *b_MasterAnaDev_targetID;   //!
TBranch        *b_MasterAnaDev_targetZ;   //!
TBranch        *b_MasterAnaDev_target_code;   //!
TBranch        *b_MasterAnaDev_vtx_module;   //!
TBranch        *b_MasterAnaDev_vtx_plane;   //!
TBranch        *b_MasterAnaDev_E_ccqe;   //!
TBranch        *b_MasterAnaDev_E_wide_window;   //!
TBranch        *b_MasterAnaDev_Q2_ccqe;   //!
TBranch        *b_MasterAnaDev_Q2_wide_window;   //!
TBranch        *b_MasterAnaDev_W_wide_window;   //!
TBranch        *b_MasterAnaDev_minos_trk_bave;   //!
TBranch        *b_MasterAnaDev_minos_trk_chi2;   //!
TBranch        *b_MasterAnaDev_minos_trk_end_u;   //!
TBranch        *b_MasterAnaDev_minos_trk_end_v;   //!
TBranch        *b_MasterAnaDev_minos_trk_end_x;   //!
TBranch        *b_MasterAnaDev_minos_trk_end_y;   //!
TBranch        *b_MasterAnaDev_minos_trk_end_z;   //!
TBranch        *b_MasterAnaDev_minos_trk_eqp;   //!
TBranch        *b_MasterAnaDev_minos_trk_eqp_qp;   //!
TBranch        *b_MasterAnaDev_minos_trk_fit_pass;   //!
TBranch        *b_MasterAnaDev_minos_trk_ndf;   //!
TBranch        *b_MasterAnaDev_minos_trk_p;   //!
TBranch        *b_MasterAnaDev_minos_trk_p_curvature;   //!
TBranch        *b_MasterAnaDev_minos_trk_p_range;   //!
TBranch        *b_MasterAnaDev_minos_trk_qp;   //!
TBranch        *b_MasterAnaDev_minos_trk_vtx_x;   //!
TBranch        *b_MasterAnaDev_minos_trk_vtx_y;   //!
TBranch        *b_MasterAnaDev_minos_trk_vtx_z;   //!
TBranch        *b_MasterAnaDev_nu_energy_recoil;   //!
TBranch        *b_MasterAnaDev_r_minos_trk_bdL;   //!
TBranch        *b_MasterAnaDev_r_minos_trk_end_dcosx;   //!
TBranch        *b_MasterAnaDev_r_minos_trk_end_dcosy;   //!
TBranch        *b_MasterAnaDev_r_minos_trk_end_dcosz;   //!
TBranch        *b_MasterAnaDev_r_minos_trk_vtx_dcosx;   //!
TBranch        *b_MasterAnaDev_r_minos_trk_vtx_dcosy;   //!
TBranch        *b_MasterAnaDev_r_minos_trk_vtx_dcosz;   //!
TBranch        *b_MasterAnaDev_recoil_E;   //!
TBranch        *b_MasterAnaDev_recoil_E_wide_window;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMinosInitProjPx;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMinosInitProjPy;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMinosInitProjPz;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMinosInitProjX;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMinosInitProjY;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMinosInitProjZ;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMnvFinalPx;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMnvFinalPy;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMnvFinalPz;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMnvFinalX;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMnvFinalY;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMnvFinalZ;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMnvInitPx;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMnvInitPy;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMnvInitPz;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMnvInitX;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMnvInitY;   //!
TBranch        *b_MasterAnaDev_t_minos_trk_primFSLepMnvInitZ;   //!
TBranch        *b_MasterAnaDev_target_dist_to_division;   //!
TBranch        *b_MasterAnaDev_target_zDist;   //!
TBranch        *b_MasterAnaDev_x_wide_window;   //!
TBranch        *b_MasterAnaDev_y_wide_window;   //!
TBranch        *b_MasterAnaDev_ref_targZ;   //!
TBranch        *b_MasterAnaDev_smeared_is_fiducial;   //!
TBranch        *b_MasterAnaDev_smeared_pass_dist_to_div;   //!
TBranch        *b_MasterAnaDev_smeared_targetID;   //!
TBranch        *b_MasterAnaDev_smeared_targetZ;   //!
TBranch        *b_MasterAnaDev_smeared_vtx_mod;   //!
TBranch        *b_MasterAnaDev_muon_vtx;   //!
TBranch        *b_MasterAnaDev_ref_dist_to_division;   //!
TBranch        *b_MasterAnaDev_ref_dist_to_target;   //!
TBranch        *b_MasterAnaDev_sys_muon_curve_energy_shift;   //!
TBranch        *b_MasterAnaDev_sys_muon_energy_shift;   //!
TBranch        *b_MasterAnaDev_sys_muon_minerva_energy_shift;   //!
TBranch        *b_MasterAnaDev_sys_muon_qSquared_shift;   //!
TBranch        *b_MasterAnaDev_sys_muon_range_energy_shift;   //!
TBranch        *b_MasterAnaDev_sys_muon_wSquared_shift;   //!
TBranch        *b_MasterAnaDev_sys_muon_xbj_shift;   //!
TBranch        *b_MasterAnaDev_sys_muon_y_shift;   //!
TBranch        *b_MasterAnaDev_sys_nu_energy_shift;   //!
TBranch        *b_MasterAnaDev_sys_recoil_energy_shift;   //!
TBranch        *b_MasterAnaDev_sys_recoil_qSquared_shift;   //!
TBranch        *b_MasterAnaDev_sys_recoil_wSquared_shift;   //!
TBranch        *b_MasterAnaDev_sys_recoil_xbj_shift;   //!
TBranch        *b_MasterAnaDev_sys_recoil_y_shift;   //!
TBranch        *b_MasterAnaDev_sys_total_qSquared_shift;   //!
TBranch        *b_MasterAnaDev_sys_total_wSquared_shift;   //!
TBranch        *b_MasterAnaDev_sys_total_xbj_shift;   //!
TBranch        *b_MasterAnaDev_sys_total_y_shift;   //!
TBranch        *b_ev_run;   //!
TBranch        *b_ev_subrun;   //!
TBranch        *b_ev_detector;   //!
TBranch        *b_ev_triggerType;   //!
TBranch        *b_ev_gate;   //!
TBranch        *b_ev_global_gate;   //!
TBranch        *b_ev_gps_time_sec;   //!
TBranch        *b_ev_gps_time_usec;   //!
TBranch        *b_mc_run;   //!
TBranch        *b_mc_subrun;   //!
TBranch        *b_mc_nInteractions;   //!
TBranch        *b_mc_MIState;   //!
TBranch        *b_mc_pot;   //!
TBranch        *b_mc_beamConfig;   //!
TBranch        *b_mc_processType;   //!
TBranch        *b_mc_nthEvtInSpill;   //!
TBranch        *b_mc_nthEvtInFile;   //!
TBranch        *b_mc_intType;   //!
TBranch        *b_mc_current;   //!
TBranch        *b_mc_charm;   //!
TBranch        *b_mc_weight;   //!
TBranch        *b_mc_XSec;   //!
TBranch        *b_mc_diffXSec;   //!
TBranch        *b_mc_incoming;   //!
TBranch        *b_mc_fluxDriverProb;   //!
TBranch        *b_mc_targetNucleus;   //!
TBranch        *b_mc_targetZ;   //!
TBranch        *b_mc_targetA;   //!
TBranch        *b_mc_targetNucleon;   //!
TBranch        *b_mc_struckQuark;   //!
TBranch        *b_mc_seaQuark;   //!
TBranch        *b_mc_resID;   //!
TBranch        *b_mc_primaryLepton;   //!
TBranch        *b_mc_incomingE;   //!
TBranch        *b_mc_Bjorkenx;   //!
TBranch        *b_mc_Bjorkeny;   //!
TBranch        *b_mc_Q2;   //!
TBranch        *b_mc_nuT;   //!
TBranch        *b_mc_w;   //!
TBranch        *b_mc_vtx;   //!
TBranch        *b_mc_incomingPartVec;   //!
TBranch        *b_mc_initNucVec;   //!
TBranch        *b_mc_primFSLepton;   //!
TBranch        *b_mc_nFSPart;   //!
TBranch        *b_mc_FSPartPx;   //!
TBranch        *b_mc_FSPartPy;   //!
TBranch        *b_mc_FSPartPz;   //!
TBranch        *b_mc_FSPartE;   //!
TBranch        *b_mc_FSPartPDG;   //!
TBranch        *b_mc_er_nPart;   //!
TBranch        *b_mc_er_ID;   //!
TBranch        *b_mc_er_status;   //!
TBranch        *b_mc_er_posInNucX;   //!
TBranch        *b_mc_er_posInNucY;   //!
TBranch        *b_mc_er_posInNucZ;   //!
TBranch        *b_mc_er_Px;   //!
TBranch        *b_mc_er_Py;   //!
TBranch        *b_mc_er_Pz;   //!
TBranch        *b_mc_er_E;   //!
TBranch        *b_mc_er_FD;   //!
TBranch        *b_mc_er_LD;   //!
TBranch        *b_mc_er_mother;   //!
TBranch        *b_mc_fr_nNuAncestorIDs;   //!
TBranch        *b_mc_fr_nuAncestorIDs;   //!
TBranch        *b_mc_fr_nuParentID;   //!
TBranch        *b_mc_fr_decMode;   //!
TBranch        *b_mc_fr_primProtonVtx;   //!
TBranch        *b_mc_fr_primProtonP;   //!
TBranch        *b_mc_fr_nuParentDecVtx;   //!
TBranch        *b_mc_fr_nuParentProdVtx;   //!
TBranch        *b_mc_fr_nuParentProdP;   //!
TBranch        *b_mc_cvweight_total;   //!
TBranch        *b_wgt;   //!
TBranch        *b_mc_cvweight_totalFlux;   //!
TBranch        *b_mc_cvweight_totalXsec;   //!
TBranch        *b_mc_ppfx1_cvweight;   //!
TBranch        *b_mc_hornCurrent_cvweight;   //!
TBranch        *b_mc_gen1_cvweight_total;   //!
TBranch        *b_gen1_wgt;   //!
TBranch        *b_mc_gen1_cvweight_totalFlux;   //!
TBranch        *b_mc_gen1_cvweight_NA49;   //!
TBranch        *b_mc_wgt_Flux_BeamFocus_sz;   //!
TBranch        *b_mc_wgt_Flux_BeamFocus;   //!
TBranch        *b_mc_wgt_gen1_Flux_Tertiary_sz;   //!
TBranch        *b_mc_wgt_gen1_Flux_Tertiary;   //!
TBranch        *b_mc_wgt_gen1_Flux_NA49_sz;   //!
TBranch        *b_mc_wgt_gen1_Flux_NA49;   //!
TBranch        *b_mc_wgt_Norm_sz;   //!
TBranch        *b_mc_wgt_Norm;   //!
TBranch        *b_mc_wgt_ppfx1_Total_sz;   //!
TBranch        *b_mc_wgt_ppfx1_Total;   //!
TBranch        *b_n_prongs;   //!
TBranch        *b_prong_nParticles;   //!
TBranch        *b_prong_part_score;   //!
TBranch        *b_prong_part_mass;   //!
TBranch        *b_prong_part_charge;   //!
TBranch        *b_prong_part_pid;   //!
TBranch        *b_prong_part_E;   //!
TBranch        *b_prong_part_pos;   //!

//fChain->SetBranchAddress("MasterAnaDev_leptonE", MasterAnaDev_leptonE, &b_MasterAnaDev_leptonE);
#endif // MasterAnaDevvars_h
