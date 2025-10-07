/*
 * LBFGSTest.cpp
 *
 *  Created on: 30 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "DimerTest.h"

#include "../../../gpr/Enums.h"

namespace gpr {
namespace tests {

DimerTest::DimerTest() : dimer(0.01)
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.
}

DimerTest::~DimerTest() { }

// TEST_F(DimerTest, RotateForce)
//{
//    gpr::Coord G01;
//    gpr::Coord orient;
//    gpr::Coord F_rot;
//    std::vector<double> F_rot_ref(6);
//    double dimer_sep = 0.01;
//
//    G01.resize(2, 6);
//    orient.resize(1, 6);
//
//    G01(0, 0) = 0.0222174753484639;
//    G01(0, 1) = -0.0222402530044595;
//    G01(0, 2) = 0.0185353722102027;
//    G01(0, 3) = -0.0173837548235882;
//    G01(0, 4) = 0.0219137885404279;
//    G01(0, 5) = -0.0258161071576406;
//    G01(1, 0) = 0.0214452456338877;
//    G01(1, 1) = -0.0331231900862804;
//    G01(1, 2) = 0.0283157675762483;
//    G01(1, 3) = -0.0158030003032378;
//    G01(1, 4) = 0.0208628464314946;
//    G01(1, 5) = -0.0383251019392285;
//
//    orient(0, 0) = -0.3810391799989235;
//    orient(0, 1) = -0.5114555099749947;
//    orient(0, 2) = 0.4411290293765586;
//    orient(0, 3) = -0.3595730639399913;
//    orient(0, 4) = -0.0734979445112717;
//    orient(0, 5) = -0.5137439516964519;
//
//    F_rot_ref = {-1.0732676278061717, 0.5286706052417522, -0.5347551668163790,
//                 -1.4747003484244667, -0.0266229840558352,
//                 0.8465087539411261};
//
//    dimer.rotateForce(G01, orient, dimer_sep, F_rot);
//
//    for (gpr::Index_t n = 0; n < F_rot.getNumCols(); ++n) {
//        EXPECT_LE(std::fabs(F_rot(0, n) - F_rot_ref[n]), threshold)
//            << "Elements of the field are not equal to the expected ones.";
//    }
//}

// TEST_F(DimerTest, RotateDimerLow)
//{
//    gpr::Coord F_rot;
//    gpr::Coord G01;
//    gpr::Coord orient;
//    double dimer_sep;
//
//    gpr::Coord R;
//    uint8_t potential = POTENTIAL_TEST2;
//    double T_anglerot = 0.;
//    double estim_Curv = 1;
//    double estim_G1 = 1;
//    gpr::Coord orient_new;
//    gpr::Coord orient_rot_new;
//    double Curv;
//    gpr::Coord G1;
//    gpr::Observation dtheta_obs;
//    gpr::Observation all_obs;
//
//    gpr::Coord orient_new_ref;
//    gpr::Coord orient_rot_new_ref;
//    gpr::Coord G1_ref;
//    gpr::Coord R1_dtheta_ref;
//    gpr::Coord G1_dtheta_ref;
//    double Curv_ref;
//    double E1_dtheta_ref;
//
//    R.resize(1, 2 * 3);
//    F_rot.resize(1, 2 * 3);
//    G01.resize(2, 2 * 3);
//    orient.resize(1, 2 * 3);
//
//    orient_new_ref.resize(1, 2 * 3);
//    orient_rot_new_ref.resize(1, 2 * 3);
//    G1_ref.resize(1, 2 * 3);
//    R1_dtheta_ref.resize(1, 2 * 3);
//    G1_dtheta_ref.resize(1, 2 * 3);
//
//    F_rot(0, 0) = -1.07243573302179e+000;
//    F_rot(0, 1) = 528.323631218742e-003;
//    F_rot(0, 2) = -534.186955024052e-003;
//    F_rot(0, 3) = -1.47328866808001e+000;
//    F_rot(0, 4) = -26.5981132823147e-003;
//    F_rot(0, 5) = 845.733468528328e-003;
//
//    G01(0, 0) = 22.2187412441741e-003;
//    G01(0, 1) = -22.2419538432585e-003;
//    G01(0, 2) = 18.5375109120299e-003;
//    G01(0, 3) = -17.3807794897030e-003;
//    G01(0, 4) = 21.9136490152291e-003;
//    G01(0, 5) = -25.8188654019664e-003;
//    G01(1, 0) = 21.4437403317027e-003;
//    G01(1, 1) = -33.1212926208074e-003;
//    G01(1, 2) = 28.3134580120646e-003;
//    G01(1, 3) = -15.8057733045760e-003;
//    G01(1, 4) = 20.8628503344299e-003;
//    G01(1, 5) = -38.3221119846374e-003;
//
//    dimer_sep = 10.0000000000000e-003;
//
//    orient(0, 0) = -381.039179998923e-003;
//    orient(0, 1) = -511.455509974995e-003;
//    orient(0, 2) = 441.129029376559e-003;
//    orient(0, 3) = -359.573063939991e-003;
//    orient(0, 4) = -73.4979445112717e-003;
//    orient(0, 5) = -513.743951696452e-003;
//
//    T_anglerot = 8.73000000000000e-003;
//
//    R(0, 0) = 8.98237316483057e+000;
//    R(0, 1) = 9.93723083577204e+000;
//    R(0, 2) = 7.89441632385049e+000;
//    R(0, 3) = 7.65248322727496e+000;
//    R(0, 4) = 9.95590549457398e+000;
//    R(0, 5) = 7.87787958998366e+000;
//
//    orient_new_ref(0, 0) = -627.706329884103e-003;
//    orient_new_ref(0, 1) = -84.1127344208758e-003;
//    orient_new_ref(0, 2) = 42.2880054326739e-003;
//    orient_new_ref(0, 3) = -770.126393092029e-003;
//    orient_new_ref(0, 4) = -51.6005574137432e-003;
//    orient_new_ref(0, 5) = 36.9358506853572e-003;
//
//    orient_rot_new_ref(0, 0) = 33.7455015336214e-003;
//    orient_rot_new_ref(0, 1) = 561.405493795794e-003;
//    orient_rot_new_ref(0, 2) = -504.801804441458e-003;
//    orient_rot_new_ref(0, 3) = -89.1352421541912e-003;
//    orient_rot_new_ref(0, 4) = 53.7874277771169e-003;
//    orient_rot_new_ref(0, 5) = 646.546264480762e-003;
//
//    Curv_ref = 34.5579478503286e-003;
//
//    G1_ref(0, 0) = 21.6381582772695e-003;
//    G1_ref(0, 1) = -22.1549166895741e-003;
//    G1_ref(0, 2) = 18.8649584801689e-003;
//    G1_ref(0, 3) = -17.6581772820319e-003;
//    G1_ref(0, 4) = 21.8289359490012e-003;
//    G1_ref(0, 5) = -25.9382925222035e-003;
//
//    R1_dtheta_ref(0, 0) = 8.97727844229888e+000;
//    R1_dtheta_ref(0, 1) = 9.93304844889349e+000;
//    R1_dtheta_ref(0, 2) = 7.89791765289864e+000;
//    R1_dtheta_ref(0, 3) = 7.64705291657440e+000;
//    R1_dtheta_ref(0, 4) = 9.95516609345974e+000;
//    R1_dtheta_ref(0, 5) = 7.87410372598883e+000;
//
//    E1_dtheta_ref = 268.422097239501e-006;
//
//    G1_dtheta_ref(0, 0) = 21.4263025463240e-003;
//    G1_dtheta_ref(0, 1) = -30.4813824740716e-003;
//    G1_dtheta_ref(0, 2) = 26.0833530773849e-003;
//    G1_dtheta_ref(0, 3) = -16.2806875633080e-003;
//    G1_dtheta_ref(0, 4) = 21.0852193621926e-003;
//    G1_dtheta_ref(0, 5) = -35.3650094447835e-003;
//
////    dimer.optimizeRotationAngle(R, orient, G01, F_rot, potential, dimer_sep,
////                          T_anglerot, estim_Curv, estim_G1, nullptr,
////                          orient_new, orient_rot_new, Curv, G1, dtheta_obs);
//
////    EXPECT_LE(std::fabs(Curv - Curv_ref), 10 * threshold)
////        << "Curvature is not equal to the expected one.";
////
////    EXPECT_LE(std::fabs(E1_dtheta_ref - E1_dtheta_ref), 10 * threshold)
////        << "Energy is not equal to the expected one.";
////
////    for (gpr::Index_t n = 0; n < orient_new_ref.getSize(); ++n) {
////        EXPECT_LE(std::fabs(orient_new[n] - orient_new_ref[n]), threshold)
////            << "Orientation field is not equal to the expected one." << n;
////    }
////    for (gpr::Index_t n = 0; n < orient_rot_new_ref.getSize(); ++n) {
////        EXPECT_LE(std::fabs(orient_rot_new[n] - orient_rot_new_ref[n]),
////                  threshold)
////            << "Rotated orientation field is not equal to the expected one."
////            << n;
////    }
////    for (gpr::Index_t n = 0; n < G1_ref.getSize(); ++n) {
////        EXPECT_LE(std::fabs(G1[n] - G1_ref[n]), threshold)
////            << "Gradient field is not equal to the expected one." << n;
////    }
////    for (gpr::Index_t n = 0; n < R1_dtheta_ref.getSize(); ++n) {
////        EXPECT_LE(std::fabs(dtheta_obs.R[n] - R1_dtheta_ref[n]), threshold)
////            << "Position field is not equal to the expected one." << n;
////    }
////    for (gpr::Index_t n = 0; n < G1_dtheta_ref.getSize(); ++n) {
////        EXPECT_LE(std::fabs(dtheta_obs.G[n] - G1_dtheta_ref[n]), threshold)
////            << "Rotated gradient field is not equal to the expected one." <<
/// n; /    }
//}

// TEST_F(DimerTest, RotateDimer)
//{
//    gpr::Coord R;
//    gpr::Coord orient;
//    gpr::Coord G01;
//    uint8_t potential = POTENTIAL_TEST1;
//    double dimer_sep;
//    double T_anglerot;
//    double estim_Curv;
//    gpr::Rotinfo rotinfo;
//    gpr::Coord orient_new;
//    double Curv;
//    gpr::Observation observation;
//    gpr::Observation all_obs;
//
//    R.resize(1, 2 * 3);
//    orient.resize(1, 2 * 3);
//    G01.resize(2, 2 * 3);
//    rotinfo.F_old.resize(1, 2 * 3);
////    rotinfo.deltaOrient.resize(1, 2 * 3);
//
//    R(0, 0) = 8.98237316483057e+000;
//    R(0, 1) = 9.93723083577204e+000;
//    R(0, 2) = 7.89441632385049e+000;
//    R(0, 3) = 7.65248322727496e+000;
//    R(0, 4) = 9.95590549457398e+000;
//    R(0, 5) = 7.87787958998366e+000;
//
//    orient(0, 0) = -627.706329884103e-003;
//    orient(0, 1) = -84.1127344208758e-003;
//    orient(0, 2) = 42.2880054326739e-003;
//    orient(0, 3) = -770.126393092029e-003;
//    orient(0, 4) = -51.6005574137432e-003;
//    orient(0, 5) = 36.9358506853572e-003;
//
//    G01(0, 0) = 22.2187412441741e-003;
//    G01(0, 1) = -22.2419538432585e-003;
//    G01(0, 2) = 18.5375109120299e-003;
//    G01(0, 3) = -17.3807794897030e-003;
//    G01(0, 4) = 21.9136490152291e-003;
//    G01(0, 5) = -25.8188654019664e-003;
//    G01(1, 0) = 21.4725532869451e-003;
//    G01(1, 1) = -21.9070976375274e-003;
//    G01(1, 2) = 18.6483116322486e-003;
//    G01(1, 3) = -17.7756480204136e-003;
//    G01(1, 4) = 21.7629297991571e-003;
//    G01(1, 5) = -25.6590672302472e-003;
//
//    dimer_sep = 10.0000000000000e-003;
//    T_anglerot = 8.73000000000000e-003;
//    estim_Curv = 0.;
//
//    rotinfo.F_old(0, 0) = -1.07243573302179e+000;
//    rotinfo.F_old(0, 1) = 528.323631218742e-003;
//    rotinfo.F_old(0, 2) = -534.186955024051e-003;
//    rotinfo.F_old(0, 3) = -1.47328866808001e+000;
//    rotinfo.F_old(0, 4) = -26.5981132823147e-003;
//    rotinfo.F_old(0, 5) = 845.733468528328e-003;
//
////    rotinfo.deltaOrient(0, 0) = -246.667149885180e-003;
////    rotinfo.deltaOrient(0, 1) = 427.342775554119e-003;
////    rotinfo.deltaOrient(0, 2) = -398.841023943885e-003;
////    rotinfo.deltaOrient(0, 3) = -410.553329152038e-003;
////    rotinfo.deltaOrient(0, 4) = 21.8973870975285e-003;
////    rotinfo.deltaOrient(0, 5) = 550.679802381809e-003;
//
//    rotinfo.num_lbfgsiter_rot = 6;
//
//    double Curv_ref = 0.;
//    double E_obs_ref = 140.981934401792e-006;
//
//    std::vector<double> R_obs_ref(6);
//    R_obs_ref = {8.97759283794573e+000, 9.93411441890023e+000,
//                 7.89488070133982e+000, 7.64432340054361e+000,
//                 9.95588249557613e+000, 7.87708117088127e+000};
//
//    std::vector<double> G_obs_ref(6);
//    G_obs_ref = {22.0169292639464e-003, -25.3799869163130e-003,
//                 20.4383773983493e-003, -17.5501833843735e-003,
//                 21.9491835909580e-003, -28.6553826196841e-003};
//
//    std::vector<double> rotinfo_F_rot_old_ref(6);
//    rotinfo_F_rot_old_ref = {53.4891599595057e-003,  -79.8015441165797e-003,
//                             -15.7096591623766e-003, -38.4990524938851e-003,
//                             22.2728497243401e-003,  -26.3255510147931e-003};
//
//    std::vector<double> orient_new_ref(6);
//    orient_new_ref = {-610.333333531775e-003, -116.295796532668e-003,
//                      43.1015651575280e-003,  -780.800289227480e-003,
//                      -45.0154141907832e-003, 20.8634991626622e-003};
//
//    std::vector<double> rotinfo_deltaR_mem_ref(12);
//    rotinfo_deltaR_mem_ref = {
//        -246.667149885180e-003, 427.342775554119e-003, -398.841023943885e-003,
//        -410.553329152038e-003, 21.8973870975285e-003,  550.679802381809e-003,
//        17.3729963523277e-003,  -32.1830621117919e-003, 813.559724854092e-006,
//        -10.6738961354504e-003, 6.58514322295999e-003,
//        -16.0723515226950e-003};
//
//    std::vector<double> rotinfo_deltaF_mem_ref(6);
//    rotinfo_deltaF_mem_ref = {1.12592489298129e+000, -608.125175335322e-003,
//                              518.477295861675e-003, 1.43478961558612e+000,
//                              48.8709630066548e-003, -872.059019543121e-003};
//
//    dimer.rotate(R, orient, G01, potential, dimer_sep, T_anglerot,
//                      estim_Curv, all_obs, nullptr, rotinfo, orient_new, Curv,
//                      observation);
//
//    EXPECT_LE(std::fabs(Curv - Curv_ref), 10 * threshold)
//        << "Curvature is not equal to the expected one.";
//
//    EXPECT_LE(std::fabs(observation.E(0, 0) - E_obs_ref), 10 * threshold)
//        << "Energy is not equal to the expected one.";
//
//    for (gpr::Index_t n = 0; n < orient_new.getSize(); ++n) {
//        EXPECT_LE(std::fabs(orient_new[n] - orient_new_ref[n]), threshold)
//            << "Orientation field is not equal to the expected one." << n;
//    }
//    for (gpr::Index_t n = 0; n < observation.R.getSize(); ++n) {
//        EXPECT_LE(std::fabs(observation.R[n] - R_obs_ref[n]), threshold)
//            << "Coordinates of the new observed location are not equal to the
//            "
//               "expected one."
//            << n;
//    }
//    for (gpr::Index_t n = 0; n < observation.G.getSize(); ++n) {
//        EXPECT_LE(std::fabs(observation.G[n] - G_obs_ref[n]), threshold)
//            << "Gradient field is not equal to the expected one." << n;
//    }
//    for (gpr::Index_t n = 0; n < rotinfo.F_old.getSize(); ++n) {
//        EXPECT_LE(std::fabs(rotinfo.F_old[n] - rotinfo_F_rot_old_ref[n]),
//                  threshold)
//            << "Rotational force field is not equal to the expected one." <<
//            n;
//    }
////    for (gpr::Index_t n = 0; n < rotinfo.deltaOrient.getSize(); ++n) {
////        EXPECT_LE(std::fabs(rotinfo.deltaOrient[n] -
/// rotinfo_deltaR_mem_ref[n]), /                  threshold) /            <<
///"Change of orientation field is not equal to the expected one." / << n; / }
////    for (gpr::Index_t n = 0; n < rotinfo.deltaF.getSize(); ++n) {
////        EXPECT_LE(std::fabs(rotinfo.deltaF[n] - rotinfo_deltaF_mem_ref[n]),
////                  threshold)
////            << "Change of rotational field is not equal to the expected
/// one." /            << n; /    }
//}

// TEST_F(DimerTest, TranslateDimer)
//{
//    gpr::Coord R;
//    gpr::Coord orient;
//    gpr::Coord F_R;
//    double curv;
//    gpr::TransInfo transinfo;
//    gpr::Coord R_new;
//    gpr::Coord R_new_ref;
//    gpr::TransitionParameters translation_param;
//
//    R_new_ref.resize(1, 3 * 2);
//    R.resize(1, 3 * 2);
//    orient.resize(1, 3 * 2);
//    F_R.resize(1, 3 * 2);
//    transinfo.F_old.resize(1, 3 * 2);
//    transinfo.deltaOrient.resize(1, 3 * 2);
//
//    R.set(0, 0, {8.982538896311098, 9.937443825779654, 7.894092264564180});
//    R.set(0, 1, {7.652281581062163, 9.955691736800297, 7.877995340784432});
//
//    orient.set(0, 0,
//               {0.674705565855243, -0.016101844764894, -0.241764934787344});
//    orient.set(0, 1,
//               {-0.651738746177718, 0.007998301613614, -0.247458037267712});
//
//    F_R.set(0, 0, {-0.021182098720882, 0.021747610563288,
//    -0.018160531471649}); F_R.set(0, 1, {0.016361342402227,
//    -0.021295086644017, 0.024996940898092});
//
//    curv = -5.156535441726833;
//
//    translation_param.step_length = 0.1;
//    translation_param.max_step_length = 0.1;
//
//    transinfo.F_old.set(
//        0, 0, {0.016573148052742, 0.021299000761832, -0.032405928631147});
//    transinfo.F_old.set(
//        0, 1, {-0.020164621279596, -0.021375777367842, 0.011575080077217});
//
//    transinfo.deltaOrient.set(
//        0, 0, {0.165731480527417, 0.212990007618319, -0.324059286311467});
//    transinfo.deltaOrient.set(
//        0, 1, {-0.201646212795956, -0.213757773678424, 0.115750800772167});
//    transinfo.deltaOrient *= 1e-3;
//
//    transinfo.num_lbfgsiter_trans = 6;
//
//    R_new_ref.set(0, 0,
//                  {8.984799577839603, 9.946773537381196, 7.883158760864689});
//    R_new_ref.set(0, 1,
//                  {7.648053020353811, 9.947018098632240, 7.884341447664946});
//
//    dimer.translate(R, orient, F_R, curv, translation_param, transinfo,
//                         R_new);
//
//    EXPECT_EQ(R_new.getSize(), 6)
//        << "Number of elements in the field is not equal to the expected
//        one.";
//
//    for (gpr::Index_t n = 0; n < R_new.getSize(); ++n) {
//        EXPECT_LE(fabs(R_new[n] - R_new_ref[n]), DBL_EPSILON * 1e3)
//            << "Elements of the field are not equal to the expected ones.";
//    }
//}

} /* namespace tests */
} /* namespace gpr */
