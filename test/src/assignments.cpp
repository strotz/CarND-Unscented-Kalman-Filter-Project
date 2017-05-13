#include <ukf_augmented.h>
#include <ukf_parts.h>
#include "gtest/gtest.h"

#include "Eigen/Dense"
#include "ukf_state.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

TEST(UKF, Eigen) {
  VectorXd a = VectorXd(1);
  a << 1;

  VectorXd b = VectorXd(1);
  b << 1;

  VectorXd c = VectorXd(1);
  c << 2;

  VectorXd d = VectorXd(2);
  d << 1, 1;

  EXPECT_EQ(a, b);
  EXPECT_NE(a, c);
  EXPECT_NE(a.size(), d.size());
}

TEST(UKF, Augmentation) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a = 0.2;

  //Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd = 0.2;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set example state
  State x = State();
  x <<   5.7441,
    1.3800,
    2.2049,
    0.5015,
    0.3528;

  //create example covariance matrix
  StateCovariance P = StateCovariance();
  P <<
    0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
    -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
    0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
    -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
    -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  //create augmented mean vector
  AugmentedState x_aug = AugmentedState(x);

  //create augmented state covariance
  AugmentedStateCovariance P_aug = AugmentedStateCovariance(P, std_a, std_yawdd);

  //create sigma point matrix
  AugmentedSpaceSigmaPoints Xsig_aug = AugmentedSpaceSigmaPoints(x_aug, P_aug, lambda);

  MatrixXd Xsig_aug_expected = MatrixXd(7, 15);
  Xsig_aug_expected <<
  5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,  5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
    1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,  1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
  2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,  2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
  0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,  0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
  0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528, 0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
       0,        0,        0,        0,        0,        0,  0.34641,        0,        0,        0,        0,        0,        0, -0.34641,        0,
       0,        0,        0,        0,        0,        0,        0,  0.34641,        0,        0,        0,        0,        0,        0, -0.34641;


  EXPECT_TRUE(((Xsig_aug - Xsig_aug_expected).norm() < 0.00001));
}

TEST(UKF, Prediction) {

  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //create example sigma point matrix
  MatrixXd Xsig_aug_data = MatrixXd(n_aug, 2 * n_aug + 1);
  Xsig_aug_data <<
    5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
    1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
    2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
    0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
    0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
    0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
    0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;


  double delta_t = 0.1; //time diff in sec

  PositionPredictor position_predictor_;

  AugmentedSpaceSigmaPoints Xsig_aug(Xsig_aug_data);

  //predict sigma points
  StateSigmaPoints Xsig_pred = position_predictor_.LoadPoints(Xsig_aug, delta_t);

  MatrixXd Xsig_pred_expected = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred_expected <<
  5.93553,  6.06251,  5.92217,   5.9415,  5.92361,  5.93516,  5.93705,  5.93553,  5.80832,  5.94481,  5.92935,  5.94553,  5.93589,  5.93401,  5.93553,
  1.48939,  1.44673,  1.66484,  1.49719,    1.508,  1.49001,  1.49022,  1.48939,   1.5308,  1.31287,  1.48182,  1.46967,  1.48876,  1.48855,  1.48939,
  2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049, 2.23954,   2.2049,  2.12566,  2.16423,  2.11398,   2.2049,   2.2049,  2.17026,   2.2049,
  0.53678, 0.473387, 0.678098, 0.554557, 0.643644, 0.543372,  0.53678, 0.538512, 0.600173, 0.395462, 0.519003, 0.429916, 0.530188,  0.53678, 0.535048,
  0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528, 0.387441, 0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528, 0.318159;

  EXPECT_TRUE(Xsig_pred.isApprox(Xsig_pred_expected, 0.0001));

}