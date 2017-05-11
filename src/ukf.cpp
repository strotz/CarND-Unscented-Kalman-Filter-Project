#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include "ukf_radar.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() :
  n_x_(5),
  n_aug_(7),
  lambda_(3 - n_aug_),
  x_(), // n_x_
  P_(),
  Xsig_pred_() {

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;


  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // TODO: Hint: one or more values initialized above might be wildly off...

  // TODO: init state covariance matrix P_

  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage &measurement) {
  if (!is_initialized_) {
    Initialize(measurement);

    is_initialized_ = true;
    return; // skip predict/update
  }

  //compute the time elapsed between the current and previous measurements
  double dt = (measurement.timestamp_ - time_us_) / 1000000.0;  //dt - expressed in seconds
  time_us_ = measurement.timestamp_;

  Prediction(dt);

  if (measurement.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
    if (use_radar_) {
      UpdateRadar(measurement);
    }
  } else if (measurement.sensor_type_ == MeasurementPackage::SensorType::LASER) {
    if (use_laser_) {
      UpdateLidar(measurement);
    }
  }
}

void UKF::Initialize(const MeasurementPackage &measurement) {

  if (measurement.sensor_type_ == MeasurementPackage::SensorType::LASER) {
    x_.set_pos_x(measurement.lidar_pos_x());
    x_.set_pos_y(measurement.lidar_pos_y());
  } else if (measurement.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
    x_.set_pos_x(measurement.radar_distance_ro() * cos(measurement.radar_angle_phi()));
    x_.set_pos_y(measurement.radar_distance_ro() * sin(measurement.radar_angle_phi()));
  } else {
    // TODO:

  }
  time_us_ = measurement.timestamp_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  //create augmented mean vector
  AugmentedState x_aug = AugmentedState(x_);

  //create augmented covariance matrix
  AugmentedStateCovariance P_aug = AugmentedStateCovariance(P_, std_a_, std_yawdd_);

  AugmentedSpaceSigmaPoints Xsig_aug = AugmentedSpaceSigmaPoints(x_aug, P_aug, lambda_);

  //predict sigma points
  Predictor predictor = Predictor(lambda_);

  Xsig_pred_ = predictor.LoadPoints(Xsig_aug, delta_t);

  //predicted state mean
  x_ = predictor.CalculateWeightedMean(Xsig_pred_);

  //predicted state covariance matrix
  P_ = predictor.CalculateCovariance(Xsig_pred_, x_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage &measurement) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage &measurement) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  RadarSpace radar(lambda_);

  //create matrix for sigma points in measurement space
  RadarSigmaPoints Zsig = radar.LoadPoints(Xsig_pred_);

  //mean predicted measurement
  RadarState z_pred = radar.CalculateWeightedMean(Zsig);

  //measurement covariance matrix S
  RadarCovariance S = radar.CalculateCovariance(Zsig, z_pred);

  //add measurement noise covariance matrix
  S.AddNoise(std_radr_, std_radphi_, std_radrd_);

//  //create example vector for incoming radar measurement
//  VectorXd z = VectorXd(n_z);
//  z <<
//    5.9214,
//    0.2187,
//    2.0062;
//

//  //create matrix for cross correlation Tc
//  MatrixXd Tc = MatrixXd(3, 3);  // TODO: 3
//
//  //calculate cross correlation matrix
//  Tc.fill(0.0);
//  for(int i=0; i < 2 * n_aug + 1; i++)
//  {
//    VectorXd x_diff = Xsig_pred_.col(i) - x;
//    VectorXd z_diff = Zsig.col(i) - z_pred;
//
//    // TODO: normalize angles z_diff(1) and x_diff(3)
//
//    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
//  }
//
//  //calculate Kalman gain K;
//  MatrixXd K = Tc * S.inverse();
//
//  //residual
//  VectorXd z_diff = z - z_pred;
//
//  //angle normalization
//  while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
//  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
//
//  //update state mean and covariance matrix
  // correction K * z_diff
  x_.ApplyCorrection(Eigen::VectorXd(5));
//  P_ = P - K * S * K.transpose();
  P_.ApplyCorrection(-Eigen::MatrixXd(5,5));
}
