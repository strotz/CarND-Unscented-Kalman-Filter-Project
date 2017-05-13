#include "ukf.h"
#include "tools.h"

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
  lambda_(3 - AugmentedSpaceDim),
  x_(),
  P_(),
  number_of_points_(SpaceBase::dimension_to_points(AugmentedSpaceDim)),
  Xsig_pred_(SpaceBase::dimension_to_points(AugmentedSpaceDim)),
  weights_(SpaceBase::dimension_to_points(AugmentedSpaceDim)),
  position_predictor_()
{

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 5; // TODO: tuning parameter

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1; // TODO: tuning parameter

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

  weights_.Initialize(lambda_, AugmentedSpaceDim);

  // init state covariance matrix P_
  P_ <<
    1, 0, 0, 0, 0,
    0, 1, 0, 0, 0,
    0, 0, 1000, 0, 0,
    0, 0, 0, 1000, 0,
    0, 0, 0, 0, 1000;

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

  StateOps ops = StateOps(x_);
  if (measurement.sensor_type_ == MeasurementPackage::SensorType::LASER) {
    ops.set_pos_x(measurement.lidar_pos_x());
    ops.set_pos_y(measurement.lidar_pos_y());
  } else if (measurement.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
    ops.set_pos_x(measurement.radar_distance_ro() * cos(measurement.radar_angle_phi()));
    ops.set_pos_y(measurement.radar_distance_ro() * sin(measurement.radar_angle_phi()));
  } else {
    cout << "Error: unsupported sensor type" << endl;
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
  Xsig_pred_ = position_predictor_.LoadPoints(Xsig_aug, delta_t);

  //predicted state mean
  x_ = position_predictor_.CalculateWeightedMean(weights_,  Xsig_pred_);

  //predicted state covariance matrix
  P_ = position_predictor_.CalculateCovariance(weights_, Xsig_pred_, x_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage &measurement) {
  /**
  TODO: implement lidar

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
  RadarSpace radar;

  //create matrix for sigma points in measurement space
  RadarSigmaPoints Zsig = radar.LoadPoints(Xsig_pred_);

  //mean predicted measurement
  RadarState z_pred = radar.CalculateWeightedMean(weights_, Zsig);

  //measurement covariance matrix S
  RadarCovariance S = radar.CalculateCovariance(weights_, Zsig, z_pred);

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(RadarSpaceDim, RadarSpaceDim);
  R << std_radr_ * std_radr_, 0, 0,
    0, std_radphi_ * std_radphi_, 0,
    0, 0, std_radrd_ * std_radrd_;

  S = S + R;

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(RadarSpaceDim);
  z << 
    measurement.radar_distance_ro(), 
    measurement.radar_angle_phi(), 
    measurement.radar_velocity_ro_dot(); 

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(SpaceDim, RadarSpaceDim);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for(int i=0; i < number_of_points_; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;

    z_diff(1) = SpaceBase::normalize_angle(z_diff(1));
    x_diff(3) = SpaceBase::normalize_angle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = SpaceBase::normalize_angle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  // calculate the radar NIS.
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
