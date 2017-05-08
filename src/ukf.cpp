#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() : n_x_(5),
             n_aug_(7),
    x_(n_x_),
    P_(n_x_),
    Xsig_pred_(n_x_, n_aug_)
  {

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  lambda_ = 3 - n_aug_;

  // TODO: init state covariance matrix

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage& measurement) {
  if (!is_initialized_) {
    Initialize(measurement);

    is_initialized_ = true;
    return; // skip predict/update
  }

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement.timestamp_ - time_us_) / 1000000.0;  //dt - expressed in seconds
  time_us_ = measurement.timestamp_;

  Prediction(time_us_);

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
  }
  else if (measurement.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
    x_.set_pos_x(measurement.radar_distance_ro() * cos(measurement.radar_angle_phi()));
    x_.set_pos_y(measurement.radar_distance_ro() * sin(measurement.radar_angle_phi()));
  }
  else {
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
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& measurement) {
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
void UKF::UpdateRadar(const MeasurementPackage& measurement) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
