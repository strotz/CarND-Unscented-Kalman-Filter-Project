#ifndef MEASUREMENT_PACKAGE_H_
#define MEASUREMENT_PACKAGE_H_

#include "Eigen/Dense"

class MeasurementPackage {
public:
  long long timestamp_;

  enum SensorType{
    LASER,
    RADAR
  } sensor_type_;

  Eigen::VectorXd raw_measurements_;

  double lidar_pos_x() const {
    return raw_measurements_(0);
  }

  double lidar_pos_y() const {
    return raw_measurements_(1);
  }

  double radar_distance_ro() const {
    return raw_measurements_(0);
  }

  double radar_angle_phi() const {
    return raw_measurements_(1);
  }

  double radar_velocity_ro_dot() const {
    return raw_measurements_(2);
  }

};

#endif /* MEASUREMENT_PACKAGE_H_ */
