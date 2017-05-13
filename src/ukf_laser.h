#ifndef UNSCENTEDKF_UKF_LASER_H
#define UNSCENTEDKF_UKF_LASER_H

#include "ukf_base.h"
#include "ukf_state.h"

const int LidarSpaceDim = 2;

class LidarState : public StateBase<LidarSpaceDim> {
public:
  explicit LidarState() :
    StateBase() {
  }

  LidarState(const Eigen::VectorXd& other) :
    StateBase() {
    StateBase::operator=(other);
  }
};

class LidarCovariance : public CovarianceBase<LidarSpaceDim> {
public:
  explicit LidarCovariance() : CovarianceBase()
  {
  }

  LidarCovariance& operator=(const Eigen::MatrixXd& other) {
    CovarianceBase::operator=(other);
    return *this;
  }
};

class LidarSigmaPoints : public SigmaPointsBase<LidarSpaceDim> {
public:
  explicit LidarSigmaPoints(int number_of_points) :
    SigmaPointsBase(number_of_points) {
  }

  Eigen::VectorXd diff_from_mean(int i, const LidarState &mean) const {
    return col(i) - mean;
  }
};

class LidarSpace : public SpaceTransformation<LidarSigmaPoints, LidarState, LidarCovariance> {
public:
  static LidarState ConvertToRadarSpace(const StateOps &point) {

    // extract values for better readibility
    double p_x = point.pos_x();
    double p_y = point.pos_y();

    // measurement model
    LidarState data;
    data << p_x, p_y;
    return data;
  }

  static LidarSigmaPoints LoadPoints(const StateSigmaPoints &state_sigma_points) {
    int number_of_points = state_sigma_points.number_of_points();
    LidarSigmaPoints result(number_of_points);
    for (int i = 0; i < number_of_points; i++) {  // iterate over sigma points
      Eigen::VectorXd t = state_sigma_points.col(i);
      StateOps ops(t);
      result.col(i) = ConvertToRadarSpace(ops);
    }
    return result;
  }
};


#endif //UNSCENTEDKF_UKF_LASER_H
