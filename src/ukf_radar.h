#ifndef UNSCENTEDKF_RADER_H
#define UNSCENTEDKF_RADER_H

#include "ukf_state.h"
#include "ukf_parts.h"

const int RadarSpaceDim = 3;


class RadarOps {
  Eigen::VectorXd& state_ref_;

public:
  RadarOps(Eigen::VectorXd &state_) : state_ref_(state_) {}

  double r() const {
    return state_ref_(0);
  }

  RadarOps& set_r(double value) {
    state_ref_(0) = value;
    return *this;
  }

  double phi() const {
    return state_ref_(1);
  }

  RadarOps& set_phi(double value) {
    state_ref_(1) = value;
    return *this;
  }

  double r_dot() const {
    return state_ref_(2);
  }

  RadarOps& set_r_dot(double value) {
    state_ref_(2) = value;
    return *this;
  }
};

class RadarState : public StateBase<RadarSpaceDim> {
public:
  explicit RadarState() :
    StateBase() {
  }

  RadarState& operator=(const Eigen::VectorXd& other) {
    StateBase::operator=(other);
    return *this;
  }
};

class RadarCovariance : public CovarianceBase<RadarSpaceDim> {
public:
  explicit RadarCovariance() : CovarianceBase()
  {
  }

  RadarCovariance& operator=(const Eigen::MatrixXd& other) {
    CovarianceBase::operator=(other);
    return *this;
  }
};

class RadarSigmaPoints : public SigmaPointsBase<RadarSpaceDim, RadarOps> {
  friend class RadarSpace;

public:
  explicit RadarSigmaPoints(int number_of_points) :
    SigmaPointsBase(number_of_points) {
  }

  Eigen::VectorXd diff_from_mean(int i, const RadarState &mean) const {
    Eigen::VectorXd diff = col(i) - mean;
    double angle = RadarOps(diff).phi();
    RadarOps(diff).set_phi(SpaceBase::normalize_angle(angle));
    return diff;
  }
};

class RadarSpace : public SpaceTransformation<RadarSigmaPoints, RadarState, RadarCovariance> {

public:
  explicit RadarSpace() :
    SpaceTransformation()
  {
  }

  static RadarState ConvertToRadarSpace(const StateOps &point) {

    // extract values for better readibility
    double p_x = point.pos_x();
    double p_y = point.pos_y();
    double v = point.velocity();
    double yaw = point.yaw_angle();

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    RadarState data;
    RadarOps(data)
      .set_r(sqrt(p_x * p_x + p_y * p_y))  //r
      .set_phi(atan2(p_y, p_x))            //phi
      .set_r_dot((p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y));   //r_dot

    return data;
  }

  static RadarSigmaPoints LoadPoints(const StateSigmaPoints &state_sigma_points) {
    int number_of_points = state_sigma_points.number_of_points();
    RadarSigmaPoints result(number_of_points);
    for (int i = 0; i < number_of_points; i++) {  // iterate over sigma points
      Eigen::VectorXd t = state_sigma_points.col(i);
      StateOps ops(t);
      result.col(i) = ConvertToRadarSpace(ops);
    }
    return result;
  }
};


#endif // UNSCENTEDKF_RADER_H