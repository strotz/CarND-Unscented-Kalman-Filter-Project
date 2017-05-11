#ifndef UNSCENTEDKF_RADER_H
#define UNSCENTEDKF_RADER_H

#include "Eigen/Dense"
#include "ukf_state.h"
#include "ukf_parts.h"

class RadarOps {
  Eigen::VectorXd &state_;

public:
  RadarOps(Eigen::VectorXd &state_) : state_(state_) {}

  double r() const {
    return state_(0);
  }

  void set_r(double value) {
    state_(0) = value;
  }

  double phi() const {
    return state_(1);
  }

  void set_phi(double value) {
    state_(1) = value;
  }

  double r_dot() const {
    return state_(2);
  }

  void set_r_dot(double value) {
    state_(2) = value;
  }
};

class RadarState : public StateBase {
public:
  RadarState() :
    StateBase(3) {
  }

  RadarState& operator=(const Eigen::VectorXd& other) {
    StateBase::operator=(other);
    return *this;
  }

};

class RadarSigmaPoints : public SigmaPointsBase<RadarOps> {
  friend class RadarSpace;

public:
  RadarSigmaPoints() :
    SigmaPointsBase(3, SpaceBase::dimension_to_points(7)) { // TODO: 3, 7
  }

  Eigen::VectorXd diff_from_mean(int i, const RadarState &mean) const {
    Eigen::VectorXd diff = col(i) - mean;
    double angle = RadarOps(diff).phi();
    RadarOps(diff).set_phi(SpaceBase::normalize_angle(angle));
    return diff;
  }
};



class RadarCovariance : public CovarianceBase {
public:
  RadarCovariance() :
    CovarianceBase(3) // TODO: 3
  {
  }

  RadarCovariance& operator=(const Eigen::MatrixXd& other) {
    Eigen::MatrixXd::operator=(other);
    return *this;
  }
};


class RadarSpace : public SpaceTransformation<RadarSigmaPoints, RadarState, RadarCovariance> {

public:
  RadarSpace(double lambda) :
    SpaceTransformation(lambda, 7) // TODO: 7
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

    RadarOps ops = RadarOps(data); // TODO: wrap it
    ops.set_r(sqrt(p_x * p_x + p_y * p_y));  //r
    ops.set_phi(atan2(p_y, p_x));             //phi
    ops.set_r_dot((p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y));   //r_dot

    return data;
  }

  static RadarSigmaPoints LoadPoints(const StateSigmaPoints &state_sigma_points) {

    RadarSigmaPoints result;

    for (int i = 0; i < state_sigma_points.number_of_points(); i++) {  // iterate over sigma points
      result.col(i) = ConvertToRadarSpace(state_sigma_points.point(i));
    }
    return result;
  }
};


#endif // UNSCENTEDKF_RADER_H