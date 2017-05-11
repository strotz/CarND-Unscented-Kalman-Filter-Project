#ifndef UNSCENTEDKF_STATE_H
#define UNSCENTEDKF_STATE_H

#include "ukf_base.h"

//
// Operations with 5 dimentional state
//
class StateOps {
public:
  StateOps(Eigen::VectorXd &state_) :
    state_ref_(state_) {
  }

  double pos_x() const {
    return state_ref_(0);
  }

  void set_pos_x(const double &value) {
    state_ref_(0) = value;
  }

  double pos_y() const {
    return state_ref_(1);
  }

  void set_pos_y(const double &value) {
    state_ref_(1) = value;
  }

  double velocity() const {
    return state_ref_(2);
  }

  void set_velocity(const double &value) {
    state_ref_(2) = value;
  }

  double yaw_angle() const {
    return state_ref_(3);
  }

  void set_yaw_angle(const double &value) {
    state_ref_(3) = value;
  }

  double yaw_rate() const {
    return state_ref_(4);
  }

  void set_yaw_rate(const double &value) {
    state_ref_(4) = value;
  }

protected:

  Eigen::VectorXd& state_ref_;
};

//
// Concrete class represents 5 dimentional state
//
class State : public StateBase {

public:
  State() : StateBase(5) { // TODO: 5
  }

  State& operator=(const Eigen::VectorXd& other) {
    StateBase::operator=(other);
    return *this;
  }
};

//
// Concrete class that holds state covariance matrix (5)
//
class StateCovariance : public CovarianceBase {
public:
  StateCovariance() : CovarianceBase(5) { // TODO: 5
  }

  StateCovariance& operator=(const Eigen::MatrixXd& other) {
    CovarianceBase::operator=(other);
    return *this;
  }
};

//
//
//
class StateSigmaPoints : public SigmaPointsBase<StateOps> {

public:

  StateSigmaPoints() : 
    SigmaPointsBase(5, SpaceBase::dimension_to_points(7)) // TODO: 5, 7
  {
  }

  StateSigmaPoints& operator=(const StateSigmaPoints& other) {
    SigmaPointsBase::operator=(other);
    return *this;
  }

  Eigen::VectorXd diff_from_mean(int i, const State& mean) const
  {
    Eigen::VectorXd diff = col(i) - mean;
    double angle = StateOps(diff).yaw_angle();
    StateOps(diff).set_yaw_angle(SpaceBase::normalize_angle(angle));
    return  diff;
  }
};

#endif // UNSCENTEDKF_STATE_H
