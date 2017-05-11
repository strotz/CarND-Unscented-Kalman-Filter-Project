#ifndef UNSCENTEDKF_STATE_H
#define UNSCENTEDKF_STATE_H

#include "ukf_base.h"

const int SpaceDim = 5;

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

  StateOps& set_pos_x(const double &value) {
    state_ref_(0) = value;
    return *this;
  }

  double pos_y() const {
    return state_ref_(1);
  }

  StateOps& set_pos_y(const double &value) {
    state_ref_(1) = value;
    return *this;
  }

  double velocity() const {
    return state_ref_(2);
  }

  StateOps& set_velocity(const double &value) {
    state_ref_(2) = value;
    return *this;
  }

  double yaw_angle() const {
    return state_ref_(3);
  }

  StateOps& set_yaw_angle(const double &value) {
    state_ref_(3) = value;
    return *this;
  }

  double yaw_rate() const {
    return state_ref_(4);
  }

  StateOps& set_yaw_rate(const double &value) {
    state_ref_(4) = value;
    return *this;
  }

protected:

  Eigen::VectorXd& state_ref_;
};

//
// Concrete class represents 5 dimentional state
//
class State : 
  public StateBase<SpaceDim>{

public:
  explicit State() : 
    StateBase()
  {
  }

  State& operator=(const Eigen::VectorXd& other) {
    StateBase::operator=(other);
    return *this;
  }
};

//
// Concrete class that holds state covariance matrix (5)
//
class StateCovariance : public CovarianceBase<SpaceDim> {
public:
  StateCovariance() : CovarianceBase() { 
  }

  StateCovariance& operator=(const Eigen::MatrixXd& other) {
    CovarianceBase::operator=(other);
    return *this;
  }
};

//
//
//
class StateSigmaPoints : public SigmaPointsBase<SpaceDim, StateOps> {

public:

  StateSigmaPoints(int number_of_points) : 
    SigmaPointsBase(number_of_points)
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
