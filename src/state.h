#ifndef UNSCENTEDKF_STATE_H
#define UNSCENTEDKF_STATE_H

#include "Eigen/Dense"

class StateTraits {
public:
  StateTraits(Eigen::VectorXd data) : state_(data) {
  }

  const Eigen::VectorXd& const_raw() const {
    return state_;
  }

  Eigen::VectorXd& raw() {
    return state_;
  }

  //
  // getters and setters
  //

  double pos_x() const {
    return state_(0);
  }

  void set_pos_x(const double& value) {
    state_(0) = value;
  }

  double pos_y() const {
    return state_(1);
  }

  void set_pos_y(const double& value) {
    state_(1) = value;
  }

  double velocity() const {
    return state_(2);
  }

  double yaw_angle() const {
    return state_(3);
  }

  double yaw_rate() const {
    return state_(4);
  }

protected:

  Eigen::VectorXd state_;
};

class State : public StateTraits
{
public:
  State(int state_dimension) : StateTraits(Eigen::VectorXd(state_dimension)),
    n_x_(state_dimension)
  {
    state_.fill(0.0);
  }

private:

  int n_x_;
};

class AugmentedState : public StateTraits
{
public:

  AugmentedState(int state_dimension, int augmented_state_dimension) :
          StateTraits(Eigen::VectorXd(augmented_state_dimension)),
          n_x_(state_dimension),
          n_aug_(augmented_state_dimension)
  {
  }

  void Load(const State& state) {
    state_.head(n_x_) = state.const_raw();
  }

private:

  int n_x_;
  int n_aug_;

};


class StateCovariance
{
public:
  StateCovariance(int state_dimension) : data_(Eigen::MatrixXd(state_dimension, state_dimension)) {
  }

  Eigen::MatrixXd& raw() {
    return data_;
  }

private:
  Eigen::MatrixXd data_;
};

class SigmaPoints
{
public:
  SigmaPoints(int state_dimension, int augmented_state_dimension) : data_(Eigen::MatrixXd(state_dimension, 2 * augmented_state_dimension + 1)) {
  }

  void Load(const AugmentedState& augmentedState) {
  }

private:
  Eigen::MatrixXd data_;
};

#endif //UNSCENTEDKF_STATE_H
