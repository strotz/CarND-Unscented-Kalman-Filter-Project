#ifndef UNSCENTEDKF_STATE_H
#define UNSCENTEDKF_STATE_H

#include "Eigen/Dense"

class State
{
public:
  State() : state_(Eigen::VectorXd(5)) {
    state_.fill(0.0);
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

private:

  Eigen::VectorXd state_;
};

#endif //UNSCENTEDKF_STATE_H
