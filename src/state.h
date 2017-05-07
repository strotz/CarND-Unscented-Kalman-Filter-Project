#ifndef UNSCENTEDKF_STATE_H
#define UNSCENTEDKF_STATE_H

#include "Eigen/Dense"

class State
{
public:
  State() : state_(Eigen::VectorXd(5)) {
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

  double& pos_x() {
    return state_(0);
  }

  double pos_y() const {
    return state_(1);
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
