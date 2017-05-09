#ifndef UNSCENTEDKF_STATE_H
#define UNSCENTEDKF_STATE_H

#include "Eigen/Dense"

class State {
public:
  State(int dimension) : n_x_(dimension),
    state_(Eigen::VectorXd(n_x_))
  {
    state_.fill(0.0);
  }

  const Eigen::VectorXd& raw() const {
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

  void LoadHead(const State& other) {

    //create augmented mean state
    state_.fill(0.0);
    state_.head(other.n_x_) = other.raw();
  }

private:

  int n_x_;
  Eigen::VectorXd state_;
};

class StateCovariance
{
public:
  StateCovariance(int dimension) : n_x_(dimension),
                                   data_(Eigen::MatrixXd(n_x_, n_x_)) {
    data_.fill(0.0);
  }

  Eigen::MatrixXd& raw() {
    return data_;
  }

  Eigen::MatrixXd sqrt() const {
    return data_.llt().matrixL();
  }

  void LoadTopLeft(const StateCovariance& other) {
    data_.topLeftCorner(5,5) = other.data_;
  }

  double& at(int x, int y) {
    return data_(x, y);
  }

private:

  int n_x_;
  Eigen::MatrixXd data_;
};

class SigmaPoints
{
public:
  SigmaPoints(int state_dimension, int augmented_state_dimension) : n_x_(state_dimension),
                                                                    n_aug_(augmented_state_dimension),
          data_(Eigen::MatrixXd(n_x_, 2 * n_aug_+ 1)) {
  }

  void Load(const State& x, const StateCovariance& P, double lambda) {

    //calculate square root of P
    Eigen::MatrixXd A = P.sqrt();

    //set first column of sigma point matrix
    data_.col(0)  = x.raw();

    //set remaining sigma points
    double t = sqrt(lambda + n_x_);
    for (int i = 0; i < n_aug_; i++)
    {
      data_.col(i + 1) = x.raw() + t * A.col(i);
      data_.col(i + 1 + n_aug_) = x.raw() - t * A.col(i);
    }
  }

private:

  int n_x_;
  int n_aug_;
  Eigen::MatrixXd data_;
};

#endif //UNSCENTEDKF_STATE_H
