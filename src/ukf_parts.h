#ifndef UNSCENTEDKF_PARTS_H
#define UNSCENTEDKF_PARTS_H

#include "Eigen/Dense"

class StateData {

public:
  StateData(int size) :
    state_(size) {
    state_.fill(0.0);
  }

  StateData(const Eigen::VectorXd& other) :
    state_(other) {
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

  void set_velocity(const double& value) {
    state_(2) = value;
  }

  double yaw_angle() const {
    return state_(3);
  }

  void set_yaw_angle(const double& value) {
    state_(3) = value;
  }

  double yaw_rate() const {
    return state_(4);
  }

  void set_yaw_rate(const double& value) {
    state_(4) = value;
  }

  double nu_a() const {
    return state_(5);
  }

  double nu_yawdd() const {
    return state_(6);
  }


private:

  Eigen::VectorXd state_;
};

class State : public StateData {
public:
  State(int dimension) :
    n_x_(dimension),
    StateData(n_x_)
  {
  }

  void LoadHead(const State& other) {
    //create augmented mean state
    raw().fill(0.0);
    raw().head(other.n_x_) = other.raw();
  }

  void set_data(const StateData& other) {
    raw() = other.raw();
  }

private:

  int n_x_;
};

class StateCovariance
{
public:
  StateCovariance(int dimension) :
    n_x_(dimension),
    data_(Eigen::MatrixXd(n_x_, n_x_))
  {
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

  void set_data(const StateCovariance& other) {
    data_ = other.data_;
  }

private:

  int n_x_;
  Eigen::MatrixXd data_;
};

class SigmaPoints
{
public:
  SigmaPoints(int state_dimension, int augmented_state_dimension) :
      n_x_(state_dimension),
      n_aug_(augmented_state_dimension),
      Xsig_(Eigen::MatrixXd(n_x_, 2 * n_aug_+ 1)) {
  }

  void Load(const State& x, const StateCovariance& P, double lambda) {

    //calculate square root of P
    Eigen::MatrixXd A = P.sqrt();

    //set first column of sigma point matrix
    Xsig_.col(0)  = x.raw();

    //set remaining sigma points
    double t = sqrt(lambda + n_x_);
    for (int i = 0; i < n_aug_; i++)
    {
      Xsig_.col(i + 1) = x.raw() + t * A.col(i);
      Xsig_.col(i + 1 + n_aug_) = x.raw() - t * A.col(i);
    }
  }

  const StateData point(int i) const {
    return StateData(Xsig_.col(i));
  }

  void set_point(int i, const StateData& data) {
    Xsig_.col(i) = data.raw();
  }

  int points() const {
    return n_aug_ * 2 + 1;
  }

protected:

  int n_x_;
  int n_aug_;

  ///
  Eigen::MatrixXd Xsig_;
};

class Transform : public SigmaPoints
{
public:
  Transform(int state_dimension, int augmented_state_dimension) :
    SigmaPoints(state_dimension, augmented_state_dimension),
    weights_(augmented_state_dimension)
  {
  }

  void InitWeights(double lambda) {
    weights_.fill(0.5 / (n_aug_ + lambda));
    weights_(0) = lambda / (lambda + n_aug_);
  }

  StateData CalculateWeightedMean() {
    StateData result(n_x_);
    for (int i = 0; i < points(); i++) {  //iterate over sigma points
      result.raw() = result.raw() + weights_(i) * Xsig_.col(i);
    }
    return result;
  }

  //angle normalization
  double normalize_angle(double angle) {
    while (angle > M_PI) angle -= 2. * M_PI;
    while (angle < -M_PI) angle += 2. * M_PI;
    return angle;
  }

  StateCovariance CalculateCovariance(const StateData& mean) {
    StateCovariance result(n_x_);
    for (int i = 0; i < points(); i++) {  //iterate over sigma points
      // state difference
      Eigen::VectorXd x_diff = Xsig_.col(i) - mean.raw();
      x_diff(3) = normalize_angle(x_diff(3));
      result.raw() = result.raw() + weights_(i) * x_diff * x_diff.transpose();
    }
    return result;
  }

private:

  ///* Weights of sigma points
  Eigen::VectorXd weights_;
};

class Process
{
public:

  StateData Predict(const StateData& current, double delta_t) {

    //extract values for better readability
    double p_x = current.pos_x();
    double p_y = current.pos_y();
    double v = current.velocity();
    double yaw = current.yaw_angle();
    double yawd = current.yaw_rate();
    double nu_a = current.nu_a();
    double nu_yawdd = current.nu_yawdd();

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    StateData predict(5); // TODO: how to pass size
    predict.set_pos_x(px_p);
    predict.set_pos_y(py_p);
    predict.set_velocity(v_p);
    predict.set_yaw_angle(yaw_p); // TODO: normalize
    predict.set_yaw_rate(yawd_p);
    return predict;
  }
};

#endif //UNSCENTEDKF_PARTS_H
