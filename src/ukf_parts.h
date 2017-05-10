#ifndef UNSCENTEDKF_PARTS_H
#define UNSCENTEDKF_PARTS_H

#include "Eigen/Dense"

//
// Base for all state classes
//
template<class Operations>
class StateBase : public Operations {

protected:

  StateBase(int size) :
    size_(size),
    state_(size),
    Operations(state_) {

    state_.fill(0.0);
  }

public: // TODO: is it posible to avoid?

  const Eigen::VectorXd& raw() const {
    return state_;
  }

protected:

  void set_raw(const Eigen::VectorXd& value) {
    state_ = value;
  }

  void set_head(int head_size, const Eigen::VectorXd& value) {
    state_.head(head_size) = value;
  }

public:

  int size() const {
    return size_;
  }

private:

  int size_;
  Eigen::VectorXd state_;
};

//
// Operations with 5 dimentional state
//
class StateOps {
public:
  StateOps(Eigen::VectorXd &state_) :
    state_(state_) {
  }

public:
  //
  // getters and setters
  //

  double pos_x() const {
    return state_(0);
  }

  void set_pos_x(const double &value) {
    state_(0) = value;
  }

  double pos_y() const {
    return state_(1);
  }

  void set_pos_y(const double &value) {
    state_(1) = value;
  }

  double velocity() const {
    return state_(2);
  }

  void set_velocity(const double &value) {
    state_(2) = value;
  }

  double yaw_angle() const {
    return state_(3);
  }

  void set_yaw_angle(const double &value) {
    state_(3) = value;
  }

  double yaw_rate() const {
    return state_(4);
  }

  void set_yaw_rate(const double &value) {
    state_(4) = value;
  }

protected:

  Eigen::VectorXd& state_;
};

//
// Concrete class represents 5 dimentional state
//
class State : public StateBase<StateOps> {
public:
  State() :
    StateBase(5) { // TODO: 5
  }

  State& operator=(const State& other) {
    set_raw(other.raw());
    return *this;
  }
};

//
// Class represents operations on augmented state (7 dim)
//
class AugmentedStateOps : public  StateOps {
public:

  AugmentedStateOps(Eigen::VectorXd &state_) : StateOps(state_) {}

  double nu_a() const {
    return state_(5);
  }

  double nu_yawdd() const {
    return state_(6);
  }
};

//
// Concrete class that represents 7 dimentional (augmented) state
//
class AugmentedState : public StateBase<AugmentedStateOps> {
private:
  AugmentedState() :
    StateBase(7) { // TODO: 7
  }

public:
  static AugmentedState Extend(const State &state) {

    AugmentedState result;

    //create augmented mean state
    result.set_head(state.size(), state.raw());

    return result;
  }
};

//
//
//

class CovarianceBase {
protected:

  CovarianceBase(int size) :
    data_(Eigen::MatrixXd(size, size)) {
    data_.fill(0.0);
  }

public:

  const Eigen::MatrixXd& raw_data() const {
    return data_;
  }

  // TODO: can it be avoided?
  Eigen::MatrixXd& raw_data() {
    return data_;
  }

protected:

  Eigen::MatrixXd data_; // TODO: how we can hide it?
};

//
// Concrete class that holds state covariance matrix (5)
//
class StateCovariance : public CovarianceBase {
public:

  StateCovariance() : CovarianceBase(5) { // TODO: 5
  }

//  Eigen::MatrixXd &raw() {
//    return data_;
//  }
//
//  const Eigen::MatrixXd &raw() const {
//    return data_;
//  }
//
//  Eigen::MatrixXd sqrt() const {
//    return data_.llt().matrixL();
//  }
//
//  double &at(int x, int y) {
//    return data_(x, y);
//  }
//
//  void set_data(const StateCovariance &other) {
//    data_ = other.data_;
//  }
//
//private:
//
//  int n_x_;

};

//
// Concrete class that holds covariance matrix for augmented state (7)
//
class AugmentedStateCovariance : public CovarianceBase {
private:

  AugmentedStateCovariance() :
    CovarianceBase(7) // TODO: 7
  {}

public:

  Eigen::MatrixXd sqrt() const {
    return data_.llt().matrixL();
  }

public:

  static AugmentedStateCovariance Extend(const StateCovariance& covariance, double std_a, double std_yawdd) {
    AugmentedStateCovariance result;
    result.data_.topLeftCorner(5, 5) = covariance.raw_data();
    result.data_(5, 5) = std_a * std_a;
    result.data_(6, 6) = std_yawdd * std_yawdd;

    return result;
  }
};

template <class Ops>
class SigmaPointsBase {
protected:

  SigmaPointsBase(int size, int number_of_points) :
    size_(size),
    number_of_points_(number_of_points),
    sigma_points_(size, number_of_points) {
  }

  SigmaPointsBase(int size) :
    SigmaPointsBase(size, SigmaPointsBase::dimension_to_points(size)) {
  }

  static int dimension_to_points(int dimension) {
    return dimension * 2 + 1;
  }

public:

  int number_of_points() const {
    return number_of_points_;
  }

  const Ops point(int i) const {
    return Ops(sigma_points_.col(i));
  }

protected:

  Eigen::VectorXd& raw_point(int i) {
    return sigma_points_.col(i);
  }

  void set_raw_point(int i, const Eigen::VectorXd& point) {
    sigma_points_.col(i) = point;
  }

protected: // TOOO: make protected getters

  int size_;
  Eigen::MatrixXd sigma_points_;

private:

  int number_of_points_;
};

//
//template<class Space>
//class SigmaPoints
//{
//
//public:
//
//  const Space point(int i) const {
//    return Space(sigma_points_.col(i));
//  }
//
//  void set_point(int i, const Space &data) {
//    sigma_points_.col(i) = data.raw();
//  }
//
//  int points() const {
//    return number_of_points_;
//  }
//
//protected:
//
//  int n_x_;
//};

class AugmentedSpaceSigmaPoints : public SigmaPointsBase<AugmentedStateOps> {
public:
  AugmentedSpaceSigmaPoints() :
    SigmaPointsBase(7) {} // TODO: 7

  static AugmentedSpaceSigmaPoints CreateSigmaPoints(
    const AugmentedState& state,
    const AugmentedStateCovariance& covariance,
    double lambda) {

    AugmentedSpaceSigmaPoints result;

    //calculate square root of P
    Eigen::MatrixXd A = covariance.sqrt();

    //set first column of sigma point matrix
    result.sigma_points_.col(0) = state.raw();

    //set remaining sigma points
    double t = sqrt(lambda + result.size_);
    for (int i = 0; i < result.size_; i++) {
      result.sigma_points_.col(i + 1) = state.raw() + t * A.col(i);
      result.sigma_points_.col(i + 1 + result.size_) = state.raw() - t * A.col(i);
    }

    return result;
  }
};

class Transform : public SigmaPointsBase<State> {
public:
  Transform() :
    SigmaPointsBase(5, SigmaPointsBase::dimension_to_points(7)), // TODO: 5, 7
    weights_(SigmaPointsBase::dimension_to_points(7)) {
  }

  void InitWeights(double lambda) {
    throw std::exception();
//    int d = lambda + n_aug_;
//    weights_.fill(0.5 / d);
//    weights_(0) = lambda / d;
  }

private:

  static State PredictNextState(const AugmentedStateOps& current, double delta_t) {

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
    State result;
    result.set_pos_x(px_p);
    result.set_pos_y(py_p);
    result.set_velocity(v_p);
    result.set_yaw_angle(yaw_p); // TODO: normalize
    result.set_yaw_rate(yawd_p);
    return result;
  }


public:

  void LoadPoints(const AugmentedSpaceSigmaPoints& Xsig_aug, double delta_t) {
    for (int i = 0; i < Xsig_aug.number_of_points(); i++) {
      State predicted = PredictNextState(Xsig_aug.point(i), delta_t);
      set_raw_point(i, predicted.raw());
    }
  }

  State CalculateWeightedMean() {
    State result;
    for (int i = 0; i < number_of_points(); i++) {  //iterate over sigma points
      result. TODO = result.raw() + weights_(i) * point(i).raw();
    }
    return result;
  }

  //angle normalization
  double normalize_angle(double angle) {
    while (angle > M_PI) angle -= 2. * M_PI;
    while (angle < -M_PI) angle += 2. * M_PI;
    return angle;
  }

  StateCovariance CalculateCovariance(const State& mean) {
    StateCovariance result;
    for (int i = 0; i < number_of_points(); i++) {  // iterate over sigma points
      // state difference
      Eigen::VectorXd x_diff = raw_point(i) - mean.raw();
      x_diff(3) = normalize_angle(x_diff(3));
      result.raw_data() = result.raw_data() + weights_(i) * x_diff * x_diff.transpose(); // TODO: function style?
    }
    return result;
  }


private:

  ///* Weights of sigma points
  Eigen::VectorXd weights_;
};

////////////////////
//
//

class RadarOps {
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

private:
  Eigen::VectorXd& state_;
};

class RadarState : public StateBase<RadarOps> {
public:
  RadarState() :
    StateBase(3) {
  }
};

class RadarSpace : public SigmaPointsBase<RadarOps> {
public:
  RadarSpace() :
    SigmaPointsBase(3, SigmaPointsBase::dimension_to_points(7)) { // TODO: 3, 7
  }

  static RadarState ConvertToRadarSpace(const StateOps& point) {

    // extract values for better readibility
    double p_x = point.pos_x();
    double p_y = point.pos_y();
    double v = point.velocity();
    double yaw = point.yaw_angle();

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    RadarState result;
    result.set_r(sqrt(p_x * p_x + p_y * p_y));  //r
    result.set_phi(atan2(p_y, p_x));             //phi
    result.set_r_dot((p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y));   //r_dot

    return result;
  }


  static RadarSpace LoadPoints(const Transform& transform) {
    RadarSpace result;

    for (int i = 0; i < transform.number_of_points(); i++) {  // iterate over sigma points

      State point = transform.point(i);
      RadarState radar = RadarSpace::ConvertToRadarSpace(point);
      result.set_raw_point(i, radar.raw());
    }

    return result;
  }
};

class StatePredictor {
public:

};


class SpaceConversion {
public:
};

#endif //UNSCENTEDKF_PARTS_H
