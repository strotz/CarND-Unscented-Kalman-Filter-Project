#ifndef UNSCENTEDKF_BASE_H
#define UNSCENTEDKF_BASE_H

#include "Eigen/Dense"

class SpaceBase
{
public:

  static int dimension_to_points(int dimension) {
    return dimension * 2 + 1;
  }

  //angle normalization
  static double normalize_angle(double angle) {
    while (angle > M_PI) angle -= 2. * M_PI;
    while (angle < -M_PI) angle += 2. * M_PI;
    return angle;
  }
};

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

public: 

  const Eigen::VectorXd& raw() const {
    return state_;
  }

  void set_raw(const Eigen::VectorXd& value) {
    state_ = value;
  }

  void set_head(int head_size, const Eigen::VectorXd& value) {
    state_.head(head_size) = value;
  }

  int size() const {
    return size_;
  }

private:

  int size_;

protected:

  Eigen::VectorXd state_;
};


//
// Base class for covariances
//

class CovarianceBase {
protected:

  CovarianceBase(int size) :
    data_(Eigen::MatrixXd(size, size)), 
    size_(size)
  {
    data_.fill(0.0);
  }

  int size_;
  Eigen::MatrixXd data_; // TODO: how we can hide it?

public:

  const Eigen::MatrixXd& data() const {
    return data_;
  }

  void set_data(const Eigen::MatrixXd& data) {
    data_ = data;
  }
};

//
// Container for Sigma Points
//

template <class Ops>
class SigmaPointsBase {
protected:

  SigmaPointsBase(int size, int number_of_points) :
    size_(size),
    number_of_points_(number_of_points),
    sigma_points_(size, number_of_points) {
  }

  SigmaPointsBase(int size) :
    SigmaPointsBase(size, SpaceBase::dimension_to_points(size)) {
  }

public:
  
  int number_of_points() const {
    return number_of_points_;
  }

  const Ops point(int i) const {
    Eigen::VectorXd c = sigma_points_.col(i);
    return Ops(c);
  }

  const Eigen::MatrixXd& raw() const {
    return sigma_points_;
  }

  void set_raw(const Eigen::MatrixXd& points) {
    sigma_points_ = points;
  }

  const Eigen::VectorXd raw_point(int i) const {
    return sigma_points_.col(i);
  }

  void set_raw_point(int i, const Eigen::VectorXd& point) {
    sigma_points_.col(i) = point;
  }

  int size() const {
    return size_;
  }

private:

  int size_;
  int number_of_points_;
  Eigen::MatrixXd sigma_points_;
};


#endif // UNSCENTEDKF_BASE_H