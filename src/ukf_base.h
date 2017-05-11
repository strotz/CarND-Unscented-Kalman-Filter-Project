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
class StateBase : public Eigen::VectorXd {

protected:

  StateBase(int size) :
    size_(size),
    Eigen::VectorXd(size)
  {
    fill(0.0);
  }

public:

  StateBase& operator=(const Eigen::VectorXd& other) {
    Eigen::VectorXd::operator=(other);
    return *this;
  }

  void set_head(int head_size, const Eigen::VectorXd& value) {
    head(head_size) = value;
  }

  int size() const {
    return size_;
  }

private:

  int size_;
};


//
// Base class for covariances
//

class CovarianceBase : public Eigen::MatrixXd {
protected:

  CovarianceBase(int size) :
    Eigen::MatrixXd(size, size),
    size_(size)
  {
    fill(0.0);
  }

public:

  CovarianceBase& operator=(const Eigen::MatrixXd& other) {
    Eigen::MatrixXd::operator=(other);
    return *this;
  }

  int size_;
};

//
// Container for Sigma Points
//

template <class Ops>
class SigmaPointsBase : public Eigen::MatrixXd {

protected:

  SigmaPointsBase(int size, int number_of_points) :
    size_(size),
    number_of_points_(number_of_points),
    Eigen::MatrixXd(size, number_of_points)
  {
  }

  SigmaPointsBase(int size) :
    SigmaPointsBase(size, SpaceBase::dimension_to_points(size)) {
  }

public:
  
  int number_of_points() const {
    return number_of_points_;
  }

  const Ops point(int i) const {
    Eigen::VectorXd c = col(i);
    return Ops(c);
  }

  int size() const {
    return size_;
  }

private:

  int size_;
  int number_of_points_;
};


#endif // UNSCENTEDKF_BASE_H