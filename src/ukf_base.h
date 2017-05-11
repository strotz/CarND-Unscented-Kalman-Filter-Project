#ifndef UNSCENTEDKF_BASE_H
#define UNSCENTEDKF_BASE_H

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
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
template<int Dim>
class StateBase : public Eigen::VectorXd {

protected:

  StateBase() :
    Eigen::VectorXd(Dim)
  {
  }

public:

  StateBase& operator=(const Eigen::VectorXd& other) {
    Eigen::VectorXd::operator=(other);
    return *this;
  }

  int size() const {
    return Dim;
  }
};


//
// Base class for covariances
//

template<int Dim>
class CovarianceBase : public Eigen::MatrixXd {

protected:

  CovarianceBase() :
    Eigen::MatrixXd(Dim, Dim)
  {
  }

public:

  CovarianceBase& operator=(const Eigen::MatrixXd& other) {
    Eigen::MatrixXd::operator=(other);
    return *this;
  }
};

//
// Container for Sigma Points
//

template <int Dim, class Ops>
class SigmaPointsBase : public Eigen::MatrixXd {

protected:

  SigmaPointsBase(int number_of_points) :
    number_of_points_(number_of_points),
    Eigen::MatrixXd(Dim, number_of_points)
  {
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
    return Dim;
  }

private:

  int number_of_points_;
};


#endif // UNSCENTEDKF_BASE_H