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

  // angle normalization
  static double normalize_angle(double angle) {
    return fmod(angle + M_PI, 2. * M_PI) - M_PI;
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
    fill(0.0);
  }

public:

  StateBase& operator=(const Eigen::VectorXd& other) {
    Eigen::VectorXd::operator=(other);
    return *this;
  }

  int dim() const {
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
    fill(0.0);
  }

public:

  CovarianceBase& operator=(const Eigen::MatrixXd& other) {
    Eigen::MatrixXd::operator=(other);
    return *this;
  }

  int dim() const {
    return Dim;
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
    fill(0.0);
  }

public:
  
  int number_of_points() const {
    return number_of_points_;
  }

  int dim() const {
    return Dim;
  }

private:

  int number_of_points_;
};


#endif // UNSCENTEDKF_BASE_H