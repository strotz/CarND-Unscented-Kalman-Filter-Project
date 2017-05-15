#ifndef UNSCENTEDKF_BASE_H
#define UNSCENTEDKF_BASE_H

#include "Eigen/Dense"

using Eigen::MatrixXd;

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

template <int Dim>
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

class Weights : public Eigen::VectorXd
{
  int number_of_points_;

public:

  Weights(int number_of_points) :
    number_of_points_(number_of_points),
    Eigen::VectorXd(number_of_points)
  {
  }

  void Initialize(double lambda, int space_dimension) {
    double d = lambda + space_dimension;
    fill(0.5 / d);
    (*this)(0) = lambda / d;
  }

  int number_of_points() const
  {
    return number_of_points_;
  }
};


template<class TargetSigmaPoints, class TargetState, class TargetStateCovariance>
class SpaceTransformation
{
public:

  TargetStateCovariance CalculateCovariance(
    const Weights& weights,
    const TargetSigmaPoints& sigma_points,
    const TargetState& mean)
  {
    TargetStateCovariance result;
    result.fill(0.0);
    int len = weights.number_of_points();
    for (int i = 0; i < len; i++) {  // iterate over sigma points
      // state difference
      Eigen::VectorXd x_diff = sigma_points.diff_from_mean(i, mean);
      result = result + weights(i) * x_diff * x_diff.transpose();
    }
    return result;
  }
};



#endif // UNSCENTEDKF_BASE_H