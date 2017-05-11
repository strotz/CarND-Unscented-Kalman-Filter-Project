#ifndef UNSCENTEDKF_AUGMENTED_H
#define UNSCENTEDKF_AUGMENTED_H

#include "ukf_state.h"

//
// Class represents operations on augmented state (7 dim)
//
class AugmentedStateOps : public  StateOps {
public:

  AugmentedStateOps(Eigen::VectorXd &state_) : StateOps(state_) {}

  double nu_a() const {
    return state_ref_(5);
  }

  double nu_yawdd() const {
    return state_ref_(6);
  }
};

//
// Concrete class that represents 7 dimentional (augmented) state
//
class AugmentedState : public StateBase {

public:
  AugmentedState(const State &state) :
    StateBase(7) // TODO: 7
  { 
    //create augmented mean state
    set_head(state.size(), state);
  }

  AugmentedState& operator=(const Eigen::VectorXd& other) {
    StateBase::operator=(other);
    return *this;
  }
};


//
// Concrete class that holds covariance matrix for augmented state (7)
//
class AugmentedStateCovariance : public CovarianceBase {

public:

  AugmentedStateCovariance(const StateCovariance& covariance, double std_a, double std_yawdd) :
    CovarianceBase(7) // TODO: 7
  {
    topLeftCorner(5, 5) = covariance;
    (*this)(5, 5) = std_a * std_a;
    (*this)(6, 6) = std_yawdd * std_yawdd;
  }

  Eigen::MatrixXd sqrt() const {
    return llt().matrixL();
  }

  AugmentedStateCovariance& operator=(const Eigen::MatrixXd& other) {
    CovarianceBase::operator=(other);
    return *this;
  }

};

//
// Concrete class that holds Sigma Points within augmented space (7)
//
class AugmentedSpaceSigmaPoints : public SigmaPointsBase<AugmentedStateOps> {
public:
  AugmentedSpaceSigmaPoints(
    const AugmentedState& state,
    const AugmentedStateCovariance& covariance,
    double lambda
    ) : SigmaPointsBase(7) 
  { // TODO: 7
    //calculate square root of P
    Eigen::MatrixXd A = covariance.sqrt();

    //set first column of sigma point matrix
    col(0) = state;

    //set remaining sigma points
    int n_aug = size();
    double t = std::sqrt(lambda + n_aug);  // TODO: 7  - n_x or n_aug

    for (int i = 0; i < n_aug; i++) {
      col(i + 1) = state + t * A.col(i);
      col(i + 1 + n_aug) = state - t * A.col(i);
    }
  }
};


#endif // UNSCENTEDKF_AUGMENTED_H
