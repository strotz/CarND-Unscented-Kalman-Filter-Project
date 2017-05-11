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
class AugmentedState : public StateBase<AugmentedStateOps> {

public:
  AugmentedState(const State &state) :
    StateBase(7) // TODO: 7
  { 
    //create augmented mean state
    set_head(state.size(), state.raw());
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
    data_.topLeftCorner(5, 5) = covariance.data();
    data_(5, 5) = std_a * std_a;
    data_(6, 6) = std_yawdd * std_yawdd;
  }

  Eigen::MatrixXd sqrt() const {
    return data_.llt().matrixL();
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
    set_raw_point(0, state.raw());

    //set remaining sigma points
    int n_aug = size();
    double t = sqrt(lambda + n_aug);  // TODO: 7  - n_x or n_aug

    for (int i = 0; i < n_aug; i++) {
      set_raw_point(i + 1, state.raw() + t * A.col(i));
      set_raw_point(i + 1 + n_aug, state.raw() - t * A.col(i));
    }
  }
};


#endif // UNSCENTEDKF_AUGMENTED_H
