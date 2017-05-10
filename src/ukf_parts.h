#ifndef UNSCENTEDKF_PARTS_H
#define UNSCENTEDKF_PARTS_H

#include "Eigen/Dense"

#include "ukf_state.h"
#include "ukf_augmented.h"

template<class TargetSigmaPoints, class TargetState, class TargetStateCovariance>
class SpaceTransformation
{
protected:

  int number_of_points_;
  Eigen::VectorXd weights_;

  SpaceTransformation(double lambda, int space_dimension) :
    number_of_points_(SpaceBase::dimension_to_points(space_dimension)),
    weights_(number_of_points_)
  {
    double d = lambda + space_dimension; // TODO: double check n_aug or n_x ?
    weights_.fill(0.5 / d);
    weights_(0) = lambda / d;
  }

public:

  TargetState CalculateWeightedMean(const TargetSigmaPoints& sigma_points) {
    TargetState result;
    for (int i = 0; i < number_of_points_; i++) {  //iterate over sigma points
      result.set_raw(result.raw() + weights_(i) * sigma_points.raw_point(i));
    }
    return result;
  }

  TargetStateCovariance CalculateCovariance(
    const TargetSigmaPoints& sigma_points,
    const TargetState& mean) 
  {
    TargetStateCovariance result;
    for (int i = 0; i < number_of_points_; i++) {  // iterate over sigma points
                                                   // state difference
      Eigen::VectorXd x_diff = sigma_points.diff_from_mean(i, mean);
      result.set_data(result.data() + weights_(i) * x_diff * x_diff.transpose()); 
    }
    return result;
  }
};


class Predictor : public SpaceTransformation<StateSigmaPoints, State, StateCovariance> {

public:

  Predictor(double lambda) :
    SpaceTransformation(lambda, 7) // TODO: 7
  {
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
    }
    else {
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

  StateSigmaPoints LoadPoints(const AugmentedSpaceSigmaPoints& Xsig_aug, double delta_t) {
    StateSigmaPoints result;
    for (int i = 0; i < Xsig_aug.number_of_points(); i++) {
      State predicted = PredictNextState(Xsig_aug.point(i), delta_t);
      result.set_raw_point(i, predicted.raw());
    }
    return result;
  }
};


#endif //UNSCENTEDKF_PARTS_H
