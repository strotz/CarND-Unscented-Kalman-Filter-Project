#ifndef UNSCENTEDKF_PARTS_H
#define UNSCENTEDKF_PARTS_H

#include "ukf_state.h"
#include "ukf_augmented.h"

class Weights : public Eigen::VectorXd
{
  int number_of_points_;

public:

  Weights(double lambda, int space_dimension) :
    number_of_points_(SpaceBase::dimension_to_points(space_dimension)),
    Eigen::VectorXd(SpaceBase::dimension_to_points(space_dimension))
  {
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
protected:

  Weights weights_;

  SpaceTransformation(const Weights& weights) :
    weights_(weights)
  {
  }

public:

  TargetState CalculateWeightedMean(const TargetSigmaPoints& sigma_points) {
    TargetState result;
    int len = weights_.number_of_points();
    for (int i = 0; i < len; i++) { //iterate over sigma points
      result = result + weights_(i) * sigma_points.col(i);
    }
    return result;
  }

  TargetStateCovariance CalculateCovariance(
    const TargetSigmaPoints& sigma_points,
    const TargetState& mean) 
  {
    TargetStateCovariance result;
    int len = weights_.number_of_points();
    for (int i = 0; i < len; i++) {  // iterate over sigma points
      // state difference
      Eigen::VectorXd x_diff = sigma_points.diff_from_mean(i, mean);
      result = result + weights_(i) * x_diff * x_diff.transpose();
    }
    return result;
  }
};


class PositionPredictor : public SpaceTransformation<StateSigmaPoints, State, StateCovariance> {

public:

  PositionPredictor(const Weights& weights) :
    SpaceTransformation(weights)
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
    State state;

    StateOps(state)
      .set_pos_x(px_p)
      .set_pos_y(py_p)
      .set_velocity(v_p)
      .set_yaw_angle(yaw_p) // TODO: normalize ?
      .set_yaw_rate(yawd_p);

    return state;
  }

public:

  StateSigmaPoints LoadPoints(const AugmentedSpaceSigmaPoints& Xsig_aug, double delta_t) {
    int number_of_points = Xsig_aug.number_of_points();
    StateSigmaPoints result(number_of_points);
    for (int i = 0; i < number_of_points; i++) {
      Eigen::VectorXd t = Xsig_aug.col(i);
      AugmentedStateOps ops(t);
      result.col(i) = PredictNextState(t, delta_t);
    }
    return result;
  }
};


#endif //UNSCENTEDKF_PARTS_H
