#ifndef UNSCENTEDKF_PARTS_H
#define UNSCENTEDKF_PARTS_H

#include "ukf_state.h"
#include "ukf_augmented.h"


class PositionPredictor : public SpaceTransformation<StateSigmaPoints, State, StateCovariance> {

public:

  explicit PositionPredictor() :
    SpaceTransformation()
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

    yaw_p = SpaceBase::normalize_angle(yaw_p);

    StateOps(state)
      .set_pos_x(px_p)
      .set_pos_y(py_p)
      .set_velocity(v_p)
      .set_yaw_angle(yaw_p)
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
