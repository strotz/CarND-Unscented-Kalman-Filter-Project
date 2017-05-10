#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse.fill(0.0);

  if (estimations.size() == 0) {
    cout << "CalculateRMSE () - Error - the estimation vector size should not be zero" << endl;
    return rmse;
  }

  if (estimations.size() != ground_truth.size()) {
    cout << "CalculateRMSE () - Error - the estimation vector size should equal ground truth vector size" << endl;
    return rmse;
  }

  //accumulate squared residuals
  for (size_t i = 0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    rmse = rmse.array() + residual.array() * residual.array();
  }

  //calculate the mean
  rmse = rmse / estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}
