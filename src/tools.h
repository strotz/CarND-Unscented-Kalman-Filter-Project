#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include <iostream>
#include "Eigen/Dense"
#include <ctime>

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

};

class dl
{
public:

  dl(const char* place) :
    place_(place),
    start_(clock())
  {
  }

  ~dl() {
    double elapsed_secs = double(clock() - start_) / CLOCKS_PER_SEC;
    if (elapsed_secs > 1.0) {
      std::cout << place_ << " " << elapsed_secs << std::endl;
    }
  }

  const char* place_;
  clock_t start_;
};

#endif /* TOOLS_H_ */
