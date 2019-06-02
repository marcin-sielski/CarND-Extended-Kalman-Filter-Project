#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

class Tools {

 public:

  /**
   * A helper method to calculate RMSE
   * @param estimations estimations
   * @param ground_truth ground truth
   * @return root mean square error
   */
  static Eigen::VectorXd CalculateRMSE(
      const std::vector<Eigen::VectorXd> &estimations,
      const std::vector<Eigen::VectorXd> &ground_truth);

  /**
   * A helper method to calculate Jacobians
   * @param x_state state vector
   * @return Jacobian
   */
  static Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

  /**
   * A helper method to convert Polar coordinates to Cartesian coordinates
   * @param x_state_polar Polar coordinates
   * @return Cartesian coordinates
   */
  static Eigen::VectorXd ConvertPolarToCartesian(
      const Eigen::VectorXd &x_state_polar);

  /**
   * A helper method to convert Cartesian coordinates to Polar coordinates
   * @param x_state Cartesian coordinates
   * @return Polar coordinates
   */
  static Eigen::VectorXd ConvertCartesionToPolar(
      const Eigen::VectorXd &x_state);

};

#endif  // TOOLS_H_
