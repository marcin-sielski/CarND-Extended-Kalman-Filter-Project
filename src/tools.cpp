#include "tools.h"
#include <iostream>


using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;


/**
 * Calculate root mean square error
 */
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    cout << __func__ <<
        ": Invalid estimation or ground_truth data, returning Zero Matrix." <<
        endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;

}


/**
 * Calculate Jacobian
 */
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3,4);
  Hj << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // pre-compute a set of terms to avoid repeated calculation
  double c1 = px*px+py*py;
  double c2 = sqrt(c1);
  double c3 = (c1*c2);

  // check division by zero
  if (fabs(c1) < 0.0001) {
    cout << __func__ << ": Division by 0, returning Zero Matrix." << endl;
    return Hj;
  }

  // compute the Jacobian matrix
  Hj <<                px/c2,                 py/c2,     0,     0,
                      -py/c1,                 px/c1,     0,     0,
       py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;

}


/*
 * Convert Polar coordinates to Cartesian coordinates
 */
VectorXd Tools::ConvertPolarToCartesian(const VectorXd &x_state_polar) {

  double rho = x_state_polar(0);
  double theta = x_state_polar(1);
  double rho_dot = x_state_polar(2);

  VectorXd x_state(4);
  x_state << rho*cos(theta), rho*sin(theta), rho_dot*cos(theta),
      rho_dot*sin(theta);

  return x_state;

}


/**
 * Convert Cartesian coordinates to polar coordinates
 */
VectorXd Tools::ConvertCartesionToPolar(const VectorXd &x_state) {

  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  VectorXd x_state_polar(3);

  double rho = sqrt(px*px+py*py);
  x_state_polar << rho, atan2(py, px),
      (fabs(rho) < 0.0001 ? 0 : (px*vx+py*vy)/rho);

  return x_state_polar;

}
