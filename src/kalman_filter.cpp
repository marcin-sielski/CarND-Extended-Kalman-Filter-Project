#include "kalman_filter.h"
#include "tools.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;


/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */
/**
 * Constructor.
 */
KalmanFilter::KalmanFilter() {}


/**
 * Destructor.
 */
KalmanFilter::~KalmanFilter() {}


/**
 * Initialize Kalman Filter
 */
void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in) {

  x_ = x_in;
  P_ = P_in;
  F_ = F_in;

}

/**
 * Returns current state vector
 */
VectorXd &KalmanFilter::GetState() {

  return x_;

}

/**
 * Returns current state covariance matrix
 */
MatrixXd &KalmanFilter::GetStateCovariance() {

  return P_;

}

/**
 * Predict state
 */
void KalmanFilter::Predict(double dt, double noise_ax, double noise_ay) {

  // predict the state

  // update the state transition matrix F according to the new elapsed time.

  F_(0, 2) = dt;
  F_(1, 3) = dt;

  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  // set the process covariance matrix Q
  Q_ = MatrixXd(4, 4);

  // update the process noise covariance matrix.
  Q_ << dt_4/4*noise_ax,               0, dt_3/2*noise_ax,               0,
                      0, dt_4/4*noise_ay,               0, dt_3/2*noise_ay,
        dt_3/2*noise_ax,               0,   dt_2*noise_ax,               0,
                      0, dt_3/2*noise_ay,               0,   dt_2*noise_ay;
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;

}


/**
 * Update state with KF step
 */
void KalmanFilter::Update(const VectorXd &z, const MatrixXd &H,
                          const MatrixXd &R) {

  // update the state by using Kalman Filter equations

  VectorXd y = z - H * x_;

  UpdateState_(y, H, R);

}


/**
 * Update state with EKF step
 */
void KalmanFilter::UpdateEKF(const VectorXd &z, const MatrixXd &H,
                             const MatrixXd &R) {

  // update the state by using Extended Kalman Filter equations

  VectorXd y = z - Tools::ConvertCartesionToPolar(x_);

  // normalize
  while(y(1) > M_PI) y(1) -= 2*M_PI;
  while(y(1) < -M_PI) y(1) += 2*M_PI;

  UpdateState_(y, H, R);

}


/**
 * Helper method which include common equations for KF and EKF update step
 */
void KalmanFilter::UpdateState_(const VectorXd &y, const MatrixXd &H,
                               const MatrixXd &R) {

  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd K = (P_ * Ht) * S.inverse();

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;

}
