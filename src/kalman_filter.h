#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_


#include "Eigen/Dense"


class KalmanFilter {

 public:

  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param x_in Initial state
   * @param P_in Initial state covariance
   * @param F_in Transition matrix
   */
  void Init(Eigen::VectorXd &x_in, Eigen::MatrixXd &P_in,
            Eigen::MatrixXd &F_in);

  /**
   * Get current state vector
   * @return current state vector
   */
  Eigen::VectorXd &GetState();

  /**
   * Get current state covariance matrix
   * @return current state covariance matrix
   */
  Eigen::MatrixXd &GetStateCovariance();

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param dt Time between k and k+1 in s
   * @param noise_ax acceleration noise x
   * @param noise_ay acceleration noise y
   */
  void Predict(double dt, double noise_ax, double noise_ay);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   * @param H Measurement matrix
   * @param R Measurement covariance matrix
   */
  void Update(const Eigen::VectorXd &z, const Eigen::MatrixXd &H,
              const Eigen::MatrixXd &R);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   * @param H Measurement matrix
   * @param R Measurement covariance matrix
   */
  void UpdateEKF(const Eigen::VectorXd &z, const Eigen::MatrixXd &H,
                 const Eigen::MatrixXd &R);

 private:

  /**
   * Helper function - common equations for KF and EKF
   * @param H Measurement matrix
   * @param R Measurement covariance matrix
   */
  void UpdateState_(const Eigen::VectorXd &y, const Eigen::MatrixXd &H,
                   const Eigen::MatrixXd &R);

  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // state transition matrix
  Eigen::MatrixXd F_;

  // process covariance matrix
  Eigen::MatrixXd Q_;

};

#endif // KALMAN_FILTER_H_
