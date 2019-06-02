#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;


/**
 * Constructor.
 */
FusionEKF::FusionEKF() {

  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ << 1.0, 0, 0, 0,
              0, 1.0, 0, 0;

  // state vector
  VectorXd x = VectorXd(4);
  x << 1, 1, 1, 1;

  // set the process and measurement noises
  MatrixXd P = MatrixXd(4, 4);
  P << 1, 0,    0,    0,
       0, 1,    0,    0,
       0, 0, 1000,    0,
       0, 0,    0, 1000;

  MatrixXd F = MatrixXd(4, 4);
  F << 1, 0, 1, 0,
       0, 1, 0, 1,
       0, 0, 1, 0,
       0, 0, 0, 1;

  ekf_.Init(x, P, F);

}


/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

/**
 * Process incoming measurement data from multiple sensors.
 */
void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /**
   * Initialization
   */
  if (!is_initialized_) {

    VectorXd &x = ekf_.GetState();

    // first measurement
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

      // convert radar from polar to Cartesian coordinates and initialize state.
      x = Tools::ConvertPolarToCartesian(measurement_pack.raw_measurements_);

    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

      // initialize state.
      x(0) = measurement_pack.raw_measurements_(0);
      x(1) = measurement_pack.raw_measurements_(1);

    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * Time is measured in seconds.
   * Use of noise_ax = 9 and noise_ay = 9 for Q matrix.
   */

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.Predict(dt, 9, 9);

  /**
   * Update
   */

  /**
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

    // Radar updates
    Hj_ = Tools::CalculateJacobian(ekf_.GetState());
    ekf_.UpdateEKF(measurement_pack.raw_measurements_, Hj_, R_radar_);

  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_, H_laser_, R_laser_);

  }

  // print the output
  cout << "x_ = " << ekf_.GetState() << endl;
  cout << "P_ = " << ekf_.GetStateCovariance() << endl;
}

/**
 * Returns current state vector.
 */
VectorXd &FusionEKF::GetState() {

  return ekf_.GetState();

}
