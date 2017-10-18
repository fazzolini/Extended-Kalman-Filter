#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
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
  // also initialize ekf_ matrices
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);


  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO: [DONE]
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  // H matrices:

  // this is the same as original H from classroom project
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // new non-linear Hj_ Jacobian
  /**
   * initialize with ones, will be updated later for each step
   * using Tools::CalculateJacobian()
   */
  Hj_ <<  1, 1, 0, 0,
          1, 1, 0, 0,
          1, 1, 1, 1;

  // Kalman Filter state transition matrix
  // same as in classroom exercise
  ekf_.F_ <<  1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;

  // Kalman Filter state covariance matrix
  // same as in clsssroom exercise
  ekf_.P_ <<  1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;

  //set the acceleration noise components
  noise_ax = 9;
  noise_ay = 9;

  // initial G_
  ekf_.G_ = MatrixXd(4,2);
  ekf_.G_ <<  1, 0,
              0, 1,
              1, 0,
              0,1;

  // initial Qv_ - will stay constant throughout
  ekf_.Qv_ = MatrixXd(2,2);
  ekf_.Qv_ << noise_ax, 0,
              0, noise_ay;

  // also: initialize our tools
  Tools tools;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO: [DONE]
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      // we need to extract all values from measurement package
      float rho = measurement_pack.raw_measurements_(0); // distance
      float phi_raw = measurement_pack.raw_measurements_(1); // angle [in radians] - can be outside range
      float phi = tools.NormalizePhi(phi_raw); // make sure phi is in range from -pi to +pi
      float delta_rho = measurement_pack.raw_measurements_(2); // speed in the direction of rho

      // put values in the state
      ekf_.x_(0) = rho * cos(phi); // distance projection on vertical axis x
      ekf_.x_(1) = rho * sin(phi); // distance projection on horizontal axis y
      ekf_.x_(2) = delta_rho * cos(phi); // speed projection on vertical axis x
      ekf_.x_(3) = delta_rho * sin(phi); // speed projection on horizontal axis y
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      // easy case: just update coordinates
      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
    }

    // update previous timestamp
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO: [DONE]
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // 1. compute the time elapsed between the current and previous measurements
  float dt;
  dt = float((measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0);  //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  // 2. modify the F matrix so that the time is integrated
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  // 3. set the process covariance matrix Q
  //    first, set G
  ekf_.G_(0,0) = dt * dt / 2;
  ekf_.G_(1,1) = dt * dt / 2;
  ekf_.G_(2,0) = dt;
  ekf_.G_(3,1) = dt;
  // secondly, calculate new Q using G and Qv (introduced in classroom example)
  ekf_.Q_ = ekf_.G_ * ekf_.Qv_ * ekf_.G_.transpose();

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO: [DONE]
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // 1. calculate Jacobian
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    // 2. update H and R matrices of the state
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    // 3. do the radar update step
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    // 1. update H and R matrices of the state
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    // 3. do the laser update step
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
