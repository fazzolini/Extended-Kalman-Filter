#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO: [DONE][OK]
    * predict the state
  */
  // the same as in in-class exercise
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO: [DONE][OK]
    * update the state by using Kalman Filter equations
  */
  // the same as in in-class exercise
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z_measured) {
  /**
  TODO: [DONE]
    * update the state by using Extended Kalman Filter equations
  */
  // NOTE: z_measured is incoming current radar measurement (rho, phi, delta_rho)
  // need to first compute (rho, phi, delta_rho) from the current state (px, py, vx, vy)
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  float rho = sqrtf(px * px + py * py);
  float phi = atanf(py / px);
  float delta_rho; // need to make sure that denominator is not zero
  if (fabsf(rho) < 0.0001) {
    delta_rho = 0;
  } else {
    delta_rho = (px * vx + py * vy) / rho;
  }

  // create a vector of predictions in polar coordinate format
  VectorXd z_predicted(3);
  z_predicted << rho, phi, delta_rho;

  // next lines are similar to the standard linear Kalman Filter
  VectorXd y = z_measured - z_predicted;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate is still the same as linear Kalman Filter
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
