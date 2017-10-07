#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  use_laser_ = true;   // if this is false, laser measurements will be ignored (except during init)
  use_radar_ = true;   // if this is false, radar measurements will be ignored (except during init)

  x_ = VectorXd(5);    // initial state vector
  P_ = MatrixXd(5, 5); // initial covariance matrix

  std_a_ = 1; // Process noise standard deviation longitudinal acceleration in m/s^2
  std_yawdd_ = 1; // Process noise standard deviation yaw acceleration in rad/s^2
  std_laspx_ = 0.15; // Laser measurement noise standard deviation position1 in m
  std_laspy_ = 0.15; // Laser measurement noise standard deviation position2 in m
  std_radr_ = 0.3; // Radar measurement noise standard deviation radius in m
  std_radphi_ = 0.03;   // Radar measurement noise standard deviation angle in rad
  std_radrd_ = 0.3; // Radar measurement noise standard deviation radius change in m/s

  // Initialize NISs
  NIS_L_ = 0.0;
  NIS_R_ = 0.0;

  /**
  Complete the initialization. See ukf.h for other member properties.
  */

  x_ << 1, 1, 0.5, 0.5, 0.5;
  P_ << 0.15, 0, 0, 0, 0,
        0, 0.15, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  is_initialized_ = false;
  use_laser_ = false;
  use_radar_ = false;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3.0 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // initialize weights
  weights_ = VectorXd(2*n_aug_+1);
  for (int i=1; i <= (2*n_aug_+1); i++)
  {
     if (i == 1)
     {
         weights_(i-1) = lambda_/(lambda_+n_aug_);
     }
     else
     {
         weights_(i-1) = 1.0/(2*(lambda_+n_aug_));
     }
  }
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package
 * The latest measurement data of either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /** switch between lidar and radar measurements.  */
  if (!is_initialized_){
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      x_(0) = meas_package.raw_measurements_(0) * cos(meas_package.raw_measurements_(1));
      x_(1) = meas_package.raw_measurements_(0) * sin(meas_package.raw_measurements_(1));
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
  }  
  // Get dt
  float dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package); //LIDAR Update
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package); //RADAR Update
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //SigmaPointsPrediction(&Xsig_pred_, delta_t);
  //PredictMeanCovariance();

  /***********************************************************
  * Generate Augumented Sigma Points
  ************************************************************/

  // Initialize augumented sigma points
  MatrixXd P_aug(n_aug_, n_aug_);
  MatrixXd Xsig_aug(n_aug_, 2*n_aug_+1);
  VectorXd x_aug(n_aug_);

  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  Xsig_aug.fill(0);
  Xsig_aug.col(0) = x_aug;

  // calculate the augmented sigma points
  double coef = sqrt(lambda_+n_aug_);
  MatrixXd A = P_aug.llt().matrixL();

  for (int i = 1; i <= n_aug_ ;i++)
  {
    Xsig_aug.col(i) = x_aug + coef * A.col(i-1);
    Xsig_aug.col(i+n_aug_) = x_aug - coef*A.col(i-1);
  }

  /***********************************************************
  * Sigma Points Prediction
  ************************************************************/
  for (int i = 0; i < (2 * n_aug_+ 1); i++)
  {
	double px_p, py_p, v_p, yaw_p, yawd_p;

    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);

    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    if (fabs(yawd) <= 0.001)
    {
        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);
    }
    else
    {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }

    v_p = v;
    yaw_p = yaw + yawd*delta_t;
    yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5* nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // cout << Xsig_pred_ << endl;
  
  /***********************************************************
  * Predicted Sigma Points Mean and Covariance
  ************************************************************/
  VectorXd x_temp(5);
  x_temp.fill(0.0);
  for (int i=0; i< 2 * n_aug_ + 1; i++){
	  x_temp = x_temp + weights_(i) * Xsig_pred_.col(i);
  }
  x_ = x_temp;

  //predict state covariance matrix
  MatrixXd P_pred(n_x_, n_x_);
  P_pred.fill(0);
  for (int i=0; i<2*n_aug_+1; i++){
      VectorXd x_diff = Xsig_pred_.col(i)-x_;
      while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

      MatrixXd p_temp = (x_diff * x_diff.transpose()) * weights_(i);
      P_pred += p_temp;
  }
  P_ = P_pred;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the lidar NIS.
  */
  /***********************************************************
  * Get measurements
  ************************************************************/

  VectorXd z(2);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1); 


  /*************************************************************
  * Predicted Sigma Points transform to measurement state: LADAR
  **************************************************************/

  MatrixXd Zsig(2, 2*n_aug_+1);
  for (int i=0; i < (2*n_aug_+1); i++)
  {
    Zsig(0,i) = Xsig_pred_(0,i);
    Zsig(1,i) = Xsig_pred_(1,i);
  }

  VectorXd z_pred(2);
  z_pred.fill(0);
  for (int i=0; i < (2*n_aug_+1); i++)
  {
    z_pred = z_pred + Zsig.col(i) * weights_(i);
  }
  
  /*************************************************************
  * Update measurement
  **************************************************************/

  MatrixXd S(2,2);
  S.fill(0);
  for (int i=0; i < (2*n_aug_+1); i++)
  {
    VectorXd diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * diff * diff.transpose();
  }

  MatrixXd R(2,2);
  R << std_laspx_*std_laspx_, 0, 0, std_laspy_*std_laspy_;
  S = S + R;

  // Update State
  int n_z = z_pred.size();
  MatrixXd T(n_x_, n_z);
  T.fill(0);
  for (int i=0; i < (2*n_aug_+1); i++)
  {
    VectorXd x_dif = Xsig_pred_.col(i) - x_;
    VectorXd z_dif = Zsig.col(i) - z_pred;
    T = T + weights_(i) * x_dif * z_dif.transpose();
  }

  MatrixXd K = T * S.inverse();
  VectorXd z_dif = z - z_pred;
  x_ = x_ + K * z_dif;
  P_ = P_ - K * S * K.transpose();

  /*************************************************************
  * NIS for Lidar
  **************************************************************/
  NIS_L_ = z_dif.transpose() * S.inverse() * z_dif;
  cout << "Lidar NIS: " << NIS_L_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  /***********************************************************
  * Get measurements
  ************************************************************/
  VectorXd z(3);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  z(2) = meas_package.raw_measurements_(2);

  /*************************************************************
  * Predicted Sigma Points transform to measurement state: RADAR
  **************************************************************/
  MatrixXd Zsig(3, 2*n_aug_+1);
  for (int i=0; i < (2*n_aug_+1); i++)
  {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double yawd = Xsig_pred_(4,i);

    double rho = sqrt(p_x*p_x + p_y*p_y);
    double phi = atan2(p_y, p_x);
    double rho_dot = ((p_x * cos(yaw) * v) + (p_y * sin(yaw) * v)) / rho;

    Zsig(0,i) = rho;
    Zsig(1,i) = phi; 
    Zsig(2,i) = rho_dot;
  }

  VectorXd z_pred(3);
  z_pred.fill(0);
  for (int i=0; i < (2*n_aug_+1); i++)
  {
    z_pred = z_pred + Zsig.col(i) * weights_(i);
  }

  /*************************************************************
  * Update measurement
  **************************************************************/
  MatrixXd S(3,3);
  S.fill(0);
  for (int i=0; i < (2*n_aug_+1); i++)
  {
    VectorXd dif = Zsig.col(i) - z_pred;
    while (dif(1)> M_PI) dif(1)-=2.*M_PI;
    while (dif(1)<-M_PI) dif(1)+=2.*M_PI;
    S = S + weights_(i) * dif * dif.transpose();
  }

  MatrixXd R(3,3);
  R << std_radr_ * std_radr_, 0, 0,
       0, std_radphi_ * std_radphi_, 0,
       0, 0, std_radrd_ * std_radrd_;
  S = S + R;

  int n_z = z_pred.size();
  MatrixXd T1(n_x_, n_z);
  T1.fill(0);

  for (int i=0; i < (2*n_aug_+1); i++){
    VectorXd x_dif = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_dif(3)> M_PI) x_dif(3)-=2.*M_PI;
    while (x_dif(3)<-M_PI) x_dif(3)+=2.*M_PI;
    
    VectorXd z_dif = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_dif(1)> M_PI) z_dif(1)-=2.*M_PI;
    while (z_dif(1)<-M_PI) z_dif(1)+=2.*M_PI;

    T1 = T1 + weights_(i) * x_dif * z_dif.transpose();
  }

  MatrixXd K = T1 * S.inverse();

  VectorXd z_dif = z - z_pred;
  // angle normalization
  while (z_dif(1)> M_PI) z_dif(1) -= 2.*M_PI;
  while (z_dif(1)<-M_PI) z_dif(1) += 2.*M_PI;

  x_ = x_ + K * z_dif;
  P_ = P_ - K * S * K.transpose();

  /*************************************************************
  * NIS for RADAR
  **************************************************************/
  NIS_R_ = z_dif.transpose() * S.inverse() * z_dif;
  cout << "Radar NIS:" << NIS_R_ << endl;
}
