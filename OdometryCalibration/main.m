#import 2d geometry utils
source "../tools/utilities/geometry_helpers_2d.m"

source "odometryCalibration.m"

params = [-1, 1, 0.3];

h = figure();

#load the calibration matrix
disp('loading the calibration matrix');
Z=load("../datasets/differential_drive_calibration_data/differential_drive_calibration.txt");

#compute the ground truth trajectory
TrueTrajectory=compute_ground_truth_trajectory(Z(:,3:5));
disp('ground truth red');
hold on;
plot(TrueTrajectory(:,1),TrueTrajectory(:,2), 'r-', 'linewidth', 2);
pause(1);

#compute the uncalibrated odometry
OdomTrajectory=compute_odometry_trajectory(params, Z(:,1:2));
disp('odometry green');
hold on;
plot(OdomTrajectory(:,1),OdomTrajectory(:,2), 'g-', 'linewidth', 2);
pause(1);

disp('computing calibration parameters');
#compute the calibration parameters
X=ls_calibrate_odometry(params, Z);
disp(X);
pause(1);

disp('computing calibrated odometry');
%COdom=apply_odometry_correction(X,Z(:,1:2));
CalTrajectory=compute_odometry_trajectory(X, Z(:,1:2));
hold on;
plot(CalTrajectory(:,1),CalTrajectory(:,2), 'b-', 'linewidth', 2);
waitfor(h);
