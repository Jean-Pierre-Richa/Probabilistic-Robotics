function A = v2t(v)
		A = [cos(v(3)), -sin(v(3)), v(1);
				 sin(v(3)),  cos(v(3)), v(2);
				 0, 				 0, 				1];
end

function v = t2v(A)
		v(1:2, 1) = A(1:2, 3);
		v(3, 1)   = atan2(A(2,1), A(1,1));
end


% If we have an encoder that gives us N clicks per revolution then the distance we travel per click (C = diameter/num_clicks)
% Where diameter_l and diameter_r are the Diameters of the left and right wheels respectively, num_clicks_l and num_clicks_r are the number
% of clicks passed on the left and right wheels respectively and Cl, and Cr are the calibrations of the left and right
% wheels respectively
% If we take the number of clicks and multiply it by the calibration of the wheel, we get the distance traveled (distance)
% So our equation becomes (distance = C*num_clicks) for distance_l, Cl, num_clicks_l and distance_r, Cr, and num_clicks_r.
% after calculating distance_l and distance_r, we can treat our model as a unicycle and
	% calculate the distance traveled by the midpoint on the rear axis (delta_s = (distance_l+distance_r)/2)
	% calculate the rotation of the robot with (delta_theta=(distance_l-distance_r)/baseline) where baseline is the baseline between the contact
	% points with the encoders on the 2 wheels
% We can neglect the above part for now, because we do not have these quantities

function h = prediction(params, data)
% As given by the slides, the right and left encoder to meters factors kl = -1, kr = 1, and the baseline = 0.3;

h = zeros(1,3);

% we calculate the left wheel motion by
Wl = params(1)*data(1); % where data(1) is the left encoder ticks
% same for the right wheel
Wr = params(2)*data(2);
% baseline
b = params(3);

% Angular displacement
delta_theta = (Wr - Wl)/b;

p_x = 1 - delta_theta^2/6 + delta_theta^4/120 - delta_theta^6/5040;
p_y = delta_theta/2 - delta_theta^3/24 + delta_theta^5/720;

% X displacement
delta_x = ((Wr + Wl)/2)*p_x;

% Y displacement
delta_y = ((Wr + Wl)/2)*p_y;

h(1) = delta_x;
h(2) = delta_y;
h(3) = delta_theta;

end

function T=compute_ground_truth_trajectory(U)
	T=zeros(size(U,1),3); % Matrix in the rows size of U(Z) and 3 columns
	P=v2t(zeros(1,3));   	 % 1x3 row vector
	for i=1:size(U,1),   	 % for i from 1 to N rows in U
		u=U(i,1:3)'; 				 % u = the ith row and the vector from the 3rd column till the end, which contains the ground truth vector
		P = P*v2t(u); 		   % here point p is being transformed into a homogeneous
											   % transformation matrix, using the rotation matrix, the
											   % values (x, y) from the 3rd ad the 4th columns and the angle theta from the 5th columns
		T(i,1:3)=t2v(P)';   % then transforming each row into a vector to get the point "P" location.
	end
end


#computes the trajectory of the robot by chaining up
#the incremental movements of the odometry vector
#U:	a Nx3 matrix, each row contains the odoemtry ux, uy utheta
#T:	a Nx3 matrix, each row contains the robot position (starting from 0,0,0)
function T=compute_odometry_trajectory(params, U)
	T=zeros(size(U,1),3);
	P=v2t(zeros(1,3));
	for i=1:size(U,1),
		u=U(i,1:2)';
		h = prediction(params, u);
		P = P*v2t(h);
		T(i,1:3)=t2v(P)';
	end
end


#computes a calibrated vector of odometry measurements
#by applying the bias term to each line of the measurements
#X: 	3x3 matrix obtained by the calibration process
#U: 	Nx3 matrix containing the odometry measurements
#C:	Nx3 matrix containing the corrected odometry measurements

%function C=apply_odometry_correction(X, U)
%	C=zeros(size(U,1),3);
%	for i=1:size(U,1),
%		u=U(i,1:3)'; %'
%		uc = X*u;
%		C(i,:)=uc;
%	end
%end


#the measurements are coming in the
#Z matrix

#this function solves the odometry calibration problem
#given a measurement matrix Z.
#Every row of the matrix contains
#z_i = [tl, tr, x, y, theta]
#Z:	The measurement matrix
#X:	the calibration matrix
#returns the bias correction matrix BIAS
function X=ls_calibrate_odometry(params, Z)
	#accumulator variables for the linear system
	H=zeros(3,3);
	b=zeros(3,1);
	#initial solution (the identity transformation)
	X=params;

	#loop through the measurements and update the
	#accumulators
	for i=1:size(Z,1),
		e=error_function(params,i,Z,X);
		A=jacobian(params,i,Z);
		H = H+A'*A;
		b = b+A'*e';
	end
	#solve the linear system
	deltaX=-H\b;
	#this reshapes the 9x1 increment vector in a 3x3 matrix
	%dX=reshape(deltaX,3,3)';%'
	#computes the cumulative solution
	X=X+deltaX';
end

#this function computes the error of the i^th measurement in Z
#given the calibration parameters
#i:	the number of the measurement
#X:	the actual calibration parameters
#Z:	the measurement matrix
#e:	the error of the ith measurement
function e=error_function(params,i,Z,X)
	uprime=Z(i,1:2);
	h=Z(i,3:5);
	u = prediction(params, uprime);
	e = h-u;
end

#derivative of the error function for the ith measurement in Z
#does not depend on the state
#i:	the measurement number
#Z:	the measurement matrix
#A:	the jacobian of the ith measurement
function A=jacobian(params,i,Z)

	noise = 0.05;
	u=Z(i,1:2);
	A=zeros(3,3);

	kl = params(1);
	kr = params(2);
	b = params(3);

	% Numerical jacobian:
	pred_kl = prediction([kl+noise, kr, b],u) - prediction([kl-noise, kr, b],u);
	pred_kr = prediction([kl, kr+noise, b],u) - prediction([kl, kr-noise, b],u);
	pred_b  = prediction([kl, kr, b+noise],u) - prediction([kl, kr, b-noise],u);
	A = -[pred_kl; pred_kr; pred_b]'/(2*noise);
	printf("A");
	disp(A);
end
