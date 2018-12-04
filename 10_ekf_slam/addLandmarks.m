function [mu, sigma, id_to_state_map, state_to_id_map] = addLandmarks(mu, sigma, measurements, id_to_state_map, state_to_id_map)

#determine how many landmarks we have seen in this step
new_measured_l = length(measurements.observation);
state_dim = size(mu,1);
new_landmarks_nb = (state_dim-3)/2;

#since I have applied the correction, I need to update my
#data with the new mu values
mu_t     = mu(1:2,:); #translational part of the robot pose
mu_theta = mu(3);     #rotation of the robot

#re precompute some quantities that come in handy later on
c   = cos(mu_theta);
s   = sin(mu_theta);
R   = [c -s; s c];  #rotation matrix
Rt  = [c,s;-s c];    #transposed rotation matrix
Rtp = [-s,c;-c,-s]; #derivative of transposed rotation matrix

#Now its time to add, if observed, the NEW landmaks, without applying any correction
for i=1:new_measured_l

  #retrieve info about the observed landmark
  measurement = measurements.observation(i);

  #fetch the position in the state vector corresponding to the actual measurement
  %spol: state_pos_of_landmark
  spol = id_to_state_map(measurement.id);

  #IF current landmark is a NEW landmark
  if(spol == -1)

    #adjust direct and reverse mappings
    %increase the number of the state vector accordingly
    new_landmarks_nb++;
    spol=new_landmarks_nb;
    state_to_id_map(new_landmarks_nb)=measurement.id;

    #landmark position in the world
    %lpw: land_pose_world = #TODO: compute the landmark position in the world (using the current robot position: mu_t, R)
    lpw = mu_t + R*[measurement.x_pose;measurement.y_pose];

    #retrieve from the index the position of the landmark block in the state
    %nlsvi = new_landmark_state_vector_index
    nlsvi=4+2*(new_landmarks_nb-1);

    #increase mu and sigma size
    mu(nlsvi:nlsvi+1,1) = lpw;

    #initial noise assigned to a new landmark
    #for simplicity we put a high value only in the diagonal.
    #A more deeper analysis on the initial noise should be made.
    initial_landmark_noise=2;
    landmark_sigma = eye(2)*initial_landmark_noise;

    #extend the structure
    sigma(nlsvi,:)   = 0;
    sigma(nlsvi+1,:) = 0;
    sigma(:,nlsvi)   = 0;
    sigma(:,nlsvi+1) = 0;

    #add the covariance block
    sigma(nlsvi:nlsvi+1, nlsvi:nlsvi+1)=landmark_sigma;
    if new_landmarks_nb>0
      new_l = landmark(state_to_id_map(1), [mu(4), mu(5)]);
      for u = 2:new_landmarks_nb
        new_l(end+1) = landmark(state_to_id_map(u), [mu(3+2*u-1), mu(3+2*u)]);
  endfor
  endif
    %printf("observed new landmark with identifier: %i \n",measurement.id);
    %fflush(stdout);
  endif

end
