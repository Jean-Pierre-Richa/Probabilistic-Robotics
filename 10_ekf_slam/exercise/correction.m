#this function computes the update (also called correction)
#step of the filter
#inputs:
#  mu: mean,
#  sigma: covariance of the robot-landmark set (x,y, theta, l_1, ..., l_N)
#
#  observations:
#            a structure containing n observations of landmarks
#            for each observation we have
#            - the index of the landmark seen
#            - the location where we have seen the landmark (x,y) w.r.t the robot
#
#  id_to_state_map:
#            mapping that given the id of the measurement, returns its position in the mu vector
#  state_to_id_map:
#            mapping that given the index of mu vector, returns the id of the respective landmark
#
#outputs:
#  [mu, sigma]: the updated mean and covariance
#  [id_to_state_map, state_to_id_map]: the updated mapping vector between landmark position in mu vector and its id

function [mu, sigma, id_to_state_map, state_to_id_map] = correction(mu, sigma, observations, id_to_state_map, state_to_id_map)

  #determine how many landmarks we have seen in this step
  num_landmarks_measured = length(observations.observation);

  #dimension of the state (robot pose + landmark positions)
  state_dim = size(mu,1);

  #dimension of the current map (how many landmarks we have already seen)
  num_landmarks = (state_dim-3)/2;

  #if I have seen no landmarks, i do nothing
  if (num_landmarks_measured == 0)
    return;
  endif

  mu_t     = mu(1:2,:); # translational part of the robot pose
  mu_theta = mu(3); # rotation of the robot

  # re precompute some quantities that come in handy later on
  c   = cos(mu_theta);
  s   = sin(mu_theta);
  R   = [c -s; s c];  #rotation matrix
  Rt  = [c,s;-s c];    #transposed rotation matrix
  Rtp = [-s,c;-c,-s]; #derivative of transposed rotation matrix

  %# for below computation, we need to count how many observations
  %# of old landmarks we have
  num_old_landmarks_measured = 0;

  # Here two cases arise, the current landmark has been already seen, i.e. REOBSERVED landmark,
  # or the current landmark is completely new, i.e. NEW landmark.
  #
  # for simplicity we can divide the problem: first analyze only the reobserved landmark
  # and work with them as in a localization procedure (of course, with full Jacobian now).
  # With this reobserved landmark we compute the correction/update of mean and covariance.
  # Then, after the correction is done, we can simply add the new landmark expanding the
  # mean and the covariance.
  #
  #
  # First of all we are interested in REOBSERVED landmark
  for i=1:num_landmarks_measured

    #retrieve info about the observed landmark
    measurement = observations.observation(i);

    #fetch the position in the state vector corresponding to the actual measurement
    state_pos_of_landmark = id_to_state_map(measurement.id);

    #compute the index (vector coordinate) in the state vector corresponding to the pose of the landmark;
    landmark_state_vector_index = 4+2*(state_pos_of_landmark-1);

    #IF current landmark is a REOBSERVED LANDMARK
    if(state_pos_of_landmark != -1)

      #increment the counter of observations originating
      #from already known landmarks
      num_old_landmarks_measured++;

      # where we see the landmark
      z_t(end+1,:) = measurement.x_pose;
      z_t(end+1,:) = measurement.y_pose;

      #fetch the position of the landmark in the state (its x and y coordinate)
      landmark_mu = mu(landmark_state_vector_index:landmark_state_vector_index+1,:);

      #compute the measurement
      measure_prediction = Rt*(landmark_mu-mu_t);

      #add prediction to prediction vector
      h_t(end+1,:) = measure_prediction(1);
      h_t(end+1,:) = measure_prediction(2);

      #compute jacobian w.r.t robot
      C_m = zeros(2, state_dim);
      C_m(1:2, 1:2) = -Rt;
      C_m(1:2, 3) = Rtp*(landmark_mu-mu_t);

      #jacobian w.r.t landmark
      C_m(:,landmark_state_vector_index:landmark_state_vector_index+1)= Rt;

      #accumulate measurement
      C_t(end+1,:) = C_m(1,:);
      C_t(end+1,:) = C_m(2,:);
    endif
  endfor

  #if I have seen again at least one landmark
  #I need to update, otherwise I jump to the new landmark case
  if ((num_old_landmarks_measured > 0))

    #observation noise
    sigma_w = 0.01; #assume constant noise
    %sigma_z = sigma_w*eye(2*num_old_landmarks_measured); #TODO: set the correct dimension of the matrix
    sigma_z = sigma_w*eye(2*num_old_landmarks_measured);

    #Kalman gain
    K = sigma*C_t'*(sigma_z+C_t*sigma*C_t');

    #state correction
    innovation = (z_t - h_t);
    correction = K*innovation;
    mu         = mu + correction;

    #covariance correction
    sigma = (eye(state_dim)-K*C_t)*sigma
  endif
end
