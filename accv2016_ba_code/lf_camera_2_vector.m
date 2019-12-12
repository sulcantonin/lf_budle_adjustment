function [vec] = lf_camera_2_vector(R,t)

quat = rotm2quat(R);
vec = [quat(:); t(:)];

end
