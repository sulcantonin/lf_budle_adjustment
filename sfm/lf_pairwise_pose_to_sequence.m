function [pose_sequence] = lf_pairwise_pose_to_sequence(pose_pairwise,H,f)
N = size(pose_pairwise,3);

pose_sequence = zeros(3,4,N);
pose_sequence(1:3,1:3,1) = eye(3);
pose_sequence(:  ,4  ,1) = [0;0;0];

for i=2:N+1
    R = pose_pairwise(1:3,1:3,i-1);
    t = pose_pairwise(:  ,4  ,i-1);
    Rt = [R, t; 0 0 0 1] * [pose_sequence(1:3,1:3,i-1),pose_sequence(:,4,i-1); 0 0 0 1];
    pose_sequence(1:3,1:3,i) = Rt(1:3,1:3);
    pose_sequence(:,4,i) = Rt(1:3,4);
end