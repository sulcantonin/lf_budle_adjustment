function [ P3D ] = lf_get_3DP_from_raysnd( rays,d,f )
%LF_GET_3DP_FROM_RAYSNM Summary of this function goes here
%   Detailed explanation goes here
P3D=cell(size(rays));
for i=1:length(rays)
    P3D{i}=zeros(3,1);
    for j=1:size(rays{i},1)
        P3D{i} = P3D{i}+lf_unproject_absolute( rays{i}(j,:), d(i), f );
    end
    P3D{i}=P3D{i}/size(rays{i},1);
end
end

