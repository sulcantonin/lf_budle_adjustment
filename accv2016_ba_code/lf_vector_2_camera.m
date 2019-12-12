function [R,t] = lf_vector_2_camera(vec)

if size(vec,2) == 7
    vec = vec';
end

assert(size(vec,1)==7);
assert(size(vec,2)==1);

R = quat2rotm(vec(1:4)');
t = vec(5:7);
end
