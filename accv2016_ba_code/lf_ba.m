function [X_n, pose_n,f_n,H_n] = lf_ba(X_0,pose_0,rays,f_0,varargin)
vec = @(i) i(:);

% parsing input parameters
argparser = inputParser;
addRequired (argparser,'X_0');
addRequired(argparser,'pose_0',@isnumeric);
addRequired(argparser,'rays',@iscell);
addRequired(argparser,'f_0',@isnumeric);
addParameter(argparser,'display','iter',@isstring);
addParameter(argparser,'niter',50,@isnumeric);

parse(argparser,X_0,pose_0,rays,f_0,varargin{:});

X_0 = argparser.Results.X_0;
pose_0 = argparser.Results.pose_0;
rays = argparser.Results.rays;
f_0 = argparser.Results.f_0;
display = argparser.Results.display;
niter = argparser.Results.niter;

assert(size(pose_0,1) == 3);
assert(size(pose_0,2) == 4);

n_cam = size(pose_0,3);
n_point = size(X_0,2);

assert(n_cam == size(rays,2));
assert(n_point == size(rays,1));
assert(length(f_0) == 1)

%% theta_0
theta_0 = lf_cameras_points_2_vector(X_0,pose_0,f_0);

% intermediate function, which passes local parameters
func = @(theta) (lf_manifold_ba(theta,rays));


opts = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','display', ...
    display,'MaxFunEvals',inf,'MaxIter',niter,'TolFun',eps,'TolX',eps);
tic
fprintf('Running %i iterations of bundle adjustment\n',niter)
theta_n = lsqnonlin(func,theta_0,[],[],opts);
fprintf('Budnle adjustment done after %f seconds\n',toc)

%% extrinsics and intrinsics from BA
[X_n,pose_n,f_n] = lf_vector_2_cameras_points(theta_n,n_cam,n_point);

end