%% Function creates a 3D mesh of points for a given depth map
% Input: 	D 	: disparity map (D \in R^{h x w})
%		H 	: intrinsic calibration matrix (H \in R^{5 x 5})
%		f 	: focal length 
% 		center 	: center view coordinates (in pixel coordiantes) (center \in R^{2})
%		(offset): offset on each side of disparity map
%		(step)  : rastering of disparity map 
% Output 	X,Y,Z 	: 3D points for the given D
%
% (c) O. Johannsen, A. Sulc, B. Goldluecke, University of Konstanz
% License: Creative Commons BY-SA 4.0,
% Please cite our paper if you use our code in your work:
% O. Johannsen, A. Sulc, B. Goldluecke: On Linear Structure from Motion for Light Field Cameras
function [ X,Y,Z] = lf_unproject_dmap( D,H,f,center,varargin)
% parsing input parameters
p = inputParser;
addRequired (p,'D',@isnumeric);
addRequired(p,'H',@isnumeric);
addRequired(p,'f',@isnumeric);
addRequired(p,'center',@isnumeric);
addParameter(p,'offset',1,@isnumeric);
addParameter(p,'step',1,@isnumeric);

parse(p,D,H,f,center,varargin{:});
% inlier threshold for SFM RANSAC
D = p.Results.D;
H = p.Results.H;
f = p.Results.f;
center = p.Results.center;
offset = p.Results.offset;
step = p.Results.step;

%%
Hsub=H([1,3],:);
[X,Y]=ndgrid(offset:step:size(D,1)-offset+1,offset:step:size(D,2)-offset+1);
sizeofX=size(X);

m8=ones(5,length(Y(:)));
m8(1,:)=center(2);
m8(2,:)=center(1);
m8(3,:)=Y(:);
m8(4,:)=X(:);

m9=ones(5,length(Y(:)));
m9(1,:)=center(2)+1;
m9(2,:)=center(1);
m9(3,:)=Y(:)-D(sub2ind(size(D),X(:),Y(:)));
m9(4,:)=X(:);

lf=Hsub*(m9-m8);

m8=ones(5,length(Y(:)));
m8(1,:)=center(2);
m8(2,:)=center(1);
m8(3,:)=Y(:);
m8(4,:)=X(:);

ray=H*m8;

Xt = zeros(3,length(X(:)));
Xt(3,:) = f ./ (1-lf(2,:)./lf(1,:));
Xt(1,:) = (ray(3,:) - ray(1,:)) .* Xt(3,:) / f + ray(1,:);
Xt(2,:) = (ray(4,:) - ray(2,:)) .* Xt(3,:) / f + ray(2,:);

X=Xt(1,:);
Y=Xt(2,:);
Z=Xt(3,:);

X = reshape(X,sizeofX);
Y = reshape(Y,sizeofX);
Z = reshape(Z,sizeofX);
end
