%% Transforms a 3D point into a 2D manifold M (checkout "On Linear Structure from Motion for Light Field Cameras")
% Input:       X : A 3D points to be transformed (they need to be in
%                       their corresponding camera coordinate frame 
%                       (X \in R^{3x n x m})
%              f : f for each camera (f \in R^{m})
% Output       M : A 2D manifold eq. (10) for each 3D point in each camera
%                       (M \in R^{2 x 5 x n x m})
%
% (c) O. Johannsen, A. Sulc, B. Goldluecke, University of Konstanz
% License: Creative Commons BY-SA 4.0,
% Please cite our paper if you use our code in your work:
% O. Johannsen, A. Sulc, B. Goldluecke: On Linear Structure from Motion for Light Field Cameras

function M = lf_M_from_X(X,f)
    assert(length(f) == size(X,3))
    ncam = length(f);
    npoints = size(X,2);

    M = zeros(2,5,npoints,ncam);
    M(1,1,:,:) = 1;
    M(2,2,:,:) = 1;
    for j=1:ncam
        x = X(1,:,j);
        y = X(2,:,j);
        z = X(3,:,j);

        M(1,3,:,j) = f(j)./z-1;
        M(2,4,:,j) = f(j)./z-1;
        
        M(1,5,:,j) =  -(f(j).*x)./z;
        M(2,5,:,j) =  -(f(j).*y)./z;
        for i=1:size(M,3)
            M(:,:,i,j) = M(:,:,i,j) ./ norm(M(:,:,i,j));
        end
    end
end

