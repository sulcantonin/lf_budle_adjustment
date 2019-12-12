%% For a list of feature matches, estimate the 4D affine subspace for the
%% projected rays for both light fields according to equation (10).
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [matches] = lf_estimate_manifolds_for_matches( matches )

    N = matches.npoints;
    for i=1:N
        M1 = lf_estimate_point_manifold( matches.rays1{i} );
        M2 = lf_estimate_point_manifold( matches.rays2{i} );
        matches.M1s{i} = M1;
        matches.M2s{i} = M2;
        
        % extract estimated disparity values from point manifolds
        % just for reference, only used to estimate errors
        M1 = M1 / M1(1,1);
        M2 = M2 / M2(1,1);
        matches.disparity1(i) = 1.0 + M1(1,3);
        matches.disparity2(i) = 1.0 + M2(1,3);

    end
    
end
