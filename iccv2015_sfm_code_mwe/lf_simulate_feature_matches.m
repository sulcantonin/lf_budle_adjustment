%% Return a list of simulated feature matches for a single light field
%% List is of the format required by lf_sfm() and related methods.
%% Also stores the ground truth data for later error estimation.
%
% Params: pre-initialized light fields,
%         number of 3D points, number of rays per light field,
%         Gaussian additive error std.deviations (s,t) plane and (u,v) plane
%         percentage of outliers (should be zero if no RANSAC implemented)
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [matches] = lf_simulate_feature_matches( LF1, LF2, N, K, sigma_st, sigma_uv, outlier_percentage )

    % correspondence data
    matches.npoints = N;
    matches.disparity1 = zeros( N,1 );
    matches.disparity2 = zeros( N,1 );
    matches.disparity1_gt = zeros( N,1 );
    matches.disparity2_gt = zeros( N,1 );
    matches.Xs = zeros( N,3 );
    [SC TC] = meshgrid( 1:LF1.S, 1:LF2.T );
    SC = SC(:); TC = TC(:);
    
    for i=1:N

        st = randperm( LF1.S * LF1.T, K );
        S = ( SC(st) - 1 ) / (LF1.S - 1.0) * 2.0 - 1.0;
        T = ( TC(st) - 1 ) / (LF1.T - 1.0) * 2.0 - 1.0;
        
        % 1. generate random uvd coordinates in light field center view
        rays1 = zeros( K,4 );
        rays1( 1,3 ) = S(1);
        rays1( 1,4 ) = T(1);
        rays1( 1,1 ) = (2.0 * rand() - 1.0);
        rays1( 1,2 ) = (2.0 * rand() - 1.0);        
        matches.disparity1_gt( i ) = 0.95 * rand() + 0.05;
        matches.Xs(i,:) = lf_unproject_relative( rays1( 1,: ), matches.disparity1_gt( i ), LF1.focal_length );

        % 2. generate K-1 correspondences in S,T range
        for k=1:K
            rays1( k,3 ) = S(k);
            rays1( k,4 ) = T(k);
            if rand() < outlier_percentage
                rays1( k,1 ) = rand()*LF1.X;
                rays1( k,2 ) = rand()*LF1.Y;
            else
                [ rays1( k,1 ) rays1( k,2 ) ] = ...
                    lf_project( matches.Xs(i,:)', rays1( k,3 ), rays1( k,4 ), ...
                                LF1.R, LF1.t, LF1.focal_length );            
            end
        end
        
        % store in ray arrays
        matches.rays1{i} = rays1;
    end

    
    % 3. transform UVST coordinates over to second light field
    for i=1:N
        Xt = LF2.R * matches.Xs(i,:)' + LF2.t;
        matches.disparity2_gt( i ) = LF2.focal_length / Xt(3);

        st = randperm( LF1.S * LF1.T, K );
        S = ( SC(st) - 1 ) / (LF1.S - 1.0) * 2.0 - 1.0;
        T = ( TC(st) - 1 ) / (LF1.T - 1.0) * 2.0 - 1.0;

        rays2 = zeros( K,4 );
        for k=1:K
            rays2( k,3 ) = S(k);
            rays2( k,4 ) = T(k);
            if rand() < outlier_percentage
                rays2( k,1 ) = rand()*LF2.X;
                rays2( k,2 ) = rand()*LF2.Y;
            else
                [ rays2( k,1 ) rays2( k,2 ) ] = ...
                    lf_project( matches.Xs(i,:)', rays2( k,3 ), rays2( k,4 ), LF2.R, LF2.t, LF2.focal_length );            
            end
        end
        
        matches.rays2{i} = rays2;
%        fprintf( '   forward warp ray %g %g %g %g\n', rays2( 1,: ) );
    end


    % 4. apply some Gaussian noise, plus a certain random chance that it is a
    % complete outlier
    sigma_X = sigma_uv;
    sigma_Y = sigma_uv;
    for i=1:N
        matches.rays1{i} = matches.rays1{i} + sigma_X * randn(K,4);
        matches.rays2{i} = matches.rays2{i} + sigma_Y * randn(K,4);
    end

    % estimate projection manifolds
    matches = lf_estimate_manifolds_for_matches( matches );
    
end
