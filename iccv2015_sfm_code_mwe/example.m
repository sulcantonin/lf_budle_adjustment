%% Demo code for our ICCV 2015 paper
%% "On Linear Structure-from-Motion for Light Field Cameras"
%% Ole Johannsen, Antonin Sulc, Bastian Goldluecke

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% to keep it simple, contains a minimal working example to demonstrate
% our method, but not the implementations of competing methods to
% generate the full table of comparisons.



%% Step 1. Setup light field parameters

% first light field setup
clear LF1;
LF1.focal_length = 0.8;
LF1.R = eye(3);
LF1.t = [0;0;0];
[LF1.Ri LF1.ti] = invert_Rt( LF1.R, LF1.t );
LF1.X = 512;
LF1.Y = 512;
LF1.S = 9;
LF1.T = 9;

% second light field setup
clear LF2;
% fixed rotation angle and translation (for testing)
a = 0.2*pi;
LF2.R = random_rotation( 30.0 );
LF2.t = random_translation( 10.0 );
LF2.focal_length = 1.0;
[LF2.Ri LF2.ti] = invert_Rt( LF2.R, LF2.t );
LF2.X = 512;
LF2.Y = 512;
LF2.S = 9;
LF2.T = 9;



%% Step 2. Generate set of 3D points with several ray correspondences

% number of 3D points to generate
npoints = 20;

% number of projections of each 3D point in each light field
projections_per_point = 10;

% noise levels
% subaperture view coordinates
sigma_st = 0.00;
% image coordinates
sigma_uv = 0.02;
% to deal with outliers, embedding into RANSAC is required
% (not implemented in this minimal code example)
outlier_percentage = 0.0;


% directions and moments, first light field
matches = ...
    lf_simulate_feature_matches( LF1, LF2, ...
                                 npoints, projections_per_point, ...
                                 sigma_st, sigma_uv, outlier_percentage );
fprintf( 'Match generation complete.\n' );




%% Step 3. Run proposed algorithm for light field pose estimation

% number of non-linear refinement iterations
% optional, can be zero, but often a small number here helps
refinement_iterations = 0;

% run main function to compute rotation R and translation t
[R t] = lf_sfm( matches, LF1.focal_length, LF2.focal_length, refinement_iterations );

% compute and display the estimation errors
errors = Rt_errors( R,t, LF2.R, LF2.t );
fprintf( 'SfM estimate complete.\n' );
fprintf( '  rotation error    angular  : %g\n', errors(1) );
fprintf( '                    absolute : %g\n', errors(2) );
fprintf( '                    relative : %g\n', errors(3) );
fprintf( '  translation error angular  : %g\n', errors(4) );
fprintf( '                    absolute : %g\n', errors(5) );
fprintf( '                    relative : %g\n', errors(6) );

