%% Main entry function for our method.
%% Computes relative rotation and translation of the second light
%% field compared to the first given a list of feature matches
%% of the form (1). Focal length must be known for both light fields.
%
% Let N be the number of matches of the form (1), then
%
% matches.rays1{1..N} contains a i x 4 matrix for each match
%                     with i uvst coordinates of rays in the first LF
% matches.rays2{1..N} contains a j x 4 matrix for each match
%                     with j uvst coordinates of rays in the second LF
%
% Focal lengths f1, f2 for both light fields must be given.
%
% The last parameter gives the number of non-linear refinement
% iterations, see equation (15), which is optional and can be zero.
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [R t] = lf_sfm( matches, f1, f2, ref_iter )

    %% Basic linear method, equation (14) and (13)
    [res R] = lf_sfm_estimate_R( matches, f1, f2 );
    [res t] = lf_sfm_estimate_t( matches, R, f1, f2 );
    %err = lf_forward_warp_error( matches, R, t, f1, f2 );
    %fprintf( 'V2: forward warp error before refining %g\n', err );

    %{
        % some example code to demonstrate how one could optimize for focal
        % length (not explained in the paper)
        % the forward warp error is minimized over a bunch of possible
        % focal length values.
        % example below optimizes only for f1 given f2=1.0, but can
        % be easily extended.
        for i=1:10

            estimator2 = @(F) lf_forward_warp_error( matches, R, t, F(1), 1.0 );
            [F opt_res] = fminsearch( estimator2, F(1), optimset('TolX',1e-8) )

            [res R E] = lf_estimate_R( matches, F(1), 1.0 );
            [res t] = lf_estimate_t( matches, R, F(1), 1.0 );
        end
    %}
    

    %% Refinement iterations (see text after equation (15))
    for i=1:ref_iter
        [res R] = lf_sfm_estimate_R_with_known_t( matches, R, t, f1, f2 );
        [res t] = lf_sfm_estimate_t( matches, R, f1, f2 );
        %err = lf_forward_warp_error( matches, R,t, f1, f2 );
        %fprintf( 'refining iteration %i forward warp error %g\n', i, err );
    end
    
end









