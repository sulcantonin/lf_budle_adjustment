%% Forward warp error. For each match and every ray in the match:
%% project ray into respective other LF, compute an error
%% measure by estimating squared norm of residual in equation (10).
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [err] = lf_forward_warp_error( matches, R, t, f1, f2 )

    npoints = length( matches.rays1 );
    assert( length( matches.rays2 ) == npoints );

    [Ri ti] = invert_Rt( R,t );
    
    err = 0;
    for i=1:npoints
        M1 = matches.M1s{i};
        M2 = matches.M2s{i};
        nprojections1 = size( matches.rays1{i}, 1 );
        nprojections2 = size( matches.rays2{i}, 1 );
        assert( size( matches.rays1{i}, 2 ) == 4 );
        assert( size( matches.rays2{i}, 2 ) == 4 );
        
        for j=1:nprojections1
            % light field 1 to light field 2
            ray1 = squeeze( matches.rays1{i} (j,:) );
            [q m] = lf_pluecker_ray_from_absolute_uvst( ray1, f1 );
            [ok ray2] = lf_pluecker_ray_to_absolute_uvst( q, m, R, t, f2 );
            if ok
                err = err + norm( M2 * [ray2;1] )^2;
            end
        end
            
        for j=1:nprojections2
            % light field 2 to light field 1
            ray2 = squeeze( matches.rays2{i}  ( j,: ) );
            [q m] = lf_pluecker_ray_from_absolute_uvst( ray2, f2 );
            [ok ray1] = lf_pluecker_ray_to_absolute_uvst( q, m, Ri, ti, f1 );
            if ok
                err = err + norm( M1 * [ray1;1] )^2;
            end
        end        
    end
    
end
