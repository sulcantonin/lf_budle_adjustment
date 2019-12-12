%% Solves for t from a list of matches given that R is already estimated
%% Generates set of equations from (13).
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [res, t] = lf_estimate_t_from_forward_warp( matches, R, f1, f2 )

    npoints = length( matches.rays1 );
    assert( length( matches.rays2 ) == npoints );

    % count total number of projections
    total_projections_1 = 0;
    total_projections_2 = 0;
    for i=1:npoints
        total_projections_1 = total_projections_1 + size( matches.rays1{i}, 1 );
        total_projections_2 = total_projections_2 + size( matches.rays2{i}, 1 );
        assert( size( matches.rays1{i}, 2 ) == 4 );
        assert( size( matches.rays2{i}, 2 ) == 4 );
    end
        
    
    % setup system of equations
    V = zeros( 2 * (total_projections_1 + total_projections_2), 3 );
    b = zeros( 2 * (total_projections_1 + total_projections_2), 1 );
    n = 1;
    for i=1:npoints

        %% First set of equations: Project rays from LF1 into LF2
        M2 = matches.M2s{i};
        P = lf_ray_projection_matrix( f2 );
        Z = M2*P;
        A = Z(:,1:3);
        B = Z(:,4:6);

        nprojections1 = size( matches.rays1{i}, 1 );
        for j=1:nprojections1
            [q, m] = lf_pluecker_ray_from_absolute_uvst( matches.rays1{i} (j,:), f1 );
            b( n:n+1, 1 ) = A*R*q + B*R*m;
            V( n:n+1, : ) = B*crossmat( R*q );
            n = n + 2;
        end
            
        %% Second set of equations: Project rays from LF2 into LF1
        M1 = matches.M1s{i};
        P = lf_ray_projection_matrix( f1 );
        Z = M1*P;
        A = Z(:,1:3);
        B = Z(:,4:6);

        nprojections2 = size( matches.rays2{i}, 1 );
        for j=1:nprojections2

            % LF 2 to 1
            [q m] = lf_pluecker_ray_from_absolute_uvst( matches.rays2{i} (j,:), f2 );
            b( n:n+1, 1 ) = A*R'*q + B*R'*m;
            V( n:n+1, : ) = - B*crossmat( R'*q )*R';
            n = n + 2;
        end


    end
    
    % Solve via least squares and compute residual
    t = V \ b;
    res = norm( V*t - b )^2;

end
