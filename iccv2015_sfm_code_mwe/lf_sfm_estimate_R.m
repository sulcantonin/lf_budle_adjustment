%% Solves for R from a list of matches using the SVD method for equation (14)
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [res Rest] = lf_sfm_estimate_R( matches, f1, f2 )

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
        
    % least squares system for first two rows of R from ray to manifold
    % correspondences
    Ae = zeros( 2 * (total_projections_1 + total_projections_2), 18 );
    Ar = zeros( 2 * (total_projections_1 + total_projections_2), 9 );
    P1 = lf_ray_projection_matrix( f1 );
    P2 = lf_ray_projection_matrix( f2 );

    n = 1;
    for i=1:npoints

        %% First set of equations: Project rays from LF1 into LF2
        M2 = matches.M2s{i};
        Z = M2*P2;
        A = Z(:,1:3);
        B = Z(:,4:6);

        nprojections1 = size( matches.rays1{i}, 1 );
        for j=1:nprojections1

            % LF 1 to 2
            [q m] = lf_pluecker_ray_from_absolute_uvst( matches.rays1{i} (j,:), f1 );

            % equation 1
            Ae( n,1:9 ) = equation_coefficients_atMb( B(1,:), q );
            Ar( n,: ) = equation_coefficients_atMb( A(1,:), q ) + equation_coefficients_atMb( B(1,:), m );
            n = n+1;
            
            % equation 2
            Ae( n,1:9 ) = equation_coefficients_atMb( B(2,:), q );
            Ar( n,: ) = equation_coefficients_atMb( A(2,:), q ) + equation_coefficients_atMb( B(2,:), m );
            n = n+1;            
            
        end
        
        
        %% Second set of equations: Project rays from LF2 into LF1
        M1 = matches.M1s{i};
        Z = M1*P1;
        A = Z(:,1:3);
        B = Z(:,4:6);

        nprojections2 = size( matches.rays2{i}, 1 );
        for j=1:nprojections2

            % LF 2 back to LF 1
            [q m] = lf_pluecker_ray_from_absolute_uvst( matches.rays2{i} (j,:), f2 );

            % equation 1
            Ae( n,10:18 ) = equation_coefficients_atMtb( B(1,:), q );
            Ar( n,: ) = equation_coefficients_atMtb( A(1,:), q ) + equation_coefficients_atMtb( B(1,:), m );
            n = n+1;
            
            % equation 2
            Ae( n,10:18 ) = equation_coefficients_atMtb( B(2,:), q );
            Ar( n,: ) = equation_coefficients_atMtb( A(2,:), q ) + equation_coefficients_atMtb( B(2,:), m );
            n = n+1;            
            
        end
    end
    N = n-1;
    
    % ignore essential matrix, just solve for R
    A = (Ae * pinv( Ae ) - eye( N )) * Ar;    
    [~,~,V] = svd( A,'econ');
    Rsolve = V(:,9);
    res = norm(A*Rsolve);
    Rest = reshape( Rsolve, 3,3 )';
    Rest = project_to_rotation( Rest );

end
