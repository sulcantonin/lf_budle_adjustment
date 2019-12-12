%% Solves for R by minimizing (15) with given t.
%% Uses SVD to get a solution with norm one, which is projected onto the
%% space of rotation matrices.
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [res Rest ] = lf_estimate_R_from_forward_warp_and_t_2( matches, R_current, t, f1, f2 )

    npoints = length( matches.rays1 );
    assert( length( matches.rays2 ) == npoints );
    [Ri ti] = invert_Rt( R_current, t );
    
    % count total number of projections
    total_projections_1 = 0;
    total_projections_2 = 0;
    for i=1:npoints
        total_projections_1 = total_projections_1 + size( matches.rays1{i}, 1 );
        total_projections_2 = total_projections_2 + size( matches.rays2{i}, 1 );
        assert( size( matches.rays1{i}, 2 ) == 4 );
        assert( size( matches.rays2{i}, 2 ) == 4 );
    end
        
    % least squares system for first two rows of R from X,Y
    % correspondences (which are independent of f)
    Ar = zeros( 2 * npoints * total_projections_1, 9 );
    n = 1;
    for i=1:npoints

        % constraints LF 1 to 2
        M2 = matches.M2s{i};
        P2 = lf_ray_projection_matrix( f2 );
        Z = M2*P2;
        A = Z(:,1:3);
        B = Z(:,4:6);

        nprojections1 = size( matches.rays1{i}, 1 );
        for j=1:nprojections1
            [q m] = lf_pluecker_ray_from_absolute_uvst( matches.rays1{i} (j,:), f1 );

            % equation 1
            Ar( n,: ) = equation_coefficients_atMb( B(1,:)*crossmat(t), q );
            Ar( n,: ) = Ar( n,: ) + equation_coefficients_atMb( A(1,:), q ) + equation_coefficients_atMb( B(1,:), m );
            n = n+1;
            
            % equation 2
            Ar( n,: ) = equation_coefficients_atMb( B(2,:)*crossmat(t), q );
            Ar( n,: ) = Ar( n,: ) + equation_coefficients_atMb( A(2,:), q ) + equation_coefficients_atMb( B(2,:), m );
            n = n+1;            
        end
        
        % use current R to set up second set of equations
        % constraints LF 2 to 1
        M1 = matches.M1s{i};
        P1 = lf_ray_projection_matrix( f1 );
        Z = M1*P1;
        A = Z(:,1:3);
        B = Z(:,4:6);

        nprojections2 = size( matches.rays2{i}, 1 );
        for j=1:nprojections2
            [q m] = lf_pluecker_ray_from_absolute_uvst( matches.rays2{i} (j,:), f2 );

            % equation 1
            Ar( n,: ) = equation_coefficients_atMtb( B(1,:)*crossmat(ti), q );
            Ar( n,: ) = Ar( n,: ) + equation_coefficients_atMtb( A(1,:), q ) + equation_coefficients_atMtb( B(1,:), m );
            n = n+1;
            
            % equation 2
            Ar( n,: ) = equation_coefficients_atMtb( B(2,:)*crossmat(ti), q );
            Ar( n,: ) = Ar( n,: ) + equation_coefficients_atMtb( A(2,:), q ) + equation_coefficients_atMtb( B(2,:), m );
            n = n+1;            
        end
        
    end
    N = n-1;
    
    % Solve for R using SVD
    [~,~,V] = svd( Ar , 'econ');
    Rsolve = V(:,9);
    res = norm(Ar*Rsolve);
    Rest = reshape( Rsolve, 3,3 )';
    Rest = project_to_rotation( Rest );

end
