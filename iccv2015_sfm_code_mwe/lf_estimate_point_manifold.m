%% Estimate the 4D affine space for a set of 4D rays according to equation (10).
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function M = lf_point_manifold( rays )

    K = size( rays,1 );
    assert( size( rays,2 ) == 4 );
    
    A = zeros( 2*K, 3 );
    b = zeros( 2*K, 1 );
    for k=1:K
        
        A( 2*k+0, 1 ) = rays( k,3 );
        A( 2*k+0, 2 ) = 1.0;
        b( 2*k+0 ) = -rays( k,1 );
        
        A( 2*k+1, 1 ) = rays( k,4 );
        A( 2*k+1, 3 ) = 1.0;
        b( 2*k+1 ) = -rays( k,2 );

    end
    
    s = A \ b;
    M = [ 1 0 s(1) 0 s(2); 0 1 0 s(1) s(3) ];
    M = M / norm(M);
    
end
