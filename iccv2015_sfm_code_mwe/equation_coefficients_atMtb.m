%% Returns the nine coefficients for estimation of a matrix M
%% from a term of the form a^T M^T b
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.

function [A] = equation_coefficients_atMtb( a, b )

    A = zeros( 1,9 );
    assert( length(a) == 3 );
    assert( length(b) == 3 );

    A(1) = a(1) * b(1);
    A(2) = a(2) * b(1);
    A(3) = a(3) * b(1);

    A(4) = a(1) * b(2);
    A(5) = a(2) * b(2);
    A(6) = a(3) * b(2);

    A(7) = a(1) * b(3);
    A(8) = a(2) * b(3);
    A(9) = a(3) * b(3);

end
