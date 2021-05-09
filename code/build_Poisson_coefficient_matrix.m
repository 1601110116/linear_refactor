function build_Poisson_coefficient_matrix(dx, nx, x_max)
% The electrostatic potential is acquired by solving a 2D Poisson equation.
%  please read DocPoisson.pdf in the 'Documentation' folder if you want 
%  to understand this code


global lap calc calc_col

tmp1 = 2*speye(nx) - diag(ones(1,nx-1), 1) - diag(ones(1, nx-1), -1);
tmp2 = -(dx^-2) * (kron(speye(nx), tmp1) + kron(tmp1, speye(nx)));
calc_col = reshape(calc(:,:,1), [], 1);
lap = tmp2(calc_col, calc_col);
