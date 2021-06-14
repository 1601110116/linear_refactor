function sphi (nx, nz)
% This function solves the Poisson equation to get phi.
%  the method is explained in DocPoisson.pdf
global w phi lap calc outside calc_col enable_parallel

if enable_parallel
	phi_par = phi(2:end-1, 2:end-1, :);
	calc_col_par = calc_col;
	lap_par = lap;
	parfor iz = 2:nz+1
		tmp_bigcol = reshape(w(2:nx+1, 2:nx+1, iz), nx*nx, 1);
		tmp_smallcol = tmp_bigcol(calc_col_par);
		tmp_smallcol = lap_par \ tmp_smallcol;
		tmp_bigcol(calc_col_par) = tmp_smallcol;
		phi_par(:, :, iz) = reshape(tmp_bigcol, nx, nx);
	end
	phi(2:end-1, 2:end-1, :) = phi_par;
else
    for iz = 2:nz+1
	    tmp_bigcol = reshape(w(2:end-1, 2:end-1, iz), nx*nx, 1);
        tmp_smallcol = tmp_bigcol(calc_col);
        tmp_smallcol = lap \ tmp_smallcol;
        tmp_bigcol(calc_col) = tmp_smallcol;
		phi(2:end-1, 2:end-1, iz) = reshape(tmp_bigcol, nx, nx);
	end
end

phi(outside) = 0;

phi = zbcs(phi);
phi = lowpass_n(phi);
