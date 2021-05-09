function f_filtered = filter_z(f, n)
% This function filters the last dimension of a 3D field using fft. 
%
% Input:
%  f: the field to filter. Its last dimension is filtered, usually z.
%     Note: f is assumed to include boundaries, and the parallel
%           boundary condition is periodic
%  n: the harmonic to keep. n=1 means the wavelength is the 
%     parallel size of the device.

f_n = fft(f(:, :, 2:end-1));
res_n = 
