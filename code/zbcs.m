function a = zbcs (a)
% This function implements z-direction boundary conditions given by zbc_mode

global zbc_mode

if zbc_mode == 1  % free bc
	a(:, :, 1) = a(:, :, 2);
	a(:, :, end) = a(:, :, end-1);
elseif zbc_mode == 2  % periodic bc
	a(:, :, 1) = a(:, :, end-1);
	a(:, :, end) = a(:, :, 2);
elseif zbc_mode == 3  % sheath bc
	error('Sheath zbc is not ready yet');
else
	error('Invalid input vairable [zbc_mode]');
end

