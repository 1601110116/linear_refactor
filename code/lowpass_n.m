function f = lowpass_n(f)

global n_up_to

persistent f_n
if n_up_to >= 0	
	f_n = fft(f(:, :, 2:end-1), [], 3);
	f_n(:, :, n_up_to+2: end-n_up_to) = 0;
	f(:, :, 2:end-1) = real(ifft(f_n, [], 3));
	f = zbcs(f);
end
