function diagnose (idiagnose, nt_per_diagnose)
% This function draws the contours fo the simulated quantities at run time
global den vi w jz ve Te phi vEx vEy visual dt xX zX denref Tref cs0 rhos0 t0 ...
	data_path diagnose_path
persistent ydiag zdiag

ydiag = length(xX) / 2;
zdiag = length(zX) / 2;
fig_position = [50, 50, 1500, 1000];
x_tick = -10: 5: 10;
font_size = 8;


datafile = fullfile(data_path, ['dat', sprintf('%4.4d', idiagnose)]);
save(datafile, 'den', 'vi', 'w', 'Te', 'phi');


if visual
	close;
	tX = idiagnose * dt * nt_per_diagnose * t0;
	figure('name', pwd, 'Visible', 'off');
	set(gcf, 'Position', fig_position);

	subplot(4,4,1);  pcolor(xX, xX, denref*den(:, :, zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$n/cm^-3$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);

	subplot(4,4,2);  pcolor(xX, xX, Tref*Te(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$T_{e}/eV$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);

	subplot(4,4,3);  pcolor(xX, xX, cs0*vi(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$v_{\parallel i}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);

	subplot(4,4,4);  pcolor(xX, xX, Tref/rhos0^2*w(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$w/\left(V*cm^{-2}\right)$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);

	subplot(4,4,5);  pcolor(xX, xX, cs0*jz(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$j_{\parallel}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);

	subplot(4,4,6);  pcolor(xX, xX, cs0*ve(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$v_{\parallel e}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);

	subplot(4,4,7);  pcolor(xX, xX, Tref*phi(:,:,zdiag));
	colormap jet;  colorbar;  shading interp;
	title('$$\phi/V$$', 'interpreter', 'latex');
	xlabel('y/cm');  ylabel('x/cm');
	set(gca, 'XTick', x_tick);
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);


	subplot(4,4,9);  pcolor(zX, xX, squeeze(denref*den(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$n/cm^-3$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);

	subplot(4,4,10);  pcolor(zX, xX, squeeze(Tref*Te(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$T_{e}/eV$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);

	subplot(4,4,11);  pcolor(zX, xX, cs0*squeeze(vi(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$v_{\parallel i}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);

	subplot(4,4,12);  pcolor(zX, xX, Tref/rhos0^2*squeeze(w(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$w/\left(V*cm^{-2}\right)$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);

	subplot(4,4,13);  pcolor(zX, xX, cs0*squeeze(jz(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$j_{\parallel}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);

	subplot(4,4,14);  pcolor(zX, xX, cs0*squeeze(ve(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$v_{\parallel e}/\left(cm/s\right)$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);

	subplot(4,4,15);  pcolor(zX, xX, Tref*squeeze(phi(:,ydiag,:)));
	colormap jet;  colorbar;  shading interp;
	title('$$\phi/V$$', 'interpreter', 'latex');
	xlabel('z/cm');  ylabel('x/cm');
	set(gca, 'YTick', x_tick);
	set(gca, 'fontSize', font_size);

	subplot(4,4,16);  xlabel(['time = ', num2str(tX), ...
		sprintf(' s\nidiagnose = '), num2str(idiagnose)]);
	set(gca, 'fontSize', font_size);
    
	fig_file = fullfile(diagnose_path, ['contours', sprintf('%4.4d', idiagnose), '.png']);
    print(gcf, '-dpng', fig_file);

end
