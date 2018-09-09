clear all; close all;

is_new = true;  % new simulation; make this false to load calculated results

if is_new
	%% Set up the common parameters.
	L0 = 1e-9;  % length unit: nm
	xrange = [-1000 1000];  % x boundaries in L0
	yrange = [-1000 1000];  % y boundaries in L0
	Npml = [10 10];  % [Nx_pml, Ny_pml]

	wvlen = 1000;
	eps_square = -100-10i;  % square permittivity (some metal)
	a = 600;  % side of square

	Ns = 100:200:900;  % list of numbers of cells to test
	n = length(Ns);
	sigmas = NaN(1,n);

	for i = 1:n
		disp(['Solving for Nx = Ny = ', num2str(Ns(i)), '.']);
		N = [Ns(i) Ns(i)];  % [Nx Ny]

		%% Expand the simulation domain to include PML.
		[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);

		%% Set up the permittivity.
		within_square = @(x, y) abs(x) < a/2 & abs(y) < a/2;  % function handle that returns true for (x,y) in square
		eps_with_sq = ones(N);
		eps_with_sq = assign_val(eps_with_sq, xrange, yrange, within_square, eps_square);

		eps_no_sq = ones(N);

		%% Assign the magnetic current source for a plane wave propagating in the y-direction.
		Mz = zeros(N);
		Mz(:, Npml(2) + 10) = 1;

		%% Solve the TE mode equation with with and without the square.
		[Hz0, Ex0, Ey0] = solveTE(L0, wvlen, xrange, yrange, eps_no_sq, Mz, Npml);
		[Hz, Ex, Ey] = solveTE(L0, wvlen, xrange, yrange, eps_with_sq, Mz, Npml);

		%% Calculate the intensity of the incident wave.
		[Sx0, Sy0] = poynting(Hz0, Ex0, Ey0);  % Poynting vector of incident wave
		I_inc = Sy0(ceil(N(1)/2), ceil(N(2)/2));  % Sx0 == 0

		%% Calculate the absorption cross section.
		[Sx, Sy] = poynting(Hz, Ex, Ey);  % Poynting vector of total wave
		offset = 20;
		ind_x = (Npml(1)+offset):(N(1)-Npml(1)-offset);
		ind_y = (Npml(2)+offset):(N(2)-Npml(2)-offset);
		P_abs = -flux_rect(Sx, Sy, ind_x, ind_y, dL);  % scattered power flux

		sigmas(i) = P_abs / I_inc;  % absorption cross section
	end
	save('conv_test', 'Ns', 'sigmas');
else  % not new simulation
	load('conv_test');
end

%% Plot the absorption cross sections.
figure;
plot(Ns, sigmas, 'o-');
title('absorption cross section');
xlabel('Nx');
ylabel('\sigma_{abs}');

%% Plote the percentage errors with respect to Nx = Ny = 900.
figure;
plot(Ns, (sigmas - sigmas(end))/sigmas(end) * 100, 'o-');
title('pencent error in absorption cross section');
xlabel('Nx');
ylabel('% error');

