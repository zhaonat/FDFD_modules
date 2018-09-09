clear all; close all; clc;

%% Set up the domain parameters.
L = [10 10 10];  % [Lx Ly Lz] in microns
N = [4 4 4];  % [Nx Ny Nz]
dL = L./N;  % [dx dy dz]

%% Create Dws for w = x,y,z and s = f,b, and then plot their sparsity patterns.
fig = figure(1);
set (fig, 'Units', 'normalized', 'Position', [0 0 1 1]);  % maximize figure window
ind_plot = 0;
for s = 'fb'
	for w = 'xyz'
		Dws = createDws_dense(w, s, dL, N);

		ind_plot = ind_plot + 1;
		subplot(2, 3, ind_plot); spy(Dws);
	end
end

Dxf = 0; Dyf = 0;
Dxb = 0; Dyb = 0;
%% Run basic tests on the created matrices.h
for w = 'xyz'
	Dwf = createDws_dense(w, 'f', dL, N);	
	
	% Check if the diagonal elements of Dxf are -1/dx.
	dw = dL('xyz' == w);
	tf = all(diag(Dwf) == -1/dw);
	disp(['- diag(D', w, 'f) == -1/d', w, '?  ', logical2str(tf), '.']);

	% Check if the sum of the elements in each row of Dxf is zero.
	tf = all(sum(Dwf, 2) == 0);
	disp(['- row-wise sum of D', w, 'f == 0?  ', logical2str(tf), '.']);

	% Check if Dxb == -Dxf^T.
	Dwb = createDws_dense(w, 'b', dL, N);
	tf = isequal(Dwb, -Dwf.');
	disp(['- D', w, 'b == -D', w, 'f^T?  ', logical2str(tf), '.']);
	
	fprintf('\n');
    if w=='x'
       Dxf = Dwf; Dxb = Dwb; 
    end
    if w == 'y' 
       Dyf = Dwf; Dyb = Dwb;
    end
    
end

