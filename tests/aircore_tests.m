%effectively, what is different from a test and an example
% a test must test an ultra specific functionality

%% Test double class
% exp = 'double';
% act = ones;
% assert(isa(act,exp))


%% specify domain
%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
wvlen = 2.0;  % wavelength in L0
xrange = [-2 2];  % x boundaries in L0
yrange = [-2 2];  % y boundaries in L0
N = [400 400];  % [Nx Ny]
Npml = 0*[10 10];  % [Nx_pml Ny_pml] need to deal with the special case where both are 
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML


%% specify aircore
air_core = air_core_structure(xrange, yrange, N, Lpml);
y_core_center = 0;
air_core_thickness = 2;
wall_thickness = 0.5;
wall_materials =[6,12];


%% grating specifications
num_cells = 2;
lattice_constant = 1;
thickness = 0.3;
epsilon_array = [12,23];
fill_factor = 0.3;

y_center = -0.5;
wall_properties = {num_cells, lattice_constant, ...
                            thickness, epsilon_array, fill_factor};
wall_properties2 = {num_cells, lattice_constant, ...
                            0.2, epsilon_array, fill_factor};
%% build the gratings first, pass the properties to the air core
grating1 = periodic_grating(xrange, yrange, N, Lpml);
grating1.add_grating_array(num_cells, lattice_constant, ...
                thickness, epsilon_array, fill_factor, y_center)
            
grating1.add_grating_array(num_cells, lattice_constant, ...
                thickness, epsilon_array, fill_factor, 0.5)
%% test the multi cor3
y_centers = [-0.5, -0.5+0.7+thickness,-0.5+2*0.7+2*thickness];

air_core. build_multi_core(y_centers, {wall_properties, wall_properties2,wall_properties})

figure();
visreal(air_core.epsilon, xrange, yrange);

%% draw outlines based on fields in the class...
for i = 1:length(air_core.grating_properties)
    structure = air_core.grating_properties{i};
    line(xrange, [structure.y_pos, structure.y_pos])
    line(xrange, [structure.wall_coords(1), structure.wall_coords(1)]);
    line(xrange, [structure.wall_coords(2), structure.wall_coords(2)]);
end


