
%% Code to debug the internal workings of the sfactor function

wrange = [-6*10^-6 6*10^-6]
s = 'f'
c = 3*10^8; c0 = c;
eps0 = 8.85*10^-12;
mu0 = 4*pi*10^-7;
Nw = 100
Nw_pml = 0;
L0 = 10^-6; wvlen = 3*L0;
omega = 2*pi*c0/(wvlen);
lnR = -16;

%% create w coordinates
wmin = wrange(1); wmax = wrange(2);
deltaw = abs(wmax-wmin)/Nw;
wcoord = linspace(wmin,wmax,Nw+1);
windex = 1:length(wcoord);

right_pml = Nw_pml; %isolate points
left_pml = length(wcoord) - Nw_pml;
%where pml will be placed
right_pml_coord = wcoord(right_pml+1);
left_pml_coord = wcoord(left_pml);
PML_thickness = abs(wrange - [right_pml_coord left_pml_coord]);


wn = [right_pml,left_pml] %coordinates where the pml starts
eta0 = sqrt(mu0/eps0);  % vacuum impedance

m = 4;  % degree of polynomial grading

lnR = -16;  % R: target reflection coefficient for normal incidence
x_pml = Nw_pml/Nw;
y_mpl = Nw;
s_max = -(m+1)*lnR/(2*eta0)/PML_thickness(1)

%% check s and create pml coordinates accordingly
if s == 'b'
    wcoord = wcoord(1:end);
elseif s == 'f'
    wcoord = (wcoord(1:end-1) + wcoord(2:end))/2; %in the forward case, we are evaluating the pml
    %impedance at the midpoint
end

%% Finally create the s-factor array
sfactor_array = ones(Nw,1);
for i = 1:wn(1)-1
    sfactor_array(i) = 1-1i*s_max*((wcoord(wn(1))-wcoord(i))/PML_thickness(1))^m/(omega*eps0);
end
for i = 0:length(wcoord)-wn(2)
   index = length(wcoord)-i
   sfactor_array(index) = 1 - 1i*s_max*((wcoord(index)-wcoord(wn(2)))/PML_thickness(1))^m/(omega*eps0);
end
figure;
plot(abs(sfactor_array))
hold on;
test = sfactor_array;

s = 'b'
%% homework answer version
eta0 = sqrt(mu0/eps0);  % vacuum impedance
m = 4;  % degree of polynomial grading
lnR = (-16);  % R: target reflection coefficient for normal incidence

w_array = linspace(wrange(1), wrange(2), Nw+1);

loc_pml = [w_array(1 + Nw_pml), w_array(end - Nw_pml)];
d_pml = abs(wrange - loc_pml);

sigma_max = -(m+1)*lnR/(2*eta0) ./ d_pml;

if s == 'b'
	ws = w_array(1:end-1);
else  % s == 'f'
	assert(s == 'f');
	ws = (w_array(1:end-1) + w_array(2:end)) / 2;
end

ind_pml = {ws < loc_pml(1), ws > loc_pml(2)};  % {} means cell array

sfactor_array = ones(1, Nw);
%for n = 1, we deal with the left side,
%for n = 2, we deal with the right side pml, the center is just ones...
for n = 1:2
	sfactor = @(l) 1 - 1i * sigma_max(n)/(omega*eps0) * (l/d_pml(n)).^m;
	sfactor_array(ind_pml{n}) = sfactor(abs(loc_pml(n) - ws(ind_pml{n})));
end

figure(2)
plot(abs(sfactor_array))


