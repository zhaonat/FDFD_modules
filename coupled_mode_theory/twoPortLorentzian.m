function [reconstructed_fit, parameters] = twoPortLorentzian(Ref, omega_scan, RT, x0)

    % input omega_scan should be in SI units rad/s
    % parameters are ordered as (gamma_tau, omega0, gamma_rad)
    % performs a coupled mode theory fit 
    % the background for the twoPort Lorentzian should be very small.
    
    if(nargin < 4)
        % initial condition
        x0 = rand(3,1);
    end
    
    if(RT == 'R')
        omega_scaled = omega_scan*1e-15; 
        %lorentzian = @(w,w0,t) (1/pi)*(t/2)/((w-w0)^2+t^2/4);

        cost_function = @(x) sum(abs(Ref - (x(3))^2./((omega_scaled-x(2)).^2+(x(1)+x(3))^2)));
        options = optimset('MaxFunEvals', 40000, 'MaxIter', 100000);
        parameters = fminsearch(cost_function,x0, options);
        reconstructed_fit = ( parameters(3))^2./((omega_scaled-parameters(2)).^2 + (parameters(1)+parameters(3))^2);
    else
        disp('here')
       omega_scaled = omega_scan*1e-15; 
        %lorentzian = @(w,w0,t) (1/pi)*(t/2)/((w-w0)^2+t^2/4);

        cost_function = @(x) sum(abs(Ref - (omega_scaled-x(2)).^2./((omega_scaled-x(2)).^2+(x(1)+x(3))^2)));
        options = optimset('MaxFunEvals', 40000, 'MaxIter', 1000000);
        parameters = fminsearch(cost_function,x0, options);
        reconstructed_fit = ( parameters(3))^2./((omega_scaled-parameters(2)).^2 + (parameters(1)+parameters(3))^2); 
    end
end