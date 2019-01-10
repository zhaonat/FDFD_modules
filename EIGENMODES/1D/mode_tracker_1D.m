
%given a cell which contains MULTIPLE eigenmode simulations, we want to
% SORT ALL THE MODES OUT BY SHAPE

function [field_tracker] = mode_tracker_1D(eigen_mode_simulations, k_eigs_store, omega_scan)

    num_sims = length(eigen_mode_simulations);
    field_tracker = cell(1);
    
    %% extract out unique fields
    for i = 1:length(eigen_mode_simulations{1})
        field_tracker{i,1} = eigen_mode_simulations{1}{i};
        field_tracker{i,2} = [];
        field_tracker{i,3} = [];

    end
    
    for i = 1:num_sims
        eigenmodes = eigen_mode_simulations{i};
        
        for j = 1:length(eigenmodes)
           field = eigenmodes{j}; 
%            plot(real(field));
%            hold on;
%            drawnow();
           %compare this field to every field in field_tracker
           d = size(field_tracker);
           mindiff = 1e8; minindex = 1;
           for k = 1:d(1)
              stored_field = field_tracker{k,1};
              if(isempty(stored_field))
                  field_tracker{k,1} = field;
                  field_tracker{k,2} = [k_eigs_store(j,i)];
                  field_tracker{k,3} = [omega_scan(i)];

              elseif(norm(abs(field)/max(abs(field))-...
                      abs(stored_field)/max(abs(stored_field))) < mindiff)
                 mindiff = norm(abs(field)/max(abs(field))-...
                      abs(stored_field)/max(abs(stored_field)));
                 minindex = k;
               
              end
               
           end
           
           field_tracker{minindex,2} = [field_tracker{minindex,2}, k_eigs_store(j,i)];
           field_tracker{minindex,3} = [field_tracker{minindex,3}, omega_scan(i)];
           
           
        end
        
    end   
 

end


% figure()
% d = size(field_tracker);
% for i = 1:d(1)
%   subplot(121)
%   plot(real(field_tracker{i,1}));
%   hold on;
%   subplot(122)
%   plot(real(field_tracker{i,2}), field_tracker{i,3}, '.')  
%   hold on;
%   plot(-abs(imag(field_tracker{i,2})), field_tracker{i,3}, '.r')  
%   hold on;
%   drawnow();
% end