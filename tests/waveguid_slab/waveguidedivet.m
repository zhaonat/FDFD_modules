%% create waveguide with divet

function eps_divet = WaveGuideDivet(eps_vac,epsilon, divetSize, width)
    N = size(eps_vac);
    xc= round(N(1)/2); yc = round(N(2)/2);
    eps_divet = eps_vac;
%     eps_divet(xc-round(width/2):xc+round(width/2),:) = epsilon;
%     eps_divet(xc+round(width/2):xc+round(width/2)+divetSize,yc:yc+divetSize) = 12;
    %% create slab
    eps_divet(:,yc-round(width/2):yc+round(width/2)) = epsilon;
    
    %% divet
    eps_divet(xc-divetSize:xc+divetSize,...
        yc+round(width/2):yc+round(width/2)+divetSize) = 12;
end