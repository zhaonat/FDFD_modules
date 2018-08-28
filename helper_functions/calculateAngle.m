function [lambda] = calculateAngle(App, Avv, Avp, Apv);
    lambda = max(abs(eig(full((App\Apv)*(Avv\Avp)))));
end