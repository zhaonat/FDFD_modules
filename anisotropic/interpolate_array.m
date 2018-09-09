function [R] = interpolate_array(epsilon, w, s)  
    sign = diff('bf' == s);  % +1 for s=='f'; -1 for s=='b'
    R = (epsilon+circshift(epsilon, -sign * ('xyz' == w)))/2;

end