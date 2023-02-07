% Function that checks if x is in the the sphere
function [c, ceq] = SphereConstraints(x)
    c = [];
    ceq = dot(x,x)-1;
end