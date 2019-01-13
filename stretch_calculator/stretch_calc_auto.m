stretch = struct;
N = numel(loc);

for i = 1:N 
    [A,t] = stretch_est(loc(1).r, loc(2).r);
    [C,T] = stretch_refine(loc(1).r,loc(2).r,A,t,2,1);
    stretch(1).C = C;
    stretch(1).T = T;
end
