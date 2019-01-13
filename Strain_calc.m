N = length(params);

for i = 1:N
    strain(i).r = params(i).pks(:,1:2);
end

for j = 3:N-1
   [A,t] = stretch_est(strain(j).r, strain(j+1).r);
   [strain(j).C,strain(j).T] = stretch_refine(strain(j).r,strain(j+1).r,A,t,2,1);
end

save 180810_stretch_cal_set2 strain params