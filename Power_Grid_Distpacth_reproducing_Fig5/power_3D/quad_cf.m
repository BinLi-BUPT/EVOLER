function [y] = quad_cf(xx)
a_end = [0.001562, 0.00194, 0.00482];
b_end = [7.92, 7.85, 7.97];
c_end = [561, 310, 78];
e_end = [300, 200, 150];
f_end = [0.0315, 0.042, 0.063];
P_i_min = [100,100,50];
ddd = length(xx);
f_end2 = 0;
for ii = 1 : ddd
    f_end2 = f_end2+c_end(ii)+b_end(ii)*xx(ii)+a_end(ii)*xx(ii)^2+abs(e_end(ii)*sin(f_end(ii)*(P_i_min(ii) - xx(ii))));
end
y = f_end2;





end
