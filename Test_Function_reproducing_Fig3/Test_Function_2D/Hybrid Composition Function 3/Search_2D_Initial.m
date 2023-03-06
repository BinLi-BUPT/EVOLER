function [Initial_index,p_max,p_min] = Search_2D_Initial(popmin,popmax,Discret_leng)

for ii = 1 : length(popmin)
    all_axis(ii,:) = linspace(popmin(ii),popmax(ii),Discret_leng);
end

for iii=1:1:Discret_leng
    for jjj=1:1:Discret_leng
        Z1(iii,jjj)=hybrid_func3([all_axis(1, iii),all_axis(2, (jjj))]);  % weierstrass: a standard test function
    end
end
%----------------------------------
[m0,n0]=find(Z1 == min(min(Z1)));                   % The center of the attention subspace, i.e., the global optimum in the representation space
Index_fro = [m0,n0];
for ii = 1 : length(popmin)
    p_max(ii) =all_axis(ii,Index_fro(ii)) +(all_axis(ii,2)-all_axis(ii,1))/2;
    p_min(ii) =all_axis(ii,Index_fro(ii)) -(all_axis(ii,2)-all_axis(ii,1))/2;
    Initial_index(ii) = all_axis(ii,Index_fro(ii));
end
% sum(abs(Initial_index - shift_op))
% Initial_min - hybrid_func3(Initial_index,shift_op)
end