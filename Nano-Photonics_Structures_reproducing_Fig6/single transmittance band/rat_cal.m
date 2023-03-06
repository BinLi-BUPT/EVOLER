function rank_est = rat_cal(R_tem,tao)
deee = double(R_tem);
[~,sdd,~] = svd(deee,'econ');
sig_val = diag(sdd);
rat_sig = [];
for ii = 1 : length(sig_val)
    rat_sig(ii) = sum(sig_val(1:ii))/sum(sig_val);
end
[~,ind1] = find(rat_sig > tao);
rank_est = ind1(1);
end