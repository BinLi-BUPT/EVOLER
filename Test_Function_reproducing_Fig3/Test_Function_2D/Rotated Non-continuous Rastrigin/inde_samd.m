function c_temp = inde_samd(R_tem3,C_sapm,rat_sig)
deee = double(R_tem3);
gap_err = [];
for ii = 1 : size(deee)
    inde_now =  deee(ii,:);
    gap_err(ii) = sum(sum(abs(inde_now - R_tem3)));
    
end
gap_err(C_sapm) = 0;
[~,ind_all] = sort(gap_err,'descend');
samp_c11 = ind_all(1:rat_sig-length(C_sapm))';
c_temp = unique([C_sapm,samp_c11]);
end