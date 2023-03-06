function c_temp = inde_samd(R_tem3,C_sapm,rat_sig)
deee = double(R_tem3);
diff_dee = var(deee,0,2);
% diff_dee2 = diff(deee,1);
% ssd = sum(abs(diff_dee2),2);


diff_dee2 = diff_dee;
diff_dee2(C_sapm) = 0;
[~,ind_all] = sort(diff_dee2,'ascend');
samp_c11 = ind_all(end-(rat_sig-length(C_sapm)-1):end)';
c_temp = unique([C_sapm,samp_c11]);
end