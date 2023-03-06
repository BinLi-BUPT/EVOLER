function [Initial_index,p_max,p_min] = Search_30D_Initial(popmax,popmin,s,Discret_leng)
D = length(popmax);
for ii = 1 : D
    all_axis(ii,:) = linspace(popmin(ii),popmax(ii),Discret_leng);
end
Index_sample =(round(Discret_leng/2));               %  index
index_temp = (round(Discret_leng/2));

%% test accuracy
for ii = 1 : Discret_leng
    for jj = 1 : Discret_leng
        for kk = 1 : Discret_leng
            for oo = 1 : Discret_leng
                for pp = 1 : Discret_leng
                    for ii1 = 1 : Discret_leng
                        for jj1 = 1 : s
                            for kk1 = 1 : s
                                for oo1 = 1 : s
                                    for pp1 = 1 : s
                                        ind_cal = [all_axis(1,ii), all_axis(2,jj), all_axis(3,kk),all_axis(4,oo), all_axis(5,pp), all_axis(6,ii1), all_axis(7,Index_sample(jj1)),all_axis(8:12,index_temp)',...
                                            all_axis(13,Index_sample(kk1)),all_axis(14:18,index_temp)',all_axis(19,Index_sample(oo1)),all_axis(20:24,index_temp)',all_axis(25,Index_sample(pp1)),all_axis(26:30,index_temp)'];
                                        A1_test(ii,jj,kk,oo,pp,ii1,jj1,kk1,oo1,pp1) = hybrid_func2(ind_cal );
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
[m_all_1]=find(tensor(A1_test) == collapse(tensor(A1_test),[1,2,3,4,5,6],@min));

%% 2
for ii = 1 : s 
for jj = 1 : Discret_leng
for kk = 1 : Discret_leng
for oo = 1 : Discret_leng
for pp = 1 : Discret_leng
for ii1 = 1 : Discret_leng
for jj1 = 1 : Discret_leng
for kk1 = 1 : s
for oo1 = 1 : s
for pp1 = 1 : s
    ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2:6,index_temp)', all_axis(7,jj), all_axis(8,kk),all_axis(9,oo), all_axis(10,pp), all_axis(11,ii1), all_axis(12,(jj1)), all_axis(13,Index_sample(kk1)),...
      all_axis(14:18,index_temp)',all_axis(19,Index_sample(oo1)),all_axis(20:24,index_temp)',all_axis(25,Index_sample(pp1)),all_axis(26:30,index_temp)'];    
    A2_test(ii,jj,kk,oo,pp,ii1,jj1,kk1,oo1,pp1) = hybrid_func2(ind_cal );
end
end
end
end
end
end
end
end
end
end

[m_all_2]=find(tensor(A2_test) == collapse(tensor(A2_test),[1,2,3,4,5,6,7],@min));
%% 3
for ii = 1 : s
    for jj = 1 : s
        for kk = 1 : Discret_leng
            for oo = 1 : Discret_leng
                for pp = 1 : Discret_leng
                    for ii1 = 1 : Discret_leng
                        for jj1 = 1 : Discret_leng
                            for kk1 = 1 : Discret_leng
                                for oo1 = 1 : s
                                    for pp1 = 1 : s
                                        ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2:6,index_temp)', all_axis(7,Index_sample(jj)),all_axis(8:12,index_temp)', all_axis(13,kk), all_axis(14,oo),all_axis(15,pp),...
                                            all_axis(16,ii1),all_axis(17,jj1),all_axis(18,kk1),all_axis(19,Index_sample(oo1)),all_axis(20:24,index_temp)',all_axis(25,Index_sample(pp1)),all_axis(26:30,index_temp)'];
                                        A3_test(ii,jj,kk,oo,pp,ii1,jj1,kk1,oo1,pp1) = hybrid_func2(ind_cal );
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
[m_all_3]=find(tensor(A3_test) == collapse(tensor(A3_test),[1,2,3,4,5,6,7,8],@min));

%% 4
for ii = 1 : s
    for jj = 1 : s
        for kk = 1 : s
            for oo = 1 : Discret_leng
                for pp = 1 : Discret_leng
                    for ii1 = 1 : Discret_leng
                        for jj1 = 1 : Discret_leng
                            for kk1 = 1 : Discret_leng
                                for oo1 = 1 : Discret_leng
                                    for pp1 = 1 : s
                                        ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2:6,index_temp)', all_axis(7,Index_sample(jj)), all_axis(8:12,index_temp)', all_axis(13,Index_sample(kk)), all_axis(14:18,index_temp)', all_axis(19,oo),all_axis(20,pp),...
                                            all_axis(21,ii1),all_axis(22,jj1),all_axis(23,kk1),all_axis(24,(oo1)),all_axis(25,Index_sample(pp1)),all_axis(26:30,index_temp)'];
                                        A4_test(ii,jj,kk,oo,pp,ii1,jj1,kk1,oo1,pp1) = hybrid_func2(ind_cal );
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
[m_all_4]=find(tensor(A4_test) == collapse(tensor(A4_test),[1,2,3,4,5,6,7,8,9],@min));

%% 54
for ii = 1 : s
    for jj = 1 : s
        for kk = 1 : s
            for oo = 1 : s
                for pp = 1 : Discret_leng
                    for ii1 = 1 : Discret_leng
                        for jj1 = 1 : Discret_leng
                            for kk1 = 1 : Discret_leng
                                for oo1 = 1 : Discret_leng
                                    for pp1 = 1 : Discret_leng
                                        ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2:6,index_temp)', all_axis(7,Index_sample(jj)), all_axis(8:12,index_temp)', all_axis(13,Index_sample(kk)), all_axis(14:18,index_temp)', all_axis(19,Index_sample(oo)), all_axis(20:24,index_temp)',all_axis(25,pp),...
                                            all_axis(26,ii1),all_axis(27,jj1),all_axis(28,kk1),all_axis(29,(oo1)),all_axis(30,(pp1))];
                                        A5_test(ii,jj,kk,oo,pp,ii1,jj1,kk1,oo1,pp1) = hybrid_func2(ind_cal );
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
[m_all_5]=find(tensor(A5_test) == collapse(tensor(A5_test),[1,2,3,4,5,6,7,8,9,10],@min));

Index_fro = [m_all_1(1,1:6) m_all_2(1,2:7) m_all_3(1,3:8) m_all_4(1,4:9) m_all_5(1,5:10)];
% Index_fro = [m_all_1(1,1:6) _all_2(1,2:7) _all_3(1,3:8) _all_4(1,4:9) _all_5(1,5:10)];
% Initial_min = double(Result_cal(A51_tensor,A52_tensor,A53_tensor,A54_tensor,A55_tensor,R_core,s,Index_fro));
for ii = 1 : D
    p_max(ii) =all_axis(ii,Index_fro(ii)) +(all_axis(ii,2)-all_axis(ii,1))/2;
    p_min(ii) =all_axis(ii,Index_fro(ii)) -(all_axis(ii,2)-all_axis(ii,1))/2;
    Initial_index(ii) = all_axis(ii,Index_fro(ii));
end


end