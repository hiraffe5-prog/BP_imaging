function [pos_list]=scatter_pos(L_vec,B_vec)
%建立
L_vec = L_vec*180/pi;
B_vec = B_vec*180/pi;



Tar_sel=[-11,-43;
         -11,-42.5;
         -11,-43.5;
         -11.5,-43;
         -11.5,-42.5;
         -11.5,-43.5;
         -10.5,-43;
         -10.5,-42.5;
         -10.5,-43.5;
         -10,-43;
         -10,-42.5;
         -10,-43.5;
         -12,-43;
         -12,-42.5;
         -12,-43.5 ];

Tar_Num = size(Tar_sel,1);
        
for ii = 1:Tar_Num
 n_ind(ii) =  round(abs(Tar_sel(ii,2) - B_vec(1))/abs(B_vec(1)-B_vec(2))) + 1;
 m_ind(ii) =  round(abs(Tar_sel(ii,1) - L_vec(1))/abs(L_vec(1)-L_vec(2))) + 1;
end

for ii = 1:Tar_Num
  pos_list(ii) = (n_ind(ii)-1)*length(L_vec)+m_ind(ii);
end

