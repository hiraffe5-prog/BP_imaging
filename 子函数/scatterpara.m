function [incidence,sigma]=scatterpara(R,pos_xyz,dem,L_vec,B_vec)
%建立
[row,col]=size(dem);
Tar_xyz=zeros(row,col,3);
for n=1:row
    for m=1:col
        %Tar((n-1)*length(L_vec)+m,:)=[L_vec(m),B_vec(n)];
        Tar_xyz(n,m,:)=(R+dem(n,m))*...
            [cos(B_vec(n)).*cos(L_vec(m)),cos(B_vec(n)).*sin(L_vec(m)),sin(B_vec(n))];
    end
end


incidence=zeros(row,col);

h1=waitbar(0,"散射系数计算中...");
for n=1:row-1
     waitbar(n/row);
    for m=1:col-1
        a=Tar_xyz(n,m,:);
        b=Tar_xyz(n,m+1,:);
        c=Tar_xyz(n+1,m,:);
        a=a(:);b=b(:);c=c(:);
        normvec=cross(a-b,a-c);
        normvec=sign(normvec(1))*normvec;%这里保证面的法向量朝外
        rlos=pos_xyz(:,1)-a;
        incidence(n,m)=acos(normvec.'*rlos/norm(rlos)/norm(normvec));
    end
end
close(h1)

% if(exist('Random_Phase.mat'))
% load Random_Phase.mat
% else
% Random_Phase=2*pi*rand(row,col);
% save Random_Phase Random_Phase
% end

disp('生成随机相位')
Random_Phase=2*pi*rand(row,col);
sigma=cos(incidence) .* exp(1j.*Random_Phase);

end