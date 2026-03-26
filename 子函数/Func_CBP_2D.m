function [amsum,amsum_new,rp_middle] = Func_CBP_2D(s1,ParaCBP,xp2,yp2,zp2,x_pos,y_pos,z_pos)

r_start = ParaCBP.Rmin;
delta_r = ParaCBP.deltaR;
Na = ParaCBP.Na;
Nr_up = ParaCBP.NrInterp;
lambda = ParaCBP.Lambda;


h = waitbar(0,'1','Name','bp处理',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)
amsum=zeros(size(xp2));%存储幅度  
rp=zeros(1,length(xp2));    %存储与雷达的距离
for i = 1:Na

    s_temp = s1(i,:);

    rp=sqrt((xp2-x_pos(i)).^2+(yp2-y_pos(i)).^2+(zp2-z_pos(i)).^2);
    n_round=round((rp-r_start)/delta_r+1);
    n_round(n_round<1) = 1; n_round(n_round>Nr_up) = Nr_up;

    amsum=amsum+s_temp(n_round).*exp(1j*4*pi*(rp)/lambda);
     
    if getappdata(h,'canceling')
        disp('程序中止')
        delete(h)
        return
    end 
    waitbar(i/Na,h,sprintf('%6.2f %%',i/Na*100)) 
end
delete(h)
% fclose(fid);


%% 相位补偿（保相）
rp_middle=(sqrt((xp2-x_pos(Na/2+1)).^2+(yp2-y_pos(Na/2+1)).^2+(zp2-z_pos(Na/2+1)).^2));
amsum_new = amsum.*exp(-1i*pi*4*rp_middle/lambda);


