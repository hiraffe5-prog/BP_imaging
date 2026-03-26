
%% 截取
slc=y_rd;
figure,imagesc(abs(slc))
[row_pos,col_pos]=find(abs(slc)==max(abs(slc(:))));
row_len=64;
col_len=64;

slc_crop=slc((row_pos-row_len):(row_pos+row_len-1),(col_pos-col_len):(col_pos+col_len-1));
figure,imagesc(abs(slc_crop))


%% 评估
upr=64;
upa=32;
range_pixel=c/2/fs;
azimuth_pixel=PRF/Na_fft;

[ rour4dB, roua4dB, Ia , Ir] = Reso_PSLR_ISLR_4dB(slc_crop,range_pixel,azimuth_pixel,upr,upa,1);
rour4dB_theory=c/2/B;
roua4dB_theory=1/Na*PRF;

