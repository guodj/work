%% load
load fdenreCIR.mat
%% creat axes
h_fig = figure;
set(h_fig,'position',[151    46   521   748]);
% 输入参数
h_N = 5; % 行数
w_N = 2; % 列数
h_gap = 0.03; % 行间间距
w_gap = 0.02;  % 列间间距
h_marg = [0.05 0.12]; % 上边缘, 下边缘
w_marg = [0.12 0.08]; % 左边缘, 右边缘
sub_creat_axes
%% draw pictures
denff=fdenreCIR;
yr=[-1.5,1.5;350,550;-20,20;-20,20;-20,20];
ytickl={-1.5:1.5:1.5;350:100:550;-20:10:20;-20:10:20;-20:10:20};
xticl=-5:5;
title_text={'SBs with CIR','SBs without CIR'};
tit_y_po=1.07;
ylabel_text={'B_{z,GSM} (nT)', 'V_{SW} (km/s)','\delta\rho at 250 km (%)',...
    '\delta\rho at 400 km (%)','\delta\rho at 550 km (%)'};
xlabel_text='Epoch Time (Days)';
rec_text={'( a )','( b )','( c )','( d )','( e )','( f )','( g )'};
lenpos=[0.4 0 0.25,0.045];
leng1='Ineffective - Effective';
leng2='Effective - Ineffective';
err=1;
ran=1;
ran_text_for='%.0f%%';
color1=[0,0.5,0];
color2=[0,0,0];
minornum=[[1,5];[1,5];[1,5];[1,5];[1,5]];
minorlenf=0.5;
figure6_draw
%% save
set(gcf,'paperpositionmode','auto')
picdir='D:\mywork\ther_dens_sect\matlabworkspace\figure6\';
picname='figure6.eps';
print (gcf,'-depsc2', [picdir,picname]) 
% saveas(gcf,[picdir,picname],'psc2');
fix_dottedline([picdir,picname])
