%% load
load fdenres.mat
%% creat axes
h_fig = figure;
set(h_fig,'position',[70   184   700   568]);
% 输入参数
h_N = 3; % 行数
w_N = 3; % 列数
h_gap = 0.02; % 行间间距
w_gap = 0.025;  % 列间间距
h_marg = [0.08 0.13]; % 上边缘, 下边缘
w_marg = [0.1 0.06]; % 左边缘, 右边缘
sub_creat_axes
%% draw pictures
denff=fdenres;
yr=[-20,20;-20,20;-20,20];
ytickl={-20:10:20;-20:10:20;-20:10:20};
xticl=-5:5;
title_text={'Solar Minimum','Solar Medium','Solar Maximum'};
tit_y_po=1.07;
ylabel_text={'\delta\rho at 250 km (%)','\delta\rho at 400 km (%)',...
    '\delta\rho at 550 km (%)'};
xlabel_text='Epoch Time (Days)';
rec_text={'( a )','( b )','( c )','( d )','( e )','( f )','( g )'};
% leng1='\textbf{N -- S SBs}';
% leng2='\textbf{S -- N SBs}';
lenpos=[0.41 0.015 0.25,0.045];
leng1='Ineffective - Effective';
leng2='Effective - Ineffective';
err=1;
ran=1;
ran_text_for='%.0f%%';
color1=[0 0.5 0];% green
color2=[0 0 0];% black
minornum=[[1,5];[1,5];[1,5]];
minorlenf=0.5;
figure5_draw
%% save
set(gcf,'paperpositionmode','auto')
picdir='D:\mywork\ther_dens_sect\matlabworkspace\figure5\';
picname='figure5.eps';
print(gcf,'-depsc2',[picdir,picname]);
fix_dottedline([picdir,picname])
