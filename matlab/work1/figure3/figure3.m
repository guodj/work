%% load
load fdenab.mat
%% creat axes
h_fig = figure;
set(h_fig,'position',[70   184   700   568]);
% 输入参数
h_N = 3; % 行数
w_N = 4; % 列数
h_gap = 0.028; % 行间间距
w_gap = 0.025;  % 列间间距
h_marg = [0.08 0.12]; % 上边缘, 下边缘
w_marg = [0.1 0.05]; % 左边缘, 右边缘
sub_creat_axes
%% draw pictures
denff=fdenab;
yr=[-6,6;-0.6,0.6;-0.06,0.06];
ytickl={-6:3:6;-0.6:0.3:0.6;-0.06:0.03:0.06};
xticl=-5:5;
err=1;
ran=1;
title_text={'March','September','June','December'};
tit_y_po=1.07;
ylabel_text={'\delta\rho at 250 km','\delta\rho at 400 km',...
    '\delta\rho at 550km'};
xlabel_text='Epoch Time (Days)';
rec_text={'( a )','( b )','( c )','( d )','( e )','( f )','( g )'};
leng1='Away - Toward';
leng2='Toward - Away';
lenpos=[0.41 0.008 0.25,0.045];
ran_text_for='%.2f';
color1=[0,0,1];
color2=[1,0,0];
minornum=[[1,3];[1,3];[1,3]];
minorlenf=0.4;
%%
% main program used to draw the figure.
figure3_draw
%--------------------------------------------%
%% save
set(gcf,'paperpositionmode','auto')
picdir='D:\mywork\ther_dens_sect\matlabworkspace\figure3\';
picname='figure3.eps';
print(gcf,'-depsc2',[picdir,picname]);
fix_dottedline([picdir,picname])
