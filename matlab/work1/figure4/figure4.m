%% load
load fdenre.mat
%% creat axes
h_fig = figure;
set(h_fig,'position',[70   184   700   568]);
% �������
h_N = 3; % ����
w_N = 4; % ����
h_gap = 0.02; % �м���
w_gap = 0.025;  % �м���
h_marg = [0.08 0.13]; % �ϱ�Ե, �±�Ե
w_marg = [0.1 0.05]; % ���Ե, �ұ�Ե
sub_creat_axes
%% draw pictures
denff=fdenre;
yr=[-20,20;-20,20;-20,20];
ytickl={-20:10:20;-20:10:20;-20:10:20};
xticl=-5:5;
title_text={'March','September','June','December'};
tit_y_po=1.07;
ylabel_text={'\delta\rho at 250 km (%)','\delta\rho at 400 km (%)',...
    '\delta\rho at 550 km (%)'};
xlabel_text='Epoch Time (Days)';
rec_text={'( a )','( b )','( c )','( d )','( e )','( f )','( g )'};
lenpos=[0.41 0.015 0.25,0.045];
leng1='Away - Toward';
leng2='Toward - Away';
err=1;
ran=1;
ran_text_for='%.0f%%';
color1=[0,0,1];
color2=[1,0,0];
minornum=[[1,5];[1,5];[1,5]];
minorlenf=0.5;
figure4_draw
%% save
set(gcf,'paperpositionmode','auto')
picdir='D:\mywork\ther_dens_sect\matlabworkspace\figure4\';
picname='figure4.eps';
print(gcf,'-depsc2',[picdir,picname]);
fix_dottedline([picdir,picname])
