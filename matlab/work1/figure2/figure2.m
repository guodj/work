% solar parameters and ap variation.
%% load
load fig2dat.mat
%% creat axes
h_fig = figure;
set(h_fig,'position',[1034  71  789 910]);
% create axes
h_N = 7; % row number
w_N = 4; % column number
h_gap = 0.02; % gap between rows
w_gap = 0.025;  % gap between columns
h_marg = [0.04 0.11]; % left and right margins
w_marg = [0.1 0.06]; % top and botome margins
sub_creat_axes
%% draw pictures
denff=fpic;
yr=[-4,4;-4,4;-1.5,1.5;-1.5,1.5;350,550;0,20;80,160];
ytickl={-4:2:4;-4:2:4;-1.5:1.5:1.5;-1.5:1.5:1.5;...
    350:100:550;0:10:20;80:40:160};
xticl=-5:5;
title_text={'March','September','June','December'};
tit_y_po=1.07;
ylabel_text={'B_{x,GSE} (nT)','B_{y,GSE} (nT)','B_{z,GSE} (nT)',...
    'B_{z,GSM} (nT)','V_{SW} (km/s)','Ap','F107'};
xlabel_text='Epoch Time (Days)';
rec_text={'( a )','( b )','( c )','( d )','( e )','( f )','( g )'};
leng1='Away - Toward';
leng2='Toward - Away';
lenpos=[0.4, 0, 0.25,0.045];
err=1;
ran=0;
% ran_text_for='%.2f';
color1=[0,0,1];
color2=[1,0,0];
minornum=[[1,5];[1,5];[1,5];[1,5];[1,5];[1,5];[1,4]];
minorlenf=0.5;
figure2_draw
%% save
set(gcf,'paperpositionmode','auto')
picdir='D:\mywork\ther_dens_sect\matlabworkspace\figure2\';
picname='figure2.eps';
print(gcf,'-depsc2',[picdir,picname]);
fix_dottedline([picdir,picname])
