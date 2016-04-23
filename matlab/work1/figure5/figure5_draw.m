% used for figure 2-6, this is only valid in this study.
% it demands the plotted values to place in a cell(figrow,figcolumn,SBtyp)
for irow=1:h_N
    for icol=1:w_N
        axes(h_axes((irow-1)*w_N+icol))
        y1=denff{irow,icol,1};
        y2=denff{irow,icol,2};
        y1_num=num2str(length(y1));
        y2_num=num2str(length(y2));
        y1m=prctile(y1,[25,50,75],1);
        y2m=prctile(y2,[25,50,75],1);
        %-------------------------------------------------
        hold on,box on
        myplot(-5:5,y1m(2,:),'color',color1);
        myplot(-5:5,y2m(2,:),'color',color2,'linestyle','-.');
        %-------------------------------------------------
        set(gca,'ylim',yr(irow,:),'ytick',ytickl{irow,1},'yticklabel',[],...
            'xtick',xticl,'xticklabel',[],'ticklength',[0.05,0.03],...
            'fontsize',11,'xgrid','on','ygrid','on')
        
        %% errorbar
        if err==1
            y1erruplen=mean(y1m(3,:)-y1m(2,:));
            y1errdolen=mean(y1m(2,:)-y1m(1,:));
            y1errpo=[3,y1m(2,3+6)];
            y2erruplen=mean(y2m(3,:)-y2m(2,:));
            y2errdolen=mean(y2m(2,:)-y2m(1,:));
            y2errpo=[4,y2m(2,4+6)];
            herr=errorbar(y1errpo(1),y1errpo(2),y1errdolen,y1erruplen,...
                'color',color1,'linewidth',1.5);
            errorbar_tick(herr,30)
            herr=errorbar(y2errpo(1),y2errpo(2),y2errdolen,y2erruplen,...
                'color',color2,'linewidth',1.5);
            errorbar_tick(herr,30)
        end
        %% range
        if ran==1
            ht=plot_text(0.43,0.85,num2str(range(y1m(2,:)),...
                ran_text_for),12);
            set(ht,'color',color1,'horizontalalign','right')
            ht=plot_text(0.43,0.7,num2str(range(y2m(2,:)),...
                ran_text_for),12);
            set(ht,'color',color2,'horizontalalign','right')
        end
        %% set picture
        %-----------------------------------------------------------------%
        if irow==1
            mytitle_3c(gca,title_text{icol},y1_num,y2_num,...
                color1,color2,tit_y_po);
        end
        if icol==1
            set(gca,'yticklabel',ytickl{irow,1})
            ylabel(ylabel_text{irow},'fontsize',10);
        end
        if irow==h_N
            set(gca,'xticklabel',xticl)
            xlabel(xlabel_text,'fontsize',10)
        end
        if icol==w_N
            plot_text(1.05,0.5,rec_text{irow})
        end
        a=legend(h_axes(1),leng1,leng2);
        set(a,'box','off','position',lenpos,'orientation','horizontal')
        minorxy(minornum(irow,:),minorlenf)
        set(gca,'xcolor','k')
    end
end
