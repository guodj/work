function z=mygriddata(x0,y0,z0,x,y)
    % A simple replacement of griddata
    % by Dongjie, 15-9-7
    ux=unique(x,'sorted');
    uy=unique(y,'sorted');
    x0=gridxy(x0,ux);
    y0=gridxy(y0,uy);
    % added on 15-11-1, by Guo
    index1=x0<ux(1) | x0>ux(end);
    index2=y0<uy(1) | y0>uy(end);
    index=index1 | index2;
    x0(index)=[];
    y0(index)=[];
    z0(index)=[];
    %-----------------------------------------
    
    [~,index]=sortrows([x0,y0],[1,2]);
    x0=x0(index);
    y0=y0(index);
    z0=z0(index);

    [~,index,~]=unique([x0,y0],'rows','stable');
    ux0=x0(index);
    uy0=y0(index);
    z=ux0*nan;

    edge=find(diff(y0));
    edge=[0;edge;length(y0)];
    for k=1:(length(edge)-1)
        z(k)=nanmean(z0(edge(k)+1:edge(k+1)));
    end

    if length(ux0)<numel(x)
        xt=reshape(x,numel(x),1);
        yt=reshape(y,numel(y),1);
        zt=xt*nan;
        index=ismember([xt,yt],[ux0,uy0],'rows');
        zt(index)=z;
        z=zt;
    end

    z=reshape(z,size(x));
end

function gridx0=gridxy(x0,x)
    % It is requied that x is in increasing order. 
    gridx0=x0;
    xi=abs(x(2)-x(1));
    % Modified on 15-11-1, by Guo
    gridx0=x(1)+round((x0-x(1))/xi)*xi; % x(i) should be a decimal with limited digits.
    % --------------------------------------------
end        
