% evaluate lagged, time-averaged values in x
% simple expressions may not be used (memory overflow)
function y = LagWtdAve( x,n,deltat,tau )
y=x;
for i=1:length(x)
    bs=(abs(i-n)+(i-n))/2+1; es=i;
    weight=exp(((bs:es)-i) *deltat/tau);
    weight=reshape(weight,size(x(bs:es)));
    y(i)=(nansum(x(bs:es).*weight))/sum(weight);
end

end