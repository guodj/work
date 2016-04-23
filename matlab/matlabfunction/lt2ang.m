function ang = lt2ang( lt )
% LT2THETA converts local time to angle in a polar coordinate
% input: lt, local time
% output: ang, the angle from x axis
ang=lt;
aa=(lt>=6 & lt<=24);
bb=(lt<=6);
ang(aa)=( (lt(aa)-6)/18 ) * (1.5*pi);
ang(bb)=( lt(bb)/6 ) * (0.5*pi) + 1.5*pi;


end

