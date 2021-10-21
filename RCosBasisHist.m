function Basis=RCosBasisHist(hist,input)
dt=1e-4;
% func=@(x,c,dc)(cos(max(-pi,min(pi,(x-c)*pi/dc/2)))+1)/2;
RCos= @ (a,c,fie,t) .5*cos(max(fie-pi,min(fie+pi,a'*log(dt*t+c)))-fie)+.5;
t=1:hist.delay;
fie=pi;
c=1;

%Basis=[RCos([2000,1000,500],c,0,t);RCos(1000,c,pi/2,t);RCos(700,c,1.2*pi,t)];
Basis=[RCos([2000],c,1*pi,t);RCos([2000],c,1.5*pi,t)];


% for f=hist.freq
%     for i=hist.delay*[0,.25,.5,.75]
%         Basis(end+1,:)=(func(t,i,f));
%     end
% end
% Basis=Basis-mean(Basis,2);
% Basis=(normalize(Basis,2));