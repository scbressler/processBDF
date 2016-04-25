function plotN100Lines(N1)
if nargin<1
    N1 = 0.1; % Auditory N1 = 0.100, Visual P1 = 0.175;
end

sLag = 0.75:0.75:3;
sLed = 0.6:0.6:3;
sCen = 0:3;


t13 = [sLag(1:end-1);sLag(1:end-1);sLag(2:end)];
t14 = [sLed(1:end-1);sLed(1:end-1);sLed(2:end)];
t1n = [sCen(1:end-1);sCen(1:end-1);sCen(2:end)];


line(N1+t13(1:2,:),ylim,'Color','b','LineStyle','-');
line(N1+t14(1:2,:),ylim,'Color','r','LineStyle','-');
line(N1+t1n(1:2,:),ylim,'Color',0.6*ones(1,3),'LineStyle','--');
end