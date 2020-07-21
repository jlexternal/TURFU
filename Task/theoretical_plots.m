% Repetition curve theoretical interpretation

xAxis = [-8:.1:8];

%equal 
par0(1) = 0;
par0(2) = 1;
p0 = par0(1)+par0(2).*xAxis;
p0 = 1./(1+exp(-p0));

%higher slope
par0_hi(1) = 0;
par0_hi(2) = 3;
p0_hi = par0_hi(1)+par0_hi(2).*xAxis;
p0_hi = 1./(1+exp(-p0_hi));

%lower slope
par0_lo(1) = 0;
par0_lo(2) = .4;
p0_lo = par0_lo(1)+par0_lo(2).*xAxis;
p0_lo = 1./(1+exp(-p0_lo));

%negative horizontal shift
par1(1) = -3;
par1(2) = 1;
p1 = par1(1)+par1(2).*xAxis;
p1 = 1./(1+exp(-p1));

%positive horizontal shift
par2(1) = 3;
par2(2) = 1;
p2 = par2(1)+par2(2).*xAxis;
p2 = 1./(1+exp(-p2));

%plots
close all;
clf;
figure(1);
plot(xAxis, p0,     'LineWidth', 2, 'Color', [0 0 0]);
hold on;
plot(xAxis, p1,     'LineWidth', 2, 'Color', [0 .7 0]);
plot(xAxis, p2,     'LineWidth', 2, 'Color', [.7 0 0]);
xline(0,':');
yline(.5,'--');
hold off;
%%
figure(2);
xAxis2 = [-4:.1:4];
plot(xAxis2, p0([find(xAxis==-4):find(xAxis==4)]),     'LineWidth', 2, 'Color', [0 0 0]+.5);
hold on;
plot(xAxis2, p0_lo([find(xAxis==-4):find(xAxis==4)]),  'LineWidth', 2, 'Color', [0 0 0]+.75);
plot(xAxis2, p0_hi([find(xAxis==-4):find(xAxis==4)]),  'LineWidth', 2, 'Color', [0 0 0]+.25);
xline(0,':');
yline(.5,'--');
hold off;



