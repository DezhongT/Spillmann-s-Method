clear all;  %#ok<*CLSCR>
close all;

dinfo = dir('*txt');

FONT = 'Arial';
FONTSIZE = 12;

pWidth = 8 ;
pHeight = 6;

h1 = figure(1);

E = 1.8e6;
h = 1.6e-3;
EI = E * pi * h^4/4;

colpos = [210 180 40;
    211 58 50;
    145 50 65;
    80 70 93;
    45 135 105;
    140 63 45;]/255; % colors

hold on

% for K = 1:length(dinfo)
%     n = 1;
%     thisfilename = dinfo(K).name;
%     data = load(thisfilename);
%     x = data(:,4);
%     x= data(:,4)/(2*pi); %R
%     
%     y = data(:,2); %F
%     
%     x = n^2*h./x;
%     y = n^2*y*h^2/EI;
%     plot(x, y, 'Color', colpos(K,:));
%     hold on;
% end

data = importdata('simDER0.80.txt');
x1 =1 - data(1:end,6);
   
y1 = data(1:end,2);

n=2;
    
x1 = n^2*h./x1;
y1 = n^2*y1*h^2/EI;
plot(x1, y1);


hold on;

Theory

hold off

xlabel('n^2 h / e', 'Fontname', FONT,'FontSize',FONTSIZE);
ylabel('n^2 F h^2 / B', 'Fontname', FONT,'FontSize',FONTSIZE);

axis([1e-3 1e-1 1e-4 1e0]);
box on
set(gca, 'ytick', [1e-3 1e-2 1e-1 1]);
set(gca,'Fontname', FONT,'FontSize',FONTSIZE);
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'XMinorTick','on','YMinorTick','on')

set(gcf, 'PaperUnits','inches', 'PaperPosition',[0 0 pWidth pHeight], ...
    'PaperSize', [pWidth pHeight]);


