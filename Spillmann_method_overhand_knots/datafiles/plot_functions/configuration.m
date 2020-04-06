clear all;
n = 2;

h = 1.6e-3;
E = 1.8e6;
h = 1.6e-3;
EI = E * pi * h^4/4;

data = importdata('simDER0.01.txt');

% data2 = importdata('simDER0.801.txt');
% data3 = importdata('simDER0.00.txt');
data_c = importdata('simDER_v0.01.txt');


x1 =data(1:end,1);
   
y1 = data(1:end,2); %F
    
% x1 = n^2*h./x1;
% y1 = n^2*y1*h^2/EI;

% 
% hold on;
% plot(data(1:end, 6), data(1:end,2), 'r');
plot(x1, y1);
% plot(data2(1:end,6), data2(1:end,2), 'b');
% plot(data3(1:end,6), data3(1:end,2), 'g');
% 
K = length(data_c)/301;
figure
% ph = plot3(data_c(1:301, 1), data_c(1:301, 2),data_c(1:301, 3));
for i = 1:K
%     ph.XData = data_c(301*(i-1)+1:301*i, 1);
%     ph.YData = data_c(301*(i-1)+1:301*i, 2);
%     ph.ZData = data_c(301*(i-1)+1:301*i, 3);
    plot3(data_c(301*(i-1)+1:301*i, 1), data_c(301*(i-1)+1:301*i, 2), data_c(301*(i-1)+1:301*i, 3));
    axis equal;
    grid on;
    drawnow
    
end









% % 
% hold on;
% plot3(data(:,1), data(:,2), data(:,3));
% plot3(data_c(:,1), data_c(:,2), data_c(:,3));
% axis equal
