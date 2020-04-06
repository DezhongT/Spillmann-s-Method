data = importdata('tdz.txt');

data1 = importdata('simDER.txt');

hold on;
plot3(data1(:,1),data1(:,2),data1(:,3), 'r');
plot3(data(:,1),data(:,2),data(:,3));
axis equal;