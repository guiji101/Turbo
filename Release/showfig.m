clear
clc
%-0.5:0.5:2
log=load('log.txt');
Npoints=log(1);
EbN0=log(2:end);
figure
for i=1:Npoints
    str=['out' num2str(i) '.txt'];
    data=load(str);
    plot(data(:,1),data(:,2),'-o','linewidth',2),hold on;
end
grid on