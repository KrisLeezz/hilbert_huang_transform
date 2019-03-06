function IMF=myexample(path)
%% 
% fs=1000; % should be > 10 times of the maximum frequency of interest
% t=0:1/fs:10-1/fs;
% 
% f1=5;
% f2=20;
% X=sin(2*pi*f1*t);
% Y=sin(2*pi*f1*t) + cos(2*pi*f2*t);
% X synchronises with Y at f1 (5Hz)
test=load(path);
data=struct2table(test);
data2=table2array(data);
x=data2(:,1);
% X=table2array(x);
X=x;
T=size(X,1);
y=data2(:,2);
% Y=table2array(y);
Y=y;
t=linspace(0,T,T);
%% Perform NA-MEMD decomposition
% No_noise_channels=1;
% 
% noise_process=zeros(length(X),No_noise_channels);
% 
% Noise_power=mean([var(X) var(Y)]);
% 
% for i=1:No_noise_channels
%     noise_process(:,i)=randn(length(X),1).*sqrt(Noise_power);
% end

memd_input=[X Y];%
disp(size(memd_input));
IMF=memd(memd_input);

%% Plot IMFs
column=2;
row=size(IMF,2);%返回列数也就是IMF的个数

figure('name','IMFs');
subplot(row+1,column,1);
plot(t,memd_input(:,1));
title('X');
ylabel('Magnitude');
axis tight;
%  ylim([-1 1]);

subplot(row+1,column,2);
plot(t,memd_input(:,2));
title('Y');
ylabel('Magnitude');
axis tight;
% ylim([-2 2]);

% subplot(row+1,column,3);plot(t,memd_input(:,3));
% title('WGN');
% ylabel('Magnitude');
% axis tight;
% ylim([-1 1]);

for i=1:column
    for j=1:row

        subplot(row+1,column,((j)*column)+i);

        if i==1
            plot(t,squeeze(IMF(1,j,:)));
        elseif i==2
            plot(t,squeeze(IMF(2,j,:)));   
        elseif i==3
            plot(t,squeeze(IMF(3,j,:)));
        end

        ylabel(num2str(j));

        if j==size(IMF,2)
            xlabel('Time (s)');
        end

        axis tight;
%         ylim([-1 1]);

    end
end
