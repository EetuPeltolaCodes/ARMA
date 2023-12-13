%% Exercise 1
clearvars
close all
clc

data=readmatrix("SeriesReport-202212160449-V.xlsx");
data=data(1:end-12,2);
% Remove last year because it isn't full data

subplot(2,1,1)
plot(data),axis tight
legend("Original data")
grid minor

Data=data(end-11*12+1:end);
% We have a large drop towards the end so we cut the data from there

subplot(2,1,2)
plot(Data),axis tight,hold on

% We can see a clear upwards trend with the data also there is clear
% seasonality because it is monthly data

X=[ones(length(Data),1) [1:length(Data)]'];
Y=Data;
b=X\Y

% Now we can calculate the trendline and see the trend.

trendline=X*b;
plot(trendline)
legend("Cut data","Trendline")
grid minor

dData=Data-trendline;
% Detrend the data

p=12;

for i=1:p
    SI(i) = mean(dData(i:p:end));
end
SI
S=repmat(SI',[11 1]);
% Repeat the plot for 11 years

figure(2)
plot(dData,"b"),hold on
plot(S,"r"),axis tight
legend("Detrended data","Seasonal index"),hold off
grid minor

I=dData-S;
% Calculate the irregular part

figure(3)
subplot(2,1,1)
plot(I),axis tight
grid minor
legend("Irregular part of the data")

h=kpsstest(I)
% Test isn't good so we difference the data
Id=diff(I);

subplot(2,1,2)
plot(dData,"b"),hold on
plot(Id,"r"),axis tight,hold off,grid minor
legend("Detrended data","Differenced irregular data")

h=kpsstest(Id)
% Now the test is good

figure(4)
subplot(2,1,1)
autocorr(Id,11)
subplot(2,1,2)
parcorr(Id,11)
% Auto correlation and partial correlation shows that points are inside the
% significant levels

% Forecasting trendline

Xf=[ones(p,1),(length(Data)+1:length(Data)+p)'];
Tf=Xf*b;

Mf=Tf+SI';

figure(5)
plot(Data,"b"),hold on
plot(length(Data):length(Data)+p-1,Mf,"r")
legend("Original data","Forecast data"),hold off,grid minor

% Arma
Models=[];
% Looking at the ACF and PACF plots I would think the smallest value is
% somewhere around 10 so lets loop from 7 to 11
for i=7:11
    for j=7:11
        M=armax(Id,[i j]);
        a=aic(M);
        Models=[Models; i j a];
    end
end
m=find(Models(:,end)==min(Models(:,end)));
Models(m,:)

ar=Models(m,1);
ma=Models(m,2);

ModelFit=armax(Id,[ar ma]);

Idf=forecast(ModelFit,Id,p);

% Forecast for irregular data
figure(6)
subplot(2,1,1)
plot(Id,"b"),hold on
plot(length(Id):length(Id)+p-1,Idf,"r"),axis tight,hold off,grid minor
legend("Differenced irregular data","Forecast for the irregular data")

Xf=(length(Data)+1:length(Data)+p)';

% Forecast for whole data
If(1)=I(end)+Idf(1);
for i=2:p
    If(i)=If(i-1)+Idf(i);
end

Forecast=Tf+SI'+If';
Forecast,format shortG

subplot(2,1,2)
plot(Data,"b"),hold on
plot(length(Data):length(Data)+p-1,Forecast,"r"),grid minor
legend("Original data","Forecast")

%% Exercise 2
clearvars
close all
clc

data=readmatrix("SeriesReport-202212190329-V.xlsx");
data=data(1:end-12,2);
% Remove the first year and the last year of data because they aren't full
% So we have 81 years of data

subplot(2,1,1)
plot(data),axis tight,grid minor
legend("Original data")

Data=data(end-13*12+1:end);
% We have a large drop towards the end so we cut the data from there

subplot(2,1,2)
plot(Data),axis tight,hold on

X=[ones(length(Data),1) [1:length(Data)]'];
Y=Data;
b=X\Y
% Now we can calculate the trendline and see the trend.

trendline=X*b;
plot(trendline)
legend("Cut data","Trendline")
grid minor

dData=Data-trendline;
% Detrend the data

p=12;
for i=1:p
    SI(i) = mean(dData(i:p:end));
end
SI
S=repmat(SI',[13 1]);
% Repeat the plot for 13 years

figure(2)
plot(dData,"b"),hold on
plot(S,"r"),axis tight
legend("Detrended data","Seasonal index"),hold off,grid minor

I=dData-S;
% Calculate the irregular part

figure(3)
subplot(2,1,1)
plot(I),axis tight,grid minor
legend("Irregular part of the data")

h=kpsstest(I)
% Test isn't good so we difference the data
Id=diff(I);

subplot(2,1,2)
plot(dData,"b"),hold on
plot(Id,"r"),axis tight,hold off,grid minor
legend("Detrended data","Differenced irregular data")

h=kpsstest(Id)
% Now the test is good

figure(4)
subplot(2,1,1)
autocorr(Id,13)
subplot(2,1,2)
parcorr(Id,13)
% Auto correlation and partial correlation aren't the greatest.

% Forecasting 
% Forecasting the trendline

Xf=[ones(p,1) (length(Data)+1:length(Data)+p)'];
Tf=Xf*b;

Mf=Tf+SI';

figure(5)
plot(data,"b"),hold on
plot(length(data):length(data)+p-1,Mf,"r")
legend("Original data","Forecast","Location","best"),hold off,grid minor

% Arma
Models=[];
for i=0:4
    for j=0:4
        M=armax(Id,[i j]);
        a=aic(M);
        Models=[Models; i j a];
    end
end
m=find(Models(:,end)==min(Models(:,end)));
Models(m,:)

ar=Models(m,1);
ma=Models(m,2);

ModelFit=armax(Id,[ar ma]);

Idf=forecast(ModelFit,Id,p);

% Forecast for irregular data
figure(6)
subplot(2,1,1)
plot(Id,"b"),hold on
plot(length(Id):length(Id)+p-1,Idf,"r"),axis tight,hold off,grid minor
legend("Differenced irregular data","Forecast for the irregular data")

Xf=(length(Data)+1:length(Data)+p)';

% Forecast for whole data
If(1)=I(end)+Idf(1);
for i=2:p
    If(i)=If(i-1)+Idf(i);
end

Forecast=Tf+SI'+If';
Forecast
format shortG

subplot(2,1,2)
plot(Data,"b"),hold on
plot(length(Data):length(Data)+p-1,Forecast,"r")
legend("Original data","Forecast","Location","best"),grid minor