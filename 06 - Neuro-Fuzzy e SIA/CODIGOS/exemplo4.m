clc
clear all


numPts=51;
x=linspace(-10,10,numPts)';

y=-2*x-x.^2;
data=[x y];
trndata=data(1:2:numPts,:);
chkdata=data(2:2:numPts,:);

%Set the number and type of membership functions:
% numMFs=5;
% mfType='gbellmf';
% % Generate the FIS-matrix and execute the ANFIS-training by 40 rounds:
% fismat=(genfis1(trndata,numMFs,mfType))
% numEpochs=40;
% 
% [fismat1,trnErr,ss,fismat2,chkErr]=anfis(trndata,fismat,numEpochs,NaN,chkdata);


numMFs=15;
mfType='gbellmf';
fismat=(genfis1(trndata,numMFs,mfType))
numEpochs=40;

[fismat1,trnErr,ss,fismat2,chkErr]=anfis(trndata,fismat,numEpochs,NaN,chkdata);