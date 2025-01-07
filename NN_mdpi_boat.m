function NN_mdpi_boat()

%% import paths
addpath('C:\Users\pitto1\Lavori CARE\Script\GIT\CompareImuQsBr');
addpath('C:\Users\pitto1\Lavori CARE\Script\GIT\Process input');
addpath('C:\Users\pitto1\Lavori CARE\Script\GIT\BoatMeasures');
addpath('C:\Users\pitto1\Lavori CARE\Script\GIT\PlotThings');
addpath('C:\Users\pitto1\Lavori CARE\Script\GIT\OBJECTS');
addpath('C:\Users\pitto1\Lavori CARE\Script\GIT\ARDUINO');

%% define the paths and the subjects
METHOD   = 'seconds';
% tested different methods, only seconds has been used


SOURCE   = 'ProcD_*';
CART     = 'D:\Lavori CARE HD\DATAxScripts\BoatMeasures\FIGURE\';
CARTSAVE = 'RandSearchBoat';


ALLDATAF = {[CART,'S01_COMFORT\ProcessedData\', SOURCE, '.mat']
            [CART,'S01_RATIO\ProcessedData\', SOURCE, '.mat']
            [CART,'S01_NORMAL\ProcessedData\', SOURCE, '.mat']

            [CART,'S02_COMFORT\ProcessedData\', SOURCE, '.mat']
            [CART,'S02_NORMAL\ProcessedData\', SOURCE, '.mat']
            
            [CART,'S03_NORMAL\ProcessedData\', SOURCE, '.mat']
            [CART,'S03_RATIO\ProcessedData\', SOURCE, '.mat']
            
            [CART,'S04_RATIO\ProcessedData\', SOURCE, '.mat']
            [CART,'S04_NORMAL\ProcessedData\', SOURCE, '.mat']
            [CART,'S04_COMFORT\ProcessedData\', SOURCE, '.mat']
            
            [CART,'S05_NORMAL\ProcessedData\', SOURCE, '.mat']
            [CART,'S05_COMFORT\ProcessedData\', SOURCE, '.mat']
            [CART,'S05_RATIO\ProcessedData\', SOURCE, '.mat']

            [CART,'S06_COMFORT\ProcessedData\', SOURCE, '.mat']
            [CART,'S06_NORMAL\ProcessedData\', SOURCE, '.mat']
            [CART,'S06_RATIO\ProcessedData\', SOURCE, '.mat']            

            [CART,'S07_COMFORT\ProcessedData\', SOURCE, '.mat']
            [CART,'S07_NORMAL\ProcessedData\', SOURCE, '.mat']
            [CART,'S07_RATIO\ProcessedData\', SOURCE, '.mat']
            
            [CART,'S08_NORMAL\ProcessedData\', SOURCE, '.mat']
            [CART,'S08_COMFORT\ProcessedData\', SOURCE, '.mat']
            [CART,'S08_RATIO\ProcessedData\', SOURCE, '.mat']

            [CART,'S09_COMFORT\ProcessedData\', SOURCE, '.mat']
            [CART,'S09_NORMAL\ProcessedData\', SOURCE, '.mat']
            [CART,'S09_RATIO\ProcessedData\', SOURCE, '.mat']

            [CART,'S10_NORMAL\ProcessedData\', SOURCE, '.mat']
            [CART,'S10_COMFORT\ProcessedData\', SOURCE, '.mat']
            [CART,'S10_RATIO\ProcessedData\', SOURCE, '.mat']

            [CART,'S11_COMFORT\ProcessedData\', SOURCE, '.mat']
            [CART,'S11_NORMAL\ProcessedData\', SOURCE, '.mat']
            
            [CART,'S12_NORMAL\ProcessedData\', SOURCE, '.mat']
            [CART,'S12_COMFORT\ProcessedData\', SOURCE, '.mat']
            [CART,'S12_RATIO\ProcessedData\', SOURCE, '.mat']};

nS = size(ALLDATAF,1);

% 1-height, 2-weight, 3-ID
SubInfo = [88   190     1
           88   190     1
           88   190     1
           80   184     2
           80   184     2
           88   194     3
           88   194     3
           82   184     4
           82   184     4
           82   184     4
           83   185     5
           83   185     5     
           83   185     5
           69   184     6
           69   184     6
           69   184     6
           82   193     7
           82   193     7
           82   193     7
           85   184     8
           85   184     8
           85   184     8
           74   180     9
           74   180     9
           74   180     9
           73   179     10
           73   179     10
           73   179     10
           72   181     11
           72   181     11
           72   185     12
           72   185     12
           72   185     12];     


%% cycle on the sessions to extract the segmented data
allData = {};
allTarg = {};
allSubj = [];
for s=1:nS
    [dataC, targC, IDsubj] = getDataCycle_article(ALLDATAF{s,1}, SubInfo(s,:));
    for c=1:length(dataC)
        allData{end+1} = dataC{c};
        allTarg{end+1} = targC{c};
        allSubj(end+1) = IDsubj(c);
    end
end


%% shuffle the cycles and divide into train and test
[dTrain, tTrain, ...
 dValid, tValid, ...
 dTest,  tTest,  IDtvt] = shuffleAndTrainTest_article(allData, allTarg, [80 20 0], 'boat', allSubj);


%% normalize by mean and SD
[XTrain, muX, sigmaX] = nomalizeMeanSD(dTrain);
[TTrain, muT, sigmaT] = nomalizeMeanSD(tTrain);
for n = 1:length(dValid)
    X        = (dValid{n} - muX) ./ sigmaX;
    T        = (tValid{n} - muT) ./ sigmaT;
    XValid{n} = X;
    TValid{n} = T;
end
for n = 1:length(dTest)
    X        = (dTest{n} - muX) ./ sigmaX;
    T        = (tTest{n} - muT) ./ sigmaT;
    XTest{n} = X;
    TTest{n} = T;
end


%% import the setup file for the architecture and training 
numChannels = size(dTrain{1},1);
numOUT      = size(tTrain{1},1);


infoNN = defineNNsetup(numChannels, numOUT);

%% loop over the different setup
nSU = length(infoNN);
 for su = 1:nSU
     % define the architechture
     layers = infoNN(su).Setup;
     disp(layers)
     % define the options
     options = trainingOptions(      "adam", ...
         MaxEpochs                 = 10000, ...
         SequencePaddingDirection  = "right", ...
         Shuffle                   = "every-epoch", ...
         Plots                     = "none",  ...
...%         Plots="training-progress",  ...
        MiniBatchSize             = 32,...
         Verbose                   = 0,...
         ValidationData            = {XValid, TValid}, ...
         LearnRateSchedule         = 'piecewise', ...
         ValidationPatience        = 20,...
         LearnRateDropPeriod       = 10000000);
     % train
     net = trainNetwork(XTrain,TTrain,layers,options);
     % evaluate the NN performance
     PerfNetTMP  = evaluateNET(net, XTest, TTest, muT, sigmaT, 1:size(TTest{1},1), 1:length(XTest), METHOD, IDtvt{3}, ['D:\Lavori CARE HD\DATAxScripts\NNvel2pwr\nets\',CARTSAVE],['Sim_',num2str(su)]);
     save(['D:\Lavori CARE HD\DATAxScripts\NNvel2pwr\nets\',CARTSAVE,'\Sim_',num2str(su),'.mat'],"PerfNetTMP","net");
     PerfNet(su) = PerfNetTMP;
 end
 save(['D:\Lavori CARE HD\DATAxScripts\NNvel2pwr\nets\',CARTSAVE,'\PerfAll.mat'],...
     "PerfNet","XTrain","TTrain","XValid","TValid","XTest","TTest","IDtvt",...
     "sigmaT","muT","METHOD","infoNN");









%% helper functions
%---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------






function [data, targ, IDsubj] = getDataCycle_article(fileName, SI)

% load the data
fileName  = dir(fileName);
tmpDATA   = load([fileName.folder,'\',fileName.name]);

% get the events from all data
EV        = load([fileName.folder,'\AllData.mat'], 'EVt');
EV        = EV.EVt;
IMU       = load([fileName.folder,'\AllData.mat'], 'IMU');
IMU       = IMU.IMU{1};
IMUbodies = {IMU.Body};
Iboat     = strcmp(IMUbodies, 'boat');
Ipelv     = strcmp(IMUbodies, 'pelvis');
Ic7       = strcmp(IMUbodies, 'trunkCerv');
Imai      = strcmp(IMUbodies, 'trunkMAI');

% get the gyro and acceleration
% get the main eigenvectors

% boat
gyro   = IMU(Iboat).Gyro.Data;
gyro   = gyro - mean(gyro);
[v,~]  = eig(cov(gyro));
g_boat = (v' * gyro')';
acc    = IMU(Iboat).Accl.Data;
acc    = acc - mean(acc);
[v,~]  = eig(cov(acc));
a_boat = (v' * acc')';
% pelvis
gyro   = IMU(Ipelv).Gyro.Data;
gyro   = gyro - mean(gyro);
[v,~]  = eig(cov(gyro));
g_pelv = (v' * gyro')';
acc    = IMU(Ipelv).Accl.Data;
acc    = acc - mean(acc);
[v,~]  = eig(cov(acc));
a_pelv = (v' * acc')';
% c7
gyro   = IMU(Ic7).Gyro.Data;
gyro   = gyro - mean(gyro);
[v,~]  = eig(cov(gyro));
g_c7   = (v' * gyro')';
acc    = IMU(Ic7).Accl.Data;
acc    = acc - mean(acc);
[v,~]  = eig(cov(acc));
a_c7   = (v' * acc')';
% mai
gyro   = IMU(Imai).Gyro.Data;
gyro   = gyro - mean(gyro);
[v,~]  = eig(cov(gyro));
g_mai  = (v' * gyro')';
acc    = IMU(Imai).Accl.Data;
acc    = acc - mean(acc);
[v,~]  = eig(cov(acc));
a_mai  = (v' * acc')';
% get the principal component
g_pelv = g_pelv(:,3);
a_pelv = a_pelv(:,3);
g_c7   = g_c7(:,3);
a_c7   = a_c7(:,3);
g_mai  = g_mai(:,3);
a_mai  = a_mai(:,3);
% compose the vector to pass to ne script
g_imu  = [g_c7 g_mai g_pelv];
a_imu  = [a_c7 a_mai a_pelv];


% - gate angles
angGate  = squeeze(mean(tmpDATA.ProcData.PCHang.Data(:,1:2,:),2));
% - gate velocities
%   - do later
% - stretcher and gate forces
%   - along X
forLegs  = -squeeze(tmpDATA.ProcData.PCHfor.Data(:, 1,:)) * 9.81;
forGate  = squeeze(sum(tmpDATA.ProcData.PCHfor.Data(:,[2 4],:),2)) * 9.81;
% - speed and acceleration of the boat
velBoat  = squeeze(tmpDATA.ProcData.PCHtel.Data(:, 1,:));
accBoat  = squeeze(tmpDATA.ProcData.PCHtel.Data(:, 3,:));
% - power
pwrPin   = squeeze(tmpDATA.ProcData.PWRc.Data(:, 2,:));
pwrArms  = squeeze(tmpDATA.ProcData.PWRc.Data(:, 4,:));
pwrTrunk = squeeze(tmpDATA.ProcData.PWRc.Data(:, 5,:));
pwrLegs  = squeeze(tmpDATA.ProcData.PWRc.Data(:, 6,:));
% - velocities wrt boat
%   - C7, MAI, pelvis
velC7    = squeeze(tmpDATA.ProcData.IMUcD.Data(:, 4,:));
velMAI   = squeeze(tmpDATA.ProcData.IMUcD.Data(:, 5,:));
velPelv  = squeeze(tmpDATA.ProcData.IMUcD.Data(:, 6,:));
% - subject info
UUU      = ones(size(angGate(:,1)));
subW     = UUU*SI(1);
subH     = UUU*SI(2);

% concatenate the cycles into data and target
nC     = size(angGate,2);
data   = {};
targ   = {};
IDsubj = [];

% interpolates the data in the new cycles FCF - 100 points
TimeC = [];
DataC = [];
TargC = [];

for c=1:nC
    % get the data normalized to the CFC cycles
    % concatenate data and time vectors
    Ttmp2 = linspaceEXACT(EV.time(c,1), EV.time(c,3), size(angGate,1))';

    % time vector
    TimeC = [TimeC; Ttmp2(1:end-1,:)];

    % target and data
    DataC = [DataC; angGate(1:end-1,c) velBoat(1:end-1,c) accBoat(1:end-1,c) velC7(1:end-1,c) velMAI(1:end-1,c) velPelv(1:end-1,c) subW(1:end-1,1) subH(1:end-1,1)];

    TargC = [TargC; forLegs(1:end-1,c) forGate(1:end-1,c) pwrPin(1:end-1,c) pwrArms(1:end-1,c) pwrTrunk(1:end-1,c) pwrLegs(1:end-1,c)];

end
% add the last frame
TimeC = [TimeC; Ttmp2(end,:)];

% target and data
DataC = [DataC; angGate(end,c) velBoat(end,c) accBoat(end,c) velC7(end,c)   velMAI(end,c)   velPelv(end,c) subW(end,1) subH(end,1)];
TargC = [TargC; forLegs(end,c) forGate(end,c) pwrPin(end,c)  pwrArms(end,c) pwrTrunk(end,c) pwrLegs(end,c)];

% interpolate the data with a constant timeframe (100Hz)
FRsamp = 100;
TimeI  = TimeC(1) : 1/FRsamp : TimeC(end) + 1/FRsamp;
DataI  = interp1(TimeC, DataC, TimeI);
DataI  = fillmissing(DataI,"previous");
TargI  = interp1(TimeC, TargC, TimeI);
TargI  = fillmissing(TargI,"previous");
imuI   = (interp1(IMU(end).Time, [a_imu g_imu], TimeI));

% compute the gate velocity based on the position (1st row) and replace it in the data
frI        = 1/mean(diff(TimeI));
tmpV       = DataI(:,1);
tmpV       = [fliplr(tmpV); tmpV; fliplr(tmpV)];
% lowpass filter
tmpV       = lowpass(tmpV, 10, 100);
tmpV       = gradient(tmpV)*frI;
tmpV       = tmpV(length(TimeI)+1:length(TimeI)*2);
DataI(:,1) = tmpV;

figure,hold all
plot(DataI(:,1))
plot(tmpV)



for c=2:nC
    % define the events
    C0 = EV.time(c-1,1);
    F0 = EV.time(c-1,2);
    C1 = EV.time(c,1);
    F1 = EV.time(c,2);
    % get the info from the previous cycle
    % average boat speed of previous cycle
    aveBV0   = mean(velBoat(1:end,c-1));
    % max boat speed of previous cycle
    maxBV0   = max(velBoat(1:end,c-1));
    % dt between max arm speed and catch
    % duration of previous drive
    tmpC0dt  = F0 - C0;
    % duration of previous recovery
    tmpF0dt  = C1 -F0;
    % duration of current drive
    tmpC1dt  = F1 - C1;
    % duration of previous cycle
    tmpR0dt  = C1 - C0;
    % index of selected frames for the cycle
    Ifr      = (TimeI >= F0) & (TimeI <= F1);

    UUU      = ones(1,sum(Ifr));
    if any(any(isnan(DataI(Ifr, :)))) || any(any(isnan(TargI(Ifr, :))))
    else
%         data{end+1}   = [DataI(Ifr, :)'; UUU*aveBV0; UUU*maxBV0; UUU*tmpC0dt; UUU*tmpF0dt; UUU*tmpC1dt; UUU*tmpR0dt;  imuI(Ifr, :)'];
        data{end+1}   = [DataI(Ifr, :)'; UUU*aveBV0; UUU*maxBV0; UUU*tmpC0dt; UUU*tmpC1dt; imuI(Ifr, :)'];
        targ{end+1}   = [TargI(Ifr, :)]';
        IDsubj(end+1) = SI(3);
    end
end













