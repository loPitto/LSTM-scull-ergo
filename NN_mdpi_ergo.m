function NN_mdpi_ergo()

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
CART     = 'D:\Lavori CARE HD\DATAxScripts\MisureErgo\TestFeedback\DATAcycles\';


ALLDATAF = {[CART,'S01.mat']
            [CART,'S02.mat']
            [CART,'S03.mat']
            [CART,'S04.mat']
            [CART,'S05.mat']
            [CART,'S06.mat']
            [CART,'S07.mat']
            [CART,'S08.mat']
            [CART,'S09.mat']
            [CART,'S10.mat']
            [CART,'S11.mat']};

nS = size(ALLDATAF,1);

% 1-height, 2-weight, 3-ID
SubInfo = [185 81 1
           182 70 2
           185 69 3
           202 90 4
           185 83 5
           178 70 6
           194 79 7
           178 77 8
           185 84 9
           180 70 10
           187 78 11];


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
 dTest,  tTest,  IDtvt] = shuffleAndTrainTest_article(allData, allTarg, [80 20 0], 'ergo', allSubj);

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


%% loop ove the different setup
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
     PerfNetTMP  = evaluateNET(net, XTest, TTest, muT, sigmaT, 1:size(TTest{1},1), 1:length(XTest), METHOD, IDtvt{3}, 'D:\Lavori CARE HD\DATAxScripts\NNvel2pwr\nets\RandSearchErgo',['Sim_',num2str(su)]);
     save(['D:\Lavori CARE HD\DATAxScripts\NNvel2pwr\nets\RandSearchErgo\Sim_',num2str(su),'.mat'],"PerfNetTMP","net");
     PerfNet(su) = PerfNetTMP;
 end

 save('D:\Lavori CARE HD\DATAxScripts\NNvel2pwr\nets\RandSearchErgo\PerfAll.mat',...
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

% list of fields, describing the conditions (cad20, cadMAX, feedback, ...)
FFF = fields(tmpDATA.d2s);

% loop on the feedback session (other data will be done differently)
III = 1;
for f=1:length(FFF)

    % If they are segmented into cycles, concatenate and interpolate (200 Hz ?)
    tmpDATAcfc = tmpDATA.d2s.(FFF{f}).CFC;
    tmpDATAcf  = tmpDATA.d2s.(FFF{f}).CF;
 
    [Pwr.Handle, Time] = concatenateInterp(tmpDATAcfc.Pwr.Handle, tmpDATAcfc.Pwr.Time);
    Pwr.Arms           = concatenateInterp(tmpDATAcfc.Pwr.Arms,   tmpDATAcfc.Pwr.Time);
    Pwr.Trunk          = concatenateInterp(tmpDATAcfc.Pwr.Trunk,  tmpDATAcfc.Pwr.Time);
    Pwr.Legs           = concatenateInterp(tmpDATAcfc.Pwr.Legs,   tmpDATAcfc.Pwr.Time);

    For.Handle         = concatenateInterp(tmpDATAcfc.For.Handle, tmpDATAcfc.Pwr.Time);
    For.Feet           = concatenateInterp(tmpDATAcfc.For.Feet,   tmpDATAcfc.Pwr.Time);

    Vel.Handle         = concatenateInterp(tmpDATAcfc.Vel.Handle, tmpDATAcfc.Pwr.Time);
    Vel.Arms           = concatenateInterp(tmpDATAcfc.Vel.Arms,   tmpDATAcfc.Pwr.Time);
    Vel.Trunk          = concatenateInterp(tmpDATAcfc.Vel.Trunk,  tmpDATAcfc.Pwr.Time);
    Vel.Legs           = concatenateInterp(tmpDATAcfc.Vel.Legs,   tmpDATAcfc.Pwr.Time);

    Pos.Handle         = concatenateInterp(tmpDATAcfc.Pos.Handle, tmpDATAcfc.Pwr.Time);
    Pos.Trunk          = concatenateInterp(tmpDATAcfc.Pos.Trunk,  tmpDATAcfc.Pwr.Time);
    Pos.Seat           = concatenateInterp(tmpDATAcfc.Pos.Seat,   tmpDATAcfc.Pwr.Time);

    [EVcfc, EVcf]      = getEVfromCycles(tmpDATAcfc.Pwr.Time, tmpDATAcf.Pwr.Time);


    if true
        close all
        figure
        subplot(4,1,1),hold all
        plot(Time,Pwr.Handle)
        plot(Time,Pwr.Arms)
        plot(Time,Pwr.Trunk)
        plot(Time,Pwr.Legs)
        plotEvents(EVcfc)
        plotEvents(EVcf)

        subplot(4,1,2),hold all
        plot(Time,For.Handle)
        plot(Time,For.Feet)
        plotEvents(EVcfc)
        plotEvents(EVcf)

        subplot(4,1,3),hold all
        plot(Time,Vel.Handle)
        plot(Time,Vel.Arms)
        plot(Time,Vel.Trunk)
        plot(Time,Vel.Legs)
        plotEvents(EVcfc)
        plotEvents(EVcf)

        subplot(4,1,4),hold all
        plot(Time,Pos.Handle)
        plot(Time,Pos.Trunk)
        plot(Time,Pos.Seat)
        plotEvents(EVcfc)
        plotEvents(EVcf)


    end

    % segment the data according t a new criterion
    % only do for the METHOD 'seconds'
    % -1- recovery -> drive
    for c = 2:size(EVcfc.time,1)-1
        % define the events
        C0 = EVcfc.time(c-1,1);
        F0 = EVcfc.time(c-1,2);
        C1 = EVcfc.time(c,1);
        F1 = EVcfc.time(c,2);
        % find the time between the finish of the previous cycle and the finish of the current cycle
        indT  = (Time>=F0) & (Time<=F1);

        % segment all the data for these frames
        tmpT  = Time(indT);

        tmpFh = For.Handle(indT);
        tmpFf = For.Feet(indT);

        tmpPh = Pwr.Handle(indT);
        tmpPa = Pwr.Arms(indT);
        tmpPt = Pwr.Trunk(indT);
        tmpPl = Pwr.Legs(indT);

        tmpVh = Vel.Handle(indT);
        tmpVa = Vel.Arms(indT);
        tmpVt = Vel.Trunk(indT);
        tmpVl = Vel.Legs(indT);

        tmpAh = gradient(tmpVh) / mean(diff(tmpT));
        tmpAa = gradient(tmpVa) / mean(diff(tmpT));
        tmpAt = gradient(tmpVt) / mean(diff(tmpT));
        tmpAl = gradient(tmpVl) / mean(diff(tmpT));

        tmpXh = Pos.Handle(indT);
        tmpXc = Pos.Trunk(indT);
        tmpXs = Pos.Seat(indT);
        
        tmpXa = tmpXh - tmpXc;
        tmpXt = tmpXc - tmpXs;

        % max handle speed of previous cycle
        tmpVhMax = max(tmpVh);
        % dt between max arm speed and catch
        tmpVhDt  = C1 - (tmpT(tmpVh==tmpVhMax));
        % duration of previous drive
        tmpC0dt  = F0 - C0;
        % duration of current drive
        tmpC1dt  = F1 - C1;

        UUU   = ones(size(tmpFh));
        % add the data to the input and target to be passed
%         data2nn{III}   = [tmpXh;            tmpXt;              tmpXs;...
%                           tmpVh;            tmpVa;              tmpVt;          tmpVl;...
%                           tmpAh;            tmpAt;              tmpAl;...
%                           UUU*SI(1);        UUU*SI(2);...
%                           UUU*tmpVhMax;     UUU*tmpVhDt;...
%                           UUU*tmpC0dt;      UUU*tmpC1dt];
        data2nn{III}   = [tmpXh;            tmpXc;          tmpXs;...
                          tmpVh;            tmpVa;          tmpVt;          tmpVl;...
                          tmpAh;            tmpAt;          tmpAl;...
                          UUU*SI(1);        UUU*SI(2);...
                          UUU*tmpVhMax;     UUU*tmpVhDt;...
                          UUU*tmpC0dt;      UUU*tmpC1dt];

        targ2nn{III}   = [tmpFf; tmpFh; tmpPh; tmpPa; tmpPt; tmpPl];
        
        % separate the MAX spm with a different ID
        if strcmp(FFF{f}, 'MAX')
            ID2nn(III) = SI(3) + 100;
        else
            ID2nn(III) = SI(3);
        end

        III = III +1;
    end













end

data   = data2nn;
targ   = targ2nn;
IDsubj = ID2nn;













%---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------

function [Dout, Tout] = concatenateInterp(Din, Tin)


FREQ = 150;% interpolate at 150Hz

nC = size(Din, 2);
tmpD = [];
tmpT = [];
for c=1:nC
    tmpD = [tmpD; Din(1:end-1,c)];
    tmpT = [tmpT; Tin(1:end-1,c)];
end
tmpD = [tmpD; Din(end,end)];
tmpT = [tmpT; Tin(end,end)];

Tout = tmpT(1):1/FREQ:tmpT(end);
Dout = interp1(tmpT, tmpD, Tout);

%---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------

function [EVcfc, EVcf] = getEVfromCycles(Tcfc, Tcf)

for q=1:size(Tcfc,2)
    EV(q,:) = [Tcfc(1,q) Tcf(end,q) Tcfc(end,q)];
end

EVcfc.time = EV;
EVcf.time =  EV(:,[2 2 3]);


%---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------

function [OUTd, OUTt] = interpolateRecoveryDrive(T,iTr,iTd,D)
if 1
    Tr   = T(iTr);
    Td   = T(iTd);
    Tri  = linspaceEXACT(Tr(1), Tr(end), 50);
    Tdi  = linspaceEXACT(Td(1), Td(end), 50);
    tmpR = interp1(Tr, D(iTr), Tri);
    tmpD = interp1(Td, D(iTd), Tdi);
    OUTd = [tmpR tmpD];
    OUTt = [Tri Tdi];
else
    Td   = T(iTd);
    Tdi  = linspaceEXACT(Td(1), Td(end), 100);
    tmpD = interp1(Td, D(iTd), Tdi);
    OUTd = tmpD;
    OUTt = Tdi;
end








