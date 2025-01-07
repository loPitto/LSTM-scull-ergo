function OUT = evaluateNET(net, dTest, tTest, muT, sigmaT, IndVAR, IndCYC, METHOD, IDtvt, CART, SIM)
%
If   = 0;
% F1 = figure('units','normalized','outerposition',[0 0 1 1]);
% F2 = figure('units','normalized','outerposition',[1 0 1 1]);

nFig = (IndCYC(end) - IndCYC(1))+1;
nRow = ceil(sqrt(nFig));
nCol = ceil(nFig/nRow);

numOUT      = size(tTest{1},1);

% initialize the differences between predictions and target
DiffYT = [];
Yall   = [];
Tall   = [];

COLOR = [0 0 0
    0 .7 0
    0 0 .7
    .7 0 0
    .7 .7 0
    0 .7 .7
    .7 0 .7];

for idx = IndCYC
    If = If + 1;
    %     for n = 1:length(dTest)
    %         X        = (dTest{n} - muX) ./ sigmaX;
    %         T        = (tTest{n} - muT) ./ sigmaT;
    %         XTest2{n} = X;%(:,1:end-1) ;
    %         TTest2{n} = T;%(:,2:end);
    %     end
    X = dTest{idx};
    T = tTest{idx};



    switch METHOD
        case {'interp','seconds'}
            net     = resetState(net);
            offset  = 1;
            [net,~] = predictAndUpdateState(net,X(:,1:offset));

            numTimeSteps           = size(X,2);
            numPredictionTimeSteps = numTimeSteps - offset;
            Y                      = zeros(numOUT,numPredictionTimeSteps);

            for t = 1:numPredictionTimeSteps
                Xt           = X(:,offset+t);
                [net,Y(:,t)] = predictAndUpdateState(net,Xt);
            end
            T = T .* sigmaT + muT;
            Y = Y .* sigmaT + muT;

            Y = [zeros(length(IndVAR),1) Y(IndVAR,:)];
            %     X = X(IndVAR,:);
            T = T(IndVAR,:);
        case 'nonRecurrent'

            Y = net.predict(X);
            T = T .* sigmaT + muT;
            Y = Y .* sigmaT + muT;
            Y = Y(IndVAR,:);
            %     X = X(IndVAR,:);
            T = T(IndVAR,:);
    end


    t0 = 1:size(T,2);
    t1 = linspaceEXACT(t0(1), t0(end), 101);
    Ti = interp1(t0,T',t1)';
    Yi = interp1(t0,Y',t1)';
    Xi = interp1(t0,X(1,:),t1);

    % append the differences
    DiffYT(:,:,end+1) = Ti-Yi;
    Yall(:,:,end+1) = Yi;
    Tall(:,:,end+1) = Ti;

    if 0
        % plot the comparison of the curves cycle by cycle
        figure(F1)
        subplot(nRow,nCol,If)
        hold all
        for v=1:length(IndVAR)
            plot(Xi(1,:),Ti(v,:),'color',COLOR(v,:)+.3,'LineWidth',1,'LineStyle','-')
            plot(Xi(1,:),Yi(v,:),'color',COLOR(v,:),'LineWidth',1,'LineStyle','-.')
            ylim([-500 2000])
            drawnow
        end
        figure(F2)
        subplot(nRow,nCol,If)
        hold all
        for v=1:length(IndVAR)
            plot(Ti(v,:),'color',COLOR(v,:)+.3,'LineWidth',1,'LineStyle','-')
            plot(Yi(v,:),'color',COLOR(v,:),'LineWidth',1,'LineStyle','-.')
            drawnow
            ylim([-500 2000])
        end
    end

end

DiffYT(:,:,1) = [];%<<<<<<<<<<<<<<<<<<<<<<<< remove first cycle (started at end+1) <<<<<<<<<<<<<<<<<
Yall(:,:,1) = [];
Tall(:,:,1) = [];

DiffYT(:,1,:) = [];%<<<<<<<<<<<<<<<<<<<<<<<< remove first sample (prediction starts at 0) <<<<<<<<<<<<<<<<<
Yall(:,1,:) = [];
Tall(:,1,:) = [];



OUT.all = computePerf(DiffYT, Yall, Tall, [CART,'\',SIM,'_AllSubj.csv']);

%% separate into different subjects
Subj = unique(IDtvt);

for s=1:length(Subj)

disp(['--- Subject ', num2str(Subj(s)),' ---'])


    IndS = find(IDtvt == Subj(s));
    OUT.("Subj_"+Subj(s)) = computePerf(DiffYT(:,:,IndS), Yall(:,:,IndS), Tall(:,:,IndS), [CART,'\',SIM,'_Subj_',num2str(Subj(s)),'.csv']);
    OUT.("Subj_"+Subj(s)).Nsubj = length(IndS);

    % compute mean and std for each target
    tmpMN_Y = mean(Yall(:,:,IndS), 3)';
    tmpMN_T = mean(Tall(:,:,IndS), 3)';

    tmpSD_Y = std(Yall(:,:,IndS), [], 3)';
    tmpSD_T = std(Tall(:,:,IndS), [], 3)';

    figure ('Position',[ 1922         476        1917         521])
    subplot(1,2,1),hold all,grid on
    plotMeanStd(1:size(tmpMN_Y,1), tmpMN_Y(:,1), tmpSD_Y(:,1), [1 0 0])
    plotMeanStd(1:size(tmpMN_T,1), tmpMN_T(:,1), tmpSD_T(:,1), [.5 0 0])
    plotMeanStd(1:size(tmpMN_Y,1), tmpMN_Y(:,2), tmpSD_Y(:,2), [0 0 1])
    plotMeanStd(1:size(tmpMN_T,1), tmpMN_T(:,2), tmpSD_T(:,2), [0 0 .5])
    ylim([-300 1300])
    ylabel('N')
    line([1 size(tmpMN_Y,1)], [0 0], 'color','k','linestyle','--')
    legend({'','Foot pred','','Foot target', '','Gate pred','','Gate target'},'Location','northwest')
    title("N cycles: " + length(IndS))

    subplot(1,2,2),hold all,grid on
    plotMeanStd(1:size(tmpMN_Y,1), tmpMN_Y(:,3), tmpSD_Y(:,3), [0 0 0])
    plotMeanStd(1:size(tmpMN_T,1), tmpMN_T(:,3), tmpSD_T(:,3), [.5 .5 .5])
    plotMeanStd(1:size(tmpMN_Y,1), tmpMN_Y(:,4), tmpSD_Y(:,4), [0 1 0])
    plotMeanStd(1:size(tmpMN_T,1), tmpMN_T(:,4), tmpSD_T(:,4), [0 .5 0])
    plotMeanStd(1:size(tmpMN_Y,1), tmpMN_Y(:,5), tmpSD_Y(:,5), [0 0 1])
    plotMeanStd(1:size(tmpMN_T,1), tmpMN_T(:,5), tmpSD_T(:,5), [0 0 .5])
    plotMeanStd(1:size(tmpMN_Y,1), tmpMN_Y(:,6), tmpSD_Y(:,6), [1 0 0])
    plotMeanStd(1:size(tmpMN_T,1), tmpMN_T(:,6), tmpSD_T(:,6), [.5 0 0])
    ylim([-300 2000])
    ylabel('W')
    line([1 size(tmpMN_Y,1)], [0 0], 'color','k','linestyle','--')
    legend({'','Gate pred','','Gate target', '','Arms pred','','Arms target', '','Trunk pred','','Trunk target', '','Legs pred','','Legs target'},'Location','northwest')


    saveas(gcf,[CART,'\',SIM,'_Subj_',num2str(Subj(s)),'.png'],'png')
    close(gcf)
end




%% helper functions
% **********************************************************************************************************************
% **********************************************************************************************************************
% **********************************************************************************************************************

%% compute the performance measures
function OUT = computePerf(DYT, Y, T, FN)


%   - MAE
MAE = mean(abs(DYT),3);
%   - RMSE
%   - ME
ME  = mean(DYT,3);

nC = size(DYT,3);
nV = size(DYT,1);

nFig = nV;
nRow = ceil(sqrt(nFig));
nCol = ceil(nFig/nRow);
% F3 = figure('units','normalized','outerposition',[1 0 1 1]);
% 
% for v=1:nV
%     for c=1:nC
%         subplot(nRow,nCol,v),hold all,grid on
%         plot(DiffYT(v,:,c),'r')
%         plot(MAE(v,:),'Color',[0 .7 0],'LineWidth',4)
%         plot(ME(v,:),'Color',[0 .7 0],'LineWidth',4,'LineStyle',':')
%         plot(Yall(v,:,c),'b')
%         plot(Tall(v,:,c),'k')
%         xlim([1 size(MAE,2)])
%         ylim([-500 2000])
%     end
% end


CC = zeros(nV,nC);

for v=1:nV
    for c=1:nC
        tmp = corrcoef(Y(v,:,c), T(v,:,c));
        CC(v,c) = tmp(2);
    end
end


% find the indexes for the effective drive phase (20-10 kg)
for q=1:size(T,3)
    FH = T(2,:,q)';
    I1 = find(FH>=20*9.81,1,"first");
    I2 = find(FH>=10*9.81,1,"last");

    Ttmp = T(:,:,q);
    Ytmp = Y(:,:,q);

    Tdrive(:,:,q) = (interp1(1:length(FH),Ttmp',linspaceEXACT(I1,I2,51)))';
    Ydrive(:,:,q) = (interp1(1:length(FH),Ytmp',linspaceEXACT(I1,I2,51)))';

    INTdrive(q,:) = [I1 I2];

end

OUT.T   = T;
OUT.Y   = Y;
OUT.DR  = INTdrive;
OUT.Tdr = Tdrive;
OUT.Ydr = Ydrive;





Ddrive = Tdrive-Ydrive;

dispMAE = mean(mean(abs(Ddrive),3),2);
dispME  = mean(abs(mean(Ddrive,3)),2);



% disp('------------------------------------------------------------------------------------------------------------')
% disp('--- NON NORMALIZED ---')
% disp('------------------------------------------------------------------------------------------------------------')
% disp('--- Drive ---')
% disp('MAE: ')
% disp(dispMAE)
% 
% disp('ME: ')
% disp(dispME)


OUT.MAE_drive_nn = dispMAE;
OUT.ME_drive_nn  = dispME;


III = 1:100;

dispMAE = mean(MAE(:,III),2);
dispME  = mean(abs(ME(:,III)),2);

% disp('--- Cycle ---')
% disp('MAE: ')
% disp(dispMAE)
% 
% disp('ME: ')
% disp(dispME)

OUT.MAE_cycle_nn = dispMAE;
OUT.ME_cycle_nn  = dispME;





% III = 50:100; 
% dispMAE = mean(MAE(:,III),2)./mean(mean(abs(T(:,III,:)),3),2)*100;
% dispME = mean(abs(ME(:,III)),2)./mean(mean(abs(T(:,III,:)),3),2)*100;
dispMAE = mean(mean(abs(Ddrive),3),2)./mean(mean(abs(Tdrive),3),2)*100;
dispME  = mean(abs(mean(Ddrive,3)),2)./mean(mean(abs(Tdrive),3),2)*100;





disp('------------------------------------------------------------------------------------------------------------')
disp('--- NORMALIZED ---')
disp('------------------------------------------------------------------------------------------------------------')
disp('--- Drive ---')
disp('MAE: ')
disp(dispMAE)

disp('ME: ')
disp(dispME)

OUT.MAE_drive_norm = dispMAE;
OUT.ME_drive_norm  = dispME;



III = 1:100;

dispMAE = mean(MAE(:,III),2)./mean(mean(abs(T(:,III,:)),3),2)*100;
dispME  = mean(abs(ME(:,III)),2)./mean(mean(abs(T(:,III,:)),3),2)*100;

disp('--- Cycle ---')
disp('MAE: ')
disp(dispMAE)

disp('ME: ')
disp(dispME)

OUT.MAE_cycle_norm = dispMAE;
OUT.ME_cycle_norm  = dispME;


disp('------------------------------------------------------------------------------------------------------------')
disp('--- CORR COEFF ---')
disp('------------------------------------------------------------------------------------------------------------')

disp('CC: ')
disp(nanmean(CC,2))


% print in a csv file
fid = fopen(FN,'w');

FFF= fieldnames(OUT);
for f = 1:length(FFF)
    fprintf(fid,FFF{f});
    fprintf(fid,';');
    d2p = OUT.(FFF{f});
    for q=1:length(d2p)
        fprintf(fid,'%f',d2p(q));
        fprintf(fid,';');
    end
    fprintf(fid,'\n');
end


fclose(fid);




















