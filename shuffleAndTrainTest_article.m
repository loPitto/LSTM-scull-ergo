function [dTrain, tTrain, dValid, tValid, dTest, tTest, IDtvt] = shuffleAndTrainTest_article(varargin)

%% parse input
data    = varargin{1}; 
target  = varargin{2}; 
percTVT = varargin{3}; 
ERGBOAT = varargin{4}; 
IDsubj  = varargin{5}; 

if length(varargin)>5
    LOOS = varargin{6}; 
else
    LOOS = [];
end

%%
if sum(percTVT)~=100
    error('perc > 100%')
end


switch ERGBOAT
    case 'boat'

        if 1
            % remove bad data (these are cycles with artifacts)
            %         III = [84 85 86 87 347];
            III = 341;
            data(III) = [];
            target(III) = [];

            IDsubj(III) = [];


            % mean of the first feature, to remove bad cycles

            for q=1:length(data)
                tmp = data{q}(1,:);
                MR1(q) = mean(tmp);
            end
            figure,plot(MR1)


            % flips the accelerations and rotations (makes sures that all acceleration and rotations have consistent directions -
            % >>>>>>>>>>> THIS IS NOT AUTOMATIC - IT WAS DONE AFTER VISUAL INSPECTIONS, like the above/below figures)<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            for IIF  = (16:19)-2
                tmp1 = data{1}(IIF,:);
                t1   = (0:length(tmp1)-1)/(length(tmp1)-1);
                for q=2:length(data)
                    tmp2 = data{q}(IIF,:);
                    t2   = (0:length(tmp2)-1)/(length(tmp2)-1);
                    tmp2 = interp1(t2,tmp2,t1);
                    cc   = corrcoef(tmp1,tmp2);
                    if cc(2)<0
                        data{q}(IIF,:) = -data{q}(IIF,:);
                    end
                end
            end

        end

        %
        if 1
            F1=figure;
            %     F2=figure;
            figure(F1)
            for q=1:18
                subplot(4,5,q),hold all
                for w=1:length(data)
                    plot(data{w}(q,:))
                end
            end
            drawnow
        end


    case 'ergo'
        if 1
            III = [29 72 115 187 489 803 1566];
            data(III) = [];
            target(III) = [];
            IDsubj(III) = [];
            for q=1:length(data)
                tmp = data{q}(4,:);
                MR1(q) = max(tmp);
            end
            figure,plot(MR1)


            for q=1:length(data)
                tmp = data{q}(8,:);
%                 if length(tmp)>101
%                     MR2(q)=tmp(102);
                if length(tmp)>397
                    MR2(q)=tmp(398);
                else
                    MR2(q) = 0;
                end
            end
            figure,plot(MR2)



        end
        %
        if 1
            F1=figure;
            %     F2=figure;
            figure(F1)
            for q=1:14
                subplot(4,5,q)
                hold all
                for w=1:length(data)
                    plot(data{w}(q,:))
                end
            end
            drawnow
        end



end


%% divide into train, valid and test

if isempty (LOOS)

    Ish       = randperm(length(data));
    % Ish       = 1:(length(data));
    data      = data(Ish);
    target    = target(Ish);
    IDsubj    = IDsubj(Ish);

    percTVT   = percTVT/100;

    numObs    = numel(data);
    idxTrain  = 1               : floor(percTVT(1)*numObs);
    idxValid  = idxTrain(end)+1 : floor((percTVT(1)+percTVT(2))*numObs);
    idxTest   = idxValid(end)+1 : numObs;
    if isempty(idxTest)
        idxTest = idxValid;
    end

else

    idxTrain = find(IDsubj ~= LOOS);
    idxValid = find(IDsubj == LOOS);
    idxTest  = find(IDsubj == LOOS);



end





dTrain    = data(idxTrain);
dValid    = data(idxValid);
dTest     = data(idxTest);

tTrain    = target(idxTrain);
tValid    = target(idxValid);
tTest     = target(idxTest);


IDtvt{1} = IDsubj(idxTrain);
IDtvt{2} = IDsubj(idxValid);
IDtvt{3} = IDsubj(idxTest);
































