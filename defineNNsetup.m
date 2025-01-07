function infoNN = defineNNsetup(numChannels, numOUT)

% define the interval of the several parameters
NlayerFC0  = [0 1];
NlayerLSTM = [1 ]% 3]; 
NlayerFC1  = [0];
NneurFC    = floor(exp(linspace(log(50), log(1000), 8)));
NneurLSTM  = round(exp(linspace(log(10), log(500), 8)));
DOperc     = 0:0.05:0.1;


III = 1;

for nLL = 1:length(NlayerLSTM)
    for nNL = 1:length(NneurLSTM)
        for nLF = 1:length(NlayerFC1)
            for nNF = 1:length(NneurFC)
                for nDO = 1:length(DOperc)


                    LSTMlayers = [];
                    for l = 1:NlayerLSTM(nLL)
                        LSTMlayers = [LSTMlayers
                                      lstmLayer(NneurLSTM(nNL),'OutputMode','sequence')];
                    end

                    FClayers = [];
                    for l = 1:NlayerFC1(nLF)
                        FClayers = [FClayers
                                    fullyConnectedLayer(NneurFC(nNF))];
                    end

                    layers = [sequenceInputLayer(numChannels)
                              fullyConnectedLayer(numChannels)
                              LSTMlayers
                              FClayers
                              dropoutLayer(DOperc(nDO))
                              fullyConnectedLayer(numOUT)
                              regressionLayer];

                    infoNN(III).Setup   = layers;


                    infoNN(III).Array = [NlayerLSTM(nLL) NneurLSTM(nNL) NlayerFC1(nLF) NneurFC(nNF) DOperc(nDO)];






                    III = III + 1;
                end
            end
        end
    end
end













