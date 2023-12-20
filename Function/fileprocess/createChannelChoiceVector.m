function [channelChoice] = createChannelChoiceVector(downchirpSync)
    rng(downchirpSync);
    channelChoice = fix(100*rand(1,1e4));