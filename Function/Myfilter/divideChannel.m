function [singalOut] = divideChannel(loraSet, singnal_in, preambleChannelChoice, iseparated)
    if iseparated == true
        [singal_out_3175] = filter_3175KHz(loraSet, singnal_in);
        [singal_out_3375] = filter_3375KHz(loraSet, singnal_in);
        [singal_out_3575] = filter_3575KHz(loraSet, singnal_in);
        [singal_out_3775] = filter_3775KHz(loraSet, singnal_in);
        singalOut = [singal_out_3175; singal_out_3375; singal_out_3575; singal_out_3775];
    else
        [singalOut_3175] = signalFrequencyShift(loraSet, singnal_in, -preambleChannelChoice(1));
        [singalOut_3175] = lowPassFilterFir(singalOut_3175, loraSet);

        [singalOut_3375] = signalFrequencyShift(loraSet, singnal_in, -preambleChannelChoice(2));
        [singalOut_3375] = lowPassFilterFir(singalOut_3375, loraSet);
        
        if(length(preambleChannelChoice) == 2)
            singalOut_3575 = zeros(1, length(singalOut_3375));
            singalOut_3775 = zeros(1, length(singalOut_3375));
        else
            [singalOut_3575] = signalFrequencyShift(loraSet, singnal_in, -preambleChannelChoice(3));
            [singalOut_3575] = lowPassFilterFir(singalOut_3575, loraSet);
            
            [singalOut_3775] = signalFrequencyShift(loraSet, singnal_in, -preambleChannelChoice(4));
            [singalOut_3775] = lowPassFilterFir(singalOut_3775, loraSet);
        end
        singalOut = [singalOut_3175; singalOut_3375; singalOut_3575; singalOut_3775];
        
    end