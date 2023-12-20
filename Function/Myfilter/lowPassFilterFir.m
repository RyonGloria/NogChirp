function [singalOut] = lowPassFilterFir(signal, loraSet)
%     b = fir1(30, 0.15, "low");
    b = fir1(30, loraSet.pass_arg, "low");
    singalOut = filter(b, 1, signal);