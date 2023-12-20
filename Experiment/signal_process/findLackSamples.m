Dir = 'D:\CHchirp_IPSN_samples\sf7\channel2\preamble4\subchirpNum16\';
file = dir(fullfile(Dir, '*.sigmf-data')); 
numArray = zeros(1, 128);
for i = 1:128
    numArray(i) = i;
end
for i = 1:length(file)
    tmpName = file(i).name;
    tmpName = tmpName(isstrprop(tmpName, 'digit'));
    tmpName = str2num(tmpName);
    numArray(tmpName) = 0;
end
numArray = numArray((numArray > 0));
fprintf(" %d,", numArray-1);