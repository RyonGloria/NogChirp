randSqeunce = zeros(1024, 1e4);
for i = 1:1024
    rng(i);
    randSqeunce(i,:) = fix(100*rand(1,1e4));
end
writematrix(randSqeunce,'E:\tmp.txt');
fclose all;