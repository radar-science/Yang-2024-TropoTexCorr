function dataInp = demInterp2(data)

[Na,Nr] = size(data);
dataLine = data(:);

index = find(~isnan(dataLine));
x = ceil(index/Na);
y = mod(index,Na);

dataLine(isnan(dataLine)) = [];

xi = 1:Nr;
yi = 1:Na;
dataInp = griddata(x,y,dataLine,xi,yi','cubic');
end
