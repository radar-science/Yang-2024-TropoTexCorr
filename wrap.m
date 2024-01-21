function [Data] = wrap(Data,w)
tmp = exp(1i*Data*(pi/w));
Data = angle(tmp)/(pi/w);
end