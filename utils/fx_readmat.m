function rawdata = fx_readmat(S)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fields = fieldnames(S);

xl = length(fields);
tl = length(S.(fields{1}));
yl = length(S.(fields{1})(:, 1));
rawdata = zeros(tl, yl, xl);
for i = 1:numel(fields)
    rawdata(:, :,  i) = S.(fields{i})';
end
end