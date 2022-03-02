function cf = fx_showCscans(inds, inam_interface, fignum)

% show the C scans for comparision
cf = figure(fignum);
row = ceil(sqrt(length(inds)));
col = ceil(length(inds) / row);

for i = 1:length(inds)
    
    ca = subplot(row, col, i);
    imagesc(ca, inam_interface(:, :, inds(i)));
    title(ca, [num2str(inds(i)) , 'th layer']);
    colorbar;
end

