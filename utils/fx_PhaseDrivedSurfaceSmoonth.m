function rear_I = fx_PhaseDrivedSurfaceSmoonth(rear_I, Inph, phase_rear, search_R)
% % increase the robustness of the back-surface location in the data by
% using phase data.

% rear_I:  unsmoothed rear surface location, 2D matrix
% Inph:    instantaneous phase, 3D data
% phase_rear_mean: expected phase of the rear surface
% error:   critical of phase errors.
 

[lx, ly, ~] = size(Inph);

% search towards the rear_I_mean
for i=1:lx
    for j=1:ly
        % find by interpolation
        % a error comes when there are duplicate points in input to 'interp1'.
        %remove duplicates using 'unique' function.
        inph_a = squeeze(Inph(i, j, :));
        [~, I_s] = min ( abs(inph_a(search_R + rear_I(i,j) ) - phase_rear) ); 
%         disp(mean(inph_a(rear_I(i,j))));
        rear_I(i, j) = rear_I(i, j) + I_s + search_R(1);
    end
end

end

