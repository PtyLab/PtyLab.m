function [R, C] = rectify_fermat(R,C)
% [R, C] = rectify_fermat(R, C);
% keeps only central square within Fermat grid, discards other points
% R: row 
% C: column
% note: the function assumes the scan grid is centered around origin

disp('--')
disp('rectify travel path')
num_points_before_rectification = length(R);
disp(['number of points before rectification: ', num2str(num_points_before_rectification)])
% get radius of circle containing original scan grid
radius = max( [abs(R); abs(C)] );
% get side length of largest square inside original circle
side_length = 2*radius / sqrt(2);
% discard points outside that square
idx  = or(abs(R) > side_length / 2, ...
          abs(C) > side_length / 2);
R(idx) = [];
C(idx) = [];
num_points_after_rectification = length(R);
disp(['number of points after rectification: ', num2str(num_points_after_rectification)])

end