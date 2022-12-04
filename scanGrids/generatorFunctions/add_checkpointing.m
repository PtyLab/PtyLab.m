function [R, C] = add_checkpointing(R, C, n)
% [R, C] = rectify_fermat(R, C);
% keeps only central square within Fermat grid, discards other points
% R: row
% C: column
% n: number of checkpoints
% note: this function assumes the first scan point to serve as checkpoint

num_scan_points = length(R);
checkpoint_spacing = round(num_scan_points / (n+1));
checkpoint_indices = (1:n)*checkpoint_spacing;

Rtemp = [ ];
Ctemp = [ ];
for k = 1:n
    if k == 1
        Rtemp = [Rtemp; R(1:checkpoint_indices(1)); R(1)];
        Ctemp = [Ctemp; C(1:checkpoint_indices(1)); R(1)];
    elseif k==n
        Rtemp = [Rtemp; R(checkpoint_indices(k-1)+1:end)];
        Ctemp = [Ctemp; C(checkpoint_indices(k-1)+1:end)];
    else
        Rtemp = [Rtemp; R(checkpoint_indices(k-1)+1:checkpoint_indices(k)); R(1)];
        Ctemp = [Ctemp; C(checkpoint_indices(k-1)+1:checkpoint_indices(k)); R(1)];
    end
end
R = Rtemp;
C = Ctemp;
end