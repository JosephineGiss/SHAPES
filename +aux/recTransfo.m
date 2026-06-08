function res = recTransfo(tab,varargin)
    % rectangular transformation for tensors
    Dim_mat = size(tab);
    % order = [1,4,5,7,2,6,8,9,3];

    % rectangular matrix (in) into tensor (out)
    if length(Dim_mat) == 2
        if all(Dim_mat == [9,1])
            res = [tab(1), tab(4), tab(5); tab(7), tab(2), tab(6); tab(8), tab(9), tab(3)];
        elseif all(Dim_mat == [3,9])
            res = zeros(3,3,3);
            for i = 1:3
                res(i,:,:) = [tab(i,1), tab(i,4), tab(i,5); tab(i,7), tab(i,2), tab(i,6); tab(i,8), tab(i,9), tab(i,3)];
            end
        elseif all(Dim_mat == [9,3])
            res = zeros(3,3,3);
            for k = 1:3
                res(:,:,k) = [tab(1,k), tab(4,k), tab(5,k); tab(7,k), tab(2,k), tab(6,k); tab(8,k), tab(9,k), tab(3,k)];
            end
        elseif all(Dim_mat == [9,9])
            res = zeros(3,3,3,3);
            kmat = [1, 4, 5; 7, 2, 6; 8, 9, 3];
            for i = 1:3
                for j = 1:3
                    k = kmat(i,j);
                    res(:,:,i,j) = [tab(1,k), tab(4,k), tab(5,k); tab(7,k), tab(2,k), tab(6,k); tab(8,k), tab(9,k), tab(3,k)];
                end
            end
    % tensor (in) into rectangular matrix (out)
        elseif all(Dim_mat == [3,3])
            res = [tab(1,1); tab(2,2); tab(3,3); tab(1,2); tab(1,3); tab(2,3); tab(2,1); tab(3,1); tab(3,2)];
        end
    elseif length(Dim_mat) == 3
        if varargin{1} == "39"
            res = zeros(3,9);
            for i = 1:3
                res(i,:) = [tab(i,1,1), tab(i,2,2), tab(i,3,3), tab(i,1,2), tab(i,1,3), tab(i,2,3), tab(i,2,1), tab(i,3,1) tab(i,3,2)];
            end
        elseif varargin{1} == "93"
            res = zeros(9,3);
            for k = 1:3
                res(:,k) = [tab(1,1,k), tab(2,2,k), tab(3,3,k), tab(1,2,k), tab(1,3,k), tab(2,3,k), tab(2,1,k), tab(3,1,k) tab(3,2,k)];
            end
        end
    elseif length(Dim_mat) == 4
        res = [tab(1,1,1,1), tab(1,1,2,2), tab(1,1,3,3), tab(1,1,1,2), tab(1,1,1,3), tab(1,1,2,3), tab(1,1,2,1), tab(1,1,3,1), tab(1,1,3,2);...
            tab(2,2,1,1), tab(2,2,2,2), tab(2,2,3,3), tab(2,2,1,2), tab(2,2,1,3), tab(2,2,2,3), tab(2,2,2,1), tab(2,2,3,1), tab(2,2,3,2);...
            tab(3,3,1,1), tab(3,3,2,2), tab(3,3,3,3), tab(3,3,1,2), tab(3,3,1,3), tab(3,3,2,3), tab(3,3,2,1), tab(3,3,3,1), tab(3,3,3,2);...
            tab(1,2,1,1), tab(1,2,2,2), tab(1,2,3,3), tab(1,2,1,2), tab(1,2,1,3), tab(1,2,2,3), tab(1,2,2,1), tab(1,2,3,1), tab(1,2,3,2);...
            tab(1,3,1,1), tab(1,3,2,2), tab(1,3,3,3), tab(1,3,1,2), tab(1,3,1,3), tab(1,3,2,3), tab(1,3,2,1), tab(1,3,3,1), tab(1,3,3,2);...
            tab(2,3,1,1), tab(2,3,2,2), tab(2,3,3,3), tab(2,3,1,2), tab(2,3,1,3), tab(2,3,2,3), tab(2,3,2,1), tab(2,3,3,1), tab(2,3,3,2);...
            tab(2,1,1,1), tab(2,1,2,2), tab(2,1,3,3), tab(2,1,1,2), tab(2,1,1,3), tab(2,1,2,3), tab(2,1,2,1), tab(2,1,3,1), tab(2,1,3,2);...
            tab(3,1,1,1), tab(3,1,2,2), tab(3,1,3,3), tab(3,1,1,2), tab(3,1,1,3), tab(3,1,2,3), tab(3,1,2,1), tab(3,1,3,1), tab(3,1,3,2);...
            tab(3,2,1,1), tab(3,2,2,2), tab(3,2,3,3), tab(3,2,1,2), tab(3,2,1,3), tab(3,2,2,3), tab(3,2,2,1), tab(3,2,3,1), tab(3,2,3,2)];
    end
end