function S_out = combineStructs(S_1, S_2)
%COMBINESTRUCTS  Combines structures "S_1","S_2" into single struct "S_out"

% Get fieldnames
f_1 = fieldnames(S_1);
f_2 = fieldnames(S_2);

%Check for clashes
for i = 1:size(f_1,1)
    f_i = f_1{i};
    if max(strcmp(f_i,f_2)) %strcmp returns logic array.
        %Pop-up msg; programme continues after clicking
        errordlg({'The following field is repeated in both structs:';'';...
            f_i;'';...
            'To avoid unintentional overwrite of this field, the program will now exit with error'},...
            'Inputs Invalid!')
        %Exit message; forces quit from programme gives a pop-up
        error('Exit due to repeated field in structure')
    end
end

% Combine  into a single output structure
for i = 1:size(f_1,1)
    S_out.(f_1{i}) = S_1.(f_1{i});
end
for i = 1:size(f_2,1)
    S_out.(f_2{i}) = S_2.(f_2{i});
end

end

