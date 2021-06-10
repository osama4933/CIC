function symbol_reformat(demod_sym_stack)
%% write symbols to txt file, to be fed to RPP0
% Each Symbol appears on a column and symbols of different packets are separated
% with -1 identifier
fileID_1 = fopen('symbols.txt','w');
demod_sym_stack(find(demod_sym_stack < 0)) = 0;
for i = 1:size(demod_sym_stack,1)
    fprintf(fileID_1,'%d\n',-1);
    fprintf(fileID_1,'%d\n',demod_sym_stack(i,:));
end
fclose(fileID_1);

end
