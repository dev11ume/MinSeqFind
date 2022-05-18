% FUNCTION TO GET seqs_mat MATRIX FROM SEQUENCES
function [seqs_mat,l_seqs]=get_seq_col_mat_only_v1(seqs,Llib,Rlib)
aa2=[1 0 2 0 0 0 3 0 0 0 0 0 0 5 0 0 0 0 0 4 0 0 0 0 0 0];
aa2(56)=5; %FOR x
if iscell(seqs)
    l_seqs=length(seqs);
    l=length(seqs{1});
    seqs_mat=[repmat(aa2(Llib-64),[l_seqs 1]) aa2(cell2mat(seqs)-64) repmat(aa2(Rlib-64),[l_seqs 1])];
else
    [l_seqs,l]=size(seqs);
    seqs_mat=[repmat(aa2(Llib-64),[l_seqs 1]) aa2(seqs-64) repmat(aa2(Rlib-64),[l_seqs 1])];
end
end
