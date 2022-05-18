% FUNCTION TO GET seqs_mat MATRIX FROM SEQUENCES
function [seqs_mat,l_seqs]=get_seq_mat_v1(seqs,Llib,Rlib)
l_seqs=length(seqs);
aa=[1 0 2 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0];
seqs_mat=[repmat(aa(Llib-64),[l_seqs 1]) aa(cell2mat(seqs)-64) repmat(aa(Rlib-64),[l_seqs 1])];
end
