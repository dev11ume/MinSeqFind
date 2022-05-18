% FUNCTION seq_col & seqs_mat MATRIX FROM SEQUENCES v2
function [seq_col,seqs_mat,l_seqs]=get_seq_col_v2(seqs,Llib,Rlib)
aa2=[1 0 2 0 0 0 3 0 0 0 0 0 0 5 0 0 0 0 0 4 0 0 0 0 0 0];
aa2(56)=5; %FOR x
nucl_to_mat=[eye(4) ; ones(1,4)]; % FOR N
if iscell(seqs)
    l_seqs=length(seqs);
    l=length(seqs{1});
    seqs_mat=[repmat(aa2(Llib-64),[l_seqs 1]) aa2(cell2mat(seqs)-64) repmat(aa2(Rlib-64),[l_seqs 1])];
else
    [l_seqs,l]=size(seqs);
    seqs_mat=[repmat(aa2(Llib-64),[l_seqs 1]) aa2(seqs-64) repmat(aa2(Rlib-64),[l_seqs 1])];
end
seq_col=reshape(nucl_to_mat(seqs_mat',:)',(l+length(Llib)+length(Rlib))*4,l_seqs)';
end
