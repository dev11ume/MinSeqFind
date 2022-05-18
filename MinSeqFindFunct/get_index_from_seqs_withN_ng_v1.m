% FUNCTION- GET SEQUENES FROM INDEXES FOR SEQUENCES WITH N WITHOUT GAPS
function seqs=get_index_from_seqs_withN_ng_v1(inx,max_nmer)
mot='NACGT';
r=rem(floor((inx-1)*power(5,1-max_nmer:0)),5)+1;
seqs=cellstr(mot(r(:,end:-1:1)));
end
