% FUNCTION- GET SEQUENES FROM INDEXES FOR SEQUENCES WITH N
function seqs2=get_index_from_seqs_withN_v1(inx,tot_l,max_half_nmer)
if max_half_nmer==-1
    seqs2=get_index_from_seqs_withN_ng_v1(inx,tot_l);
else
    seqs2=get_index_from_seqs_withN_gp_v1(inx,tot_l,max_half_nmer);
end
end
