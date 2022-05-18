% FUNCTION TO GET REVERSE COMPLEMENT INDEXES 
function inx_rev=get_rev_inx_v2(inx,max_half_nmer,tot_l)
if max_half_nmer==-1
    inx_rev=get_rev_inx_ng_v2(inx,tot_l);
else
    inx_rev=get_rev_inx_gp_v2(inx,max_half_nmer);
end
end
