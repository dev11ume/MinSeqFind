% FUNCTION TO GET PALINDROME INDEXES FOR K-MERS WITHOUT GAPS 
function inx_pal=get_palindrome_inx_ng_v1(inx,tot_l)
inx_rev=get_rev_inx_ng_v1(inx,tot_l);
inx_pal=inx(inx==inx_rev);
end
