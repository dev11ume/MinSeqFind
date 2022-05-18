% FUNCTION TO GET PALINDROME INDEXES FOR MINSEQS WITH GAPS 
function inx_pal=get_palindrome_inx_gp_v1(inx,max_half_nmer)
inx_rev=get_rev_inx_gp_v2(inx,max_half_nmer);
inx_pal=inx(inx==inx_rev);
end
