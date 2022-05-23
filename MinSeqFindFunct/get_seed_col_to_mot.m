% FUNCTION - GET seed_mot FROM seed_col MATRIX
function seed_mot=get_seed_col_to_mot(seed_col)
nucl='ACMGRSVTWYHKDBN';
l=length(seed_col)/4;
seed_mot=zeros(1,l);
xx=[power(2,[3:-1:0])]';
for i=1:l
    pos=seed_col(4*i:-1:4*i-3)*xx;
    seed_mot(i)=nucl(pos);
end
seed_mot=char(seed_mot);
end
