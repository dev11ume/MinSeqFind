% FUNCTION - GET seed_mot FROM seed_col MATRIX
function seed_mot=get_seed_col_to_mot(seed_col)

nucl='ACMGRSVTWYHKDBN';
l=length(seed_col)/4;
seed_mot=zeros(1,l);
for i=1:l
    pos=bin2dec(num2str(seed_col(4*i:-1:4*i-3)));
    seed_mot(i)=nucl(pos);
end
seed_mot=char(seed_mot);
end
