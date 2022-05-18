% FUNCTION- GET seed_col MATRIX FROM seed_mot
function seed_col=get_mot_to_seed_col(seed_mot)

nucl='ACMGRSVTWYHKDBN';

l=length(seed_mot);
seed_col=zeros(1,4*l);
for i=1:l
    pos=find(nucl==seed_mot(i));
    seed_col(4*i:-1:4*i-3)=bin2dec(dec2bin(pos,4)');
end

end
