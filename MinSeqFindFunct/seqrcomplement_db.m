% FUNCTION TO GET SEQUENCE REVERSE COMPLEMENT 
function y=seqrcomplement_db(x)
  y=flip(x);
  y=lower(y);
  y=strrep(y,'a','T');
  y=strrep(y,'c','G');
  y=strrep(y,'g','C');
  y=strrep(y,'t','A');
  y=strrep(y,'n','N');
  y=strrep(y,'k','M');
  y=strrep(y,'m','K');
  y=strrep(y,'r','Y');
  y=strrep(y,'y','R');
  y=strrep(y,'s','W');
  y=strrep(y,'w','S');
  y=strrep(y,'b','V');
  y=strrep(y,'v','B');
  y=strrep(y,'h','D');
  y=strrep(y,'d','H');
end
