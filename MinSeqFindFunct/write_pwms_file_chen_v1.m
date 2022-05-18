% FUNCTION TO WRITE PWMs IN CHEN FORMAT 
function write_pwms_file_chen_v1(dir_nm,pathx)
if exist([dir_nm '/' pathx '/resf.mat'])
    load ([dir_nm '/' pathx '/resf.mat'],'pwm_all_arr');
    for j=1:length(pwm_all_arr)
        if ~isempty(pwm_all_arr{j})
            nmx=[dir_nm '/CHEN/' pathx '_' num2str(j)];
            for i=1:length(pwm_all_arr{j})
                if ~isempty(pwm_all_arr{j}{i})
                    nm=[nmx '_' num2str(i)];
                    pwm1=reshape(pwm_all_arr{j}{i},4,length(pwm_all_arr{j}{i})/4);
                    pwm1r=reshape(fliplr(pwm_all_arr{j}{i}')',4,length(pwm_all_arr{j}{i})/4);
                    
                    fid3=fopen([nm '.chen'],'w');
                    fprintf(fid3,'>PWM1\n');
                    for pos=1:size(pwm1,2)
                        for k=1:4
                            fprintf(fid3,'%f\t',pwm1(k,pos)./sum(pwm1(:,pos)));
                        end
                        fprintf(fid3,'\n');
                    end
                    fclose(fid3);
                    
                    fid3=fopen([nm 'r.chen'],'w');
                    fprintf(fid3,'>PWM1\n');
                    for pos=1:size(pwm1r,2)
                        for k=1:4
                            fprintf(fid3,'%f\t',pwm1r(k,pos)./sum(pwm1r(:,pos)));
                        end
                        fprintf(fid3,'\n');
                    end
                    fclose(fid3);
                end
            end
        end
    end
end
