clear all
close all
%extract data from PIVlab batch processing
%iterate through all the root folders, find the folders for each timepoint,
%extract the data into a single file, save as csv
extractroots = ["Jan28_WT-C1TA_OD1\median_blur","Feb2_WT-C1TA_OD01\median_blur","Feb3_WT-C1TA_OD001\median_blur","Feb4_WT_OD1\median_blur","Feb7_WT-C1TA_OD0001","Feb8_WT-C1TA_OD001\median_blur","Feb22_WT-C1TA_OD0001","Feb23_WT-C1TA_OD1","Feb24_WT-C1TA_OD001","Mar1_WT-C1TA_OD01","Mar3_WT-C1TA_OD1","Mar7_WT-C1TA_OD01","Mar9_WT-C1TA_OD0001"];


for bigloop = 1:size(extractroots,2);
rootFold = strcat("D:\Sean\SurfaceColonyPIV\",extractroots(bigloop));

disp(rootFold);
d = dir(rootFold);
folds = d([d(:).isdir]) ;
folds = folds(~ismember({folds(:).name},{'.','..'}));
folds = {folds.name};

title = extractroots(bigloop);
if contains(extractroots(bigloop),"\median_blur")
    title = extractBefore(extractroots(bigloop),"\median_blur");
end;
write_title = strcat(title,"_new_v-mags",".csv");
%writematrix([0 0 0 0 0], write_title); %initialize the destination file (drop the zeros in R)

for fold = 1:size(folds,2);
    disp(folds(fold));
    disp(strcat(string(fold)," of ",string(size(folds,2))));
   directory = strcat(rootFold,filesep,char(folds(fold)));
   hold_dir = directory;
   %find the .mat file in this directory and load it
   load(strcat(directory, filesep, dir(fullfile(directory, "*.mat")).name))
   directory = hold_dir; %in case loading the .mat file fucks up directory's value
   %apply the correlation filter. Throws out vectors from coordinates that
   %had poor correlation between image pairs (i.e. noise)
   %since correlation is a bit noisy between frames, take the maximum
   %projection and use that instead
   correlation_med_proj = imgaussfilt(median(cat(3, correlation_map{:}),3),2);
   med_map = correlation_map; %set up an array of the same size
   for i=1:size(correlation_map,1)
    med_map{i} = correlation_med_proj;
   end
   	for PIVresult=1:size(x,1)
	
        %[u_filt{PIVresult,1},v_filt{PIVresult,1}] = PIVlab_correlation_filter(u{PIVresult,1},v{PIVresult,1},0.75,correlation_map{PIVresult,1});
        [u_filt{PIVresult,1},v_filt{PIVresult,1}] = PIVlab_correlation_filter(u{PIVresult,1},v{PIVresult,1},0.5,med_map{PIVresult,1});
    end;
   
   % put the two vectors beside each other
   %vectors = [cellfun(@(x) fillmissing(x,'constant',0), u_filt,'UniformOutput',false) cellfun(@(x) fillmissing(x,'constant',0), v_filt,'UniformOutput',false)];
    vectors = [u v];
    vectors_filt = [u_filt v_filt];
   
   %get the vector magnitude
   mag = cell(size(vectors,1),1);
   mag_filt = cell(size(vectors,1),1);
   mean_u = mean(cat(3,u_filt{:}),3,'omitnan');
   mean_v = mean(cat(3,v_filt{:}),3,'omitnan');
   for j = 1:size(vectors,1);
        mag{j} = (sqrt(vectors{j,1}.^2+vectors{j,2}.^2));

   end;
   for j = 1:size(vectors_filt,1);
        mag_filt{j} = (sqrt(vectors_filt{j,1}.^2+vectors_filt{j,2}.^2));

   end;
   val = mean(cat(3,mag{:}),3,'omitnan'); %magnitude first, then mean : mean_vector-magnitude
   val2 = sqrt(mean_u.^2+mean_v.^2); %mean first, then magnitude : mean-vector_magnitude
   %imwrite(uint16(val*10000),directory+"\mean_vector-magnitude.png");
   %imwrite(uint16(val2*10000),directory+"\mean-vector_magnitude.png");
   %imwrite(uint16(correlation_med_proj*10000),directory+"\median_correlation_map.png");
   bfsave(val2,strcat(convertStringsToChars(directory),'\mean-vector_magnitude.tif'),'bigtiff',true);
   bfsave(correlation_med_proj,strcat(convertStringsToChars(directory),'\median_correlation_map.tif'),'bigtiff',true);
   %val = 0;
   %for a = 1:size(dat,1);
    %val = val + dat{a,1};
   %end;
   %val = val/size(dat,1);
   
   [x y] = meshgrid(1:size(val,2),1:size(val,1));
   output = [repmat(str2num(folds{fold}),size(x(:),1),1) x(:) y(:) val(:) val2(:)];
   writematrix(output, write_title, 'WriteMode','append');
end;




end;