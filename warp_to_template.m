function warp_to_template (u_file, in_img)

addpath(genpath('/project/6007967/sudesnac/HCP/Tools/spm12'));

matlabbatch{1}.spm.tools.dartel.crt_warped.flowfields = {u_file};
matlabbatch{1}.spm.tools.dartel.crt_warped.images = {{in_img}};
matlabbatch{1}.spm.tools.dartel.crt_warped.jactransf = 0;
matlabbatch{1}.spm.tools.dartel.crt_warped.K = 6;
matlabbatch{1}.spm.tools.dartel.crt_warped.interp = 1;
spm_jobman('run', matlabbatch);





end
