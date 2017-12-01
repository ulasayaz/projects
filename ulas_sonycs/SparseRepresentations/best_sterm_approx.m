function [sel,sr]=best_sterm_approx(w, layout, fig)
% returns figure handle

if nargin <1
    w=Curvelet2_clab();
end
if nargin <2
    layout=1;
end
if nargin <3
    fig=1;
end

sel={'einstein','barbara','cameraman','peppers'};
selfn={'einstein.tiff','barbara.pgm','cameraman.bmp','peppers512x512.tiff'};
imgfilefunc=@(fn) fullfile(get_projectpath(),'testdata','images',fn);

sr=SparseRep(w); 
sr.imgfiles=cellfun(imgfilefunc,selfn,'UniformOutput',false);
sr.do_ssim=true;
sr.fig=fig;
sr.w.set_wcoarsestlev(w.wcL);
sr.figopt.pixsizeX=1000;  % for ipython notebook figures width should be 1000 pix max.
suppl_output=layout;  % output of each image into separate figure
sr.show_quality_images({1,2,{0,1,2}},suppl_output); 


end

