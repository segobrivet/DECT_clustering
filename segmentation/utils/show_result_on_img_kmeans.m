
function show_result_on_img_kmeans(T, img_4D, focus, vars)

for i=1:min(vars.max_slic_subplot,length(vars.slic_show))
    
    slic_r = img_4D(:,:,vars.slic_show(i),1);
    slic_g = slic_r; slic_b = slic_r;
    
    subplot(vars.len_sp,4,2*i);
    for cl_id=1:vars.K
        cc = find(T(:,:,i) == cl_id);
        if cl_id == vars.max_cl
            slic_r(cc) = 1; slic_g(cc) = 0; slic_b(cc) = 0;  % red
        else
            slic_r(cc) = vars.clr(cl_id,1); slic_g(cc) = vars.clr(cl_id,2); slic_b(cc) = vars.clr(cl_id,3);
        end
        slic_rgb = cat(3,slic_r, slic_g, slic_b);
    end
    imshow(slic_rgb.*focus(:,:,vars.slic_show(i)))
    hold on
    plot_tumor_contour(vars.tumor_contour_list{i}, [vars.rmin, vars.cmin], [0.99,0.99,0.99]);
end