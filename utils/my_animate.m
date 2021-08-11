function my_animate(gcf,fig_pth,whos_data)
    % rotates a 3D plot through all possible angles and saves as a gif
    axis tight
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    set(gca,'zticklabel',[]);
    az = 0;
    el = 90;
    view([az,el]);
    degStep = 1;
    detlaT = 0.1;
    fCount = 71;
    f = getframe(gcf);
    [im,map] = rgb2ind(f.cdata,256,'nodither');
    im(1,1,1,fCount) = 0;
    k = 1;
    % spin 45°
    for i = 0:-degStep:-45
      az = i;
      view([az,el]);
      f = getframe(gcf);
      im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
      k = k + 1;
    end
    % tilt down
    for i = 90:-degStep:15
      el = i;
      view([az,el]);
      f = getframe(gcf);
      im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
      k = k + 1;
    end
    % spin left
    for i = az:-degStep:-90
      az = i;
      view([az,el]);
      f = getframe(gcf);
      im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
      k = k + 1;
    end
    % spin right
    for i = az:degStep:0
      az = i;
      view([az,el]);
      f = getframe(gcf);
      im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
      k = k + 1;
    end
    % tilt up to original
    for i = el:degStep:90
      el = i;
      view([az,el]);
      f = getframe(gcf);
      im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
      k = k + 1;
    end
    fig_name = fullfile(fig_pth,whos_data,'state_space.gif');
    imwrite(im,map,fig_name,'DelayTime',detlaT,'LoopCount',inf)
end % my_animate