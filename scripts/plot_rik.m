function [] = plot_rik(band, velocity, iteration, field);


if (field ~= "slip" & field ~= "rupture" & field ~= "rise" & field ~= "strike" & ...
    field ~= "dip"  & field ~= "rake" & field ~= "rough")
    error('3rd argument is invalid. Expected values: slip, rupture, rise, strike, dip, rake, rough');
end

switch field
   case "slip"
      skip = 0;
   case "rise"
      skip = 3;
   case "rupture"
      skip = 6;
   case "strike"
      skip = 9;
   case "dip"
      skip = 10;
   case "rake"
      skip = 11;
   case "rough"
      skip = 14;
end

% plot all fault plane(s), including their geometry, slip 


% files are named as node_<plane>_<velocity>_<iteration>.bin
FolderInfoSRT  = dir(['node_band' num2str(band) '*_vel' num2str(velocity) '_iter' num2str(iteration) '.bin']);

nplanes = size(FolderInfoSRT, 1);

disp(['Number of planes... ' num2str(nplanes)])

fig = figure; hold on;  

for ip = 1:nplanes

   %FolderInfoSRT(ip)

   disp(['Working on file... ' FolderInfoSRT(ip).name])

   fid = fopen(FolderInfoSRT(ip).name,'r');

   % total number of patches (sum over all refinements)
   npatches = fread(fid, 1, 'int32');

   nrefins = fread(fid, 1, 'int32');

   if (npatches > 2000e+00)
      style = 'None';
   else
      %style = '-';
      style = 'None';
   end

   x = zeros(3, npatches); 
   y = zeros(3, npatches);
   z = zeros(3, npatches); 

   if (field == "strike" | field == "dip")
      values = zeros(1, npatches);
   else
      values = zeros(3, npatches);
   end

   np = 0;

   for ref = 1:nrefins

      nvtr = fread(fid, 1, 'int32');
      nutr = fread(fid, 1, 'int32');

      if (ref == 1); UPPER_EDGE = zeros(3, nutr + 1); end;
      if (ref == nrefins); LOWER_EDGE = zeros(3, nutr + 1 - mod(nvtr,2)); end;

      LEFT_EDGE = zeros(3, round(nvtr/2 + 0.5));
      RIGHT_EDGE = zeros(3, round(nvtr/2 + 0.5));

      for j = 1:nvtr
         for i = 1:2*nutr - 1

            np = np + 1;

            y(:, np) = fread(fid, 3, 'real*4');  
            x(:, np) = fread(fid, 3, 'real*4');  
            z(:, np) = fread(fid, 3, 'real*4');   

	    fread(fid, skip, 'real*4');
            
            values(:,np) = fread(fid, min(size(values)), 'real*4');
            
            fread(fid, 17-skip-min(size(values)), 'real*4');

            pt = round(i / 2 + 0.5);
            vt = round(j / 2 + 0.5);

            if (ref == 1 & j == 1 & mod(i,2) == 1) 
               UPPER_EDGE(:, pt) = [x(1,np) y(1,np) z(1,np)];
               if (i == 2*nutr - 1); UPPER_EDGE(:, pt + 1) = [x(3,np) y(3,np) z(3,np)]; end;
            end
            
            if (ref == nrefins & j == nvtr & mod(i,2) == 1)  
               if (mod(nvtr,2) == 0) 
                  LOWER_EDGE(:, pt) = [x(1,np) y(1,np) z(1,np)];
	          if (i == 2*nutr - 1); LOWER_EDGE(:, pt + 1) = [x(3,np) y(3,np) z(3,np)]; end;
               else
                  LOWER_EDGE(:, pt) = [x(2,np) y(2,np) z(2,np)];
               end
            end

            if (mod(j,2) == 1)
               if (i == 1); LEFT_EDGE(:, vt) = [x(1,np) y(1,np) z(1,np)]; end;
               if (i == 2*nutr - 1); RIGHT_EDGE(:, vt) = [x(3,np) y(3,np) z(3,np)]; end;
            end

            if (j == nvtr)
               if (i == 1); LEFT_EDGE(:, vt) = [x(1,np) y(1,np) z(1,np)]; end;
	       if (i == 2*nutr - 1); RIGHT_EDGE(:, vt) = [x(3,np) y(3,np) z(3,np)]; end;
            end

         end
      end

      figure(fig);
      plot3(LEFT_EDGE(1,:), LEFT_EDGE(2,:), LEFT_EDGE(3,:), 'k', 'LineWidth', 1);
      plot3(RIGHT_EDGE(1,:), RIGHT_EDGE(2,:), RIGHT_EDGE(3,:), 'k', 'LineWidth', 1);      

   end   

   scaling = fread(fid, 1, 'real*4');

   % rise(slip <= 0)    = NaN;
   % rupture(slip <= 0) = NaN;
   % slip(slip <= 0)    = NaN;

   if (field == "slip")
      values = values * scaling; 
   elseif (field == "strike" | field == "dip" | field == "rake")
      values = values * 180/pi;
   end

   figure(fig);

   patch(x, y, z, values, 'LineStyle', style); 

   plot3(UPPER_EDGE(1,:), UPPER_EDGE(2,:), UPPER_EDGE(3,:), 'r', 'LineWidth', 3); 
   plot3(LOWER_EDGE(1,:), LOWER_EDGE(2,:), LOWER_EDGE(3,:), 'k', 'LineWidth', 1);

   fclose(fid);

end

figure(fig); set(gca, 'ZDir', 'Reverse'); axis equal; grid on; colormap jet;
title(field), xlabel('EW'), ylabel('NS'), zlabel('Z');


%figure(fig4); colorbar; c = get(gca, 'colorbar');

%for i = 1:length(c.Ticks)
%   if (c.Ticks(i) < 0); c.TickLabels{i} = num2str(360 + c.Ticks(i)); end;
%end




end

