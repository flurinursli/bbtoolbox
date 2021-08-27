function [] = plot_rough(normals);

% input strike ranges between 0 and 360
% input dip ranges between 0 and 180
% input rake ranges between -180 and 180
 
FolderInfoPatch = dir('ROUGH_PATCH*.BBT');
FolderInfoFSDR  = dir('ROUGH_FSDR*.BBT');
FolderInfoNRM   = dir('ROUGH_NRM*.BBT');

nplanes = size(FolderInfoPatch,1);

fig1 = figure(1); hold on;
fig2 = figure(2); hold on;
fig3 = figure(3); hold on;

for ip = 1:nplanes

   FolderInfoPatch(ip)

   fid = fopen(FolderInfoPatch(ip).name,'r');
   tmp = fread(fid,9,'real*4');
   npatch = tmp(1);
   tmp = fread(fid,3*3*npatch,'real*4');
   tmp = reshape(tmp,[3 3 npatch]);
   fclose(fid);

   EAST  = squeeze(tmp(:,1,:));
   NORTH = squeeze(tmp(:,2,:));
   VERT  = squeeze(tmp(:,3,:));

   fid = fopen(FolderInfoFSDR(ip).name,'r');
   tmp = fread(fid,5,'real*4');

   STRIKE_REF = tmp(2);
   DIP_REF    = tmp(3);

   STRIKE_REF = STRIKE_REF * ones(npatch,1);
   DIP_REF    = DIP_REF * ones(npatch,1);

   tmp = fread(fid,5*npatch,'real*4');
   tmp = reshape(tmp,[5 npatch]);

   FRACT    = squeeze(tmp(1,:));
   STRIKE   = squeeze(tmp(2,:));
   DIP      = squeeze(tmp(3,:));
   RAKE     = squeeze(tmp(4,:));

   RAKE_REF = squeeze(tmp(5,:));

   fclose(fid);


   fid = fopen(FolderInfoNRM(ip).name,'r');
   tmp = fread(fid,6*npatch,'real*4');
   tmp = reshape(tmp,[6 npatch]);
   fclose(fid);

   EAST_C  = squeeze(tmp(1,:));
   NORTH_C = squeeze(tmp(2,:));
   VERT_C  = squeeze(tmp(3,:));
   NORM    = squeeze(tmp(4:6,:));

   clear tmp;

   %{
   fid = fopen('ROUGH_SV.BBT','r');

   tmp = fread(fid,6*npatch,'real*4');

   tmp = reshape(tmp,[6 npatch]);

   fclose(fid);

   SV = squeeze(tmp(4:6,:));

   clear tmp;
   %}

   if (npatch > 1200)
      style = 'None';
   else
      style = '-';
   end

   % compute deviations from smooth case
   for i = 1:npatch

      DSTRIKE(i) = STRIKE(i) - STRIKE_REF(i);
      if (abs(DSTRIKE(i)) > 180) 
         DSTRIKE(i) = DSTRIKE(i) - sign(DSTRIKE(i))*360; 
         STRIKE(i) = STRIKE(i) + sign(DSTRIKE(i))*360;
      end
   
      DDIP(i) = DIP(i) - DIP_REF(i);

      %if (DIP(i) > 90); DIP(i) = DIP(i) - 180; end;
 
      DRAKE(i) = RAKE(i) - RAKE_REF(i);
      if (abs(DRAKE(i)) > 180)
         DRAKE(i) = DRAKE(i) - sign(DRAKE(i))*360;
         RAKE(i) = RAKE(i) + sign(DRAKE(i))*360;
      end      

   end
 
   'STRIKE ', min(STRIKE_REF),max(STRIKE_REF),min(STRIKE),max(STRIKE)
   'DIP ', min(DIP_REF),max(DIP_REF),min(DIP),max(DIP)
   'RAKE ', min(RAKE_REF),max(RAKE_REF),min(RAKE),max(RAKE)

   % set azimuth elevation for camera
   if (ip == 1) 
      AZ = max(STRIKE_REF(:)) + 50;
      EL = 35;
   end


   figure(fig1)
   subplot(1,2,1)
   patch(EAST,NORTH,VERT,DSTRIKE,'LineStyle',style); colormap jet; hold on;
   if (normals == 'y')
      quiver3(EAST_C,NORTH_C,VERT_C,NORM(1,:),NORM(2,:),NORM(3,:),'r')
   end
   view([AZ EL]);
   xlabel('EAST'); ylabel('NORTH'); zlabel('Z')
   axis equal;
   set(gca,'ZDir','Reverse');
   grid on;
   colorbar;

   subplot(1,2,2)
   patch(EAST,NORTH,VERT,STRIKE,'LineStyle',style); colormap jet; hold on;
   view([AZ EL]);
   axis equal;
   set(gca,'ZDir','Reverse');
   xlabel('EAST'); ylabel('NORTH'); zlabel('Z')
   grid on;
   colorbar;
   avg_strike(ip) = mean(STRIKE);
   orig_strike(ip) = mean(STRIKE_REF);


   figure(fig2)
   subplot(1,2,1)
   patch(EAST,NORTH,VERT,DDIP,'LineStyle',style); colormap jet; hold on;
   view([AZ EL]);
   axis equal;
   set(gca,'ZDir','Reverse');
   xlabel('EAST'); ylabel('NORTH'); zlabel('Z')
   grid on;
   colorbar;
   avg_dip(ip) = mean(abs(DIP));

   subplot(1,2,2)
   patch(EAST,NORTH,VERT,DIP,'LineStyle',style); colormap jet; hold on;
   view([AZ EL]);
   axis equal;
   set(gca,'ZDir','Reverse');
   xlabel('EAST'); ylabel('NORTH'); zlabel('Z')
   grid on;
   colorbar;
   avg_dip(ip) = mean(abs(DIP));
   orig_dip(ip) = mean(abs(DIP_REF));


   figure(fig3)
   subplot(1,2,1)
   patch(EAST,NORTH,VERT,DRAKE,'LineStyle',style); colormap jet; hold on;
   view([AZ EL]);
   axis equal;
   set(gca,'ZDir','Reverse');
   xlabel('EAST'); ylabel('NORTH'); zlabel('Z')
   grid on;
   colorbar;

   subplot(1,2,2)
   patch(EAST,NORTH,VERT,RAKE,'LineStyle',style); colormap jet; hold on;
   view([AZ EL]);
   axis equal;
   set(gca,'ZDir','Reverse');
   xlabel('EAST'); ylabel('NORTH'); zlabel('Z')
   grid on;
   colorbar;
   avg_rake(ip) = mean(RAKE);
   orig_rake(ip) = mean(RAKE_REF);

   clear EAST NORTH VERT RAKE DRAKE DIP DDIP STRIKE DSTRIKE
   clear RAKE_REF DIP_REF STRIKE_REF
 
end

text = 'Strike delta - ';
for ip = 1:nplanes
   text = [text,'s',num2str(ip),'=',num2str(orig_strike(ip),'%5.1f'),' '];
end
figure(fig1); subplot(1,2,1),title(text);

text = 'Final Strike - ';
for ip = 1:nplanes
   text = [text,'s',num2str(ip),'=',num2str(avg_strike(ip),'%5.1f'),' '];
end
figure(fig1); subplot(1,2,2),title(text)


text = 'Dip delta - ';
for ip = 1:nplanes
   text = [text,'d',num2str(ip),'=',num2str(orig_dip(ip),'%5.1f'),' '];
end
figure(fig2); subplot(1,2,1),title(text);

text = 'Final Dip - ';
for ip = 1:nplanes
   text = [text,'d',num2str(ip),'=',num2str(avg_dip(ip),'%5.1f'),' '];
end
figure(fig2); subplot(1,2,2),title(text);


text = 'Rake delta - ';
for ip = 1:nplanes
   text = [text,'r',num2str(ip),'=',num2str(orig_rake(ip),'%5.1f'),' '];
end
figure(fig3); subplot(1,2,1),title(text);

text = 'Final Rake - ';
for ip = 1:nplanes
   text = [text,'r',num2str(ip),'=',num2str(avg_rake(ip),'%5.1f'),' '];
end
figure(fig3); subplot(1,2,2),title(text);

%figure(fig2); subplot(1,2,1),title(['Dip delta - avg=',num2str(avg_strike/nplanes)]);
%figure(fig3); subplot(1,2,1),title(['Rake delta - avg=',num2str(avg_strike/nplanes)]);
