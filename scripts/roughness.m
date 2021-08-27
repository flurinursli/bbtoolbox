function [] = roughness(plane, iteration);

% files are named as node_<refinement>_<plane>_<iteration>.bin
FolderInfoSRT  = dir(['roughness_*' num2str(plane) '_' num2str(iteration) '.bin']);

nrefins = size(FolderInfoSRT, 1);

disp(['Number of refinements... ' num2str(nrefins)])

fig1 = figure; hold on;   

for ref = 1:nrefins

   %FolderInfoSRT(ip)

   disp(['Working on file... ' FolderInfoSRT(ref).name])

   fid = fopen(FolderInfoSRT(ref).name,'r');

   nutr = fread(fid, 1, 'int32');
   nvtr = fread(fid, 1, 'int32');

   npatches = nvtr * (2*nutr - 1);

   if (npatches > 3000e+00)
      style = 'None';
   else
      style = '-';
   end

   x = zeros(3, npatches); 
   y = zeros(3, npatches);
   c = zeros(3, npatches); 

   LOWER_EDGE = zeros(2, nutr + 1 - mod(nvtr,2)); 

   np = 0;

   for j = 1:nvtr
      for i = 1:2*nutr - 1

         np = np + 1;

	 for icr = 1:3
            x(icr, np) = fread(fid, 1, 'real*4');  
            y(icr, np) = fread(fid, 1, 'real*4');  
            c(icr, np) = fread(fid, 1, 'real*4');
         end

         pt = round(i / 2 + 0.5);
         vt = round(j / 2 + 0.5);

         if (j == nvtr & mod(i,2) == 1)  
            if (mod(nvtr,2) == 0) 
               LOWER_EDGE(:, pt) = [x(1,np) y(1,np)];
               if (i == 2*nutr - 1); LOWER_EDGE(:, pt + 1) = [x(3,np) y(3,np)]; end;
            else
               LOWER_EDGE(:, pt) = [x(2,np) y(2,np)];
            end
         end

      end
   end

   figure(fig1)
   %patch(y, x, mean(c,1), 'LineStyle', style); 
   patch(x, y, c, 'LineStyle', style); 

   plot(LOWER_EDGE(1,:), LOWER_EDGE(2,:), 'k', 'LineWidth', 1);

   sigma = fread(fid, 1, 'real*4');

   fclose(fid);

end

figure(fig1); set(gca, 'YDir', 'Reverse'); grid on; axis equal; colormap jet;

% interpolate over a regular grid of points
F = scatteredInterpolant(reshape(x, [3*np,1]), reshape(y, [3*np,1]), reshape(c, [3*np,1]), 'linear');

dh = min([x(3,1)-x(1,1), y(2,1)-y(1,1)])

xo = [mean(x(:,1)):dh:mean(x(:,2*nutr - 1))];
yo = [mean(y(1:2,1)):dh:mean(y(1:2,npatches))];

[XO,YO] = meshgrid(xo,yo);

ZO = F(XO, YO);

size(ZO)

figure, pcolor(XO,YO,ZO); shading interp; set(gca, 'YDir', 'Reverse'); axis equal; grid on;

% at this point, matrix "ZO" is DIPxSTRIKE  

%figure, pcolor(XO(:,1:2),YO(:,1:2),ZO(:,1:2)); set(gca, 'YDir', 'Reverse'); axis equal; grid on;

for d = 1:2

   if (d == 1) 

      text = 'PSD along strike';

      corr = zeros(size(ZO'));

      for i = 1:size(corr,2)
         corr(:,i) = autocorr(ZO(i,:));      % AC along strike
      end 

   elseif (d == 2) 

      text = 'PSD along dip'; 

      corr = zeros(size(ZO));

      for i = 1:size(corr,2)
         corr(:,i) = autocorr(ZO(:,i));      % AC along dip
      end 

   end

   n = size(corr, 1);

   psd = zeros(n, size(corr, 2));

   for i = 1:size(corr, 2)
      psd(:, i) = pwr(corr(:, i));
   end

   % theoretical 1D PSD using discrete resolution ('PWR' double number of points)
   ks = 1 / dh;
   k  = 2 * pi * ks * (0:n-1) / (2*n - 1);
   dk = 2 * pi / (2*n - 1) / dh;

   fun = sigma * k.^-3; 
   %fun = fun / dh / n;

   dpsd = mean(psd, 2);

   figure; hold on;

   title(text);

   for i = 1:size(psd, 2)
      plot(k, psd(:, i), 'Color', [0.75 0.75 0.75], 'LineWidth', 0.25);
   end

   h1 = plot(k, dpsd, 'b', 'LineWidth', 2);
   %h2 = plot(k, fun, 'r', 'LineWidth', 2);

   set(gca, 'Xscale', 'log', 'Yscale', 'log');

   grid on;

end

end

%-------------------------------------------------------------------------------

function x = autocorr(y);
% one-sided autocorrelation of input sequence

    n = length(y);

    z = xcorr(y, 'biased');

    x = z(n:end);

end

%-------------------------------------------------------------------------------

function X = pwr(x);
% PSD from input autocorrelation

    % length of one-sided acr
    n = length(x);

    % length of two sided: we will have always odd number of points
    npts = 2*n - 1;

    y = zeros(npts, 1);

    % negative lags
    y(1:n - 1) = x(end:-1:2);
    % positive lags
    y(n:end)   = x;

    Y = fft(y, npts);

    % one-sided amplitude spectrum: since npts is always odd, we never
    % cover nyquist wavenumber
    X        = abs(Y(1:n)) / npts;
    X(2:end) = 2 * X(2:end);

end
