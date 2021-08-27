function [] = plot_sw4(fname, fcut, poles, pass, der, intg, plotspec, samefig);

%e.g. plot_sw4('REF.txt', [8 0], 2, 2, 0, 0, 1, 1)

fontsize = 14;
linewidth = 1;

mytext{1} = 'NS';
mytext{2} = 'EW';
mytext{3} = 'UD';

startRow = 14;

%formatSpec = '%12f%25f%25f%f%[^\n\r]';
formatSpec = '%f%f%f%f%[\n\r]';

fileID = fopen(fname,'r');

%dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'HeaderLines' ,startRow-1, 'EndOfLine', '\r\n');

fclose(fileID);

data = [dataArray{1:end-1}];

clear filename startRow formatSpec fileID dataArray ans;

dt = data(2,1) - data(1,1);
fs = 1 / dt;

% derive
for i = 1:der
   for j = 2:4
      data(:,j) = [0 diff(data(:,j))'] / dt;
   end
end

% integrate
for i = 1:intg
   for j = 2:4
      data(:,j) = cumtrapz(data(:,1), data(:,j));
   end
end
	 	

%data(:, 2:4) = data(:, 2:4).^2;

% LP butterworth filter
if fcut(2) > 0 & fcut(1) == 0
   
   [b, a] = butter(poles, fcut(2) / (fs / 2), 'low');
   %[b, a] = cheby1(poles, 3.486, fcut(1) / (fs / 2), 'low');
      
   for j = 2:4
      if pass == 2	
         data(:,j) = filtfilt(b, a, data(:,j));
      else
         data(:,j) = filter(b, a, data(:,j));
      end
   end

% HP butterworth filter   
elseif fcut(2) == 0 & fcut(1) > 0
  
   [b, a] = butter(poles, fcut(1) / (fs / 2), 'high');
   %[b, a] = cheby1(poles, 3.486, fcut(2) / (fs / 2), 'high');
   
   for j = 2:4
      if pass == 2	
         data(:,j) = filtfilt(b, a, data(:,j));
      else
         data(:,j) = filter(b, a, data(:,j));
      end
   end

% BANDPASS butterworth filter
elseif fcut(1) ~= 0 & fcut(2) ~= 0

   [b, a] = butter(poles, [fcut(1) fcut(2)] / (fs / 2));

   for j = 2:4
      if pass == 2
         data(:,j) = filtfilt(b, a, data(:,j));
      else
         data(:,j) = filter(b, a, data(:,j));
      end
   end

end

if nargin == 6; plotspec = 0; end;
if nargin == 6 | nargin == 7; samefig = 0; end;

if samefig == 0
   figure;
   line = '-';
else  
   figure(gcf);
   line = '-';
end    
    
if plotspec == 1
   
   for j = 2:4

      c = 2*j - 3;
      
      subplot(3, 2, c), plot(data(:,1), data(:,j), 'LineStyle', line, 'LineWidth', linewidth, 'LineStyle', '-'); grid on;
      set(gca, 'LineWidth', 0.1);
      set(gca, 'fontsize', fontsize);
      
      hold on;
      
      title(mytext{j-1});		
     
      if j == 4
        xlabel('Time (s)', 'fontsize', fontsize + 2);
      end  
      
      % length of signal
      l = length(data(:,1));

      n = 2^nextpow2(l);
      n = l;

      Y = fft(data(:,j), n);
     
      %psd = abs(Y(1:n/2 + 1)).^2;
      %psd = psd / (fs * n);

      psd = abs(Y(1:n/2 + 1)) / n;
      
      %mean(abs(Y)/n)
      
      % one sided (exclude 0 and nyquist)
      psd(2:end-1) = 2 * abs(psd(2:end-1));    
   
      mean(psd)
   
      % frequency (Hz) vector
      f = fs * (0:(n/2))/n;

      subplot(3, 2, c + 1), loglog(f, psd, 'LineStyle', line, 'LineWidth', linewidth, 'LineStyle', '-'); grid on;
      set(gca, 'LineWidth', 0.1);
      set(gca, 'fontsize', fontsize);
      
      if j == 4
        xlabel('Frequency (Hz)', 'fontsize', fontsize + 2);
      end  
      
      hold on;
      
   end	

else

   for j = 2:4

      subplot(3, 1, j-1), plot(data(:,1), data(:,j), 'LineWidth', linewidth, 'LineStyle', '-'); grid on;
      set(gca, 'LineWidth', 0.1);
      set(gca, 'fontsize', fontsize);      
      
      title(mytext{j-1});

   end	

end





