function plotStressStrain()
%% -----------------------------------------------------------------------------
% Define global variables, arrays, and structures
% ------------------------------------------------------------------------------
global rArr sArr uArr vari
global plotWhat mat

%% -----------------------------------------------------------------------------
% Define rim origional centers and radii
% ------------------------------------------------------------------------------
numRims = length(plotWhat.rims);
oc = zeros(numRims,2);
or = plotWhat.rims;

%% -----------------------------------------------------------------------------
% Create a custom plot
% ------------------------------------------------------------------------------
% To create additional custom plots copy and paste this section to create
% as many custom plots as desired.

if strcmp(plotWhat.custom1, 'yes')
  nr = rArr ./ min(rArr);
  
  srRad = sArr(3,:,1) ./ mat.stren{1}(3);
  srHoop = sArr(1,:,1) ./ mat.stren{1}(1);
  
%   srRadInner = sArr(3,1:30,1) ./ mat.stren{1}(3);
%   srRadOuter = sArr(3,31:end,1) ./ mat.stren{2}(3);
%   srRad(1:30) = srRadInner;
%   srRad(31:60) = srRadOuter;
%   
%   srHoopInner = sArr(1,1:30,1) ./ mat.stren{1}(1);
%   srHoopOuter = sArr(1,31:end,1) ./ mat.stren{2}(1);
%   srHoop(1:30) = srHoopInner;
%   srHoop(31:60) = srHoopOuter;
  
  radStr = figure();
  haRadData = csvread('Ha99_GFRP_optimized_radialStress.csv');
  hahoopData = csvread('Ha99_GFRP_optimized_hoopStress.csv');
  hold on
  plot(haRadData(:,1),haRadData(:,2),'kv-', 'MarkerFaceColor', 'k')
  plot(hahoopData(:,1),hahoopData(:,2),'k^-', 'MarkerFaceColor', 'k')
  plot(nr,srRad,'b--s', 'LineWidth', 1)
  plot(nr,srHoop, 'r--o', 'LineWidth', 1)
  
%   axialStr = figure();
%   plot(rArr,sArr(2,:,1),'bo-');
% %   axis([0.5, 1, -1, 3]);
%   grid on
%   xlabel('r/r_{min}')
%   ylabel('Normalized Stress')
%   set(gca, 'FontSize', 12)
%   legend('Ha 1999 Radial', 'Ha 1999 Circumfrential','Model Radial','Model Circumfrential', 'Location', 'SouthEast')
%   
  fprintf('Custom plot 1: Complete\n')
end

%% -----------------------------------------------------------------------------
% Plot displacement
% ------------------------------------------------------------------------------
% Radial displacement gif
if strcmp(plotWhat.disGif,'yes')
  prog = waitbar(0,'0','Name','Radial Displacement .gif', 'CreateCancelBtn',...
    'setappdata(gcbf,''Canceling'',1)');
  setappdata(prog,'Canceling',0);
  
  disFig = figure('Visible', 'off');
  
  [r,t] = meshgrid(rArr, linspace(0,2*pi,length(rArr)));
  x = r .* cos(t);
  y = r .* sin(t);
  z = zeros(length(x),length(y));
  c = ones(length(x),length(y)) .* uArr(1,:) * 1000;

  surf(x,y,z,c, 'EdgeColor', 'none');
  colormap('jet');
  c = colorbar('Eastoutside', 'FontSize', 12);
  c.Label.String = 'Displacement [mm]';
  
  viscircles(oc, or, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);
  
  view(0,90)
  xlabel('Length [m]')
  ylabel('Length [m]')
  
  currT = 0;
  str = ['Time: ', num2str(currT)];
  loc = [0.613, .621, .3,.3];
  annotation('textbox',loc,'String',str, 'FitBoxToText', 'on'); 
  set(disFig, 'nextplot', 'replacechildren')
  
  F = getframe(disFig);
  [frames(:,:,1,1),map] = rgb2ind(F.cdata,256,'nodither');
  for b = 1:vari
    c = ones(length(x),length(y)) .* uArr(b,:) * 1000;

    surf(x,y,z,c, 'EdgeColor', 'none');
    c = colorbar('Eastoutside', 'FontSize', 12);
    c.Label.String = 'Displacement [mm]';
    
    if mod(b,plotWhat.interval) == 0
      currT = b / plotWhat.interval;
    end
    str = ['Time: ', num2str(currT)];
    annotation('textbox',loc,'String',str, 'FitBoxToText', 'on');
    
    viscircles(oc, or, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);
    
    view(0,90)
    xlabel('Length [m]')
    ylabel('Length [m]')
    set(disFig, 'nextplot', 'replacechildren')
    
    F = getframe(disFig);
    frames(:,:,1,b) = rgb2ind(F.cdata,map,'nodither');
    
    perc = (b / vari);
    waitbar(perc,prog,sprintf('%1.0f %%',perc*100)) 
    if getappdata(prog,'Canceling')
      delete(prog)
      return
    end
    
  end
  imwrite(frames, map, plotWhat.disGifName, 'DelayTime', plotWhat.delay, 'LoopCount', inf)
  
  delete(prog)
  fprintf('Radial Displacement gif: Complete\n')
end

% Radial displacement plot
if strcmp(plotWhat.radDis, 'yes')
  radDis = figure('Visible','on');
  
%   TzengInitial = csvread('Tzeng2001RadialDispInitial.csv');
%   TzengMid = csvread('Tzeng2001RadialDispMid.csv');
%   TzengInf = csvread('Tzeng2001RadialDispInf.csv');
  
  hold on
%   plot(TzengInitial(:,1),TzengInitial(:,2),'k^-','MarkerFaceColor','k')
%   plot(TzengMid(:,1),TzengMid(:,2),'kv-','MarkerFaceColor','k')
%   plot(TzengInf(:,1),TzengInf(:,2),'k>-','MarkerFaceColor','k')
%   
  plot(rArr,uArr(1,:),'b-o')
  plot(rArr,uArr(2,:),'g-o')
  plot(rArr,uArr(end,:),'r-o')  
%   plot(rArr*39.3701,uArr(1,:)*39.3701,'b-o')
%   plot(rArr*39.3701,uArr(2,:)*39.3701,'g-o')
%   plot(rArr*39.3701,uArr(end,:)*39.3701,'r-o')  
  
  plot(rArr,uArr(1,:)*1e-6,'b-o')
  plot(rArr,uArr(2,:)*1e-6,'r-s')
  plot(rArr,uArr(end,:)*1e-6,'k-^')
  
  xlabel('Radius [mm]')
  ylabel('Radial Displacement [mm]')
  legend('Initial', '1 Year', '10 Years')
%   legend('Tzeng Initial', 'Tzeng 10 years', 'Tzeng Infinite', 'Initial','10 Years', 'Infinite')
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Radial Displacement Plot: Complete\n')
end

%% -----------------------------------------------------------------------------
% Plot radial stress
% ------------------------------------------------------------------------------
% Radial Stress Gif
if strcmp(plotWhat.radGif, 'yes')
  prog = waitbar(0,'0','Name','Radial Stress .gif', 'CreateCancelBtn',...
    'setappdata(gcbf,''Canceling'',1)');
  setappdata(prog,'Canceling',0);
  
  radFig = figure('Visible','off');
    
  [r,t] = meshgrid(rArr, linspace(0,2*pi,length(rArr)));
  x = r .* cos(t);
  y = r .* sin(t);
  z = zeros(length(x),length(y));
  c = ones(length(x),length(y)) .* sArr(3,:,1) * 1e-6;

  surf(x,y,z,c, 'EdgeColor', 'none');
  colormap('jet');
  c = colorbar('Eastoutside', 'FontSize', 12);
  c.Label.String = 'Radial Stress [MPa]';
  
  viscircles(oc, or, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);
  
  view(0,90)
  xlabel('Length [m]')
  ylabel('Length [m]')
  
  currT = 0;
  str = ['Time: ', num2str(currT)];
  loc = [0.613, .621, .3,.3];
  annotation('textbox',loc,'String',str, 'FitBoxToText', 'on'); 
  set(radFig, 'nextplot', 'replacechildren')
  
  F = getframe(radFig);
  [frames(:,:,1,1),map] = rgb2ind(F.cdata,256,'nodither');
  for b = 1:vari
    c = ones(length(x),length(y)) .* sArr(3,:,b) * 1e-6;

    surf(x,y,z,c, 'EdgeColor', 'none');
    c = colorbar('Eastoutside', 'FontSize', 12);
    c.Label.String = 'Radial Stress [MPa]';
    
    if mod(b,plotWhat.interval) == 0
      currT = b / plotWhat.interval;
    end
    str = ['Time: ', num2str(currT)];
    annotation('textbox',loc,'String',str, 'FitBoxToText', 'on');
    
    viscircles(oc, or, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);
    
    view(0,90)
    xlabel('Length [m]')
    ylabel('Length [m]')
    set(radFig, 'nextplot', 'replacechildren')
    F = getframe(radFig);
    frames(:,:,1,b) = rgb2ind(F.cdata,map,'nodither');
    
    perc = (b / vari);
    waitbar(perc,prog,sprintf('%1.0f %%',perc*100)) 
    if getappdata(prog,'Canceling')
      delete(prog)
      return
    end
  end
  imwrite(frames, map, plotWhat.radialGifName, 'DelayTime', plotWhat.delay, 'LoopCount', inf)
  
  delete(prog)
  fprintf('Radial Sress gif: Complete\n')
end

% Radial Stress Plot
if strcmp(plotWhat.radStr, 'yes')
  radStr = figure('Visible','on');
%   
%   TzengInitial = csvread('Tzeng2001RadialStrInitial.csv');
%   TzengMid = csvread('Tzeng2001RadialStrMid.csv');
%   TzengInf = csvread('Tzeng2001RadialStrInf.csv');
  
  hold on
%   plot(TzengInitial(:,1),TzengInitial(:,2),'k^-','MarkerFaceColor','k')
%   plot(TzengMid(:,1),TzengMid(:,2),'kv-','MarkerFaceColor','k')
%   plot(TzengInf(:,1),TzengInf(:,2),'k>-','MarkerFaceColor','k')
  plot(rArr,sArr(3,:,1)*1e-6,'b-o')
  plot(rArr,sArr(3,:,2)*1e-6,'r-s')
  plot(rArr,sArr(3,:,end)*1e-6,'k-^')
%   plot(rArr*39.3701,sArr(3,:,2)*0.000145038,'g-o')
%   plot(rArr*39.3701,sArr(3,:,end)*0.000145038,'r-o')
  
  xlabel('Radius [m]')
  ylabel('Radial Stress [MPa]')
  legend('Initial', '1 Year', '10 Years')
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Radial Sress Plot: Complete\n')
end

%% -----------------------------------------------------------------------------
% Plot hoop stress
% ------------------------------------------------------------------------------
% Hoop stress gif
if strcmp(plotWhat.hoopGif, 'yes')
  prog = waitbar(0,'0','Name','Hoop Stress .gif', 'CreateCancelBtn',...
    'setappdata(gcbf,''Canceling'',1)');
  setappdata(prog,'Canceling',0);
  hoopFig = figure('Visible','off');  
    
  [r,t] = meshgrid(rArr, linspace(0,2*pi,length(rArr)));
  x = r .* cos(t);
  y = r .* sin(t);
  z = zeros(length(x),length(y));
  c = ones(length(x),length(y)) .* sArr(1,:,1) * 1e-6;

  surf(x,y,z,c, 'EdgeColor', 'none');
  colormap('jet');
  c = colorbar('Eastoutside', 'FontSize', 12);
  c.Label.String = 'Hoop Stress [MPa]';
  
  viscircles(oc, or, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);
  
  view(0,90)
  xlabel('Length [m]')
  ylabel('Length [m]')
  
  currT = 0;
  str = ['Time: ', num2str(currT)];
  loc = [0.613, .621, .3,.3];
  annotation('textbox',loc,'String',str, 'FitBoxToText', 'on'); 
  set(hoopFig, 'nextplot', 'replacechildren')
  
  F = getframe(hoopFig);
  [frames(:,:,1,1),map] = rgb2ind(F.cdata,256,'nodither');
  for b = 1:vari
    c = ones(length(x),length(y)) .* sArr(1,:,b) * 1e-6;

    surf(x,y,z,c, 'EdgeColor', 'none');
    c = colorbar('Eastoutside', 'FontSize', 12);
    c.Label.String = 'Hoop Stress [MPa]';
    
    if mod(b,plotWhat.interval) == 0
      currT = b / plotWhat.interval;
    end
    str = ['Time: ', num2str(currT)];
    annotation('textbox',loc,'String',str, 'FitBoxToText', 'on');
    
    viscircles(oc, or, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);
    
    view(0,90)
    xlabel('Length [m]')
    ylabel('Length [m]')
    set(hoopFig, 'nextplot', 'replacechildren')
    F = getframe(hoopFig);
    frames(:,:,1,b) = rgb2ind(F.cdata,map,'nodither');
    
    perc = (b / vari);
    waitbar(perc,prog,sprintf('%1.0f %%',perc*100)) 
    if getappdata(prog,'Canceling')
      delete(prog)
      return
    end
  end
  imwrite(frames, map, plotWhat.hoopGifName, 'DelayTime', plotWhat.delay, 'LoopCount', inf)
  
  delete(prog)
  fprintf('Hoop Stress gif: Complete\n')
end

% Hoop Stress plot
if strcmp(plotWhat.hoopStr, 'yes')
  hoopStr = figure('Visible','on'); %#ok<*NASGU>
%   TzengInitial = csvread('Tzeng2001HoopInitial.csv');
%   TzengMid = csvread('Tzeng2001HoopMid.csv');
%   TzengInf = csvread('Tzeng2001HoopInf.csv');
  hold on
%   plot(TzengInitial(:,1),TzengInitial(:,2),'k^-','MarkerFaceColor','k')
%   plot(TzengMid(:,1),TzengMid(:,2),'kv-','MarkerFaceColor','k')
%   plot(TzengInf(:,1),TzengInf(:,2),'k>-','MarkerFaceColor','k')
  plot(rArr,sArr(1,:,1)*1e-6,'b-o')
  plot(rArr,sArr(1,:,2)*1e-6,'r-s')
  plot(rArr,sArr(1,:,end)*1e-6,'k-^')
%   plot(rArr*39.3701,sArr(1,:,2)*0.000145038,'g-o')
%   plot(rArr*39.3701,sArr(1,:,end)*0.000145038,'r-*')
  xlabel('Radius [m]')
  ylabel('Circumferential Stress [MPa]')
  legend('Initial', '1 Year', '10 Years')
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Circumferential Stress Plot: Complete\n')
end

%% Axial Stress Plots
if strcmp(plotWhat.axialStr, 'yes')
  axialStr = figure('Visible','on'); %#ok<*NASGU>

  hold on
  
%   HaAxial = csvread('Ha2001_axialStress.csv');
%   plot(HaAxial(:,1),HaAxial(:,2),'kv-','MarkerFaceColor','k')
  plot(rArr,sArr(2,:,1)*1e-6,'b-o')
  plot(rArr,sArr(2,:,2)*1e-6,'r-s')
  plot(rArr,sArr(2,:,3)*1e-6,'k-^')
%   try
%     plot(rArr,sArr(2,:,2),'g-o')
%     plot(rArr,sArr(2,:,end),'r-*')
%   catch
    
  end
  xlabel('Radius [m]')
  ylabel('Axial Stress [MPa]')
  legend('Initial', '1 Year', '10 Years')
%   legend('Tzeng Initial', 'Tzeng 10 years', 'Tzeng Infinite', 'Initial','10 Years', 'Infinite')
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Axial Stress Plot: Complete\n')
end
