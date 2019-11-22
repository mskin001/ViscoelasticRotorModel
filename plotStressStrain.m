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
    
end

%% ------------------------------------------------------------------------
%  Plots
%  ------------------------------------------------------------------------

%------ Radial Displacement------------------------------------------------
if strcmp(plotWhat.radDis, 'yes')
  radDis = figure('Visible','on');
  
  hold on  
  plot(rArr*1000,uArr(1,:)*1000,'b-o')  
  
  xlabel('Radius [mm]')
  ylabel('Radial Displacement [mm]')

  set(gca, 'FontSize', 12)
  grid on
  fprintf('Radial Displacement Plot: Complete\n')
end

%------ Radial Stress -----------------------------------------------------
if strcmp(plotWhat.radStr, 'yes')
  radStr = figure('Visible','on');
  
  hold on

  plot(rArr*1000,sArr(3,:,1)*10^-6,'b-o')
  
  xlabel('Radius [mm]')
  ylabel('Radial Stress [MPa]')

  set(gca, 'FontSize', 12)
  grid on
  fprintf('Radial Sress Plot: Complete\n')
end

%------ Hoop Stress -------------------------------------------------------
if strcmp(plotWhat.hoopStr, 'yes')
  hoopStr = figure('Visible','on'); %#ok<*NASGU>

  hold on
  plot(rArr*1000,sArr(1,:,1)*10^-6,'b-o')
  xlabel('Radius [mm]')
  ylabel('Hoop Stress [MPa]')
  set(gca, 'FontSize', 12)
  grid on
  fprintf('Hoop Stress Plot: Complete\n')
end

%------ 1/R - Strength Ratio ----------------------------------------------


%% ------------------------------------------------------------------------
% Gifs
% -------------------------------------------------------------------------

%------ Radial Displacement -----------------------------------------------
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

%------ Radial Stress -----------------------------------------------------
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

%------ Hoop Stress -------------------------------------------------------
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

