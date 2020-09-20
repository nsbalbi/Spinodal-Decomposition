function spinodal_decomp(D,gamma,options)
% SPINODAL_DECOMP Generates and records Cahn-Hilliard spinodal
% decomposition models using Euler's method.
% 
%    SPINODAL_DECOMP() generates a Cahn-Hilliad spindoal decomposition
%    model with sample coefficients. The result is saved as an AVI named
%    'spindoal_decomposition'.
%
%    SPINODAL_DECOMP(D,gamma) generates a Cahn-Hilliad spindoal 
%    decomposition model with diffusion coefficient D and square length of
%    transitional regions gamma. The result is saved as an AVI named 
%    'spindoal_decomposition'. See "Inputs" for details.
%    
%    SPINODAL_DECOMP(...,PARAM1,VAL1,PARAM2,VAL2,...) are additional
%    name-value pairs that can be used to change default function
%    parameters. See "Parameters" for details.
%
% Inputs
%    D: double. Diffusion coefficient. Default is 10
%    gamma: double. Sqaure length of transitional regions between domains.
%       Default is 5
% 
% Parameters
%    'dt': double. Time increment between iterations of Euler's method. 
%       (Default = 0.005)
%    'GridSize': int. Edge length for the square grid used for the
%       model. Note generation time increases exponentially with grid size.
%       (Default = 200)
%    'NumIterations': int. Total number of model iterations. Each iteration
%       represents a step of Euler's method. (Default = 10000)
%    'FrameSpacing': int. Number of iterations between captured frames, 
%       only relevant if capture mode is standard. (Default = 10)
%    'CaptureMode': char. Method of video cature. Possible inputs below.
%         'standard' - Constant num of iterations between frames. (Default)
%      'incremental' - Num iterations between frames increases over time
%    'ImgStyle': char. Method of frame generation. Possible inputs below.
%         'binary' - Concentrations are binarized to major species. (Default)
%           'true' - Concentrations are mapped to the colormap.
%    'Colormap': char. Colormap used for visualization. Supports all 
%       default MATLAB colormaps (Default = 'pink')
%    'FileName': char. Name of video (Default = 'spinodal_decomposition')
%       
% Examples
%    Example 1: Generate and record default model.
%       spinodal_decomp();
%    Example 2: Generate and record model with custom constants.
%       spinodal_decomp(20,3);
%    Example 3: Generate and record model with custom constants and capture
%    mode.
%       spinodal_decomp(10,10,'CaptureMode','incremental');
%    Example 4: Generate and record model with custom constants and
%    multiple custom parameters.
%       spinodal_decomp(10,10,...
%                      'CaptureMode','incremental',...
%                      'Colormap','jet',...
%                      'ImgStyle','true',...
%                      'NumIterations',25000);   

arguments
   D double = 10
   gamma double = 5
   options.dt double = 0.005
   options.GridSize double = 200
   options.NumIterations double = 10000
   options.FrameSpacing double = 10
   options.CaptureMode char = 'standard'
   options.ImgStyle char = 'binary'
   options.Colormap char = 'pink'
   options.FileName char = 'spinodal_decomposition'
end

% Random Starting Concentrations (-1 and 1 represent different species)
u = 2*randi(2,options.GridSize)-3;

writer = VideoWriter(options.FileName);
open(writer);
figure(1)
colormap(options.Colormap)
image(u,'CDataMapping','scaled');
frame = getframe(1);
writeVideo(writer,frame);

% Variables for incremental capture
count = 1;
frameStep = 1;

for i = 1:options.NumIterations
   u = iterate(u,D,gamma,options.dt);
   % Incremental video mode
   if (strcmp(options.CaptureMode,'incremental'))
       if (count == frameStep)
           if (strcmp(options.ImgStyle,'true'))
               image(u,'CDataMapping','scaled');
           else
               uMod = round((u+1)/2); % Binarizes image
               image(uMod,'CDataMapping','scaled');
           end
           frame = getframe(1);
           writeVideo(writer,frame);
           count = 0;
           frameStep = frameStep + 1;
       end
       count = count+1;
   % Standard video mode
   else
       if (mod(i,options.FrameSpacing) == 0)
           if (strcmp(options.ImgStyle,'true'))
               image(u,'CDataMapping','scaled');
           else
               uMod = round((u+1)/2);
               image(uMod,'CDataMapping','scaled');
           end
           frame = getframe(1);
           writeVideo(writer,frame);
       end
   end
end

close(writer);

fprintf('Done!\n');

% Forward Euler Method iteration of model
function uOut = iterate(u,D,gamma,dt) 
    % Calculates laplacian of concentration field
    uLaplace = laplacian(u);
    % Calculates chemical potentials
    uMu = u.^3 - u - gamma*uLaplace;
    % Laplacian of chemical potentials
    muLaplace = laplacian(uMu);
    % Cahn-Hilliard Equation
    duT = D*muLaplace;
    % Foreward Euler Method
    uOut = u + dt*duT;
end 

end
