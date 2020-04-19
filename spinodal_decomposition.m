function spinodal_decomposition(D,gamma,dt,gridSize,numIterations,captureMode)
% Generates Cahn-Hilliard Phase Separation Model, Exports Video
% D: diffusion coefficient. Standard values is 10
% gamma: sqrt(gamma) is the domain transition length. Standard value of 5
% dt: time increment between iterations. Standard value of 0.005
% gridSize: grid size. Standard values are 100, 200, 400
% numIterations: number of iterations. Standard values between 5000, 25000
% captureMode: "standard" - 25 iterations between frames
%              "incremental" - 1 iteration between frames, then 2, then 3...

len = gridSize;
% Random Starting Concentrations
u = 2*randi(2,len)-3;

colormap pink

writer = VideoWriter('spinodal_decomposition.avi');
open(writer);

figure(1)
image(u,'CDataMapping','scaled');
frame = getframe(1);
writeVideo(writer,frame);

% Variables for incremental capture
count = 1;
frameStep = 1;

for i = 1:numIterations
   u = iterate(u,D,gamma,dt);
   % Incremental video mode
   if (captureMode == "incremental")
       if (count == frameStep)
           uMod = round((u+1)/2);
           image(uMod,'CDataMapping','scaled');
           frame = getframe(1);
           writeVideo(writer,frame);
           count = 0;
           frameStep = frameStep + 1;
       end
       count = count+1;
   % Standard video mode
   else
       if (mod(i,10) == 0)
       uMod = round((u+1)/2);
       image(uMod,'CDataMapping','scaled');
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