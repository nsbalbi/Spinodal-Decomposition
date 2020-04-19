function arrOut = laplacian(arrIn)
% Calculates laplacian for a square matrix using 3x3 convolution with edge
% wrapping
len = size(arrIn,1);
arrOut = zeros(len,len);
for i = 1:len
   for j = 1:len
      lInd = i-1;
      if lInd < 1
         lInd = len; 
      end
      rInd = i+1;
      if rInd > len
         rInd = 1; 
      end
      uInd = j-1;
      if uInd < 1
         uInd = len; 
      end
      dInd = j+1;
      if dInd > len
         dInd = 1; 
      end
      arrOut(i,j) = -1*arrIn(i,j)+0.2*(arrIn(lInd,j)+arrIn(rInd,j)+arrIn(i,uInd)+arrIn(i,dInd))+0.05*(arrIn(lInd,uInd)+arrIn(rInd,uInd)+arrIn(lInd,dInd)+arrIn(rInd,dInd));
   end
end
end