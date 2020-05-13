function arrOut = laplacian(arrIn)
% Calculates laplacian for a square matrix using 3x3 convolution with edge
% wrapping
len = size(arrIn,1);
arrOut = zeros(len,len);
for i = 1:len
   for j = 1:len
      lInd = i-1;
      rInd = i+1;
      if lInd < 1
         lInd = len; 
      elseif rInd > len
         rInd = 1; 
      end
      uInd = j-1;
      dInd = j+1;
      if uInd < 1
         uInd = len; 
      elseif dInd > len
         dInd = 1; 
      end
      arrOut(i,j) = -1*arrIn(i,j)+0.25*(arrIn(lInd,j)+arrIn(rInd,j)+arrIn(i,uInd)+arrIn(i,dInd));
   end
end
end
