function Img = func_DrawLine(Img, X0, Y0, X1, Y1, nG)
% Connect two pixels in an image with the desired graylevel
%
% Command line
% ------------
% result = func_DrawLine(Img, X1, Y1, X2, Y2)
% input:    Img : the original image.
%           (X1, Y1), (X2, Y2) : points to connect.
%           nG : the gray level of the line.
% output:   result
% result = func_DrawLine(zeros(5, 10), 2, 1, 5, 10, 1)
% result =
%      0     0     0     0     0     0     0     0     0     0
%      1     1     1     0     0     0     0     0     0     0
%      0     0     0     1     1     1     0     0     0     0
%      0     0     0     0     0     0     1     1     1     0
%      0     0     0     0     0     0     0     0     0     1

Img(X0, Y0) = nG;
Img(X1, Y1) = nG;
if abs(X1 - X0) <= abs(Y1 - Y0)
   if Y1 < Y0
      k = X1; X1 = X0; X0 = k;
      k = Y1; Y1 = Y0; Y0 = k;
   end
   if (X1 >= X0) & (Y1 >= Y0)
      dy = Y1-Y0; dx = X1-X0;
      p = 2*dx; n = 2*dy - 2*dx; tn = dy;
      while (Y0 < Y1)
         if tn >= 0
            tn = tn - p;
         else
            tn = tn + n; X0 = X0 + 1;
         end
         Y0 = Y0 + 1; Img(X0, Y0) = nG;
      end
   else
      dy = Y1 - Y0; dx = X1 - X0;
      p = -2*dx; n = 2*dy + 2*dx; tn = dy;
      while (Y0 <= Y1)
         if tn >= 0
            tn = tn - p;
         else
            tn = tn + n; X0 = X0 - 1;
         end
         Y0 = Y0 + 1; Img(X0, Y0) = nG;
      end
   end
else if X1 < X0
      k = X1; X1 = X0; X0 = k;
      k = Y1; Y1 = Y0; Y0 = k;
   end
   if (X1 >= X0) & (Y1 >= Y0)
      dy = Y1 - Y0; dx = X1 - X0;
      p = 2*dy; n = 2*dx-2*dy; tn = dx;
      while (X0 < X1)
         if tn >= 0
            tn = tn - p;
         else
            tn = tn + n; Y0 = Y0 + 1;
         end
         X0 = X0 + 1; Img(X0, Y0) = nG;
      end
   else
      dy = Y1 - Y0; dx = X1 - X0;
      p = -2*dy; n = 2*dy + 2*dx; tn = dx;
      while (X0 < X1)
         if tn >= 0
            tn = tn - p;
         else
            tn = tn + n; Y0 = Y0 - 1;
         end
         X0 = X0 + 1; Img(X0, Y0) = nG;
      end
   end
end