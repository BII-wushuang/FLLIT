function [centre_of_mass,theta] = findCoM(locX,locY)
    numPoints = length(locX);
    centroid_x = sum(locX)/numPoints;
    centroid_y = sum(locY)/numPoints;
    centre_of_mass = [centroid_x centroid_y];

    % Defining the Rotation Matrices
    coeff = pca([locX locY]);
    if (coeff(1,1) > 0 && coeff(2,1) > 0)
        theta = acosd(coeff(1,1));
    elseif (coeff(1,1) < 0 && coeff(2,1) > 0)
        theta = 180 - asind(coeff(2,1));
    elseif (coeff(1,1) > 0 && coeff(2,1) < 0)
        theta = 360 - acosd(coeff(1,1));
    elseif (coeff(1,1) < 0 && coeff(2,1) < 0)
        theta = 180 + abs(asind(coeff(2,1)));    
    end

    theta = mod(theta - 90,360);
end
