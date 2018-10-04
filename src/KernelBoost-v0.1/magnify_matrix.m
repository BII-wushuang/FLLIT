function [m_mat] = magnify_matrix(in_mat, factor)

m_mat = zeros(size(in_mat, 1)*factor, size(in_mat, 2)*factor);

for r = 1 : size(in_mat, 1)
    for c = 1 : size(in_mat, 2)
        m_mat((r-1)*factor + 1 : r*factor, (c-1)*factor + 1 : c*factor) = in_mat(r, c);
    end
end 

end
