function y = mrf(m,n,count)

field = 10 * rand(m,n);

for t = 1:count
    row = randi([1, m], 1,1);
    column = randi([1,n], 1,1);
    new = 10 * rand;
    if (row > 1 && row < m) && (column >1 && column < n)
        new = (field(row-1,column-1)+field(row-1,column)+field(row-1,column+1)+field(row,column-1)+field(row,column+1)+field(row+1,column-1)+field(row+1,column)+field(row+1,column+1))/8;
        p1 = exp(abs((field(row-1,column-1)+field(row-1,column)+field(row-1,column+1)+field(row,column-1)+field(row,column+1)+field(row+1,column-1)+field(row+1,column)+field(row+1,column+1))/8-field(row,column)));
        p2 = exp(abs((field(row-1,column-1)+field(row-1,column)+field(row-1,column+1)+field(row,column-1)+field(row,column+1)+field(row+1,column-1)+field(row+1,column)+field(row+1,column+1))/8-new));
    end
    if row == 1 && (column > 1 && column < n)
        new = (field(row,column-1)+field(row,column+1)+field(row+1,column-1)+field(row+1,column)+field(row+1,column+1))/5;
        p1 = exp(abs((field(row,column-1)+field(row,column+1)+field(row+1,column-1)+field(row+1,column)+field(row+1,column+1))/5-field(row,column)));
        p2 = exp(abs((field(row,column-1)+field(row,column+1)+field(row+1,column-1)+field(row+1,column)+field(row+1,column+1))/5-new));
    end
    if row == m && (column > 1 && column < n)
        new = (field(row-1,column-1)+field(row-1,column)+field(row-1,column+1)+field(row,column-1)+field(row,column+1))/5;
        p1 = exp(abs((field(row-1,column-1)+field(row-1,column)+field(row-1,column+1)+field(row,column-1)+field(row,column+1))/5-field(row,column)));
        p2 = exp(abs((field(row-1,column-1)+field(row-1,column)+field(row-1,column+1)+field(row,column-1)+field(row,column+1))/5-new));
    end
    if (row > 1 && row < m) && column == 1
        new = (field(row-1,column)+field(row-1,column+1)+field(row,column+1)+field(row+1,column)+field(row+1,column+1))/5;
        p1 = exp(abs((field(row-1,column)+field(row-1,column+1)+field(row,column+1)+field(row+1,column)+field(row+1,column+1))/5-field(row,column)));
        p2 = exp(abs((field(row-1,column)+field(row-1,column+1)+field(row,column+1)+field(row+1,column)+field(row+1,column+1))/5-new));
    end
    if (row > 1 && row < m) && column == n
        new = (field(row-1,column-1)+field(row-1,column)+field(row,column-1)+field(row+1,column-1)+field(row+1,column))/5;
        p1 = exp(abs((field(row-1,column-1)+field(row-1,column)+field(row,column-1)+field(row+1,column-1)+field(row+1,column))/5-field(row,column)));
        p2 = exp(abs((field(row-1,column-1)+field(row-1,column)+field(row,column-1)+field(row+1,column-1)+field(row+1,column))/5-new));
    end
    if row == 1 && column == 1
        new = (field(row,column+1)+field(row+1,column)+field(row+1,column+1))/3;
        p1 = exp(abs((field(row,column+1)+field(row+1,column)+field(row+1,column+1))/3-field(row,column)));
        p2 = exp(abs((field(row,column+1)+field(row+1,column)+field(row+1,column+1))/3-new));
    end
    if row == 1 && column == n
        new = (field(row,column-1)+field(row+1,column-1)+field(row+1,column))/3;
        p1 = exp(abs((field(row,column-1)+field(row+1,column-1)+field(row+1,column))/3-field(row,column)));
        p2 = exp(abs((field(row,column-1)+field(row+1,column-1)+field(row+1,column))/3-new));
    end
    if row == m && column == 1
        new = (field(row-1,column)+field(row-1,column+1)+field(row,column+1))/3;
        p1 = exp(abs((field(row-1,column)+field(row-1,column+1)+field(row,column+1))/3-field(row,column)));
        p2 = exp(abs((field(row-1,column)+field(row-1,column+1)+field(row,column+1))/3-new));
    end
    if row == m && column == n
        new = (field(row-1,column-1)+field(row-1,column)+field(row,column-1))/3;
        p1 = exp(abs((field(row-1,column-1)+field(row-1,column)+field(row,column-1))/3-field(row,column)));
        p2 = exp(abs((field(row-1,column-1)+field(row-1,column)+field(row,column-1))/3-new));
    end
    
    r = p1/p2;
    if r >= 1
        field(row,column) = new;
%     else
%        u = rand;
%        if u <= r
%            field(row, column) = new;
%        end
    end
end
y = field;