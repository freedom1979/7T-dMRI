function NB = region_nn(p, mask)
NN(:,1) = [p(1)+1
    p(1)+1
    p(1)+1
    p(1)
    p(1)
    p(1)
    p(1)-1
    p(1)-1
    p(1)-1
    p(1)+1
    p(1)+1
    p(1)+1
    p(1)
    p(1)
    p(1)-1
    p(1)-1
    p(1)-1
    p(1)+1
    p(1)+1
    p(1)+1
    p(1)
    p(1)
    p(1)
    p(1)-1
    p(1)-1
    p(1)-1];

NN(:,2) =[p(2)-1
    p(2)
    p(2)+1
    p(2)-1
    p(2)
    p(2)+1
    p(2)-1
    p(2)
    p(2)+1
    p(2)-1
    p(2)
    p(2)+1
    p(2)-1
    p(2)+1
    p(2)-1
    p(2)
    p(2)+1
    p(2)-1
    p(2)
    p(2)+1
    p(2)-1
    p(2)
    p(2)+1
    p(2)-1
    p(2)
    p(2)+1];

NN(:,3) =[p(3)+1
    p(3)+1
    p(3)+1
    p(3)+1
    p(3)+1
    p(3)+1
    p(3)+1
    p(3)+1
    p(3)+1
    p(3)
    p(3)
    p(3)
    p(3)
    p(3)
    p(3)
    p(3)
    p(3)
    p(3)-1
    p(3)-1
    p(3)-1
    p(3)-1
    p(3)-1
    p(3)-1
    p(3)-1
    p(3)-1
    p(3)-1];

    NB = [];
    for i = 1:26
        if mask(NN(i,1),NN(i,2),NN(i,3)) == 1
            NB = [NB;[NN(i,1) NN(i,2) NN(i,3)]];
        end
    end
end