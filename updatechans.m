function chanlist = updatechans(app)
    chanvec = zeros(32,1);
    for chan = 1:32
        chanfield = ['CheckBox_' num2str(chan)];
        if app.(chanfield).Value
            chanvec(chan) = 1;
        end
    end
    chanlist = find(chanvec);
end
