function mvec = map(vec, lb, ub)
    
mvec = vec;
br = ub - lb;

mvec(mvec > ub) = mod(mvec(mvec > ub), br);
mvec(mvec < lb) = mod(mvec(mvec < lb), -br);
mvec(mvec > ub) = mvec(mvec > ub) - br;
mvec(mvec < lb) = br + mvec(mvec < lb);

end
