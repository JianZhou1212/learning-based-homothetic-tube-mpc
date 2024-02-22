function samples = polytope_sample(poly, N_sam)
    V = poly.V;
    ver_x = V(:, 1);
    ver_y = V(:, 2);
    min_x = min(ver_x);
    max_x = max(ver_x);
    min_y = min(ver_y);
    max_y = max(ver_y);
    samples = zeros(2, N_sam);
    i = 1;
    while i <= N_sam
        x = (max_x-min_x).*rand(1) + min_x;
        y = (max_y-min_y).*rand(1) + min_y;
        if poly.A*[x; y] <= poly.b
            samples(:, i) = [x; y];
            i = i + 1;
        else
            continue;
        end
    end
end