function xx = ABinv(x) % f=AB^{-1}
    global K B
    xx = K*(B(x));
end
