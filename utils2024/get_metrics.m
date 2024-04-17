function metrics = get_metrics(img,roi_inc,roi_back,method,freq)
    metrics.mean_inc = mean(img(roi_inc));
    metrics.mean_back =  mean(img(roi_back));
    metrics.std_inc = std(img(roi_inc));
    metrics.std_back = std(img(roi_back));
    metrics.cnr = abs(metrics.mean_inc - metrics.mean_back)/...
        sqrt(metrics.std_inc.^2 + metrics.std_back.^2);
    metrics.method = method;
    metrics.freq = freq;
end