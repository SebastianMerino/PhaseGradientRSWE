function metrics = get_metrics_v2(img,roi_inc,roi_back,method,freq, gt_inc, gt_back)
    metrics.mean_inc = mean(img(roi_inc));
    metrics.mean_back =  mean(img(roi_back));
    metrics.std_inc = std(img(roi_inc));
    metrics.std_back = std(img(roi_back));
    metrics.cnr = abs(metrics.mean_inc - metrics.mean_back)/...
        (sqrt(metrics.std_inc.^2 + metrics.std_back.^2) + 1e-5);
    metrics.method = method;
    metrics.freq = freq;
    metrics.bias_inc = (metrics.mean_inc - gt_inc)./ gt_inc * 100;
    metrics.bias_back = (metrics.mean_back - gt_back)./ gt_back * 100;
end