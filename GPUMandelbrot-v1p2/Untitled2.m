
        try
            d = gpuDevice();
        catch err
            error( 'mandelbrotViewer:NoGPU', 'mandelbrotViewer requires a GPU and none appear to be availble. Type "gpuDevice" for more information.' );
        end
        if ~d.DeviceSupported
            error( 'mandelbrotViewer:GPUNotSupported', 'The selected GPU is not supported. Type "gpuDevice" for more information.' );
        end