function logCount = computeMandelbrotCUDAKernel2( xlim, numx, ylim, numy, maxIters )
        % Use pre-existing CUDA/C++ code.
        % The final way in which MATLAB can use the GPU is by calling some
        % hand-written CUDA code. The "CUDAKernel" interface allows the
        % function to be specified along with the number of threads and
        % blocks to use. This requires some knowledge of how GPUs work, but
        % does allow you to easily use existing CUDA kernels with MATLAB
        % data.
        data = createData();
        gui = createGUI();
        % Create the input arrays
        escapeRadius = 20;
        x = parallel.gpu.GPUArray.linspace( xlim(1),  xlim(2), numx );
        y = parallel.gpu.GPUArray.linspace( ylim(1),  ylim(2), numy );
        [x0,y0] = meshgrid(x, y);
        
        % Make sure we have sufficient blocks to cover the whole array
        numElements = numel( x0 );
        data.Kernel.ThreadBlockSize = [data.Kernel.MaxThreadsPerBlock,1,1];
        data.Kernel.GridSize = [ceil(numElements/data.Kernel.MaxThreadsPerBlock),1];
        
        % Call the kernel
        logCount = parallel.gpu.GPUArray.zeros( size( x0 ) );
        logCount = feval( data.Kernel, logCount, ...
            x0, y0, ...
            escapeRadius, maxIters, numElements );
        logCount = gather( logCount );
    end % computeMandelbrotCUDAKernel


