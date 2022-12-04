classdef ptyLab < matlab.mixin.Copyable
    % note: matlab.mixin.Copyable allows to make a hard copy of an object
    
    % required properties
    properties
        
        % operation
        operationMode       % 'FPM' or 'CPM': defines operation mode (FP/CP: Fourier/conventional ptychography)
        
        % physical properties
        wavelength          % (operational) wavelength, scalar quantity
        spectralDensity     % spectral density S = S(lambda), vectorial quantity
        % note: spectral density is required for polychromatic operation.
        % In this case, lambda is still scalar and determines the lateral
        % pixel size of the meshgridgrid that all other wavelengths are 
        % interpolated onto. 
        
        % (entrance) pupil / probe sampling
        dxp                 % pixel size (entrance pupil plane)
        Np                  % number of pixel (entrance pupil plane)
        xp                  % 1D coordinates (entrance pupil plane)
        Xp                  % 2D meshgrid in x-direction (entrance pupil plane)
        Yp                  % 2D meshgrid in y-direction (entrance pupil plane)
        Lp                  % field of view (entrance pupil plane)
        zp                  % distance to next plane of interest
        
        % object sampling
        dxo                 % pixel size (object plane)
        No                  % number of pixel (object plane)
        xo                  % 1D coordinates (object plane)
        Xo                  % 2D meshgrid in x-direction (object plane)
        Yo                  % 2D meshgrid in y-direction (object plane)
        Lo                  % field of view (object plane)
        zo                  % distance to next plane of interest
        
        % detector sampling
        dxd                 % pixel size (detector plane)
        Nd                  % number of pixel (detector plane)
        xd                  % 1D coordinates (detector plane)
        Xd                  % 2D meshgrid in x-direction (detector plane)
        Yd                  % 2D meshgrid in y-direction (detector plane)
        Ld                  % field of view (detector plane)
        
        % measured intensities
        ptychogram      % intensities [Nd, Nd, numPos]
        numFrames       % number of measurements (positions (CPM) / LED tilts (FPM))
        background      % background
        binningFactor   % binning factor that was applied to raw data
        
        % measured positions
        positions0      % initial positions in pixel units (real space for CPM, Fourier space for FPM)
        positions       % estimated positions in pixel units (real space for CPM, Fourier space for FPM)
        % note: Positions are given in row-column order and refer to the
        % pixel in the upper left corner of the respective data matrix;
        % -1st example: suppose the 2nd row of positions0 is [3, 4] and the
        % operation mode is 'CPM'. That implies that the second intensity
        % in the spectrogram updates an object patch that has
        % its left uppper corner pixel at the pixel coordinates [3, 4]
        % -2nd example: suppose the 2nd row of positions0 is [3, 4] and the
        % operation mode is 'FPM'. That implies that the second intensity
        % in the spectrogram is updates a patch which has pixel coordinates
        % [3,4] in the high-resolution Fourier transform
        
    end
    
    % optional properties
    properties
        params          % contains algorithmic reconstruction params
        monitor         % contains axes and colormaps
        % note: see initialParams function for a list of parameters
        
        propagator      % contains all information required to propagate from object/object spectrum to detector plane
        export          % contains relevant paths and export options
        entrancePupilDiameter   % entrance pupil diameter
    end
    
    % estimated quantities
    properties
        
        % reconstructed object
        object      % [numPixels, numPixels, nosm]
        
        % reconstructed illumination
        probe       % [numPixels, numPixels, npsm]
        
    end
    
    % constructor
    methods
        
        function obj = ptyLab( varargin )
            disp('create object')
            
            % parse optional inputs
            % note that by default 'CPM' mode is used
            p = inputParser;
            p.addParameter('operationMode','CPM')
            p.parse( varargin{:} )
            
            obj.operationMode = p.Results.operationMode;
        end
        
    end
    
end

