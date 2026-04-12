clc; clear; close all;

%% ===================== FILE LIST =====================
files = { ...
    'input.dat', ...
    'diano_8modes.dat', ...
    'diano_16modes.dat', ...
    'diano_32modes.dat' ...
    };

labels = { ...
    'Input', ...
    '8 modes', ...
    '16 modes', ...
    '32 modes' ...
    };

nFiles = numel(files);

%% ===================== READ DATA =====================
fields = cell(nFiles,1);

for i = 1:nFiles
    fields{i} = readTecplotScalarField(files{i});
    fprintf('Loaded %s of size %d x %d\n', files{i}, size(fields{i},1), size(fields{i},2));
end

%% ===================== FIELD PLOTS =====================
figure('Name','Field Comparison','Color','w');
tiledlayout(1, nFiles, 'Padding','compact', 'TileSpacing','compact');

for i = 1:nFiles
    nexttile;
    imagesc(fields{i});
    axis equal tight;
    set(gca,'YDir','normal');
    colorbar;
    title(labels{i});
end

%% ===================== ENERGY SPECTRUM =====================
kList = cell(nFiles,1);
EList = cell(nFiles,1);

for i = 1:nFiles
    [kList{i}, EList{i}] = radialSpectrum2D(fields{i});
end

%% ===================== RAW SPECTRUM =====================
figure('Name','Energy Spectrum Comparison','Color','w');
hold on;

for i = 1:nFiles
    k = kList{i};
    E = EList{i};

    valid = (k > 0) & (E > 0);
    loglog(k(valid), E(valid), 'LineWidth', 1.8);
end

set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('Wavenumber k');
ylabel('E(k)');
legend(labels, 'Location','best');
title('Spectrum Comparison: Different Number of Modes');

%% ===================== NORMALIZED SPECTRUM =====================
figure('Name','Normalized Spectrum','Color','w');
hold on;

for i = 1:nFiles
    k = kList{i};
    E = EList{i};

    E = E / sum(E + eps);

    valid = (k > 0) & (E > 0);
    loglog(k(valid), E(valid), 'LineWidth', 1.8);
end

set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('Wavenumber k');
ylabel('Normalized E(k)');
legend(labels, 'Location','best');
title('Normalized Spectrum: Different Number of Modes');

%% ===================== -2/3 REFERENCE LINE =====================
figure('Name','Spectrum with Reference Slope','Color','w');
hold on;

for i = 1:nFiles
    k = kList{i};
    E = EList{i};

    valid = (k > 0) & (E > 0);
    loglog(k(valid), E(valid), 'LineWidth', 1.8);
end

k_ref = 10;
k_ref = max(k_ref, min(kList{1}));
k_ref = min(k_ref, max(kList{1}));

E_ref = interp1(kList{1}, EList{1}, k_ref, 'linear', 'extrap');
k_line = [k_ref, max(kList{1})];
E_line = E_ref * (k_line / k_ref).^(-3);

loglog(k_line, E_line, 'k--', 'LineWidth', 2.2);

set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('Wavenumber k');
ylabel('E(k)');
legend([labels, {'-2/3 reference'}], 'Location','best');
title('Spectrum with Reference -2/3 Slope');

%% ===================== WRITE TECPLOT FILES =====================
writeTecplotSpectrumFile('modes_spectrum.dat', labels, kList, EList, false);
writeTecplotSpectrumFile('modes_spectrum_normalized.dat', labels, kList, EList, true);

fprintf('Tecplot spectrum files written successfully.\n');

%% ============================================================
%% ===================== LOCAL FUNCTIONS ======================
%% ============================================================

function field = readTecplotScalarField(filename)
    raw = fileread(filename);
    lines = splitlines(raw);

    numericRows = {};
    for i = 1:numel(lines)
        line = strtrim(lines{i});
        if isempty(line)
            continue;
        end

        if ~isempty(regexp(line, '^[\s\+\-\d\.Ee,]+$', 'once'))
            numericRows{end+1,1} = line; %#ok<AGROW>
        end
    end

    if isempty(numericRows)
        error('No numeric data found in %s', filename);
    end

    data = [];
    for i = 1:numel(numericRows)
        vals = sscanf(strrep(numericRows{i}, ',', ' '), '%f');
        if ~isempty(vals)
            data = [data; vals(:).']; %#ok<AGROW>
        end
    end

    if isempty(data)
        error('Could not parse numeric data in %s', filename);
    end

    nCols = size(data,2);

    if nCols == 1
        scalar = data(:,1);
        N = round(sqrt(numel(scalar)));
        if N*N ~= numel(scalar)
            error('Scalar data in %s is not square. Please reshape manually.', filename);
        end
        field = reshape(scalar, [N, N]);

    elseif nCols >= 3
        x = data(:,1);
        y = data(:,2);
        scalar = data(:,end);

        ux = unique(x);
        uy = unique(y);

        nx = numel(ux);
        ny = numel(uy);

        if nx * ny ~= numel(scalar)
            N = round(sqrt(numel(scalar)));
            if N*N ~= numel(scalar)
                error('Could not infer grid dimensions from %s. Check file format.', filename);
            end
            field = reshape(scalar, [N, N]);
        else
            field = nan(ny, nx);
            for j = 1:numel(scalar)
                ix = find(ux == x(j), 1);
                iy = find(uy == y(j), 1);
                field(iy, ix) = scalar(j);
            end

            if any(isnan(field(:)))
                error('Some grid points could not be mapped in %s.', filename);
            end
        end
    else
        error('Unsupported number of columns in %s.', filename);
    end
end

function [kVals, E_k] = radialSpectrum2D(field)
    field = field - mean(field(:));

    F = fftshift(fft2(field));
    E2 = abs(F).^2;

    [Ny, Nx] = size(field);

    kx = (-floor(Nx/2)):(ceil(Nx/2)-1);
    ky = (-floor(Ny/2)):(ceil(Ny/2)-1);
    [KX, KY] = meshgrid(kx, ky);
    K = sqrt(KX.^2 + KY.^2);

    kBins = round(K);
    kMax = max(kBins(:));

    E_k = zeros(kMax+1,1);
    count = zeros(kMax+1,1);

    for kk = 0:kMax
        mask = (kBins == kk);
        count(kk+1) = sum(mask(:));
        if count(kk+1) > 0
            E_k(kk+1) = sum(E2(mask)) / count(kk+1);
        end
    end

    kVals = (0:kMax).';

    valid = (kVals > 0) & (E_k > 0);
    kVals = kVals(valid);
    E_k = E_k(valid);
end

function writeTecplotSpectrumFile(outFile, labels, kList, EList, normalizeFlag)
    fid = fopen(outFile, 'w');
    if fid == -1
        error('Could not open %s for writing.', outFile);
    end

    fprintf(fid, 'TITLE = "Radial Energy Spectrum Comparison"\n');
    fprintf(fid, 'VARIABLES = "k" "E(k)"\n');

    nFiles = numel(labels);

    for i = 1:nFiles
        k = kList{i};
        E = EList{i};

        valid = (k > 0) & (E > 0);
        k = k(valid);
        E = E(valid);

        if normalizeFlag
            E = E / sum(E + eps);
        end

        nPts = numel(k);
        fprintf(fid, 'ZONE T="%s", I=%d, F=POINT\n', labels{i}, nPts);

        for j = 1:nPts
            fprintf(fid, '%.8e %.8e\n', k(j), E(j));
        end

        fprintf(fid, '\n');
    end

    fclose(fid);
end