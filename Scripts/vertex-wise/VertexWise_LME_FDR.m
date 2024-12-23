clear
%% 1) set up
hemi = {'lh', 'rh'};
measure = {'volume','thickness','area'};

base_path = '/CUSTOM_THE_PATH/OpenAccess-UShape/'; %adapt base_path

%% Hemisphere loading
for i = 1:length(hemi)
    currentHemi = hemi{i};
    disp(['Current hemisphere: ' currentHemi]);

    % 1. Read the spherical surface and cortex label
    sphere = fs_read_surf(['/usr/local/freesurfer/subjects/fsaverage/surf/' currentHemi '.sphere']); %change path if needed
    cortex = fs_read_label(['/usr/local/freesurfer/subjects/fsaverage/label/' currentHemi '.cortex.label']); %change path if needed

    for j = 1:length(measure)
        currentMeasure = measure{j};
        disp(['Current Measure: ' currentMeasure]);

        % 2. Load cortical thickness/volume/area
        [Y, mri] = fs_read_Y([base_path 'Data/vertex-wise/mri_preproc_response/' currentHemi '.' currentMeasure '.EAEnMHdifT.10.allses.mgh']);
        
        % 3. Building the design matrix
        Qdec = fReadQdec([base_path 'Data/vertex-wise/qdec_table/qdec.table.EAEnMHdifT.Final.dat']);
        Qdec = rmQdecCol(Qdec,1); %(removes first column)
        sID = Qdec(2:end,1); %(grabs the subjects' IDs)
        Qdec = rmQdecCol(Qdec,1);  %(removes the subjects'ID column)
        M = Qdec2num(Qdec);  %(converts to a numeric matrix)
        [M,Y,ni] = sortData(M,1,Y,sID);  %(sorts the data)

        %M,Y,ni are now: the ordered assessment variables, the ordered data
        %and a vector with the number of repeated measures for each subject
        %respectively.

        X = [ones(length(M),1) M M(:,1).*M(:,3) M(:,1).*M(:,4) M(:,2).*M(:,3) M(:,2).*M(:,4)]; %we add the interaction terms M(:,3).*M(:,6) M(:,3).*M(:,8)

        %In our design matrix, covariables go at the end, after the
        %interactions

        covariables = X(:, 6:9); % Extract column 6 (etiv) and 9 (time interval 1-2)
        X(:, 6:9) = []; % Remove column 6 and 9 from X
        X = [X covariables]; % Append covariables as the last columns (10 to 13)

        %we avoid 0 values in the interaction
        startColumn = 6; endColumn = 9; % Column range

        for col = startColumn:endColumn
            nonZeroValues = X(:, col);
            meanValue = mean(nonZeroValues(~isnan(nonZeroValues) & nonZeroValues ~= 0));
            X(nonZeroValues ~= 0, col) = X(nonZeroValues ~= 0, col) - meanValue;
        end

        %% 2) Parameter estimation - spatiotemporal model

        stats = lme_mass_fit_vw(X, [1], Y, ni, cortex); %fitting the spatiotemporal model

        % Save stats separately for lh and rh and measure
        if strcmp(currentHemi, 'lh')
            if strcmp(currentMeasure, 'volume')
                lh_data.volume = stats;
            elseif strcmp(currentMeasure, 'area')
                lh_data.area = stats;
            elseif strcmp(currentMeasure, 'thickness')
                lh_data.thickness = stats;
            end
            lh_data.cortex = cortex;

        elseif strcmp(currentHemi, 'rh')
            if strcmp(currentMeasure, 'volume')
                rh_data.volume = stats;
            elseif strcmp(currentMeasure, 'area')
                rh_data.area = stats;
            elseif strcmp(currentMeasure, 'thickness')
                rh_data.thickness = stats;
            end
            rh_data.cortex = cortex;
        end

        disp('Iteration completed.');
        disp(' ');
    end
end

%% 3) Inference
c_name = {'T1GM', 'T1nGM', 'T2GM', 'T2nGM', 'T1nGM-GM', 'T2nGM-GM'};
c_num = 1:length(c_name); %change contrats number
%C1 = linear effect on gestational mothers.
%C2 = linear effect on non-gestational mothers.
%C3 = quadratic effect on gestational mothers.
%C4 = quadratic effect on nongestational mothers.
%C5 = differences between gest and non-gest mothers in the linear term.
%C6 = differences between gest and non-gest mothers in the quadratic term.

contrast_matrix = [...
    0 0 0 0 0 -1 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 -1 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 1 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 1 0 0 0 0;...
    0 0 0 0 0 1 -1 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 1 -1 0 0 0 0];

for c = 1:length(c_num)
    currentContrast = c_num(c);
    CM.C = contrast_matrix(currentContrast,:);
    allcontrasts = [];
    all_fvalues_lh = [];
    all_fvalues_rh = [];
    all_petasq_lh = [];
    all_petasq_rh = [];
    all_pvalues_lh =[];
    all_pvalues_rh =[];
    all_F_lhstats = cell(1, length(measure));
    all_F_rhstats = cell(1, length(measure));
    counter = 1;

    stats_lh = {lh_data.volume, lh_data.thickness, lh_data.area};
    stats_rh = {rh_data.volume, rh_data.thickness, rh_data.area};

    for m = 1:length(measure)
        currentMeasure = measure{m};

        F_lhstats = lme_mass_F(stats_lh{m},CM); %fstatistic
        F_rhstats = lme_mass_F(stats_rh{m},CM); %fstatistic

        P = [F_lhstats.pval(lh_data.cortex) F_rhstats.pval(rh_data.cortex)];

        fvalues_lh = F_lhstats.F;
        fvalues_rh = F_rhstats.F;

        petasq_lh = ((fvalues_lh).*F_lhstats.df(1,:))./(((fvalues_lh).*F_lhstats.df(1,:) + F_lhstats.df(2,:))); %unsigned
        petasq_rh = ((fvalues_rh).*F_rhstats.df(1,:))./(((fvalues_rh).*F_rhstats.df(1,:) + F_rhstats.df(2,:))); %unsigned

        % export significance maps
        mri1 = mri; mri1.volsz(4) = 1; %This is needed for fs_write_Y

        fs_write_fstats(F_lhstats,mri,[base_path '/Data/vertex-wise/lme_outputs/sigmaps/lh/sig_' c_name{currentContrast} '_' currentMeasure '_EAEnMHdifT_lh_FDR.mgz'],'sig');
        fs_write_fstats(F_rhstats,mri,[base_path '/Data/vertex-wise/lme_outputs/sigmaps/rh/sig_' c_name{currentContrast} '_' currentMeasure '_EAEnMHdifT_rh_FDR.mgz'],'sig');

        all_fvalues_lh = [all_fvalues_lh; fvalues_lh];
        all_fvalues_rh = [all_fvalues_rh; fvalues_rh];

        all_pvalues_lh = [all_pvalues_lh; F_lhstats.pval];
        all_pvalues_rh = [all_pvalues_rh; F_rhstats.pval];

        all_petasq_lh = [all_petasq_lh; petasq_lh];
        all_petasq_rh = [all_petasq_rh; petasq_rh];

        all_F_lhstats{counter} = F_lhstats;
        all_F_rhstats{counter} = F_rhstats;
        counter = counter + 1;

        allcontrasts = [allcontrasts P];
    end

    % allcontrasts_vec = allcontrasts(:)';
    dvtx = lme_mass_FDR(allcontrasts,0.05); %multiple corrections
    pcor = -log10(dvtx);

    all_petasq_lh_thr=all_petasq_lh;
    all_petasq_rh_thr=all_petasq_rh;

    all_petasq_lh_thr(all_pvalues_lh > dvtx) = NaN; %eveything below the threshold is now a 0
    all_petasq_rh_thr(all_pvalues_rh > dvtx) = NaN; %eveything below the threshold is now a 0

    % export pcor for each contrast
    writematrix(pcor,[base_path 'Data/vertex-wise/lme_outputs/pcor/pcor' c_name{currentContrast} '_EAEnMHdifT_FDR.txt']);
    
    percetnsig_bymeasure = [];

    % export PARTIAL ETA SQ
    for m1 = 1:length(measure)
        currentMeasure1 = measure{m1};
        F_lhstats1 = all_F_lhstats{m1};
        F_rhstats1 = all_F_rhstats{m1};
        F_lhstats1.pval = all_petasq_lh_thr(m1,:);
        F_rhstats1.pval = all_petasq_rh_thr(m1,:);
        mri1 = mri; mri1.volsz(4) = 1; %This is needed for fs_write_Y

        significantvoxels_lh = sum(~isnan(F_lhstats1.pval));
        significantvoxels_rh = sum(~isnan(F_rhstats1.pval));

        percentageofsign = ((significantvoxels_lh + significantvoxels_rh)/(length(cortex)*2))*100;
        percetnsig_bymeasure = [percetnsig_bymeasure percentageofsign];

        fs_write_fstats(F_lhstats1,mri1,[base_path 'Data/vertex-wise/lme_outputs/petamaps/lh/petasq_' c_name{currentContrast} '_' currentMeasure1 '_EAEnMHdifT_lh_FDR.mgz'],'pval');
        fs_write_fstats(F_rhstats1,mri1,[base_path 'Data/vertex-wise/lme_outputs/petamaps/rh/petasq_' c_name{currentContrast} '_' currentMeasure1 '_EAEnMHdifT_rh_FDR.mgz'],'pval');
        
        F_lhstats1.pval = all_petasq_lh(m1,:);
        F_rhstats1.pval = all_petasq_rh(m1,:);

        fs_write_Y(F_lhstats1.pval,mri1,[base_path 'Data/vertex-wise/lme_outputs/petamaps/lh/unthresholded_petasq_' c_name{currentContrast} '_' currentMeasure1 '_EAEnMHdifT_lh_FDR.mgz'])
        fs_write_Y(F_rhstats1.pval,mri1,[base_path 'Data/vertex-wise/lme_outputs/petamaps/rh/unthresholded_petasq_' c_name{currentContrast} '_' currentMeasure1 '_EAEnMHdifT_rh_FDR.mgz'])
        
    end
    percetnsig_bymeasure_table = array2table(percetnsig_bymeasure,"VariableNames",measure);
    writetable(percetnsig_bymeasure_table, [base_path 'Data/vertex-wise/lme_outputs/percent_signif/psignif_' c_name{currentContrast} '_EAEnMHdifT_FDR.txt'], 'Delimiter', '\t')
end