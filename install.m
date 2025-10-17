function install(option)
    %INSTALL  Setup paths for the MVector toolbox (src & test).
    %
    % Usage
    %   install                % add paths for current session
    %   install('save')        % add paths and attempt to save permanently
    %   install('reset')       % remove previously-added src/test paths
    %
    % Layout (expected)
    %   <repo root>/
    %       install.m
    %       src/    % contains MVector.m (and any subfolders)
    %       test/   % contains Test.m (and any subfolders)
    %
    % Notes
    %   - We add: <root>, <root>/src/**, <root>/test/**
    %   - 'reset' removes only src/test paths (keeps <root> so install remains callable)
    %   - Extra arguments (strings/chars) are accepted case-insensitively.

    if nargin < 1
        option = '';
    end
    if isstring(option); option = char(option); end
    option = lower(strtrim(option));

    % --- Locate repo root (directory where this file resides)
    rootDir = fileparts(mfilename('fullpath'));
    srcDir  = fullfile(rootDir, 'src');
    testDir = fullfile(rootDir, 'test');

    % --- Build path lists: root (single), src recursive, test recursive
    rootPaths = {rootDir}; % keep root on path so 'install' is callable
    srcPaths  = {};
    testPaths = {};

    if exist(srcDir,'dir')
        srcPaths = strsplit(genpath(srcDir), pathsep);
    end
    if exist(testDir,'dir')
        testPaths = strsplit(genpath(testDir), pathsep);
    end

    % --- Combine (and keep order: root, then src, then test)
    allPaths = [rootPaths, srcPaths, testPaths];

    % --- Clean empties and duplicates
    allPaths = allPaths(~cellfun(@isempty, allPaths));
    allPaths = unique(allPaths, 'stable');

    % --- Exclusions (ignore common non-user dirs if present)
    skipPatterns = { [filesep '.git'], [filesep '.github'], [filesep '.svn'], ...
        [filesep '.hg'],  [filesep '.vscode'], 'node_modules', ...
        [filesep 'Trash'], [filesep 'codegen'], '__pycache__' };

    keepMask = true(1, numel(allPaths));
    for i = 1:numel(allPaths)
        p = allPaths{i};
        for k = 1:numel(skipPatterns)
            if contains(p, skipPatterns{k})
                keepMask(i) = false;
                break;
            end
        end
    end
    allPaths = allPaths(keepMask);

    % --- Separate the ones we add vs. the ones we remove on 'reset'
    % We always add root/src/test; on 'reset' we remove only src/test.
    addPaths   = allPaths;                 % root + src + test
    removeMask = ~strcmp(allPaths, rootDir); % don't remove root on reset
    rmPaths    = allPaths(removeMask);     % src + test only

    switch option
        case ''   % ----- temporary add
            added = 0;
            for i = 1:numel(addPaths)
                addpath(addPaths{i});
                added = added + 1;
            end
            fprintf('[MVector] Added %d folder(s) to MATLAB path.\n', added);

            % quick sanity check
            if exist('MVector','class') == 8
                fprintf('[MVector] Class detected: MVector is ready.\n');
            else
                warning('[MVector] MVector class not detected on path. Ensure src/MVector.m exists.');
            end

        case 'save'   % ----- add and persist
            added = 0;
            for i = 1:numel(addPaths)
                addpath(addPaths{i});
                added = added + 1;
            end
            fprintf('[MVector] Added %d folder(s) to MATLAB path.\n', added);

            status = savepath;
            if status == 0
                fprintf('[MVector] Path saved successfully. Available in future sessions.\n');
            else
                warning('[MVector] Unable to save the MATLAB path automatically. Try running MATLAB with elevated/administrator privileges, or use "Set Path" UI to save.');
            end

            if exist('MVector','class') == 8
                fprintf('[MVector] Class detected: MVector is ready.\n');
            else
                warning('[MVector] MVector class not detected on path. Ensure src/MVector.m exists.');
            end

        case 'reset'  % ----- remove src/test only
            removed = 0;
            for i = 1:numel(rmPaths)
                if contains(path, rmPaths{i})
                    rmpath(rmPaths{i});
                    removed = removed + 1;
                end
            end
            fprintf('[MVector] Removed %d folder(s) from MATLAB path (src/test only). Root kept.\n', removed);

        otherwise
            error('MVector:install:Option','Unknown option "%s". Use: install, install(''save''), or install(''reset'').', option);
    end
end
