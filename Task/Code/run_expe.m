function [expe,aborted,errmsg] = run_expe(expe,eyetrack)
%  RUN_EXPE  Run ACTOBS experiment
%
%  Usage: [expe,aborted,errmsg] = RUN_EXPE(expe,eyetrack)
%
%  where the experiment structure expe must contain a block substructure blck
%  containing the information regarding the experiment blocks, and a header
%  substructure hdr containing the subject number.
%
%  The function returns the updated experiment structure expe with an updated
%  block substructure blck, and two additional fields:
%    * the event substructure evnt containing the timestamps of all events of
%      every experiment block (in relative time from block start),
%    * the results substructure rslt containing the responses and latencies of
%      all choices of every experiment block.
%
%  The function also returns the flag aborted indicating if the experiment was
%  aborted prematurely, and an optional error message structure errmsg which
%  contains caughts errors during the experiment if any (empty otherwise). If
%  errmsg is not asked as output, the function will re-throw any caught error
%  after shutting down Psychtoolbox and additional hardware interfaces.
%
%  Valentin Wyart <valentin.wyart@ens.fr> - 09/2015

% check input arguments
if nargin < 2
    % assuming not CENIR MEG setup
    eyetrack = false;
end
if nargin < 1
    error('Missing experiment structure!');
end
if ~all(isfield(expe,{'hdr','blck'}))
    error('Incomplete experiment structure!');
end

% check existence of fitting functions
if ~exist('fminsearch') || ~exist('fminbnd')
    error('Missing fitting functions!');
end

% add toolboxes to path
addpath('./Toolboxes/Rand');
addpath('./Toolboxes/IO');
addpath('./Toolboxes/Stimuli/Visual');

% set screen parameters
if eyetrack
    % eye-tracker configuration
    iscr     = 1;
    res      = [1920,1080];
    fps      = 60;
    ppd      = 40;
    syncflip = true;
else
    % default configuration for flat panel displays
    iscr     = 1;
    res      = 'max';
    fps      = [];
    ppd      = 40;
    syncflip = false;
end

% set stimulation parameters
lumibg   = 128/255; % background luminance
fxtndmtr = deg2pix(24/60,2); % fixation point diameter
probdmtr = deg2pix(40/60,2); % response probe diameter
probwdth = deg2pix(4/60,1); % response probe width
printmsg = true; % print messages in command window?
%eyeoffset = deg2pix(3.30,1); % volatility banner position offset
voloffset = deg2pix(6.60,1); % volatility banner position offset
% set file saving parameters
filepath = fileparts(mfilename('fullpath'));
datapath = [filepath,filesep,'..',filesep,'Data'];

aborted = false;
errmsg  = [];
video   = [];

try
    
    % hide cursor and stop spilling key presses into MATLAB windows
    HideCursor;
    FlushEvents;
    ListenChar(2);
    
    % set response keys/buttons
    KbName('UnifyKeyNames');
    keywait = KbName('space');
    keyquit = KbName('ESCAPE');
    if mod(expe.hdr.suj,2) == 1 % 1:index 2:middle
        keyresp = KbName({'s','a','j','i'});
    else % 1:middle 2:index
        keyresp = KbName({'a','s','i','j'});
    end
    
    % open main window
    % set screen resolution and refresh rate if required
    if ~isempty(res)
        r = Screen('Resolutions',iscr);
        if strcmp(lower(res),'max')
            res = [max([r.width]),max([r.height])];
        elseif numel(res) == 2
            i = find([r.width] == res(1) & [r.height] == res(2));
            if isempty(i)
                error('Cannot set screen to %d x %d.',res(1),res(2));
            end
        else
            error('Invalid screen resolution.');
        end
        if ~isempty(fps)
            if ~any([r(i).hz] == fps)
                error('Cannot set screen to %d x %d at %d Hz.',res(1),res(2),fps);
            end
            Screen('Resolution',iscr,res(1),res(2),fps);
        else
            Screen('Resolution',iscr,res(1),res(2));
        end
    end
    % set screen synchronization properties
    % see 'help SyncTrouble',
    %     'help BeampositionQueries' or
    %     'help ConserveVRAMSettings' for more information
    if syncflip
        if ispc
            % soften synchronization test requirements
            Screen('Preference','SyncTestSettings',[],[],0.2,10);
            % enforce beamposition workaround for missing VBL interval
            Screen('Preference','ConserveVRAM',bitor(4096,Screen('Preference','ConserveVRAM')));
        end
        Screen('Preference','VisualDebuglevel',3);
        Screen('Preference','SuppressAllWarnings',0);
    else
        % skip synchronization tests altogether
        Screen('Preference','SkipSyncTests',1);
        Screen('Preference','VisualDebuglevel',0);
        Screen('Preference','SuppressAllWarnings',1);
    end
    % set font properties
    if ismac
        txtfnt = 'Helvetica';
        txtsiz = round(1.0*ppd);
    elseif ispc
        txtfnt = 'Arial'; % closest to Helvetica
        txtsiz = round(2/3*ppd); % text size is ~2/3 smaller in Windows than MacOSX
    end
    Screen('Preference','TextAlphaBlending',1);
    Screen('Preference','DefaultFontName',txtfnt);
    Screen('Preference','DefaultFontSize',txtsiz);
    Screen('Preference','DefaultFontStyle',0);
    % prepare configuration and open main window
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask','General','UseFastOffscreenWindows');
    PsychImaging('AddTask','General','NormalizedHighresColorRange');
    video.i = iscr;
    video.res = Screen('Resolution',video.i);
    video.h = PsychImaging('OpenWindow',video.i,0);
    [video.x,video.y] = Screen('WindowSize',video.h);
    video.ifi = Screen('GetFlipInterval',video.h,100,50e-6,10);
    Screen('BlendFunction',video.h,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    Priority(MaxPriority(video.h));
    Screen('ColorRange',video.h,1);
    Screen('FillRect',video.h,lumibg);
    Screen('Flip',video.h);
    % check screen refresh rate
    if ~isempty(fps) && fps > 0 && round(1/video.ifi) ~= fps
        error('Screen refresh rate not equal to expected %d Hz.',fps);
    end

    if eyetrack
        % initialize eye-tracker connection
        if EyelinkInit() ~= 1
            error('could not initialize eye-tracker connection!');
        end
        [~,ev] = Eyelink('GetTrackerVersion');
        fprintf('Connection to %s eye-tracker initialized successfully.\n',ev);
        % setup eye-tracker
        el = EyelinkInitDefaults(video.h);
        setup_eyelink;
    end
    
    % create card frame texture
    cfg = [];
    cfg.ppd = ppd;
    cfg.cardtype = 0; % empty frame
    img = make_card(cfg);
    framtex = Screen('MakeTexture',video.h,img,[],[],2);
    clear('img');
    cardrec = CenterRectOnPoint(Screen('Rect',framtex),video.x/2,video.y/2);
    
    % create fixation point texture
    img = cat(3,zeros(fxtndmtr),CreateCircularAperture(fxtndmtr));
    fxtntex = Screen('MakeTexture',video.h,img,[],[],2);
    clear('img');
    fxtnrec = CenterRectOnPoint(Screen('Rect',fxtntex),video.x/2,video.y/2);
    
    % create response probe texture
    img = cat(3,zeros(probdmtr),CreateCircle(probdmtr,probwdth));
    probtex = Screen('MakeTexture',video.h,img,[],[],2);
    clear('img');
    probrec = CenterRectOnPoint(Screen('Rect',probtex),video.x/2,video.y/2);
    
    % create instruction texture
    if mod(expe.hdr.suj,2) == 1 % 1:index 2:middle
        k = 1;
    else % 1:middle 2:index
        k = 2;
    end
    instrtex = zeros(2,2);
    for i = 1:2 % taskid => 1:observer or 2:actor
        for j = 1:2 % epimap => 1:pink=good|left or 2:blue=good|left
            img = double(imread(fullfile(filepath,'img',sprintf('instr%d%d%d.png',i,j,k))))/255;
            instrtex(i,j) = Screen('MakeTexture',video.h,img,lumibg,[],2);
        end
    end
    clear('img');
    instrrec = Screen('Rect',instrtex(1,1));
    Screen('FillRect',video.h,lumibg); 

    % rescale instruction textures on-the-fly
    s = 800/RectWidth(instrrec); % scaling factor
    instrrec = CenterRectOnPoint(ScaleRect(instrrec,s,s),video.x/2,video.y/2);%+eyeoffset);
    
    
    % create volatility banners
    volnrtex = zeros(2,1);
    volnrimg = double(imread(fullfile(filepath,'img','epivol1.png')))/255;
    volnrtex(1) = Screen('MakeTexture',video.h,volnrimg,[],[],2);
    volnrrec(1,:) = CenterRectOnPoint(Screen('Rect',volnrtex(1)),video.x/2,video.y/2+voloffset);
    volnrimg = double(imread(fullfile(filepath,'img','epivol2.png')))/255;
    volnrtex(2) = Screen('MakeTexture',video.h,volnrimg,[],[],2);
    volnrrec(2,:) = CenterRectOnPoint(Screen('Rect',volnrtex(2)),video.x/2,video.y/2+voloffset);
    
    % first flip
    tflip = Screen('Flip',video.h);

    nblck = length(expe.blck);
    for iblck = 1:nblck
        
        blck = expe.blck(iblck);
        nseq = length(blck.seqang);
        iseq = 0;
        
        seqcat = zeros(1,nseq); % sequence category => 1:orange or 2:blue
        seqang = cell(1,nseq);  % sequence sample angles (rad)
        catang = zeros(1,nseq); % categorization angle (rad) => will *not* be saved!
        resp   = zeros(1,nseq+1); % response => 1:left or 2:rght
        conf   = zeros(1,nseq+1); % confidence => 1:low or 2:high
        rt     = zeros(1,nseq+1); % response time (s)
        
        ivol = expe.blck(iblck).condtn;
        
        % initialize event structure
        evnt      = [];
        evnt.iseq = [];
        evnt.time = [];
        evnt.type = {};
        
        if eyetrack
            % calibrate eye-tracker
            EyelinkDoTrackerSetup(el,el.ENTER_KEY);
            tflip = Screen('Flip',video.h);
        end
        
        if printmsg
            fprintf('\n\nBLOCK %d: TASKID = %d, CONDTN = %d, EPIMAP = %d\n\n', ...
                iblck,blck.taskid,blck.condtn,blck.epimap);
        end

        % draw instruction screen
        Screen('DrawTexture',video.h,instrtex(blck.taskid,blck.epimap),[],instrrec);
        if ivol < 3
        Screen('DrawTexture',video.h,volnrtex(ivol),[],volnrrec(ivol,:));
        end
        Screen('DrawingFinished',video.h);
        Screen('Flip',video.h,tflip+roundfp(2.000,0.500));
        WaitKeyPress(keywait);

        % draw card frame
        Screen('DrawTexture',video.h,framtex,[],cardrec);
        Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
        Screen('DrawingFinished',video.h);
        tflip = Screen('Flip',video.h);
        tstart = tflip;
        if eyetrack
            tonset = tflip+roundfp(10.000,0.500);
        else
            tonset = tflip+roundfp(2.000,0.500);
        end
        
        if eyetrack
            % start eye-tracker recording
            Eyelink('OpenFile','ACTOBS');
            Eyelink('StartRecording');
            eye_used = Eyelink('EyeAvailable');
        end
        
        % draw response probe
        Screen('DrawTexture',video.h,framtex,[],cardrec);
        Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
        Screen('DrawTexture',video.h,probtex,[],probrec);
        Screen('DrawingFinished',video.h);
        tflip = Screen('Flip',video.h,tonset);
        % wait for 1st response
        tgo = tflip;
        key = 0;
        while ~key
            [key,tkey] = CheckKeyPress(keyresp);
        end
        resp(1) = ceil(key/2);
        conf(1) = mod(key-1,2)+1;
        rt(1) = tkey-tgo;
        flag_evnt(tkey,'RESP');
        Screen('DrawTexture',video.h,framtex,[],cardrec);
        Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
        Screen('DrawingFinished',video.h);
        tflip = Screen('Flip',video.h);
        if printmsg
            fprintf('[sequence 000] repdir = %d, reptim = %.3f s\n',resp(1),rt(1));
            fprintf('               repcon = %d\n\n',conf(1));
        end
        
        for iseq = 1:nseq
            
            % check if abort key is pressed
            if CheckKeyPress(keyquit)
                aborted = true;
                break
            end
            
            % set inter-trial interval
            tonset = tflip+roundfp(1.250,0.250);
            
            % set sequence category and sample angles
            seqcat(iseq) = get_seqcat(blck.taskid,blck.epimap,blck.seqdir(iseq),resp(iseq));
            [seqang{iseq},catang(iseq)] = get_seqang(blck.seqtlt{iseq},blck.catmap(iseq),seqcat(iseq));
            if printmsg
                fprintf('[sequence %03d] seqdir = %d, seqllr = %+05.1f, seqcat = %d\n', ...
                    iseq,blck.seqdir(iseq),blck.seqllr(iseq),seqcat(iseq));
            end
            
            % make card textures
            cfg  = [];
            cfg.ppd  = ppd;
            cfg.axis = catang(iseq);
            cardtex = cell(1,2);
            cfg.cardtype = 1;
            img = make_card(cfg);
            cardtex{1} = Screen('MakeTexture',video.h,img,[],[],2);
            ncards = blck.seqlen(iseq);
            cardtex{2} = zeros(ncards,1);
            for i = 1:ncards
                cfg.cardtype = 3;
                cfg.angl = seqang{iseq}(i);
                img = make_card(cfg);
                cardtex{2}(i) = Screen('MakeTexture',video.h,img,[],[],2);
            end
            
            % check timing
            if printmsg && GetSecs > tonset
                fprintf('               missed deadline during sequence generation!\n');
                tonset = 0;
            end
            
            % draw mapping cue
            Screen('DrawTexture',video.h,cardtex{1},[],cardrec);
            Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
            Screen('DrawingFinished',video.h);
            tflip = Screen('Flip',video.h,tonset);
            flag_evnt(tonset,'CUE');
            if eyetrack
                Eyelink('Message',sprintf('SEQ%03d_CUE',iseq));
            end
            tonset = tflip+roundfp(1.000,0.050);
            
            % draw cards
            for i = 1:ncards
                Screen('DrawTexture',video.h,cardtex{2}(i),[],cardrec);
                Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
                Screen('DrawingFinished',video.h);
                while GetSecs < tonset
                    % check for early key presses
                    if CheckKeyPress(keyresp) > 0
                        flag_evnt(GetSecs,'RESP_EARLY');
                        break
                    end
                end
                tflip = Screen('Flip',video.h,tonset);
                flag_evnt(tonset,'CARD');
                if eyetrack
                    Eyelink('Message',sprintf('SEQ%03d_CARD',iseq));
                end
                toffset = tflip+roundfp(0.100);
                tonset = tflip+roundfp(0.500,0.050);
                
                Screen('DrawTexture',video.h,cardtex{1},[],cardrec);
                Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
                Screen('DrawingFinished',video.h);
                while GetSecs < toffset
                    % check for early key presses
                    if CheckKeyPress(keyresp) > 0
                        flag_evnt(GetSecs,'RESP_EARLY');
                        break
                    end
                end
                Screen('Flip',video.h,toffset);
                
            end
            tonset = tflip+roundfp(0.500,0.050);

            % draw response probe
            Screen('DrawTexture',video.h,framtex,[],cardrec);
            Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
            Screen('DrawTexture',video.h,probtex,[],probrec);
            Screen('DrawingFinished',video.h);
            tflip = Screen('Flip',video.h,tonset);
            % wait for response
            tgo = tflip;
            key = 0;
            while ~key
                [key,tkey] = CheckKeyPress(keyresp);
            end
            resp(iseq+1) = ceil(key/2);
            conf(iseq+1) = mod(key-1,2)+1;
            rt(iseq+1) = tkey-tgo;
            flag_evnt(tkey,'RESP');
            if eyetrack
                Eyelink('Message',sprintf('SEQ%03d_RESP',iseq));
            end
            Screen('DrawTexture',video.h,framtex,[],cardrec);
            Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
            Screen('DrawingFinished',video.h);
            tflip = Screen('Flip',video.h);
            tonset = tflip+roundfp(0.250);
            if printmsg
                fprintf('               repdir = %d, reptim = %.3f s\n',resp(iseq+1),rt(iseq+1));
                fprintf('               repcon = %d\n\n',conf(iseq+1));
            end
            
            % draw empty card frame
            Screen('DrawTexture',video.h,framtex,[],cardrec);
            Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
            Screen('DrawingFinished',video.h);
            tflip = Screen('Flip',video.h,tonset);
            
            % close textures
            Screen('Close',cardtex{1});
            Screen('Close',cardtex{2});
            
        end
        
        % update block substructure
        blck.seqcat = seqcat; % sequence category => 1:orange or 2:blue
        blck.seqang = seqang; % list of sample angles (rad)
        
        % initialize results substructure
        rslt      = [];
        rslt.resp = resp;
        rslt.rt   = rt;
        rslt.conf = conf;
        
        % update experiment structure
        expe.blck(iblck) = blck;
        expe.evnt(iblck) = evnt;
        expe.rslt(iblck) = rslt;
        
        if eyetrack
            % stop recording
            Eyelink('StopRecording');
            Eyelink('CloseFile');
            % retrieve eye-tracker datafile
            if ~aborted
                fpath = datapath;
                fname = sprintf('ACTOBS_C_S%02d_b%d_%s',expe.hdr.suj,iblck,datestr(now,'yyyymmdd-HHMM'));
                fname = fullfile(fpath,fname);
                for i = 1:10
                    status = Eyelink('ReceiveFile',[],[fname,'.edf']);
                    if status >= 0 % should be > not >=
                        break
                    end
                end
                if status < 0 % should be <= not <
                    warning('Could not retrieve eye-tracker datafile %s!',fname);
                else
                    save([fname,'.mat'],'expe');
                end
            end
        end
        
        if aborted
            break
        end
        
        if printmsg
            % compute p(correct)
            pcor_sub = mean(resp(2:end) == blck.seqdir);
            fprintf('p(correct) = %.1f %%\n',pcor_sub*100);
            fprintf('average rt = %.3f s, %d response(s) > 2 s\n\n',mean(rt(2:end)),nnz(rt(2:end) > 2));
            if blck.condtn < 3
                % compute logistic regression
                xreg = [blck.seqllr',3-2*resp(1:nseq)',ones(nseq,1)];
                yreg = resp(2:end)' == 1;
                breg_sub = logreg(xreg,yreg,'probit');
                % run model
                [pcor_ref,breg_ref] = run_model(blck,0.5);
                fprintf('p(correct) = %.1f %% vs. %.1f %% (subject vs. reference)\n',pcor_sub*100,pcor_ref*100);
                fprintf('b[seqllr]  = %+.2f  vs. %+.2f  (subject vs. reference)\n',breg_sub(1),breg_ref(1));
                fprintf('b[resp-1]  = %+.2f  vs. %+.2f  (subject vs. reference)\n\n',breg_sub(2),breg_ref(2));
            end
        end
        
    end
    
    if eyetrack
        % reset eye-tracker to default configuration
        reset_eyelink;
        % close eye-tracker connection
        Eyelink('Shutdown');
    end
    
    % close Psychtoolbox
    Priority(0);
    Screen('CloseAll');
    FlushEvents;
    ListenChar(0);
    ShowCursor;
    
catch
    
    if eyetrack
        if exist('el','var')
            % stop recording
            Eyelink('StopRecording');
            % reset EyeLink system to default configuration
            reset_eyelink;
            % close eye-tracker connection
            Eyelink('Shutdown');
        end
    end
    
    % close Psychtoolbox
    Priority(0);
    Screen('CloseAll');
    FlushEvents;
    ListenChar(0);
    ShowCursor;
    
    % handle thrown error
    if nargout > 2
        errmsg = lasterror;
        errmsg = rmfield(errmsg,'stack');
    else
        rethrow(lasterror);
    end
    
end

    function [t] = roundfp(t,dt)
        % apply duration rounding policy for video flips
        % where t  - desired (input)/rounded (output) duration
        %       dt - desired uniform jitter on duration (default: none)
        n = round(t/video.ifi);
        % apply uniform jitter
        if nargin > 1 && dt > 0
            m = round(dt/video.ifi);
            n = n+ceil((m*2+1)*rand)-(m+1);
        end
        % convert frames to duration
        t = (n-0.5)*video.ifi;
    end

    function [p] = deg2pix(d,b)
        % convert degrees of visual angle into screen pixels
        % where d - degrees of visual angle
        %       b - rounding factor (default: none)
        p = d*ppd;
        if nargin > 1 && b > 0
            p = max(round(p/b),1)*b;
        end
    end

    function [x] = bit2dec(b)
        % convert binary into decimal values
        % where b - binary array (in ascending order of powers of 2)
        x = sum(b(:).*2.^(0:numel(b)-1)');
    end

    function flag_evnt(etime,etype)
        % flag event with timestamp and type
        % where etime - event timestamp (from Psychtoolbox clock)
        %       etype - event type (string)
        evnt.iseq(end+1) = iseq;
        evnt.time(end+1) = etime-tstart;
        evnt.type{end+1} = etype;
    end

    function setup_eyelink
        % setup eye-tracker to custom configuration
        % disable sounds during eye-tracker calibration
        el.targetbeep = 0;
        el.feedbackbeep = 0;
        EyelinkUpdateDefaults(el);
        % setup camera
        Eyelink('Command','active_eye = RIGHT'); % LEFT or RIGHT
        Eyelink('Command','binocular_enabled = NO'); % YES:binocular or NO:monocular
        Eyelink('Command','simulation_screen_distance = 560'); % in mm
        Eyelink('Command','use_ellipse_fitter = NO'); % YES:ellipse or NO:centroid
        Eyelink('Command','pupil_size_diameter = YES'); % YES:diameter or NO:area
        Eyelink('Command','sample_rate = 1000'); % 1000 or 500 or 250
        Eyelink('Command','elcl_tt_power = 2'); % 1:100% or 2:75% or 3:50%
        % setup calibration/validation
        Eyelink('Command','calibration_type = HV5'); % HV5 or HV9
        Eyelink('Command','generate_default_targets = NO'); % YES:default or NO:custom
        cnt = [video.x/2,video.y/2];
        off = ppd*8; % distance from center
        pnt = zeros(5,2);
        pnt(1,:) = cnt;         % 0:center
        pnt(2,:) = cnt-[0,off]; % 1:top
        pnt(3,:) = cnt+[0,off]; % 2:bottom
        pnt(4,:) = cnt-[off,0]; % 3:left
        pnt(5,:) = cnt+[off,0]; % 4:right
        pnt = num2cell(reshape(pnt',[],1));
        Eyelink('Command','calibration_targets = %d,%d %d,%d %d,%d %d,%d %d,%d',pnt{:});
        Eyelink('Command','calibration_samples = 6');
        Eyelink('Command','calibration_sequence = 0,1,2,3,4,0');
        Eyelink('Command','randomize_calibration_order = YES');
        Eyelink('Command','validation_targets = %d,%d %d,%d %d,%d %d,%d %d,%d',pnt{:});
        Eyelink('Command','validation_samples = 5');
        Eyelink('Command','validation_sequence = 0,1,2,3,4');
        Eyelink('Command','randomize_validation_order = YES');
        % setup file output
        Eyelink('Command','file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE');
        Eyelink('Command','file_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS');
    end

    function reset_eyelink
        % reset eye-tracker to default configuration
        Eyelink('Command','sample_rate = 1000');
        Eyelink('Command','calibration_type = HV9');
        Eyelink('Command','generate_default_targets = YES');
        Eyelink('Command','file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('Command','file_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,BUTTON,INPUT');
    end

end