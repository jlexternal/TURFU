function run_pres(pres)

% add toolboxes to path
addpath('./Toolboxes/Rand');
addpath('./Toolboxes/IO');
addpath('./Toolboxes/Stimuli/Visual');

% default configuration for flat panel displays
iscr     = 'max';
res      = [1920,1080];
fps      = 60;
ppd      = 40;
syncflip = true;

% set stimulation parameters
lumibg   = 128/255; % background luminance
fxtndmtr = deg2pix(24/60,2); % fixation point diameter
probdmtr = deg2pix(40/60,2); % response probe diameter
probwdth = deg2pix(4/60,1); % response probe width

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
    keyresp = KbName({'E','P'});
    
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
    
    scat = {'orange','bleu'};
    labeltxt = sprintf('paquet %s',scat{pres.seqcat(1)});
    labelrec = CenterRectOnPoint(Screen('TextBounds',video.h,labeltxt),video.x/2,video.y/2-5*ppd);
    
    nseq = length(pres.seqlen);
    seqang = cell(1,nseq);
    catang = nan(1,nseq);
    
    % wait for key to start
    Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
    Screen('DrawTexture',video.h,framtex,[],cardrec);
    Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
    Screen('DrawTexture',video.h,probtex,[],probrec);
    Screen('DrawingFinished',video.h);
    Screen('Flip',video.h);
    WaitKeyPress(keywait);
    Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
    Screen('DrawTexture',video.h,framtex,[],cardrec);
    Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
    Screen('DrawingFinished',video.h);
    tflip = Screen('Flip',video.h);
    
    for iseq = 1:nseq
        
        % check if abort key is pressed
        if CheckKeyPress(keyquit)
            aborted = true;
            break
        end
        
        % set inter-trial interval
        tonset = tflip+roundfp(1.250,0.250);
        
        % set sample angles
        [seqang{iseq},catang(iseq)] = get_seqang(pres.seqtlt{iseq},pres.catmap(iseq),pres.seqcat(iseq));
        
        % make card textures
        cfg  = [];
        cfg.ppd  = ppd;
        cfg.axis = catang(iseq);
        cardtex = cell(1,2);
        cfg.cardtype = 1;
        img = make_card(cfg);
        cardtex{1} = Screen('MakeTexture',video.h,img,[],[],2);
        ncards = pres.seqlen(iseq);
        cardtex{2} = zeros(ncards,1);
        for i = 1:ncards
            cfg.cardtype = 3;
            cfg.angl = seqang{iseq}(i);
            img = make_card(cfg);
            cardtex{2}(i) = Screen('MakeTexture',video.h,img,[],[],2);
        end
        
        % draw mapping cue
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        Screen('DrawTexture',video.h,cardtex{1},[],cardrec);
        Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
        Screen('DrawingFinished',video.h);
        tflip = Screen('Flip',video.h,tonset);
        tonset = tflip+roundfp(1.000,0.050);
        
        % draw cards
        for i = 1:ncards
            Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
            Screen('DrawTexture',video.h,cardtex{2}(i),[],cardrec);
            Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
            Screen('DrawingFinished',video.h);
            tflip = Screen('Flip',video.h,tonset);
            toffset = tflip+roundfp(0.100);
            tonset = tflip+roundfp(0.500,0.050);
            Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
            Screen('DrawTexture',video.h,cardtex{1},[],cardrec);
            Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
            Screen('DrawingFinished',video.h);
            Screen('Flip',video.h,toffset);
        end
        tonset = tflip+roundfp(1.000,0.050);
        
        % draw response probe
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        Screen('DrawTexture',video.h,framtex,[],cardrec);
        Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
        Screen('DrawTexture',video.h,probtex,[],probrec);
        Screen('DrawingFinished',video.h);
        Screen('Flip',video.h,tonset);
        WaitKeyPress(keywait);
        Screen('DrawText',video.h,labeltxt,labelrec(1),labelrec(2),0);
        Screen('DrawTexture',video.h,framtex,[],cardrec);
        Screen('DrawTexture',video.h,fxtntex,[],fxtnrec);
        Screen('DrawingFinished',video.h);
        tflip = Screen('Flip',video.h);
        
    end
    
    % close Psychtoolbox
    Priority(0);
    Screen('CloseAll');
    FlushEvents;
    ListenChar(0);
    ShowCursor;
    
catch
    
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

end