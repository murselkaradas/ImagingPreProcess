classdef ExampleUserClass2 < ws.UserClass

    % This is a very simple user class.  It writes to the console when
    % things like a sweep start/end happen.
    
    % Information that you want to stick around between calls to the
    % functions below, and want to be settable/gettable from outside the
    % object.
    properties
        Greeting = 'Hello, there!'
        TimeAtStartOfLastRunAsString_ = ''  
          % TimeAtStartOfLastRunAsString_ should only be accessed from 
          % the methods below, but making it protected is a pain.
        DMD
    end

    methods        
        function self = ExampleUserClass2()
            % creates the "user object"

        end
        
        function wake(self, rootModel)  %#ok<INUSD>
            % creates the "user object"
            fprintf('%s  Waking an instance of ExampleUserClass.\n', ...
                    self.Greeting);
                
            self.DMD.dll_name = 'alpV42basic';
            dll_header = 'alpbasic';
            dll = [self.DMD.dll_name, '.dll'];
            head = [dll_header, '.h'];
            if ~libisloaded(dll)
               loadlibrary(dll,head) %Load the alp42basic library (X64 bit)
            end
            return_lib = libisloaded(self.DMD.dll_name);
            if return_lib == 1
                return_lib = 'Library is loaded';
                disp(return_lib)
            else 
                return_lib = 'Error: Library was not loaded';
            end
            if libisloaded('alpV42basic') 
                %%Allocate the DMD
                deviceid = uint32(0);
                hdevice = uint32(1);
                hdeviceptr = libpointer('longPtr', hdevice); %make an outpointer to write
                disp('Allocate DMD')
                [return_allocate, self.DMD.hdevice] = calllib(self.DMD.dll_name, 'AlpbDevAlloc', deviceid, hdeviceptr);
                
                disp(return_allocate)
                self.DMD.reset_mode = int32(4); %For global
                self.DMD.reset_address = int32(0); %For global
                disp('Reset DMD')
                [return_reset] = calllib(self.DMD.dll_name, 'AlpbDevReset', self.DMD.hdevice, self.DMD.reset_mode, self.DMD.reset_address);               
                disp(return_reset)

                %% Clear the DMD
                first_block = int32(0);
                last_block = int32(15);
                disp('Clear DMD')
                [return_clear] = calllib(self.DMD.dll_name, 'AlpbDevClear', self.DMD.hdevice, first_block, last_block);
                disp(return_clear)
                
                query = int32(4); % Device Serial
                return_query = int32(0);
                return_queryptr = libpointer('int32Ptr', return_query);
                [return_inquiry, return_query] = calllib(self.DMD.dll_name, 'AlpbDevInquire', self.DMD.hdevice, query , return_queryptr);
                
                [return_reset] = calllib(self.DMD.dll_name, 'AlpbDevReset', self.DMD.hdevice, self.DMD.reset_mode, self.DMD.reset_address); 
                disp(return_reset)
                
                self.DMD.Nsweep = 0;
                fileList = struct2cell(dir(fullfile('Samples','*.png')));
                Nimage = size(fileList ,2);
                self.DMD.imgseq = [];
                self.DMD.Nimage_seq = 0;
                for i = 1: Nimage
                    image = (imread(fullfile(fileList{2,i},fileList{1,i})));
                    %% prepare the image to follow 8 bit format, MSB should be 0/1
                    if max(max(image)) == 1
                        image = image.*255;
                    elseif max(max(image)) == 255
                        image = image;
                    else disp('Please check, Image is not in the right format')
                    end
                    if isequal(size(image),[768 1024])
                        self.DMD.Nimage_seq = self.DMD.Nimage_seq + 1;
                        self.DMD.imgseq(:,:,self.DMD.Nimage_seq) = image;
                    end
                end
            end   
        end
        
        function delete(self)
            % Called when there are no more references to the object, just
            % prior to its memory being freed.
            fprintf('%s  An instance of ExampleUserClass is being deleted.\n', ...
                    self.Greeting);
        end
        
        function willSaveToProtocolFile(self, wsModel)  %#ok<INUSD>
            fprintf('%s  Saving to protocol file in ExampleUserClass.\n', ...
                    self.Greeting);
        end        
        
        % These methods are called in the frontend process
        function startingRun(self, wsModel)  %#ok<INUSD>
            % Called just before each set of sweeps (a.k.a. each
            % "run")
            self.TimeAtStartOfLastRunAsString_ = datestr( clock() ) ;
            fprintf('%s  About to start a run.  Current time: %s\n', ...
                    self.Greeting,self.TimeAtStartOfLastRunAsString_);

        end
        
        function completingRun(self, wsModel)  %#ok<INUSD>
            % Called just after each set of sweeps (a.k.a. each
            % "run")
            fprintf('%s  Completed a run.  Time at start of run: %s\n', ...
                    self.Greeting,self.TimeAtStartOfLastRunAsString_);

        end
        
        function stoppingRun(self, wsModel)  %#ok<INUSD>
            % Called if a sweep is manually stopped
            fprintf('%s  User stopped a run.  Time at start of run: %s\n', ...
                    self.Greeting,self.TimeAtStartOfLastRunAsString_);
           [return_free] = calllib(self.DMD.dll_name, 'AlpbDevFree', self.DMD.hdevice);
           unloadlibrary(self.DMD.dll_name); 
        end        
        
        function abortingRun(self, wsModel)  %#ok<INUSD>
            % Called if a run goes wrong, after the call to
            % abortingSweep()
            fprintf('%s  Oh noes!  A run aborted.  Time at start of run: %s\n', ...
                    self.Greeting,self.TimeAtStartOfLastRunAsString_);
			[return_free] = calllib(self.DMD.dll_name, 'AlpbDevFree', self.DMD.hdevice);
            unloadlibrary(self.DMD.dll_name);
        end
        
        function startingSweep(self, wsModel)  %#ok<INUSD>
            % Called just before each sweep
            fprintf('%s  About to start a sweep.  Time at start of run: %s\n', ...
                    self.Greeting,self.TimeAtStartOfLastRunAsString_);
            tic   
            first_row = int32(0);
            last_row = int32(767);
            self.DMD.Nsweep = self.DMD.Nsweep + 1;
            self.DMD.image   = squeeze(self.DMD.imgseq(:,:,mod(self.DMD.Nsweep,self.DMD.Nimage_seq+1)));
            imageptr = libpointer('uint8Ptr', self.DMD.image');
            tic
            [return_reset] = calllib(self.DMD.dll_name, 'AlpbDevReset', self.DMD.hdevice, self.DMD.reset_mode, self.DMD.reset_address);
            [return_load, image2] = calllib(self.DMD.dll_name, 'AlpbDevLoadRows', self.DMD.hdevice, imageptr, first_row, last_row);
            [return_reset] = calllib(self.DMD.dll_name, 'AlpbDevReset', self.DMD.hdevice, self.DMD.reset_mode, self.DMD.reset_address);
%                 return_check(return_load)
            disp(return_load)
            toc
        end
        
        function completingSweep(self, wsModel)  %#ok<INUSD>
            % Called after each sweep completes
            fprintf('%s  Completed a sweep.  Time at start of run: %s\n', ...
                    self.Greeting,self.TimeAtStartOfLastRunAsString_);
        end
        
        function stoppingSweep(self, wsModel)  %#ok<INUSD>
            % Called if a sweep goes wrong
            fprintf('%s  User stopped a sweep.  Time at start of run: %s\n', ...
                    self.Greeting,self.TimeAtStartOfLastRunAsString_);
        end        
        
        function abortingSweep(self, wsModel)  %#ok<INUSD>
            % Called if a sweep goes wrong
            fprintf('%s  Oh noes!  A sweep aborted.  Time at start of run: %s\n', ...
                    self.Greeting,self.TimeAtStartOfLastRunAsString_);
        end        
        
        function dataAvailable(self, wsModel)
            % Called each time a "chunk" of data (typically 100 ms worth) 
            % has been accumulated from the looper.
            analogData = wsModel.getLatestAIData() ;
            digitalData = wsModel.getLatestDIData() ; 
            nAIScans = size(analogData,1) ;
            nDIScans = size(digitalData,1) ;
            fprintf('%s  Just read %d scans of analog data and %d scans of digital data.\n', self.Greeting, nAIScans, nDIScans) ;
        end
        
        function startingEpisode(self, refiller)  %#ok<INUSD>
            % Called just before each episode
            fprintf('%s  About to start an episode.\n',self.Greeting);
        end
        
        function completingEpisode(self, refiller)  %#ok<INUSD>
            % Called after each episode completes
            fprintf('%s  Completed an episode.\n',self.Greeting);
        end
        
        function stoppingEpisode(self, refiller)  %#ok<INUSD>
            % Called if a episode goes wrong
            fprintf('%s  User stopped an episode.\n',self.Greeting);
        end        
        
        function abortingEpisode(self, refiller)  %#ok<INUSD>
            % Called if a episode goes wrong
            fprintf('%s  Oh noes!  An episode aborted.\n',self.Greeting);
        end
    end  % methods
    
    methods 
        % Allows access to private and protected variables for encoding.
        function out = getPropertyValue_(self, name)
            out = self.(name);
        end
        
        % Allows access to protected and protected variables for encoding.
        function setPropertyValue_(self, name, value)
            self.(name) = value;
        end        
    end  % protected methods block
    
    methods
        function mimic(self, other)
            ws.mimicBang(self, other) ;
        end
    end    
    
    methods
        % These are intended for getting/setting *public* properties.
        % I.e. they are for general use, not restricted to special cases like
        % encoding or ugly hacks.
        function result = get(self, propertyName) 
            result = self.(propertyName) ;
        end
        
        function set(self, propertyName, newValue)
            self.(propertyName) = newValue ;
        end           
    end  % public methods block            
    
end  % classdef

