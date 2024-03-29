%DM3 reader for matlab that reads all of the tags and data objects.
%All data types (except packed complex) are read properly
%Output is a struct containing each object (image, spectra, SI, or data cube)
%The tags are optionally saved in a text file with the original image name.
%Set the writeTags flag to 1 in order to write out the tags to a file.
%
%Based on DM3_reader.java for ImageJ
%
%Version 4

%History:
%Written by Peter Ercius 2008/09/05
% - Properly retrieves the sizes, scale, and origin of SI images
% - Retrieves all data objects (i.e. multiple spectra in one DM3 file)
%
% Modified by Huolin Xin 2008/11/19 to enable auto recognization of datatype
%  - energy axis (if any) is passed out
%  - multiple lines in dm can be read
%  - corrected the global variable problem (matlab is evil in global!!!)
%  - no matter what ur file's datatype is
%  (int,uint,ushort,float,double...). The image or SI can be extracted
%  properly
%  - added pick file user interface is no filename is passed in
%  - rgb, complex8 and complex16 data are handled properly now
%  - packed complex8 is not handled yet
%
%Modified P Ercius 2008/11/20
%	- Uses struct and cell
%	- Struct output for spectra and SI have energy axis struct.en()
%	- Distinguishes image, spectra, SI, and other 3D data cubes
% Modified P Ercius 2011/05/20
%	- Modified the array reading function to actually read in some array types.
%		The advantage is that strings are now read properly for dimensions.
%	- Added some parsing code to catch interesting unit names for images and spectra
% Modified P Ercius 2012/08
%   - Script can handle missing "scale" data now without error for other unexpected Data types
%	- add writeTags flag for controlling writing out tags to disk
% Modified P Ercius 2013/01
%   - Fixed bug with global variable writeTags which stopped tags from writing out
% Modified P Ercius 2013/04
%   - Moved global variables to the start of the function to fix conflict
%   between local and global writeTags variable
%
% Need:
%   - add unit check for X-ray data: "eV"
%%
%Main Function
function outputData = dm3Reader(inName,writeTags)
    
    %Clear globals
	clear global curGroupLevel MAXDEPTH curGroupAtLevelX curGroupNameAtLevelX;
    clear global curTagAtLevelX curTagName FIDout writeTagsGlobal;
    clear global dataSize dataOffset xsize ysize zsize scale scaleUnit origin numObjects dataType;
    %Define globals
	global curGroupLevel MAXDEPTH curGroupAtLevelX curGroupNameAtLevelX;
	global curTagAtLevelX curTagName FIDout writeTagsGlobal;
	global dataOffset xsize ysize zsize scale scaleUnit origin numObjects dataType;
    
    if nargin == 0
        [fName,pathName] = uigetfile('*.dm3');
        writeTagsGlobal = 0;
        if fName == 0
            error('No file opened.');
        end
        inName = fullfile(pathName,fName);
    elseif nargin == 1
        fName = inName;
    	writeTagsGlobal = 0; % do not write out the tags if not requested
    elseif nargin == 2
    	writeTagsGlobal = writeTags; %specify the writeTags flag as a global variable
    else
        disp('Too many inputs')
        outputData = [];
        return;
    end
    
	xsize = 0;
	ysize = 0;
	zsize = 1;
	numObjects = 0;
	dataType = 0;
	
	curGroupLevel = 0; %track how deep we currently are in a group
	MAXDEPTH = 64; %maximum number of group levels allowed
	curGroupAtLevelX = zeros(1,MAXDEPTH,'int16'); %track group at current level
	
	curTagAtLevelX = zeros(1,MAXDEPTH); %track tag number at current level
	curTagName = ''; %name of the current tag data item
    
    %Open input file
    
	[FID, FIDmessage] = fopen(inName,'rb');
    if FID == -1
		error(['Error opening file: ', FIDmessage])
    end
    
    %Open output file if necessary
	if(writeTagsGlobal)
		[FIDout, FIDoutmessage] = fopen([inName '_tags.txt'],'wt');
        if(FIDout == -1)
            error(['Issue opening output file: ',FIDoutmessage]);
        end
	end
	
	%Check file for version and data format
	if validDM3(FID) ~= 1
		disp('This does not seem to be a DM3 file...quitting')
		outputData = [];
        return;
	end
	
	%DM3 files have an unnamed root group which contains everything in
	%the file
	curGroupNameAtLevelX = ''; %set the name of the root group
	
	%Read the root tag group (and all of its subgroups)
	readTagGroup(FID);
	
	%Read in the data for each object now that all of the tags and relevant
	%information has been gathered from the tags
	%The first object is not usefule so start the $for$ loop at 2
	fseek(FID,0,-1); %seek to the start of the file
	jj = 0; %counter for the useful number of objects (xsize > 0) Start at 1 because the first object is always a thumbnail

	%Use a struct to hold the data
    %scienceData
	%scienceData = struct('scale',{},'scaleUnit',{},'filename','');
    outputData.filename = fName;
	for ii = 2:numObjects
		if( (xsize(ii) > 0))
            jj = jj+1; %counter for useful objects (xsize > 0)
            fseek(FID,dataOffset(ii),-1); %seek to the start of the data
            stype = getprecisionstr(dataType(ii)); %retrieve correct data type
			if(length(scaleUnit{jj}) >= 0)
                outputData(jj).scaleUnit = scaleUnit;
            end
            if strcmp(stype,'rgb')
                temp = readData(FID,[xsize(ii)*4 ysize(ii) zsize(ii)],'uint8');
                imr = temp( 1:4:(xsize(ii)-1)*4+1,:,:); %red channels
                img = temp( (1:4:(xsize(ii)-1)*4+1)+1,:,:); %green channels
                imb = temp( (1:4:(xsize(ii)-1)*4+1)+2,:,:); %blue channels
                outputData(jj).image = cat(3,imr,img,imb); %combine the channels into RGB image
            elseif strcmp(stype,'complex8')
                temp = readData(FID,[xsize(ii)*2 ysize(ii) zsize(ii)],'float32');
                imr = temp(1:2:(xsize(ii)-1)*2+1,:,:); %real part
                imi = temp((1:2:(xsize(ii)-1)*2+1)+1,:,:); %imaginary part
                outputData(jj).image = imr+imi*sqrt(-1); %create complex array
            elseif strcmp(stype,'complex16')
                temp = readData(FID,[xsize(ii)*2 ysize(ii) zsize(ii)],'double');
                imr = temp(1:2:(xsize(ii)-1)*2+1,:,:);
                imi = temp((1:2:(xsize(ii)-1)*2+1)+1,:,:);
                outputData(jj).image = imr+imi*sqrt(-1);
			else
				if ysize(ii) == 0 || ysize(ii)==1 %data is a spectra
					outputData(jj).sig = fread(FID,[xsize(ii) 1],stype);
					outputData(jj).en = (-origin(ii)*scale(ii) + (0:(xsize(ii)-1))*scale(ii))';
				elseif zsize == 1 %data is an image
					outputData(jj).image = readData(FID,[xsize(ii) ysize(ii) zsize(ii)],stype);
                    if max(size(scale)) > 0
                        outputData(jj).scale = scale(ii);
                    end
				elseif zsize(ii)>1 && max(strcmp(scaleUnit,'eV')) %data is a 3D spectrum image
                    %disp('Reading spectrum image...')
                    outputData(jj).SI = readData(FID,[xsize(ii) ysize(ii) zsize(ii)],stype);
					outputData(jj).en = (-origin(jj)*scale(jj) + (0:(zsize(ii)-1))*scale(jj))';
				else %data is a 3D array
                    %disp('Reading 3D data array...')
					outputData(jj).cube = readData(FID,[xsize(ii) ysize(ii) zsize(ii)],stype);
                    outputData(jj).scale = scale;
				end
			end
		end
	end
    
	fclose(FID);
	
	if(writeTagsGlobal)
		fclose(FIDout);
	end
    
	clear global curGroupLevel MAXDEPTH curGroupAtLevelX curGroupNameAtLevelX;
    clear global curTagAtLevelX curTagName FIDout writeTagsGlobal;
    clear global dataSize dataOffset xsize ysize zsize scale origin numObjects dataType;
end

function str = getprecisionstr(datatype)
	switch datatype
	    case 6
	        str = 'uint8';
	    case 10
	        str = 'uint16';
	    case 11
	        str = 'uint32';
	    case 9
	        str = 'int8';
	    case 1
	        str = 'int16';
	    case 7
	        str = 'int32';
	    case 2
	        str = 'float32';
	    case 12
	        str = 'double';
	    case 14
	        str = 'bit1';
	    case 23
	        str = 'rgb';
	    case 3
	        str = 'complex8';
	    case 13
	        str = 'complex16';
	    case 5
	        str = 'packedcomplex8';
	end
end

%Determine if the file is a valid DM4
function output = validDM3(FID)

	output = 1; %output will stay == 1 if the file is a true DM3 file

	%The first 4 4byte integers equal the file size
	%It is written in Big-endian format.
	head = fread(FID,3,'uint32','ieee-be'); 
	%FileSize = head(2);
	
	%Gatan file version. This must == 3 to continue
	if head(1) ~= 3
		disp('File does not seem to be version DM3')
		output = 0;
	end
	
	%Test for binary endian type. PC data is written as little endian.
	if head(3) ~= 1
		disp('Data not written in little endian (PC) format.');
		output = 0;
	end
	
end

% Read in a three dimensional binary data array from a file
% Based on read3d.m
%	fid: File ID of the file to be read. Opened with ex. fopen(fname,'rb')
%	count: [m,n,p] matrix denoting the size of the matrix
%	precision: type of numbers to read from the binary data (string)
function out = readData(fid,count,precision)
    out = reshape(fread(fid,prod(count),[precision,'=>',precision]),count);    
end

%%
%Tag functions
function readTagGroup(FID)
	global curGroupLevel curGroupAtLevelX curTagAtLevelX curGroupNameAtLevelX;

	curGroupLevel = curGroupLevel + 1; %go down a level; this is a new group
	curGroupAtLevelX(curGroupLevel) = curGroupAtLevelX(curGroupLevel) + 1;
	%Set the # of current tags at this level to -1 since the readTagEntry
	%routine pre-increments. ie the first tag will be labeled 0
	curTagAtLevelX(curGroupLevel) = 0;
	
	fread(FID,1,'uchar'); %isSorted
	fread(FID,1,'uchar'); %isOpen
	nTags = fread(FID,1,'uint32','ieee-be');

	%Iterate over the number of tag entries in this group
	oldTotalTag = curGroupNameAtLevelX;
	for ii = 1:nTags
		readTagEntry(FID);		
	end
	
	%Go back up a level now that we are finished reading this group
	curGroupLevel=1;
	curGroupNameAtLevelX = oldTotalTag;
end

function readTagEntry(FID)
	global curGroupLevel curTagAtLevelX curTagName curGroupNameAtLevelX;
	
	isData = fread(FID,1,'uchar','ieee-be'); %use big endian
	
	%Record that we have found a tag at this level
	curTagAtLevelX(curGroupLevel) = curTagAtLevelX(curGroupLevel)+1;
	
	%Get the tag if one exists
	lenTagLabel=fread(FID,1,'uint16','ieee-be');
    if lenTagLabel ~= 0
        tagLabel = fscanf(FID,'%c',lenTagLabel); %read in the string
    else
        tagLabel = num2str(curTagAtLevelX(curGroupLevel)); %unlabeled tag. Dont need to save this anymore
    end
	
%    disp(['readTagEntry, isData:' num2str(isData) ',tagLabel:' tagLabel])
    
    oldGroupName = curGroupNameAtLevelX;

	if isData == 21
		%This tag entry is data
		
		%Name the piece of data
		%curTagName = [makeGroupString '.' tagLabel];
		%curGroupNameAtLevelX = [curGroupNameAtLevelX '.' tagLabel];
		curTagName = tagLabel;
		
		%Now get the tag data
		readTagType(FID);
	else
		%This tag entry is another tag group
		%Awkward that this cant be done in readTagGroup() but:
		
		%Store the name of the group at the new level
		%curGroupNameAtLevelX(curGroupLevel+1) = {tagLabel};
		curGroupNameAtLevelX = [curGroupNameAtLevelX '.' tagLabel];
		readTagGroup(FID);
	end	
	curGroupNameAtLevelX = oldGroupName;
end

function readTagType(FID)
	Delim = fscanf(FID,'%c',4); %should always be '%%%%'

	if ~strcmp(Delim,'%%%%')
		disp('Tag type delimiter is not "%%%%"')
	end
	fread(FID,1,'uint32','ieee-be'); %nInTag: unnecessary redundant info
	
	readAnyData(FID);
	
end

%%
%Higher level function which dispatches to functions handling specific data
%types
function readAnyData(FID)
	global curTagName;

	%This specifies what kind of type we are dealing with: short, long, struct, array, etc.
	%global encodedType; %by Huolin Xin
    encodedType = fread(FID,1,'uint32','ieee-be');
	
	%Find size of encoded type
	etSize = encodedTypeSize(encodedType);
	
%disp(['readAnyData, etSize:' num2str(etSize) ' ,Type:' num2str(encodedType)] )

	if etSize > 0
		%must be a regular data type, so read it and store a tag for it
		storeTag(curTagName, readNativeData(FID,encodedType));
	elseif encodedType == 18 %it is a String
		stringSize = fread(FID,1,'uint32','ieee-be');
		out = readStringData(FID,stringSize);
		storeTag(curTagName,out)
	elseif encodedType == 15 %it is a STRUCT
		%Stores fields in curly braces but does not store field names. In
		%fact the code will break for non-zero field names
		structTypes = readStructTypes(FID);
		out = readStructData(FID,structTypes);
		storeTag(curTagName,out);
	elseif encodedType == 20 %it is an ARRAY
		%Not read. Only stores a tag which it defined to indicate the size of the
		%data chunks that are skipped
		arrayTypes = readArrayTypes(FID);
		arrData = readArrayData(FID,arrayTypes);
		storeTag(curTagName,arrData);
	end
end

%Match the encodedtype to the Gatan type and get its size.
%Returns the size in bytes of the data type
function width = encodedTypeSize(encodedType)
	%Setup constants for the different encoded *data* types used in DM3 files
	%Structs, arrays, Strings, etc are handled elsewhere
    % 	VAL_SHORT   = 2;
    % 	VAL_LONG    = 3;
    % 	VAL_USHORT  = 4;
    % 	VAL_ULONG   = 5;
    % 	VAL_FLOAT   = 6;
    % 	VAL_DOUBLE  = 7;
    % 	VAL_BOOLEAN = 8;
    % 	VAL_CHAR    = 9;
    % 	VAL_OCTET   = 10;
	%   -1 will signal an unlisted type
	
    switch encodedType
		case {0} %just in case...
			width = 0;
		case {8, 9, 10}		%{VAL_BOOLEAN, VAL_CHAR, VAL_OCTET}
			width = 1; %data is 1 byte each
		case {2, 4}			%{VAL_SHORT, VAL_USHORT}
			width = 2;
		case {3, 5, 6}		%{VAL_LONG, VAL_ULONG, VAL_FLOAT}
			width = 4;
		case {7} %{VAL_DOUBLE}
			width = 8;
		otherwise
			%disp('encodedTypeSize: Type Unknown')
			width=-1;
    end
end

function rString = readStringData(FID, stringSize)
	
	%reads string data
 	if ( stringSize <= 0 )
		rString = '';
	else	
		rString = readString(FID, stringSize);
	end
end

%Function to read ordinary data types
function val = readNativeData(FID,encodedType)

	%reads ordinary data types
% 	VAL_SHORT   = 2;
% 	VAL_LONG    = 3;
% 	VAL_USHORT  = 4;
% 	VAL_ULONG   = 5;
% 	VAL_FLOAT   = 6;
% 	VAL_DOUBLE  = 7;
% 	VAL_BOOLEAN = 8;
% 	VAL_CHAR    = 9;
% 	VAL_OCTET   = 10;
	
%     VAL_UINT64 = 11;
    
	%These need to be read as little endian!
	if ( encodedType == 2 )
		val = fread(FID,1,'short');
	elseif ( encodedType == 3 )
		val = fread(FID,1,'int32');
	elseif ( encodedType == 4 )
		val = fread(FID,1,'uint16');
	elseif ( encodedType == 5 )
		val = fread(FID,1,'uint32');
	elseif ( encodedType == 6 )
		val = fread(FID,1,'float');
	elseif ( encodedType == 7 )
		val = fread(FID,1,'double');
	elseif ( encodedType == 8 )
		val = fread(FID,1,'uchar');
	elseif ( encodedType == 9 )
		val = fread(FID,1,'*char');
        %val = fscanf(FID,'%c',1);
	elseif ( encodedType == 10)
		val = fread(FID,1,'*char');
        %val = fscanf(FID,'%c',1);   % difference with char???
%     elseif(encodedType == 11)
%         val = fread(FID,1,'uint64');
	else
		error(['Unknown data type ' num2str(encodedType)])
	end

end

%%
%STRUCT and ARRAY Functions

%Analyzes the data types in a struct
function fieldTypes = readStructTypes(FID)

	fread(FID,1,'uint32','ieee-be'); %structNameLength
	nFields = fread(FID,1,'uint32','ieee-be');

	if ( nFields > 100 )
		error('Too many fields');
	end
		
	fieldTypes = zeros(1,nFields);
	for i=1:nFields
		nameLength = fread(FID,1,'uint32','ieee-be'); %nameLength
		fieldType = fread(FID,1,'uint32','ieee-be');
		fieldTypes(i) = fieldType;
	end

end

%Reads struct data based on type info in structType
function struct = readStructData(FID,structTypes)

	struct = zeros(1,length(structTypes));

	for i=1:length(structTypes)
		encodedType = structTypes(i);
		etSize = encodedTypeSize(encodedType);

		%disp(['Tag Type = ' num2str(encodedType) ', Tag Size = ' num2str(etSize)])

		%get data
		struct(i) = readNativeData(FID, encodedType);
	end
end

%Determines the data types in an array data type
function itemTypes = readArrayTypes(FID)
	
	arrayType = fread(FID,1,'uint32','ieee-be');
	
	itemTypes=[];
	
	if ( arrayType == 15 ) %STRUCT
		itemTypes = readStructTypes(FID);
	elseif ( arrayType == 20 ) %ARRAYS
		itemTypes = readArrayTypes(FID);
	else
		s = length(itemTypes);
		itemTypes(s+1) = arrayType;
	end

end

%Reads array data
function arrOut=readArrayData(FID,arrayTypes)
	
	global curTagName scale scaleUnit origin scale_temp origin_temp numObjects curGroupNameAtLevelX;

	arraySize = fread(FID,1,'uint32','ieee-be');

	itemSize = 0;
	encodedType = 0;
	
	for i=1:length(arrayTypes)
		encodedType = arrayTypes(i);
		etSize = encodedTypeSize(encodedType);
		itemSize = itemSize + etSize;
	end

	bufSize = arraySize * itemSize;
	loc = ftell(FID);
    
	if(encodedType == 4)
        storeTag([ curTagName '.Size'], bufSize);
		storeTag([ curTagName '.Offset'], loc);
        %Read in the character array
        stringData = fread(FID,bufSize);

        %Return only the characters as a character string
        arrOut = char(stringData(stringData>0))';
        
        %Catch unit names
        fullTagName = [curGroupNameAtLevelX '.' curTagName];
        if(~isempty(strfind(fullTagName,'Dimension')) && ~isempty(strfind(fullTagName,'Units')) && numObjects > 0 )
            tt = length(scale); %the number of dimension scales and units already saved
            scale(tt+1) = scale_temp;
            scaleUnit{tt+1} = arrOut;
            origin(tt+1) = origin_temp;
        end
 
	else
		%treat as binary data and skip array without reading it
		storeTag([ curTagName '.Size'], bufSize);
		storeTag([ curTagName '.Offset'], loc);
		fseek(FID, bufSize, 0);
		arrOut = ['Array data unread. Encoded type = ' num2str(encodedType)];
    end
end

%%
%Function to write out the tags to a TXT file and find interesting tags
%(Data size, offset, etc.)
function storeTag(curTagName,curTagValue)
	global FIDout dataSize dataOffset curGroupNameAtLevelX numObjects;
	global xsize ysize zsize dataType scale_temp origin_temp writeTagsGlobal;
	
	if ischar(curTagValue)
		totalTag = [curGroupNameAtLevelX '.' curTagName '=' curTagValue];
    else
		totalTag = [curGroupNameAtLevelX '.' curTagName ' = ' num2str(curTagValue)];
    end
    
	if strfind(curTagName,'Data.Size')
		numObjects = numObjects + 1; %add 1 to the number of objects
		dataSize(numObjects) = curTagValue;
	elseif strfind(curTagName,'Data.Offset')
		dataOffset(numObjects) = curTagValue;
	elseif strfind(curTagName,'DataType')
		dataType(numObjects) = curTagValue;
	elseif strfind(totalTag,'Dimensions.1')
		xsize(numObjects) = curTagValue;
        ysize(numObjects) = 1;
		zsize(numObjects) = 1;
	elseif strfind(totalTag,'Dimensions.2')
		ysize(numObjects) = curTagValue;
	elseif strfind(totalTag,'Dimensions.3')
		zsize(numObjects) = curTagValue;
	elseif (strfind(totalTag,'Dimension.') & strfind(totalTag,'.Scale'))
		scale_temp = curTagValue;
        %scale(numObjects) = curTagValue;
	elseif (strfind(totalTag,'Dimension.') & strfind(totalTag,'.Origin'))
		origin_temp = curTagValue;
        %origin(numObjects) = curTagValue;
    end
	
	if(writeTagsGlobal)
		fwrite(FIDout,totalTag);
		fprintf(FIDout,'\n');
	end
	
end