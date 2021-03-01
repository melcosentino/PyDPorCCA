##PORCC - STEPS INCLUDE
# 1) Getting the necessary values for PorCC
# 2) Accessing individual clicks in.pgdf files
# 3) Estimating probabilities and all other parameters
# 4) Saving if classified as LQ or HQ
# 5) Remove clicks with the same date keeping the best one

##################################
#### PAMGuard functions
#################################
def findBinaryFiles(folder, containsName):
    # FINDBINARYFILES finds all binary menufile paths within a folder and subfolders
    #  [FILENAMES] = FINDBINARYFILES(FOLDER) finds all.pgdf files within a
    #    folder and it's sub folders. The def returns a cell array of
    #    FILENAMES.
    #    [FILENAMES] = FINDBINARYFILES(FOLDER, CONTAINSNAME) finds all .pgdf files
    #    within a folder that all also contain a string CONTAINSNAME e.g.this
    #    could be 'clicks' which would find all click binary files.

    # if nargin & lt 2
    containsName = []
    # end

    filenames = []
    subFiles = dir(folder)
    for i in subFiles:
        if (strcmp(subFiles(i).name, '.') == 1 or strcmp(subFiles(i).name, '..') == 1)
            continue
        # end

        if (subFiles(i).isdir == 1)
            subFolderName = [folder, '\', subFiles(i).name]
                             filenames = cat(2, filenames, findBinaryFiles(subFolderName, containsName))
            else
            binaryFileName = [folder, '\',subFiles(i).name]
                              binaryChar = char(binaryFileName)
            fileEnd = binaryChar(length(binaryChar) - 3: length(binaryChar))
            # add to binary menufile list
            if (strcmp(fileEnd, 'pgdf') == 1):
                if
            (isempty(containsName))
            # if no contansName specified then add to list
            filenames = cat(2, filenames, binaryFileName(1,:))
            elif (~isempty(strfind(subFiles(i).name, containsName))):
            # if a contansName is specified only add menufile if it contains
            # the specified string
            filenames = cat(2, filenames, binaryFileName(1,:))
            # end
            # end if
            # end

            # end for
            return filenames

    # end

    def header =

    readFileHeader(file, readExtra):
    if nargin < 2
        readExtra = 0
    # end
    header.length = fread(file, 1, 'int32')
    header.identifier = fread(file, 1, 'int32')
    header.fileFormat = fread(file, 1, 'int32')
    header.pamguard = char(fread(file, 12, 'uchar')
    ')
    header.version = readJavaUTFString(file)
    header.branch = readJavaUTFString(file)
    header.dataDate = millisToDateNum(fread(file, 1, 'int64'))
    header.analysisDate = millisToDateNum(fread(file, 1, 'int64'))
    header.startSample = fread(file, 1, 'int64')
    header.moduleType = readJavaUTFString(file)
    header.moduleName = readJavaUTFString(file)
    header.streamName = readJavaUTFString(file)
    header.extraInfoLen = fread(file, 1, 'int32')
    if readExtra:
        header.extraInfo = fread(file, header.extraInfoLen, 'int8')
    else:
        fseek(file, header.extraInfoLen, 'cof')
        # end
        # end


def readJavaUTFString(file):
    # def[str len] = readJavaUTFString(menufile)
    # read a string written in Java UTF-8 format from a menufile.
    # The first two bytes are the length of the string, then
    # it's the string.
    len2 = fread(file, 1, 'int16')
    str2 = char(fread(file, len2, 'uchar')
    ')
    return (str2, len2)


# end

def millisToDateNum(millis):
    datenum = double(millis) / 86400000.0 + 719529
    return datenum


# end

def readStdModuleHeader(file):
    # reads the module header information common to all modules.Differs from
    # the legacy code in that it does not read in or skip any information
    # specific to a module.
    header.length = fread(file, 1, 'int32')
    header.identifier = fread(file, 1, 'int32')
    header.version = fread(file, 1, 'int32')
    header.binaryLength = fread(file, 1, 'int32')
    return header


# end

def [data


hasAnnotation] = readPamData(fid, fileInfo):
# Reads in the object data that is common to all modules.This reads up to
# (but not including) the object binary length, and then calls a def
# to read the module-specific data.
# Inputs: fid = menufile identifier
# fileInfo = structure holding the menufile header, module header, a handle
# to the function to read module - specific data, etc.
# Output:
# data = structure containing data from a single object
# set constants to match flag bitmap constants in class
# DataUnitBaseData.java. The following contstants match header version 4.
TIMEMILLIS = hex2dec('1')
TIMENANOS = hex2dec('2')
CHANNELMAP = hex2dec('4')
UID = hex2dec('8')
STARTSAMPLE = hex2dec('10')
SAMPLEDURATION = hex2dec('20')
FREQUENCYLIMITS = hex2dec('40')
MILLISDURATION = hex2dec('80')
TIMEDELAYSSECS = hex2dec('100')
HASBINARYANNOTATIONS = hex2dec('200')

# initialize a new variable to hold the data
data = []
data.flagBitmap = 0
hasAnnotation = 0

# calculate where the next object starts, in case there is an error trying
# to read this one
objectLen = fread(fid, 1, 'int32')
curObj = ftell(fid)
nextObj = curObj + objectLen

# first thing to check is that this is really the type of object we think
# it should be, based on the menufile header. If not, warn the user, move the
# pointer to the next object, and exit
data.identifier = fread(fid, 1, 'int32')
if (any(data.identifier == fileInfo.objectType)):
# do nothing here - couldn 't figure out a clean way of checking if
# number wasn't in array
else
    disp(['Error - Object Identifier does not match ' fileInfo.fileHeader.moduleType ' type.  Aborting data read.'])
    fseek(fid, nextObj, 'bof')
    return
# end

# read the data, starting with the standard data that every data unit has
# version = fileInfo.fileHeader.fileFormat
try
    data.millis = fread(fid, 1, 'int64')
    if (version == 2 or (bitand(data.flagBitmap, TIMENANOS)~=0)):
        data.timeNanos = fread(fid, 1, 'int64')
    # end
    if (bitand(data.flagBitmap, UID) == UID):
        data.UID = fread(fid, 1, 'int64')
    # end
    if (bitand(data.flagBitmap, STARTSAMPLE)~=0):
        data.startSample = fread(fid, 1, 'int64')
    # end
    if (bitand(data.flagBitmap, MILLISDURATION)~=0):
        data.millisDuration = fread(fid, 1, 'float')
    # end
    if (bitand(data.flagBitmap, TIMEDELAYSSECS)~=0):
        data.numTimeDelays = fread(fid, 1, 'int16')
        td = zeros(1, data.numTimeDelays)
        for i = 1:data.numTimeDelays:
        td(i) = fread(fid, 1, 'float')
# end
data.timeDelays = td
# end try

# set date, to maintain backwards compatibility
data.date = millisToDateNum(data.millis)
# now read the module - specific data
if (isa(fileInfo.readModuleData, 'def_handle')):
    [data, error] = fileInfo.readModuleData(fid, fileInfo, data)
    if (error):
        disp(['Error - cannot retrieve ' fileInfo.fileHeader.moduleType ' data properly.'])
        fseek(fid, nextObj, 'bof')
        return
    # end
# end
# now check to see if there are standard annotations to the main data.
if (bitand(data.flagBitmap, HASBINARYANNOTATIONS)~=0):
    hasAnnotation = 1
    anStart = ftell(fid)
    anTotLength = fread(fid, 1, 'int16')
    nAn = fread(fid, 1, 'int16')
    for i in nAn:
        filePos = ftell(fid)
        anLength = fread(fid, 1, 'int16') - 2  # this length does no include itself !
        anId = readJavaUTFString(fid)
        anVersion = fread(fid, 1, 'int16')
        switch(anId)
        case
        'Beer'
        data.annotations.beamAngles = readBeamFormerAnnotation(fid, anId, anLength, fileInfo, anVersion)
    case
    'TMAN'
    data.annotations.targetMotion = readTMAnnotation(fid, anId, anLength, fileInfo, anVersion)
case
'TDBL'
data.annotations.toadAngles = readTDBLAnnotation(fid, anId, anLength, fileInfo, anVersion)
otherwise
fprintf('Unknown anotation type "%s" length %d version %d in file\n', ...
anId, anLength, anVersion)
fseek(fid, filePos + anLength, 'bof')
# end switch
endPos = ftell(fid)
if (endPos ~= filePos+anLength):
    disp('Possible annotation read size error in file')
    fseek(fid, filePos + anLength, 'bof')
    endPos = ftell(fid)
# end if
# end for
if (endPos ~= anStart+anTotLength):
    fseek(fid, anStart + anTotLength, 'bof')
# end
else:
data.annotations = []
# end

catch
mError
disp('Error loading object data')
disp(data)
disp(getReport(mError))
fseek(fid, nextObj, 'bof')
# end try
# end function

def checkArrayAllocation(array, reqLength, sampleObject):
    # global array
    if isempty(array):
        currentLength = 0
        array = []
    else:
        currentLength = numel(array)
    # end
    if (currentLength >= reqLength):
        return
    # end
    allocStep = round(sqrt(reqLength))
    allocStep = max(10, min(allocStep, 10000))
    array(reqLength + allocStep) = sampleObject
    return
    return array


# end function

def readClickData(fid, fileInfo, data)
    error = false
    try
        # read click detector specific data
        dataLength = fread(fid, 1, 'int32')
        if (dataLength == 0)
            return
        # end

        if (fileInfo.moduleHeader.version <= 3)
            data.startSample = fread(fid, 1, 'int64')
            data.channelMap = fread(fid, 1, 'int32')
        # end

        data.triggerMap = fread(fid, 1, 'int32')
        data.type = fread(fid, 1, 'int16')
        if (fileInfo.moduleHeader.version >= 2)
            data.flags = fread(fid, 1, 'int32')
        else
            data.flags = 0
        # end

        if (fileInfo.moduleHeader.version <= 3)
            nDelays = fread(fid, 1, 'int16')
            if (nDelays)
                data.delays = fread(fid, nDelays, 'float')
            # end
        # end

        nAngles = fread(fid, 1, 'int16')
        if (nAngles)
            data.angles = fread(fid, nAngles, 'float')
        # end

        if fileInfo.moduleHeader.version >= 3:
            nAngleErrors = fread(fid, 1, 'int16')
            data.angleErrors = fread(fid, nAngleErrors, 'float')
        else:
            data.angleErrors = []
        # end

        if fileInfo.moduleHeader.version <= 3:
            data.duration = fread(fid, 1, 'int16')
        else
            data.duration = data.sampleDuration  # duplicate
            # the value to maintain backwards compatibility
        # end

        data.nChan = countChannels(data.channelMap)
        maxVal = fread(fid, 1, 'float')
        data.wave = fread(fid, [data.duration, data.nChan], 'int8') * maxVal / 127

    catch
    mError
    disp(['Error reading ' fileInfo.fileHeader.moduleType '  data object.  Data read:'])
    disp(data)
    disp(getReport(mError))
    error = true
    return (data, error)
    # end try


# end function

def countChannels(channelMap)
    # count the numebr of set bits in the channel map
    nC = 0
    j = 1
    for i in 32:
        if (bitand(channelMap, j))
            nC = nC + 1
        # end
        j = j * 2
    # end
    return nC


# end

def readClickFooter(file):
    # reads module footer information for the Click Detector module. Note that
    # sometimes there is no additional footer information, so check first
    # whether or not the binaryLength variable is 0.
    footer = readStdModuleFooter(file)
    if (footer.binaryLength ~= 0)
        footer.typesCountLength = fread(file, 1, 'int16')
        footer.typesCount = fread(file, footer.typesCountLength, 'int32')
    # end
    return footer


# end function

def readStdModuleFooter(file):
    # reads the module footer information common to all modules.Differs from
    # the legacy code in that it does not read in or skip any information
    # specific to a module.
    footer.length = fread(file, 1, 'int32')
    footer.identifier = fread(file, 1, 'int32')
    footer.binaryLength = fread(file, 1, 'int32')
    return footer


# end

def readFileFooterInfo(fid, version)
    footer.length = fread(fid, 1, 'int32')
    footer.identifier = fread(fid, 1, 'int32')
    footer.nObjects = fread(fid, 1, 'int32')
    footer.dataDate = millisToDateNum(fread(fid, 1, 'int64'))
    footer.analysisDate = millisToDateNum(fread(fid, 1, 'int64'))
    footer.endSample = fread(fid, 1, 'int64')
    if version >= 3:
        footer.lowestUID = fread(fid, 1, 'int64')
        footer.highestUID = fread(fid, 1, 'int64')
    # end if version
    footer.fileLength = fread(fid, 1, 'int64')
    footer.endReason = fread(fid, 1, 'int32')
    return footer
    # end if
# end function
