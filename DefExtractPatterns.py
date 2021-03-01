## EXTRACT PATTERNS IN  CT
import pandas as pd

CPar = pd.read_csv('C:/Mel/CPODvsDPorCCA/CorrecteddB/Detected_Clicks_110815_230428.csv')
# FileName = SelectedFolder + "/ClickParameters.csv"
# CPar = pd.read_csv(FileName)

# Parameters
global MinLength
global Fs
MinLength = 16
SamFq = 576000


def NewICI(myTable, fs):
    StartSample = myTable["start_sample"]
    myTable['ICI'] = StartSample.diff() / (fs / 1000)
    myTable["CPS"] = 1000 / myTable["ICI"]
    myTable.at[0, 'CPS'] = 0.0
    myTable.at[0, 'ICI'] = 0.0
    return myTable


def FindCT(self, CP, Fs):
    # add rown numbers to see if what we removed is another click train
    self.CP["RowN"] = range(0, len(CP))
    # Save in another variable
    CPOrig = self.CP
    CP = NewICI(CP, Fs)

    # Remove outliers
    CP["AmpDiff"] = CP.amplitude.diff()
    CP = CP.drop(CP[(CP.ICI < 2) & (CP.AmpDiff < -5)].index)
    CP.reset_index(inplace=True, drop=True)
    CP = NewICI(CP, Fs)
    # remove local min
    S1 = len(CP)
    S2 = 1
    while S1 != S2:
        S1 = len(CP)
        CP["AmpDiff"] = CP.amplitude.diff()
        CP = CP.drop(CP[(CP.amplitude.shift(1) > CP.amplitude) & (CP.amplitude.shift(-1) > CP.amplitude) & (
                    CP.AmpDiff < -5)].index)
        print(len(CP))
        CP = NewICI(CP, Fs)
        S2 = len(CP)

    CP = CP.drop(CP[(CP.CPS.diff() > 140.0)].index)
    CP = NewICI(CP, Fs)

    ## Find click trains
    # Bit by bit
    CP.reset_index(drop=True, inplace=True)
    CPSDiff = CP.CPS.diff()
    CP.CPS.at[0] = 0.1
    CPSDiff.at[0] = 0.1
    PercChangeS = (CPSDiff / CP.CPS) * 100
    PercChangeS = abs(PercChangeS[0:-1])
    PercCPSDiff = pd.DataFrame(columns=["MeanAvg"])
    PercCPSDiff["MeanAvg"] = PercChangeS.rolling(window=4).mean()
    PercCPSDiff["Row"] = PercCPSDiff.index
    StartRowTemp = PercCPSDiff["Row"][PercCPSDiff["MeanAvg"] < 5]
    StartRow = StartRowTemp
    DiffStartRow = StartRow.diff()

    # go into the click trains
    RowsToKeep = []
    Keep = 1
    if len(StartRow) == 0:
        RowN = 1
    else:
        RowN = StartRow(Here(1))  # Low variability in CPS (for the next 4)
        RowsToKeep[1, 1] = RowN
        FirstRowN = RowN

        ## Find gaps in StartRow
        DiffStartRow = StartRow.diff()
        Gaps = find(DiffStartRow > 25)
        CTInfo = pd.DataFrame()
        EndAll = []
        StartAll = []
        for g = 1:len(Gaps, 1)
        RowsToKeep = []
        if g == 1:
            Start = 1
            End = StartRow[Gaps[g] - 1] + 30
        else:
            Start = StartRow[Gaps[g - 1] + 1 - 30
            End = StartRow[Gaps[g]] + 30
            EndAll[end + 1, 1] = End
            StartAll[end + 1, 1] = Start

            CTTemp = CP[Start:End, :]
            if len(CTTemp, 1) > 6900 or len(CTTemp, 1) < 10:
            # do nothing
                else:
            # Look backwards
            RowN = StartRow[Gaps[g]] - Start - 30
        while RowN > 10
            ClickCPS = CTTemp.CPS[RowN, 1]
            ClickAmp = CTTemp.Amp[RowN, 1]
            ClickStartSample = CTTemp.startSample[RowN, 1]
            ICIs = abs(ClickStartSample - CTTemp.startSample[RowN - 9:RowN - 1, 1]] / (Fs / 1000)
            CPSs = 1000. / ICIs
            Amps = CTTemp.Amp[RowN - 9:RowN - 1, 1]
            Amps = abs(ClickAmp - Amps)
            DiffCPSs = abs(ClickCPS - CPSs)
            [~, ixCPS] = min(DiffCPSs)
            if Amps[ixCPS] < 5:
                DiffCPSs[ixCPS] = 1000  # high arbitrary numnber
                [~, ixCPS] = min(DiffCPSs)
                RowN = RowN - ixCPS
            else:
                RowN = RowN - ixCPS
        Keep = Keep + 1
        RowsToKeep[Keep, 1] = RowN

    # Look forwards
    RowN = FirstRowN
    while RowN < len(CTTemp, 1) - 10
        ClickCPS = CTTemp.CPS[RowN, 1]
        ClickAmp = CTTemp.Amp[RowN, 1]
        ClickStartSample = CTTemp.startSample[RowN, 1]
        ICIs = abs(CTTemp.startSample[RowN + 1:RowN + 9, 1] - ClickStartSample) / (Fs / 1000)
        CPSs = 1000. / ICIs
        Amps = CTTemp.Amp[RowN + 1:RowN + 9, 1]
        Amps = abs(Amps - ClickAmp)
        DiffCPSs = abs(CPSs - ClickCPS)
        [~, ixCPS] = min(DiffCPSs)
        if Amps[ixCPS] < 5:
            DiffCPSs[ixCPS] = 1000
            [~, ixCPS] = min[DiffCPSs)
            RowN = RowN + ixCPS
        else:
            RowN = RowN + ixCPS
    Keep = Keep + 1
    RowsToKeep[Keep, 1] = RowN


RowsToKeep[end + 1:end + 20, 1] = [1:10, len(CTTemp, 1) - 9: len(CTTemp, 1)]
RowsToKeep = sort(unique(RowsToKeep))
RowsToKeep[RowsToKeep == 0] = []

# Delete loose clicks
L1 = len(CTTemp, 1)
L2 = 1
while L2 != L1
    L1 = len(CTTemp, 1)
    CTTemp.ICI[1, 1] = 0
    CTTemp.ICI[2:end, 1] = diff(CTTemp.startSample) / (Fs / 1000)
    CTTemp.CPS[1:len(CTTemp, 1), 1) = 1000. / CTTemp.ICI[1:end, 1)
    CTTemp.CPS[1, 1] = 0
    LooseClicks = find(CTTemp.ICI > 250)
    Positions = find(diff(LooseClicks) == 1)
    RowsToDelete = LooseClicks(Positions)
    CTTemp[RowsToDelete, :] = []
    L2 = len(CTTemp, 1)

if len(RowsToKeep, 1) > 1:
    RowsToKeep[RowsToKeep > len(CTTemp, 1)) = []
    clear
    NewCT
    NewCT = CTTemp[RowsToKeep, :]
    NewCT.CT[1:end, 1] = g
    NewCT = NewICI(NewCT, Fs)
    if g == 1:
        ClickTrains = NewCT
    else:
        ClickTrains[end + 1:end + len(NewCT, 1), :] = NewCT

    CTNum = NewCT.CT[1, 1]
    CheckEmpty = isfield(CTInfo, 'CTNum')
    if CheckEmpty == 0:
        RowNCT = 1
    else:
        RowNCT = RowNCT + 1

    # Store in CTInfo
    CTInfo[RowNCT].CTNum = CTNum  # new CT to separate
    CTInfo[RowNCT].Date = datestr(NewCT.Date[1, 1])  # Start date
    # Location
    CTInfo[RowNCT].InFile = NewCT.InFile[1, 1]
    # Length
    CTInfo[RowNCT].Length = len(NewCT, 1)
    # day/night
    if NewCT.Date[1, 1] >= srise and NewCT.Date[1, 1] < sset:
        CTInfo[RowNCT].DayNight = 'Day'
    else:  # if DATE >= sset and DATE <= srise
        CTInfo[RowNCT].DayNight = 'Night'

    # Type
    CTInfo[RowNCT].Species = '-'
    CTInfo[RowNCT].Behav = '-'
    CTInfo[RowNCT].Calf = '-'
    CTInfo[RowNCT].Saved = 0
    CTInfo[RowNCT].Notes = ' '
else:
# nothing

## SAVE
ClickParNames = strcat(thisFolder, '\ClickParametersTested.mat')
save(ClickParNames, 'CP')
CTInfoNames = strcat(thisFolder, '\CTInfo.mat')
save(CTInfoNames, 'CTInfo')
CTInfoNames = strcat(thisFolder, '\ClickTrains.mat')
save(CTInfoNames, 'ClickTrains')
