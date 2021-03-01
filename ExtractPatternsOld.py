import math

import numpy as np
import pandas as pd
from pyporcc import porcc


def CTInfoMaker(CTInfo, CTTemp):  # srise, sset
    # Store in CTInfo
    CTNum = CTTemp.CT[0]
    # Location
    # CTInfo.InFile[RowN] = CTTemp.InFile[0]
    # day/night
    # LAT = 55.576197
    # LONG = 9.848555
    # TODAY = CTTemp.Date[0]
    # [srise, sset, ~] = sunrise(LAT, LONG, 0, [], TODAY)
    #
    # if CTTemp.Date[0] >= srise and CTTemp.Date[0] < sset:
    #     CTInfo.DayNight[RowN] = 'Day'
    # else: # if DATE >= sset and DATE <= srise
    #     CTInfo.DayNight[RowN] = 'Night'
    # # end
    # Type
    Type = Species(CTTemp)
    CTInfo = CTInfo.append({'CTNum': CTNum, 'Date': CTTemp.datetime[0], 'Length': len(CTTemp), 'Species': Type,
                            'Behaviour': '-', 'Calf': '-', 'Notes': ' '}, ignore_index=True)
    return CTInfo


# end


def Species(CTTemp):
    CFDiff = CTTemp.CF.diff()
    PercChangeCF = (CFDiff / CTTemp.CF[0:-1 - 1]) * 100
    MedianPercChangeCF = abs(PercChangeCF).median()
    SDCF = CTTemp.CF.std()
    CPSDiff = CTTemp.CPS.diff()
    CPSDiff = CPSDiff.drop([0])
    PercChange = (CPSDiff / CTTemp.CPS[1:-1 - 1]) * 100
    MedianPercChangeCPS = abs(PercChange[1:-1]).median()
    if len(CTTemp) < 10:
        Type = 'Non-NBHF'
    elif MedianPercChangeCF < 0.5 and MedianPercChangeCPS < 0.05 and SDCF < 300:
        Type = 'Non-NBHF'
    elif MedianPercChangeCPS > 70 or MedianPercChangeCF > 4:
        Type = 'Non-NBHF'
    elif MedianPercChangeCPS < 30 or (MedianPercChangeCPS < 30 and MedianPercChangeCF > 2.65):
        Type = 'NBHF'
    else:
        Type = 'LQ-NBHF'
    # end
    return Type


# end


# def Behaviour(CPS, MeanPF, MeanPercChangeCF):
#     MedianCPS = CPS.median()
#     L = len(CPS)
#     SortedCPS = CPS.values.copy()
#     SortedCPS.sort()
#     CPS90Perc1 = SortedCPS[0:math.floor(0.90*L)]
#     CPS20 = CPS[0:math.floor(L*0.2)].mean()
#     CPS50 = CPS[math.floor(0.2*L):math.floor(0.8*L)].mean()
#     CPS90 = CPS[math.floor(0.8*L):-1].mean()
#     if MeanPF > 140000 and 8.5 > MedianCPS > 7.1 and MeanPercChangeCF < 0.5:
#         Behav = 'Sonar'
#     elif CPS90Perc1 < 100:
#         Behav = 'Orientation'
#     else:
#         CPS90Perc2 = SortedCPS[math.floor(0.10*L):-1]
#         if CPS90Perc2 > 100:
#             Behav = 'Socialising'
#         else:
#             BreakingPoint = find(CPS > 100)
#             DiffBP = BreakingPoint.diff()
#             BP = find(DiffBP == 1)
#             if len(BP) > 0:
#                 Pos = BreakingPoint(BP(1))
#                 if len(BP) > 0 and Pos > 5 and len(CPS) > Pos+5:
#                     Before = CPS[Pos-5:Pos].mean()
#                     After = CPS[Pos:Pos+5].mean()
#                 else:
#                     Before = 0
#                     After = 0
#                 # end
#                 if Before < 100 and After > 100:
#                     Behav = 'Foraging'
#                 else:
#                     if CPS90Perc1.mean() > 120:
#                         Behav = 'Socialising'
#                     elif CPS20 > 100 and CPS50 < 100 and CPS90 < CPS50:
#                         Behav = 'Socialising'
#                     else:
#                         Behav = 'Unknown'
#                     # end
#                 # end
#             else
#                 if mean(CPS90Perc1) > 150:
#                     Behav = 'Socialising'
#                 elif CPS20 > 100 and CPS50 < 100 and CPS90 < 100:
#                     Behav = 'Socialising'
#                 else:
#                     Behav = 'Unknown'
#                 # end
#             # end
#         # end
#     # end
# # end
#     return Behav


def NewICI(myTable, fs):
    start_sample = myTable["start_sample"]
    myTable['ICI'] = start_sample.diff() / (fs / 1000)
    myTable["CPS"] = 1000 / myTable.ICI
    myTable.at[0, 'CPS'] = 0.0
    myTable.at[0, 'ICI'] = 0.0
    return myTable


def ExtractPatterns(CP, Fs):
    CTNum = 0
    Keep = 0
    ColNames = ['CTNum', 'Date', 'Day_Night', 'Length', 'Species', 'Behaviour', 'Calf', 'Notes']
    CTInfo = pd.DataFrame(data=None, columns=ColNames)
    # add rown numbers to see if what we removed is another click train
    CP['RowN'] = range(0, len(CP))
    # Good click trains that only have echoes
    CP = NewICI(CP, Fs)
    # remove local min
    S1 = len(CP)
    S2 = 1
    while S1 != S2:
        S1 = len(CP)
        CP["AmpDiff"] = CP.amplitude.diff()
        CP = CP.drop(CP[(CP.amplitude.shift(1) > CP.amplitude) & (CP.amplitude.shift(-1) > CP.amplitude) & (
                CP.AmpDiff < -5)].index)
        CP.reset_index(inplace=True, drop=True)
        print(len(CP))
        CP = NewICI(CP, Fs)
        S2 = len(CP)

    CP = CP.drop(CP[(CP.CPS.diff() > 140.0)].index)
    CP.reset_index(inplace=True, drop=True)
    CP = NewICI(CP, Fs)
    ColNames = list(CP.columns)
    CTrains = pd.DataFrame(data=None, columns=ColNames)
    # # Find click trains
    TimeGaps = CP[(CP.ICI > 700.0) | (CP.ICI < 0.0)].index.to_series()
    TimeGaps.reset_index(inplace=True, drop=True)
    DiffTimeGaps = TimeGaps.diff()
    DiffTimeGaps.at[0] = TimeGaps[0] - 0
    # Find very long CT and reduce them to CT with fewer than 1000 clicks
    LongCTs = DiffTimeGaps[DiffTimeGaps > 1000].index

    if len(LongCTs) > 0:
        for m in range(0, len(LongCTs)):
            Length = DiffTimeGaps(LongCTs(m))
            CTsInside = math.floor(Length / 1000) + 1  # integer less than
            # Add Positions to TimeGaps
            PosInCP = TimeGaps(LongCTs(m))
            NextPos = TimeGaps(LongCTs(m) + 1)
            NewPos = 0  # [PosInCP:(Length/(CTsInside):NextPos
            TimeGaps[-1 + 1:-1 + CTsInside - 1] = round(NewPos[2:-1 - 1])
        # end
        TimeGaps = TimeGaps.sort()
        DiffTimeGaps = TimeGaps.diff()
        CTs = DiffTimeGaps[DiffTimeGaps > 9].index
        CTs.reset_index(inplace=True, drop=True)
    else:
        CTs = DiffTimeGaps[DiffTimeGaps > 9].index
        # CTs.reset_index(inplace=True, drop=True)
    # end

    for j in CTs:  # j runs through all the CT c
        print(j)
        if j == 0:
            Start = 0
            End = TimeGaps[j]
        else:
            Start = TimeGaps[j - 1]
            End = TimeGaps[j]
        CTTemp = CP[Start:End]
        CTTemp.reset_index(inplace=True, drop=True)
        CTTemp = NewICI(CTTemp, Fs)
        CTNum = CTNum + 1
        CTTemp['CT'] = CTNum
        # Delete loose clicks
        L1 = len(CTTemp)
        L2 = 1
        while L2 is not L1:
            L1 = len(CTTemp)
            CTTemp = NewICI(CTTemp, Fs)
            LooseClicks = CTTemp[CTTemp.ICI > 400].index.to_series()
            # Positions = find(diff(LooseClicks) == 1)
            RowsToDelete = LooseClicks[LooseClicks.diff() == 1]
            CTTemp = CTTemp.drop(RowsToDelete)
            L2 = len(CTTemp)
        # end
        CTTemp.reset_index(inplace=True, drop=True)
        # A large difference indicates echoes
        SortedCPS = CTTemp.CPS.values.copy()
        SortedCPS.sort()
        DiffSorted = pd.DataFrame(SortedCPS).diff()
        DiffSorted = DiffSorted.drop([0])
        MaxDiffSorted = DiffSorted.max().values
        # ExploreCTTemp(CTTemp)

        if MaxDiffSorted <= 40:
            NewCT = CTTemp
        else:
            CTTemp = CTTemp.drop(CTTemp[CTTemp.CPS >= int(MaxDiffSorted) - 30].index)
            CTTemp = NewICI(CTTemp, Fs)
            SortedCPS = CTTemp.CPS.values.copy()
            SortedCPS.sort()
            DiffSorted = pd.DataFrame(SortedCPS).diff()
            DiffSorted = DiffSorted.drop([0])
            MaxDiffSorted = DiffSorted.max().values
            # ExploreCTTemp(CTTemp)
            if MaxDiffSorted <= 50:
                NewCT = CTTemp.copy()
            elif len(CTTemp) > 20:
                # Finding the stable areas
                CPSDiff = CTTemp.CPS.diff()
                PercChangeS = (CPSDiff / CTTemp.CPS[1:-1 - 1]) * 100
                PercChangeS = abs(PercChangeS[2:-1])
                window_size = 3
                i = 0
                PercCPSDiff = []
                while i < len(PercChangeS) - window_size + 1:
                    this_window = PercChangeS[i: i + window_size]
                    # get current window
                    window_average = sum(this_window) / window_size
                    PercCPSDiff.append(window_average)
                    i += 1
                PercCPSDiff = pd.DataFrame(PercCPSDiff)
                PercCPSDiff.reset_index(inplace=True, drop=True)
                print(PercCPSDiff)
                StartRow = PercCPSDiff[PercCPSDiff[0] < 20.0].index.to_series()
                StartRow.reset_index(inplace=True, drop=True)
                DiffStartRow = StartRow.diff()
                Here = DiffStartRow[DiffStartRow == 1].index.to_series()
                Here.reset_index(inplace=True, drop=True)
                if len(StartRow) < 2:
                    NewCT = CTTemp
                    NewCT = NewICI(NewCT, Fs)
                else:  # go into the CT
                    RowN = StartRow[0]  # Low variability in CPS(for the next 4)
                    RowsToKeep = pd.DataFrame(data=None, columns=None)  # # ok<FNDSB>
                    RowsToKeep[0] = Here
                    FirstRowN = RowN

                    # Look backwards
                    while RowN > 5:
                        ClickCPS = CTTemp.CPS.iloc[RowN]
                        ClickAmp = CTTemp.amplitude.iloc[RowN]
                        Clickstart_sample = CTTemp.start_sample.iloc[RowN]
                        ICIs = abs(Clickstart_sample - CTTemp.start_sample[RowN - 5:RowN - 1]) / (Fs / 1000)
                        CPSs = 1000 / ICIs
                        Amps = CTTemp.amplitude[RowN - 5:RowN - 1]
                        Amps = abs(ClickAmp - Amps)
                        DiffCPSs = abs(ClickCPS - CPSs)
                        DiffCPSs = pd.DataFrame(DiffCPSs)
                        MinVal = DiffCPSs.start_sample.min()
                        ixCPS = DiffCPSs[DiffCPSs.start_sample == MinVal].index
                        if Amps[ixCPS[0]] < 5:
                            DiffCPSs[ixCPS[0]] = 1000  # high arbitrary numnber
                            MinVal = DiffCPSs.start_sample.min()
                            ixCPS = DiffCPSs[DiffCPSs.start_sample == MinVal].index
                            RowN = RowN - ixCPS[0]
                        else:
                            RowN = RowN - ixCPS[0]
                        # end
                        Keep = Keep + 1
                        RowsToKeep[Keep] = RowN
                    # end

                # Look forwards
                RowN = FirstRowN
                while RowN < len(CTTemp) - 10:
                    ClickCPS = CTTemp.CPS.iloc[RowN]
                    ClickAmp = CTTemp.amplitude.iloc[RowN]
                    Clickstart_sample = CTTemp.start_sample.iloc[RowN]
                    ICIs = abs(CTTemp.start_sample[RowN + 1:RowN + 9] - Clickstart_sample) / (Fs / 1000)
                    CPSs = 1000 / ICIs
                    Amps = CTTemp.amplitude[RowN + 1:RowN + 9]
                    Amps = abs(Amps - ClickAmp)
                    DiffCPSs = abs(ClickCPS - CPSs)
                    DiffCPSs = pd.DataFrame(DiffCPSs)
                    MinVal = DiffCPSs.start_sample.min()
                    ixCPS = DiffCPSs[DiffCPSs.start_sample == MinVal].index
                    if Amps[ixCPS[0]] < 5:
                        DiffCPSs[ixCPS[0]] = 1000  # high arbitrary numnber
                        MinVal = DiffCPSs.start_sample.min()
                        ixCPS = DiffCPSs[DiffCPSs.start_sample == MinVal].index
                        RowN = RowN + ixCPS[0]
                    else:
                        RowN = RowN + ixCPS[0]
                    # end
                    Keep = Keep + 1
                    RowsToKeep[Keep] = RowN
                    # end
                    Keep = Keep + 1
                    RowsToKeep[Keep] = RowN
                # end
                # adding the first 10 and the last 10 rows (although might not be needed)

                # RowsToKeep.append(CTTemp[0:10].index)
                # RowsToKeep.append(CTTemp[-1-9:-1].index)
                RowsToKeep = np.unique(RowsToKeep)
                RowsToKeep = np.delete(RowsToKeep, np.where(RowsToKeep <= 0))
                RowsToKeep = np.delete(RowsToKeep, np.where(RowsToKeep > len(CTTemp)))
                print(RowsToKeep)
                NewCT = CTTemp[RowsToKeep]
                NewCT = NewICI(NewCT, Fs)
                # end
            else:
                NewCT = CTTemp
                NewCT = NewICI(NewCT, Fs)
            # end
        # end
        CTInfo = CTInfoMaker(CTInfo, CTTemp)
        print(CTInfo)
        # Put result into a final file
        if j == 0:
            CTrains = NewCT
        else:
            CTrains.append(NewCT)
        # end
    # end
    return CTrains, CTInfo, CP


LQ = 0.6
HQ = 0.999999
Fs = 576000
classifier = porcc.PorCC(load_type='manual', config_file='default')
classifier.th1 = HQ
classifier.th2 = LQ

Clicks = pd.read_pickle('C:/Mel/CPODvsDPorCCA/CorrecteddB/Detected_Clips_110815_230428.pkl')
Clicks = classifier.classify_matrix(Clicks)
Clicks = Clicks.drop(Clicks[Clicks.pyPorCC == 3][:].index)
Clicks = Clicks.drop(['wave'], axis=1)
Clicks = Clicks.drop(Clicks[Clicks.duration_us > 500][:].index)
CP = Clicks.copy()
CTrains, CTInfo, CP = ExtractPatterns(CP, Fs)
CTInfo.to_csv('C:/Mel/CPODvsDPorCCA/CorrecteddB/CTInfo.csv')
CTrains.to_csv('C:/Mel/CPODvsDPorCCA/CorrecteddB/CTrains.csv')
CP.to_csv('C:/Mel/CPODvsDPorCCA/CorrecteddB/CP.csv')

# end
