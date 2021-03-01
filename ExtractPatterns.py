import math
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def forceRange(v, max):
    # force v to be >= 0 and < max
    if v < 0:
        return v + max
    elif v >= max:
        return v - max
    return v


def getSunriseTime(day, month, year, longitude, latitude):
    return calcSunTime(day, month, year, longitude, latitude, True)


def getSunsetTime(day, month, year, longitude, latitude):
    return calcSunTime(day, month, year, longitude, latitude, False)


def calcSunTime(day, month, year, longitude, latitude, isRiseTime, zenith=90.8):
    TO_RAD = math.pi / 180

    # 1. first calculate the day of the year
    N1 = math.floor(275 * month / 9)
    N2 = math.floor((month + 9) / 12)
    N3 = (1 + math.floor((year - 4 * math.floor(year / 4) + 2) / 3))
    N = N1 - (N2 * N3) + day - 30

    # 2. convert the longitude to hour value and calculate an approximate time
    lngHour = longitude / 15

    if isRiseTime:
        t = N + ((6 - lngHour) / 24)
    else:  # sunset
        t = N + ((18 - lngHour) / 24)

    # 3. calculate the Sun's mean anomaly
    M = (0.9856 * t) - 3.289

    # 4. calculate the Sun's true longitude
    L = M + (1.916 * math.sin(TO_RAD * M)) + (0.020 * math.sin(TO_RAD * 2 * M)) + 282.634
    L = forceRange(L, 360)  # NOTE: L adjusted into the range [0,360)

    # 5a. calculate the Sun's right ascension

    RA = (1 / TO_RAD) * math.atan(0.91764 * math.tan(TO_RAD * L))
    RA = forceRange(RA, 360)  # NOTE: RA adjusted into the range [0,360)

    # 5b. right ascension value needs to be in the same quadrant as L
    Lquadrant = math.floor(L / 90) * 90
    RAquadrant = math.floor(RA / 90) * 90
    RA = RA + (Lquadrant - RAquadrant)

    # 5c. right ascension value needs to be converted into hours
    RA = RA / 15

    # 6. calculate the Sun's declination
    sinDec = 0.39782 * math.sin(TO_RAD * L)
    cosDec = math.cos(math.asin(sinDec))

    # 7a. calculate the Sun's local hour angle
    cosH = (math.cos(TO_RAD * zenith) - (sinDec * math.sin(TO_RAD * latitude))) / (cosDec * math.cos(TO_RAD * latitude))

    if cosH > 1:
        return {'status': False, 'msg': 'the sun never rises on this location (on the specified date)'}

    if cosH < -1:
        return {'status': False, 'msg': 'the sun never sets on this location (on the specified date)'}

    # 7b. finish calculating H and convert into hours

    if isRiseTime:
        H = 360 - (1 / TO_RAD) * math.acos(cosH)
    else:  # setting
        H = (1 / TO_RAD) * math.acos(cosH)

    H = H / 15

    # 8. calculate local mean time of rising/setting
    T = H + RA - (0.06571 * t) - 6.622

    # 9. adjust back to UTC
    UT = T - lngHour
    UT = forceRange(UT, 24)  # UTC time in decimal format (e.g. 23.23)

    # 10. Return
    hr = forceRange(int(UT), 24)
    Min = round((UT - int(UT)) * 60, 0)
    return hr, Min


def CTInfoMaker(CTInfo, CTTemp):  # srise, sset
    # Store in CTInfo
    CTNum = CTTemp.CT[0]
    # Location
    # CTInfo.InFile[RowN] = CTTemp.InFile[0]
    # day/night
    latitude = 55.576197  # Will need to get it from the menu
    longitude = 9.848555  # Will need to get it from the menu
    TODAY = CTTemp.datetime.iloc[0]
    day = int(TODAY[8:10])
    month = int(TODAY[5:7])
    year = int(TODAY[0:4])
    sriseH, sriseM = getSunriseTime(day, month, year, longitude, latitude)
    ssetH, ssetM = getSunsetTime(day, month, year, longitude, latitude)
    # I don't know which format time is returned here, need to correct when I do
    HH = TODAY[11:13]
    MM = TODAY[14:16]

    if (int(HH) > ssetH and int(MM) > ssetM) and (int(HH) < sriseH and int(MM) < sriseM):
        DayNight = 'Day'
    else:  # if DATE >= sset and DATE <= srise
        DayNight = 'Night'
    # end
    # Type
    Type = Species(CTTemp)
    if Type == 'Non-NBHF':
        Behav = '-'
    else:
        Behav = Behaviour(CTTemp)
    CTInfo = CTInfo.append({'CTNum': CTNum, 'Date': CTTemp.datetime[0], 'Length': len(CTTemp), 'Species': Type,
                            'DayNight': DayNight, 'Behaviour': Behav, 'Calf': '-', 'Notes': ' '}, ignore_index=True)
    CTInfo.to_csv('C:/Mel/DPorCCA_FROM_SCRATCH/PGDF_SoundTrapNewTest/20150812/CTInfo.csv', index=False)

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


def Behaviour(CTTemp):
    CFDiff = CTTemp.CF.diff()
    PercChangeCF = (CFDiff / CTTemp.CF[0:-1 - 1]) * 100
    MeanPF = CTTemp.CF.mean()
    MeanPercChangeCF = abs(PercChangeCF).mean()
    CPS = CTTemp.CPS
    MedianCPS = CPS.median()
    L = len(CPS)
    SortedCPS = CPS.values.copy()
    SortedCPS.sort()
    CPS90Perc1 = SortedCPS[0:math.floor(0.90 * L)]
    CPS20 = CPS[0:math.floor(L * 0.2)].mean()
    CPS50 = CPS[math.floor(0.2 * L):math.floor(0.8 * L)].mean()
    CPS90 = CPS[math.floor(0.8 * L):-1].mean()
    if MeanPF > 140000 and 8.5 > MedianCPS > 7.1 and MeanPercChangeCF < 0.5:
        Behav = 'Sonar'
    elif all(CPS90Perc1 < 100):
        Behav = 'Orientation'
    else:
        CPS90Perc2 = SortedCPS[math.floor(0.10 * L):-1]
        if all(CPS90Perc2 > 100):
            Behav = 'Socialising'
        else:
            BreakingPoint = CPS[CPS > 100].index.to_series()
            DiffBP = BreakingPoint.diff()
            BP = DiffBP[DiffBP == 1].index
            if len(BP) > 0:
                Pos = BreakingPoint[BP[0]]
                if len(BP) > 0 and Pos > 5 and len(CPS) > Pos + 5:
                    Before = CPS[Pos - 5:Pos].mean()
                    After = CPS[Pos:Pos + 5].mean()
                else:
                    Before = 0
                    After = 0
                # end
                if Before < 100 and After > 100:
                    Behav = 'Foraging'
                else:
                    if CPS90Perc1.mean() > 120:
                        Behav = 'Socialising'
                    elif CPS20 > 100 and 100 > CPS50 > CPS90:
                        Behav = 'Socialising'
                    else:
                        Behav = 'Unknown'
                    # end
                # end
            else:
                if CPS90Perc1.mean() > 150:
                    Behav = 'Socialising'
                elif CPS20 > 100 and CPS50 < 100 and CPS90 < 100:
                    Behav = 'Socialising'
                else:
                    Behav = 'Unknown'
                # end
            # end
        # end
    # end
    # end
    return Behav


def NewICI(myTable, fs):
    start_sample = myTable["start_sample"]
    myTable = myTable.assign(ICI=start_sample.diff() / (fs / 1000))
    myTable = myTable.assign(CPS=1000 / myTable.ICI)
    myTable.at[0, 'CPS'] = 0.0
    myTable.at[0, 'ICI'] = 0.0
    return myTable


def ExtractPatterns(myCP, myFs):
    CTNum = 0
    Keep = 0
    ColNames = ['CTNum', 'Date', 'DayNight', 'Length', 'Species', 'Behaviour', 'Calf', 'Notes']
    CTInfo = pd.DataFrame(data=None, columns=ColNames)
    # add rown numbers to see if what we removed is another click train
    myCP.reset_index(inplace=True, drop=True)
    # Good click trains that only have echoes
    myCP = NewICI(myCP, myFs)
    # remove local min
    S1 = len(myCP)
    S2 = 1
    while S1 != S2:
        S1 = len(myCP)
        myCP["AmpDiff"] = myCP.amplitude.diff()
        myCP = myCP.drop(
            myCP[(myCP.amplitude.shift(1) > myCP.amplitude) & (myCP.amplitude.shift(-1) > myCP.amplitude) & (
                    myCP.AmpDiff < -5)].index)
        myCP.reset_index(inplace=True, drop=True)
        # print(len(myCP))
        myCP = NewICI(myCP, myFs)
        S2 = len(myCP)

    myCP = myCP.drop(myCP[(myCP.myCPS.diff() > 140.0)].index)
    myCP.reset_index(inplace=True, drop=True)
    myCP = NewICI(myCP, Fs)
    ColNames = list(myCP.columns)
    CTrains = pd.DataFrame(data=None, columns=ColNames)
    # # Find click trains
    TimeGaps = myCP[(myCP.ICI > 700.0) | (myCP.ICI < 0.0)].index.to_series()
    TimeGaps.reset_index(inplace=True, drop=True)
    DiffTimeGaps = TimeGaps.diff()
    DiffTimeGaps.at[0] = TimeGaps[0] - 0
    # Find very long CT and reduce them to CT with fewer than 1000 clicks
    LongCTs = DiffTimeGaps[DiffTimeGaps > 900].index

    if len(LongCTs) > 0:
        for m in range(0, len(LongCTs)):
            # Length = int(DiffTimeGaps[LongCTs[m]])
            # CTsInside = math.floor(Length / 600) + 1  # integer less than
            # Add Positions to TimeGaps
            PosInCP = TimeGaps[LongCTs[m]]
            NextPos = TimeGaps[LongCTs[m] + 1]
            TempCP = myCP.loc[PosInCP:NextPos]
            NewPos = TempCP.loc[(TempCP['ICI'] > 500)].index.to_series()  # probar con values
            TimeGaps.append(NewPos)
        # end
        # TimeGaps = TimeGaps.sort()
        DiffTimeGaps = TimeGaps.diff()
        CTs = DiffTimeGaps[DiffTimeGaps > 9].index
        # CTs.reset_index(inplace=True, drop=True)
    else:
        CTs = DiffTimeGaps[DiffTimeGaps > 9].index
        # CTs.reset_index(inplace=True, drop=True)
    # end

    for j in range(0, len(CTs)):  # j runs through all the CT c
        print(j, 'of', len(CTs))
        if j == 0:
            Start = 0
            End = TimeGaps[CTs[j]]
        else:
            Start = TimeGaps[CTs[j] - 1]
            End = TimeGaps[CTs[j]]
        CTTemp = myCP[Start:End]
        CTTemp.reset_index(inplace=True, drop=True)
        CTTemp = NewICI(CTTemp, myFs)
        CTNum = CTNum + 1
        CTTemp = CTTemp.assign(CT=CTNum)
        # Delete loose clicks
        L1 = len(CTTemp)
        L2 = 1
        while L2 != L1:
            L1 = len(CTTemp)
            CTTemp = NewICI(CTTemp, myFs)
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
            NewCT = CTTemp.copy()
        else:
            CTTemp = CTTemp.drop(CTTemp[CTTemp.CPS >= int(MaxDiffSorted) - 30].index)
            CTTemp = NewICI(CTTemp, myFs)
            CTTemp.reset_index(inplace=True, drop=True)
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
                # print(PercCPSDiff)
                StartRow = PercCPSDiff[PercCPSDiff[0] < 20.0].index.to_series()
                StartRow.reset_index(inplace=True, drop=True)
                DiffStartRow = StartRow.diff()
                Here = DiffStartRow[DiffStartRow == 1].index.to_series()
                Here.reset_index(inplace=True, drop=True)
                if len(StartRow) < 2:
                    NewCT = CTTemp.copy()
                    NewCT = NewICI(NewCT, myFs)
                else:  # go into the CT
                    RowN = StartRow[0]  # Low variability in CPS(for the next 4)
                    RowsToKeep = np.array(Here)  # # ok<FNDSB>
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
                        RowsToKeep = np.append(RowsToKeep, RowN)
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
                            RowsToKeep = np.append(RowsToKeep, RowN)
                    # end
                    RowsToKeep = np.append(RowsToKeep, RowN)
                    # end
                    # adding the first 10 and the last 10 rows (although might not be needed)

                    # RowsToKeep.append(CTTemp[0:10].index)
                    # RowsToKeep.append(CTTemp[-1-9:-1].index)
                    RowsToKeep = np.unique(RowsToKeep)
                    RowsToKeep = np.delete(RowsToKeep, np.where(RowsToKeep <= 0))
                    RowsToKeep = np.delete(RowsToKeep, np.where(RowsToKeep > len(CTTemp)))
                    # print(RowsToKeep)
                    NewCT = CTTemp.loc[RowsToKeep]
                    NewCT = NewICI(NewCT, myFs)
                    # end
            else:
                NewCT = CTTemp.copy()
                NewCT = NewICI(NewCT, myFs)
            # end
        # end
        CTInfo = CTInfoMaker(CTInfo, CTTemp)
        # print(CTInfo)
        # Put result into a final file
        if j == 0:
            CTrains = NewCT.copy()
        else:
            CTrains = CTrains.append(NewCT.copy())
        # end
    # end
    return CTrains, CTInfo, myCP


CP = pd.read_csv('C:/Mel/DPorCCA_FROM_SCRATCH/PGDF_SoundTrapNewTest/20150812/CPar.csv')
Fs = 576000
CTrains, CTInfo, CP = ExtractPatterns(CP, Fs)
CTInfo.to_csv('C:/Mel/DPorCCA_FROM_SCRATCH/PGDF_SoundTrapNewTest/20150812/CTInfo.csv', index=False)
CTrains.to_csv('C:/Mel/DPorCCA_FROM_SCRATCH/PGDF_SoundTrapNewTest/20150812/CTrains.csv', index=False)

thisFolder = 'C:/Mel/CPODvsDPorCCA/Simulations/NoiseSamples/JustClicks'
# thisFolder = FolderNames[j]
FilesInFolder = os.listdir(thisFolder)
if any("." in s for s in FilesInFolder):
    DetClicks = [s for s in FilesInFolder if "Detected_Clicks" in s]
    CP = pd.DataFrame()
    for thisFile in DetClicks:
        CPFile = thisFolder + '/' + thisFile
        TempCP = pd.read_csv(CPFile, index_col=0)
        TempCP['id'] = TempCP.index.values
        CP = CP.append(TempCP, ignore_index=True)
        CPFileName = thisFolder + '/CP.csv'
        CP.to_csv(CPFileName, index=False, sep=';')

Fs = 576000
CTrains, CTInfo, CP = ExtractPatterns(CP, Fs)
# CTrains.CT = int(CTrains.CT)
CTInfoFileName = thisFolder + '/CTInfo.csv'
CTrainsFileName = thisFolder + '/CTrains.csv'
CTInfo.to_csv(CTInfoFileName, index=False)
CTrains.to_csv(CTrainsFileName, index=False)

CP = pd.read_pickle('C:/Mel/CPODvsDPorCCA/Simulations/NoiseSamples/JustClicks/Detected_Clips_010120_000100.pkl')
CP1 = pd.read_csv('C:/Mel/CPODvsDPorCCA/Simulations/NoiseSamples/JustClicks/Detected_Clicks_010120_000100.csv')

fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].scatter(CP.start_sample, CP.amplitude, s=2.0)
ax[1].scatter(CP.start_sample, CP.prob_lq, s=2.0)
ax[0].set_title('CP detections')
ax[1].set_title('pyporcc detections')
ax[0].set_ylabel('Amplitude rms [db]')
ax[1].set_ylabel('Amplitude rms [db]')
ax[1].set_xlabel('Samples')
plt.show()

fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].scatter(CP1.start_sample, CP1.amplitude, s=2.0)
ax[1].scatter(CP1.start_sample, CP1.prob_lq, s=2.0)
ax[0].set_title('CP1 detections')
ax[1].set_title('pyporcc detections')
ax[0].set_ylabel('Duration [samples]')
ax[1].set_ylabel('Duration [samples]')
ax[1].set_xlabel('Samples')
plt.show()
