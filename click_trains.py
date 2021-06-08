import math
import warnings

import numpy as np
import pandas as pd
from tqdm import tqdm

import isoutlier
import sunrise


def ExtractPatterns(myCP, myFs, lat, long):
    """
     THIS FUNCTION IS CURRENTLY UNDER DEVELOPMENT. THE MATLAB VERSION SHOULD BE USED INSTEAD

     Locates acoustic events and identifies underlying patterns by keeping consecutive clicks with regular variations
     of inter-click intervals (or repetition rates) and amplitude.

     Parameters:
        myCP = table containing parameters of each click identified in the data, which were already classified as either
          high- or low-quality porpoise clicks.
        myFs = sampling frequency
        lat = latitude of the location of the device. Used to estimate whether the acoustic event occurred during the
          day or the night
        long = longitude of the location of the device. Used to estimate whether the acoustic event occurred during the
          day or the night

     Variables used:

         CTTemp: temporary click train. Using a minimum separation time
               period of 700 ms, clicks are separated into click trains.
               Click trains that have more than 1000 clicks are divided
               into smaller ones of the same size.
         FinalCT: final version of the click train. It is stored in another
               file.

    """
    ColNames = ['CTNum', 'Date', 'DayNight', 'Length', 'CTType', 'Behav', 'Calf', 'Notes']
    myCTInfo = pd.DataFrame(data=None, columns=ColNames)
    myCP = NewICI(myCP, myFs)
    # remove local min
    S1 = len(myCP)
    S2 = 1
    while S1 != S2:
        S1 = len(myCP)
        myCP["AmpDiff"] = myCP.amplitude.diff()
        myCP = myCP.drop(
            myCP[(myCP.amplitude.shift(1) > myCP.amplitude) & (myCP.amplitude.shift(-1) > myCP.amplitude) & (
                    myCP.AmpDiff < -6)].index)
        myCP.reset_index(inplace=True, drop=True)
        myCP = NewICI(myCP, myFs)
        S2 = len(myCP)
    myCP["AmpDiff"] = myCP.amplitude.diff()
    myCP = myCP.drop(myCP[(myCP.CPS.diff() > 80.0)].index)
    myCP.reset_index(inplace=True, drop=True)
    myCP = NewICI(myCP, myFs)
    myCP.reset_index(inplace=True, drop=True)
    ColNames = list(myCP.columns)
    Clicks = pd.DataFrame(data=None, columns=ColNames)
    # # Find click trains
    TimeGaps = myCP.loc[(myCP['ICI'] > 700.0) | (myCP['ICI'] < 0.0)].index.to_series()
    TimeGaps.reset_index(inplace=True, drop=True)
    if len(TimeGaps) == 0:
        Vals = myCP.sort_values('ICI', ascending=False).head(int(myCP.shape[0] * .008)).copy()
        Th = int(Vals.ICI.iloc[-1])
        TimeGaps = myCP.loc[(myCP['ICI'] > Th) | (myCP['ICI'] < 0.0)].index.to_series()
        TimeGaps.reset_index(inplace=True, drop=True)

    DiffTimeGaps = TimeGaps.diff()
    # Find very long CT and reduce them to CT with fewer than 1000 clicks
    LongCTs = DiffTimeGaps[DiffTimeGaps > 1000.0].index

    if len(LongCTs) > 0:
        TimeGapsAdd = pd.Series()
        for m in range(0, len(LongCTs)):
            LongLength = int(DiffTimeGaps[LongCTs[m]])  # Original number of clicks
            NNewCT = math.floor(LongLength / 1000) + 1  # Number of divisions
            NewLength = int(LongLength / NNewCT)  # New length of each
            PosInCP = int(TimeGaps[LongCTs[m] - 1])
            for i in np.arange(1, NNewCT):
                TimeGapsAdd.loc[i] = PosInCP + i * NewLength + 1

        TimeGaps = TimeGaps.append(TimeGapsAdd)
        TimeGaps = TimeGaps.sort_values()

        TimeGaps.reset_index(inplace=True, drop=True)
        DiffTimeGaps = TimeGaps.diff()
        CTs = DiffTimeGaps[DiffTimeGaps > 9].index
    else:
        CTs = DiffTimeGaps[DiffTimeGaps > 9].index
    print(TimeGaps)
    print(CTs)
    for j in tqdm(range(0, len(CTs) + 1), total=len(CTs) + 1):  # j runs through all the CT c
        if j == 0:
            Start = TimeGaps[0]
            End = TimeGaps[CTs[0]] - 1
        elif j == (len(CTs) + 1):
            Start = TimeGaps[CTs[-1]]
            End = len(myCP)
        else:
            i = j - 1
            Start = TimeGaps[CTs[i] - 1]
            End = TimeGaps[CTs[i]] - 1
        CT = myCP[Start:End]
        CT.reset_index(inplace=True, drop=True)
        CT = NewICI(CT, myFs)
        CTNum = j + 1
        CT = CT.assign(CT=CTNum)

        # A large difference indicates echoes
        SortedCPS = CT.CPS.values.copy()
        SortedCPS.sort()
        DiffSorted = pd.DataFrame(SortedCPS).diff()
        MaxDiffSorted = DiffSorted.max().values

        if MaxDiffSorted <= 40:
            FinalCT = CT.copy()
        else:

            Outlier = isoutlier.isoutliers(CT[['CPS']])
            HighCPS = CT.CPS > 100
            CT.drop(index=CT[HighCPS & Outlier].index)
            CT.reset_index(inplace=True, drop=True)
            CT = NewICI(CT, myFs)
            SortedCPS = CT.CPS.values.copy()
            SortedCPS.sort()
            DiffSorted = pd.DataFrame(SortedCPS).diff()
            MaxDiffSorted = DiffSorted.max().values

            if MaxDiffSorted <= 50:
                FinalCT = CT.copy()
            elif len(CT) > 20:
                # Finding stable areas
                CPSDiff = CT.CPS.diff()
                PercChangeS = (CPSDiff / CT.CPS[::-1]) * 100
                PercChangeS = abs(PercChangeS[1::])
                window_size = 3
                i = 0
                PercCPSDiff = []
                # moving average
                while i < len(PercChangeS) - window_size + 1:
                    this_window = PercChangeS[i: i + window_size]
                    # get current window
                    window_average = sum(this_window) / window_size
                    PercCPSDiff.append(window_average)
                    i += 1
                PercCPSDiff = pd.DataFrame(PercCPSDiff)
                PercCPSDiff.reset_index(inplace=True, drop=True)
                StartRow = PercCPSDiff[PercCPSDiff[0] < 20.0].index.to_series()
                StartRow.reset_index(inplace=True, drop=True)
                DiffStartRow = StartRow.diff()
                Here = DiffStartRow[DiffStartRow == 1].index.to_series()
                Here.reset_index(inplace=True, drop=True)
                if len(StartRow) < 2:
                    FinalCT = CT.copy()
                    FinalCT = NewICI(FinalCT, myFs)
                else:  # go into the CT
                    RowN = StartRow[0]  # Low variability in CPS (in the next 4 clicks)
                    RowsToKeep = np.array(Here)
                    FirstRowN = RowN

                    # Look backwards
                    while RowN > 5:
                        ClickCPS = CT.CPS.iloc[RowN]
                        ClickAmp = CT.amplitude.iloc[RowN]
                        Clickstart_sample = CT.start_sample.iloc[RowN]
                        ICIs = abs(Clickstart_sample - CT.start_sample[RowN - 5:RowN - 1]) / (myFs / 1000)
                        CPSs = 1000 / ICIs
                        Amps = CT.amplitude[RowN - 5:RowN - 1]
                        Amps = abs(ClickAmp - Amps)
                        DiffCPSs = abs(ClickCPS - CPSs)
                        DiffCPSs = pd.DataFrame(DiffCPSs)
                        MinVal = DiffCPSs.start_sample.min()
                        ixCPS = DiffCPSs[DiffCPSs.start_sample == MinVal].index
                        if Amps[ixCPS[0]] < 5:
                            DiffCPSs[ixCPS[0]] = 1000  # high arbitrary number
                            MinVal = DiffCPSs.start_sample.min()
                            ixCPS = DiffCPSs[DiffCPSs.start_sample == MinVal].index
                            RowN = RowN - ixCPS[0]
                        else:
                            RowN = RowN - ixCPS[0]

                        RowsToKeep = np.append(RowsToKeep, RowN)
                    # Look forwards
                    RowN = FirstRowN
                    while RowN < len(CT) - 10:
                        ClickCPS = CT.CPS.iloc[RowN]
                        ClickAmp = CT.amplitude.iloc[RowN]
                        Clickstart_sample = CT.start_sample.iloc[RowN]
                        ICIs = abs(CT.start_sample[RowN + 1:RowN + 9] - Clickstart_sample) / (myFs / 1000)
                        CPSs = 1000 / ICIs
                        Amps = CT.amplitude[RowN + 1:RowN + 9]
                        Amps = abs(Amps - ClickAmp)
                        DiffCPSs = abs(ClickCPS - CPSs)
                        DiffCPSs = pd.DataFrame(DiffCPSs)
                        MinVal = DiffCPSs.start_sample.min()
                        ixCPS = DiffCPSs[DiffCPSs.start_sample == MinVal].index
                        if Amps[ixCPS[0]] < 6:
                            DiffCPSs[ixCPS[0]] = 1000  # high arbitrary number
                            MinVal = DiffCPSs.start_sample.min()
                            ixCPS = DiffCPSs[DiffCPSs.start_sample == MinVal].index
                            RowN = RowN + ixCPS[0]
                        else:
                            RowN = RowN + ixCPS[0]

                        RowsToKeep = np.append(RowsToKeep, RowN)

                    RowsToKeep = np.append(RowsToKeep, RowN)
                    RowsToKeep = np.unique(RowsToKeep)
                    RowsToKeep = np.delete(RowsToKeep, np.where(RowsToKeep <= 0))
                    RowsToKeep = np.delete(RowsToKeep, np.where(RowsToKeep > len(CT) - 1))
                    if len(RowsToKeep) > 0:
                        FinalCT = CT.iloc[RowsToKeep]
                        FinalCT.reset_index(inplace=True, drop=True)
                        FinalCT = NewICI(FinalCT, myFs)
                    else:
                        FinalCT = []
            else:
                FinalCT = CT.copy()
                FinalCT = NewICI(FinalCT, myFs)

        if len(FinalCT) > 1:
            FinalCT = NewICI(FinalCT, myFs)
            myCTInfo = CTInfoMaker(myCTInfo, FinalCT, lat, long)
            # Put result into a final file
            if j == 0:
                Clicks = FinalCT.copy()
            else:
                Clicks = Clicks.append(FinalCT.copy())
        else:
            warnings.warn("This click train was empty, the number will be skipped")
    return Clicks, myCTInfo, myCP


def NewICI(myTable, fs):
    """
        Calculates inter-click intervals (ICI) and repetition rates (clicks per second - CPS)
    :parameters:
        myTable: pandas dataframe with at least the following column
            start_sample: indicates the sample number where each clicks begins in the wav file
        fs: sampling frequency

    :return:
        myTable updated
    """
    start_sample = myTable["start_sample"]
    myTable = myTable.assign(ICI=start_sample.diff() / (fs / 1000))
    myTable = myTable.assign(CPS=1000 / myTable.ICI)
    myTable.at[0, 'CPS'] = 0.0
    myTable.at[0, 'ICI'] = 0.0
    return myTable


def CTInfoMaker(myCTInfo, myCTTemp, myLat, myLong):
    """
        CTInfoMaker: for each identified click train (myCTTemp), this function estimates a series of summary parameters
            and generates a table called CTInfo. These parameters are:
                CTNum = click train number
                datetime: date and time
                DayNight = whether it was detected during the day or at night
                Length = length (in number of clicks)
                CTType = type of click train (NBHF, LQ-NBHF, Noise, Sonar)
        Parameters:
            myCTInfo: pandas dataframe to store summary data of click trains
            myCTTemp: click train
            myLat: latitude of the location of the recording device
            myLong: longitude of the location of the recording device
    """
    # Store in CTInfo
    CTNum = myCTTemp.CT[0]
    # day/night
    today = myCTTemp.datetime.iloc[0]
    day = int(today[8:10])
    month = int(today[5:7])
    year = int(today[0:4])
    sriseH, sriseM = sunrise.getSunriseTime(day, month, year, myLong, myLat)
    ssetH, ssetM = sunrise.getSunsetTime(day, month, year, myLong, myLat)
    # I don't know which format time is returned here, need to correct when I do
    HH = today[11:13]
    MM = today[14:16]

    if int(sriseH) < int(HH) < int(ssetH):
        DayNight = 'Day'
    elif int(sriseH) == int(HH):
        if int(MM) >= int(sriseM):
            DayNight = 'Day'
        else:
            DayNight = 'Night'
    elif int(ssetH) == int(HH):
        if int(MM) >= int(ssetM):
            DayNight = 'Night'
        else:
            DayNight = 'Day'
    elif int(sriseH) > int(HH) > int(ssetH):
        DayNight = 'Night'
    else:
        DayNight = 'Night'
    # Type
    Type = CTType(myCTTemp)
    if Type == 'Noise':
        Behav = '-'
    else:
        Behav = Behaviour(myCTTemp)
    myCTInfo = myCTInfo.append({'CTNum': CTNum, 'Date': myCTTemp.datetime[0], 'Length': len(myCTTemp), 'CTType': Type,
                                'DayNight': DayNight, 'Behav': Behav, 'Calf': '-', 'Notes': ' '}, ignore_index=True)
    return myCTInfo


def CTType(CT):
    """
        Estimates a series of parameters for each click train and classify them into either of three categories:
            NBHF: narrow-band, high-frequency click trains, with a high probability of being produced by harbour
                porpoises or species that emit similar signals.
            LQ-NBHF: low-quality NBHF click trains, with higher false alarms levels.
            Noise: high-frequency background noise.
        Parameters:
            CT: click train
        Returns:
            Type of click train
    """
    CFDiff = CT.CF.diff()
    PercChangeCF = (CFDiff / CT.CF[::-1]) * 100
    MedianPercChangeCF = abs(PercChangeCF[1::]).median()
    SDCF = CT.CF.std()
    CPSDiff = CT.CPS.diff()
    CT['CPS'] = CT.CPS.replace(0.0, np.nan)
    PercChange = (CPSDiff / CT.CPS[::-1]) * 100
    MedianPercChangeCPS = abs(PercChange[1::]).median()
    if len(CT) < 10:
        Type = 'Noise'
    elif MedianPercChangeCF < 0.5 and MedianPercChangeCPS < 0.05 and SDCF < 300:
        Type = 'Noise'
    elif MedianPercChangeCPS > 70 or MedianPercChangeCF > 4:
        Type = 'Noise'
    elif MedianPercChangeCPS < 30 or (MedianPercChangeCPS < 30 and MedianPercChangeCF > 2.65):
        Type = 'NBHF'
    else:
        Type = 'LQ-NBHF'
    return Type


def Behaviour(CT):
    """
        THIS FUNCTION IS STILL UNDER DEVELOPMENT

        Based on the patterns in repetition rates in NBHF and LQ-NBHF click trains, this functions classify these
        patterns into either of 5 categories:
            Orientation: low repetition rates (below 100 clicks per second), indicating the animal is exploring the
                environment without focusing on any object in particular
            Foraging: increase in click production to well over 100 clicks per second, reaching up to over 600
            Socialising: repetition rates over 100 clicks per second, or decreasing pattern
            Unknown: neither of the patterns described above
            Sonar: fixed frequency and repetition rate. These vary depending on the area.

        Parameters:
            CT = click train
    """
    CFDiff = CT.CF.diff()
    PercChangeCF = (CFDiff / CT.CF[0:-1]) * 100
    MeanPF = CT.CF.mean()
    MeanPercChangeCF = abs(PercChangeCF).mean()
    CPS = CT.CPS
    MedianCPS = CPS.median()
    L = len(CPS)
    SortedCPS = CPS.values.copy()
    SortedCPS.sort()
    CPS90Perc1 = SortedCPS[0:math.floor(0.90 * L)]
    CPS20 = CPS[0:math.floor(L * 0.2)].mean()
    CPS50 = CPS[math.floor(0.2 * L):math.floor(0.8 * L)].mean()
    CPS90 = CPS[math.floor(0.8 * L)::].mean()
    if MeanPF > 140000 and 8.5 > MedianCPS > 7.1 and MeanPercChangeCF < 0.5:
        Behav = 'Sonar'
    elif all(CPS90Perc1 < 100):
        Behav = 'Orientation'
    else:
        CPS90Perc2 = SortedCPS[math.floor(0.10 * L)::]
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
                if Before < 100 and After > 100:
                    Behav = 'Foraging'
                else:
                    if CPS90Perc1.mean() > 120:
                        Behav = 'Socialising'
                    elif CPS20 > 100 and 100 > CPS50 > CPS90:
                        Behav = 'Socialising'
                    else:
                        Behav = 'Unknown'
            else:
                if CPS90Perc1.mean() > 150:
                    Behav = 'Socialising'
                elif CPS20 > 100 and CPS50 < 100 and CPS90 < 100:
                    Behav = 'Socialising'
                else:
                    Behav = 'Unknown'
    return Behav
