import math
import warnings

import pandas as pd
import numpy as np
from tqdm import tqdm

import isoutlier
import sunrise


def ct_steps(time_gaps, ct_long, time_gap_diff, old_time_gaps, cp):
    while len(ct_long) > 0:
        for m in range(0, len(ct_long)):
            length_ct = time_gap_diff[ct_long[m]].astype(int)  # Original number of clicks
            start = old_time_gaps[ct_long[m] - 1].astype(int)  # Location of the beginning of the click train
            stop = start + length_ct - 1  # Location of the last click in the click train
            long_ct = cp[start:stop]  # Click train to split
            values = long_ct.sort_values('ICI', ascending=False).head(int(long_ct.shape[0] * .01)).copy()
            values.reset_index(inplace=True, drop=True)
            if len(values) > 0:
                th = int(values.ICI.min())
                new_pos = long_ct.loc[(cp['ICI'] > th) | (long_ct['ICI'] < 0.0)].index.to_series()  # Split positions
                new_pos_series = pd.Series(new_pos)
                time_gaps = time_gaps.append(new_pos_series)
        time_gaps = time_gaps.sort_values()
        time_gaps = time_gaps.unique()
        time_gaps = pd.Series(time_gaps)
        time_gaps.reset_index(inplace=True, drop=True)
        new_time_gap_diff = time_gaps.diff()
        new_time_gap_diff.reset_index(inplace=True, drop=True)
        ct_long = new_time_gap_diff[new_time_gap_diff > 1000].index
        old_time_gaps = time_gaps.copy()
        time_gap_diff = new_time_gap_diff.copy()

    click_trains = new_time_gap_diff[new_time_gap_diff > 9].index
    return click_trains, time_gaps


def extract_patterns(myCP, lat, long):
    """
     Locates acoustic events and identifies underlying patterns by keeping consecutive clicks with regular variations
     of inter-click intervals (or repetition rates) and amplitude.

     Parameters:
        myCP = table containing parameters of each click identified in the data, which were already classified as either
          high- or low-quality porpoise clicks.
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
    myCP.datetime = pd.to_datetime(myCP.datetime)
    myCP = new_ici(myCP)
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
        myCP = new_ici(myCP)
        S2 = len(myCP)
    myCP["AmpDiff"] = myCP.amplitude.diff()
    myCP = myCP.drop(myCP[(myCP.CPS.diff() > 80.0)].index)
    myCP.reset_index(inplace=True, drop=True)
    myCP = new_ici(myCP)
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

    # Compute the number of clicks per "gap" chosen
    DiffTimeGaps = TimeGaps.diff()
    DiffTimeGaps[0] = TimeGaps[0]  # (first one from 0)
    DiffTimeGaps[len(DiffTimeGaps)] = len(myCP) - TimeGaps.iloc[-1]  # Add last CT from last time gap until the end

    # Find very long CT and reduce them to CT with fewer than 1000 clicks
    LongCTs = DiffTimeGaps[DiffTimeGaps > 1000.0].index

    # Compute CT, which is a pandas index with all the indexes of TimeGaps where a CT starts.
    # TimeGaps then refers to the index of myCP (dataframe with all clicks) where the CT starts
    if len(LongCTs) > 0:
        new_time_gaps = TimeGaps.copy()
        CTs, TimeGaps = ct_steps(new_time_gaps, LongCTs, DiffTimeGaps, TimeGaps, myCP)
    else:
        CTs = DiffTimeGaps[DiffTimeGaps > 9].index

    if len(CTs) > 0:
        # j runs through all the CT
        # Then selects Start and Stop of click train (index of myCP where the CT starts and stops)
        for j, gap_idx in tqdm(enumerate(CTs), total=len(CTs)):
            if gap_idx == 0:
                # First one goes from 0 to the first time gap
                Start = 0
                End = TimeGaps[gap_idx]
            elif gap_idx == len(TimeGaps):
                # Last one goes from last time gap until the last click of the dataframe
                Start = TimeGaps.iloc[-1]
                End = len(myCP)
            else:
                Start = TimeGaps[CTs[j] - 1]
                End = TimeGaps[CTs[j]]

            # Create CT, a dataframe with all the clicks of the Click Train
            CT = myCP[Start:End]
            CT.reset_index(inplace=True, drop=True)
            CT = new_ici(CT)
            CT = CT.assign(CT=j)

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
                CT = new_ici(CT)
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
                        FinalCT = new_ici(FinalCT)
                    else:  # go into the CT
                        RowN = StartRow[0]  # Low variability in CPS (in the next 4 clicks)
                        RowsToKeep = np.array(Here)
                        FirstRowN = RowN

                        # Look backwards
                        while RowN > 5:
                            ClickCPS = CT.CPS.iloc[RowN]
                            ClickAmp = CT.amplitude.iloc[RowN]
                            Clickstart_datetime = CT.datetime.iloc[RowN]
                            ICIs = abs((Clickstart_datetime - CT.datetime[RowN - 5:RowN - 1]).dt.microseconds / 1000)
                            CPSs = 1000 / ICIs
                            Amps = abs(ClickAmp - CT.amplitude[RowN - 5:RowN - 1])
                            DiffCPSs = abs(ClickCPS - CPSs)
                            ixCPS = DiffCPSs.idxmin()
                            if Amps[ixCPS] < 5:
                                DiffCPSs[ixCPS] = 1000  # high arbitrary number
                                ixCPS = DiffCPSs.idxmin()
                                RowN = RowN - ixCPS
                            else:
                                RowN = RowN - ixCPS

                            RowsToKeep = np.append(RowsToKeep, RowN)
                        # Look forwards
                        RowN = FirstRowN
                        while RowN < len(CT) - 10:
                            ClickCPS = CT.CPS.iloc[RowN]
                            ClickAmp = CT.amplitude.iloc[RowN]
                            Clickstart_datetime = CT.datetime.iloc[RowN]
                            ICIs = abs((CT.datetime[RowN + 1:RowN + 9] - Clickstart_datetime).dt.microseconds / 1000)
                            CPSs = 1000 / ICIs
                            Amps = abs(CT.amplitude[RowN + 1:RowN + 9] - ClickAmp)
                            DiffCPSs = abs(ClickCPS - CPSs)
                            ixCPS = DiffCPSs.idxmin()
                            if Amps[ixCPS] < 6:
                                DiffCPSs[ixCPS] = 1000  # high arbitrary number
                                ixCPS = DiffCPSs.idxmin()
                                RowN = RowN + ixCPS
                            else:
                                RowN = RowN + ixCPS

                            RowsToKeep = np.append(RowsToKeep, RowN)

                        RowsToKeep = np.append(RowsToKeep, RowN)
                        RowsToKeep = np.unique(RowsToKeep)
                        RowsToKeep = np.delete(RowsToKeep, np.where(RowsToKeep <= 0))
                        RowsToKeep = np.delete(RowsToKeep, np.where(RowsToKeep > len(CT) - 1))
                        if len(RowsToKeep) > 0:
                            FinalCT = CT.iloc[RowsToKeep]
                            FinalCT.reset_index(inplace=True, drop=True)
                            FinalCT = new_ici(FinalCT)
                        else:
                            FinalCT = []
                else:
                    FinalCT = CT.copy()
                    FinalCT = new_ici(FinalCT)

            if len(FinalCT) > 1:
                FinalCT = new_ici(FinalCT)
                myCTInfo = ct_info_maker(myCTInfo, FinalCT, lat, long)
                # Put result into a final file
                if j == 0:
                    Clicks = FinalCT.copy()
                else:
                    Clicks = Clicks.append(FinalCT.copy())
            else:
                warnings.warn("This click train was empty, the number will be skipped")
    else:
        warnings.warn("There are no click trains in this file.")
    return Clicks, myCTInfo, myCP


def new_ici(myTable):
    """
        Calculates inter-click intervals (ICI) and repetition rates (clicks per second - CPS)
    :parameters:
        myTable: pandas dataframe with at least the datetime (should be a datetime in pandas, not a string!)

    :return:
        myTable updated
    """
    fs = 576000
    myTable['ICI'] = myTable.start_sample.diff()/(fs/1000)
                    #  myTable..datetime.diff().dt.total_seconds() * 1000
    myTable = myTable.assign(CPS=1000 / myTable.ICI)
    myTable.at[0, 'CPS'] = 0.0
    myTable.at[0, 'ICI'] = 0.0
    return myTable


def ct_info_maker(myCTInfo, myCTTemp, myLat, myLong):
    """
        ct_info_maker: for each identified click train (myCTTemp), this function estimates a series of summary parameters
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
    day = today.day
    month = today.month
    year = today.year
    sriseH, sriseM = sunrise.getSunriseTime(day, month, year, myLong, myLat)
    ssetH, ssetM = sunrise.getSunsetTime(day, month, year, myLong, myLat)
    # I don't know which format time is returned here, need to correct when I do
    HH = today.hour
    MM = today.minute

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
    Type = ct_type(myCTTemp)
    if Type == 'Noise':
        Behav = '-'
    else:
        Behav = Behaviour(myCTTemp)
    myCTInfo = myCTInfo.append({'CTNum': CTNum, 'Date': myCTTemp.datetime[0], 'Length': len(myCTTemp), 'CTType': Type,
                                'DayNight': DayNight, 'Behav': Behav, 'Calf': '-', 'Notes': ' '}, ignore_index=True)
    return myCTInfo


def ct_type(CT):
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
