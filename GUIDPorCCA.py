"""
Module : GUIDPorCCA.py
Authors : Mel Cosentino
Institution : University of Strathclyde (UK), Aarhus University (Denmark)
Last Accessed : 11/10/2020
"""

__author__ = "Mel Cosentino"
__version__ = "0"
__credits__ = "Mel Cosentino"
__email__ = "orcinus.orca.1758@gmail.com"
__status__ = "Development"

import math
import os
import tkinter as tk
from datetime import datetime
from datetime import time
from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyhydrophone as pyhy
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import soundfile
from PyQt5 import QtCore, QtGui
from PyQt5 import QtWidgets
from pyporcc import click_detector
from pyporcc import porcc
from scipy import signal

# from PyQt5.QtWidgets import QTableWidgetItem

pd.options.mode.chained_assignment = None
# import PAMGuardFunc
# TODO import Clea's algorithms
# import PyHydrophones # need to install first

global BrowseSelectedFolder, CTInfo, CTTemp, NHyd, CP, topLevelFolderMetrics
global thisFolder, sset, srise, Fs, PosPorMin, SummaryTable, numberOfFoldersMetrics


def NewICI(myTable, fs):
    start_sample = myTable["start_sample"]
    myTable = myTable.assign(ICI=start_sample.diff() / (fs / 1000))
    myTable = myTable.assign(CPS=1000 / myTable.ICI)
    myTable.at[0, 'CPS'] = 0.0
    myTable.at[0, 'ICI'] = 0.0
    return myTable


def FromOrdinal(x):
    ix = int(x)
    dt = datetime.fromordinal(ix)
    remainder = float(x) - ix
    hour, remainder = divmod(24 * remainder, 1)
    minute, remainder = divmod(60 * remainder, 1)
    second, remainder = divmod(60 * remainder, 1)
    microsecond = int(1e6 * remainder)
    if microsecond < 10:
        microsecond = 0  # compensate for rounding errors
    # for some strange reason it is 1 year over the actual date!!
    dt = datetime(dt.year - 1, dt.month, dt.day, int(hour), int(minute),
                  int(second), microsecond)

    # if microsecond > 999990:  # compensate for rounding errors
    # dt += timedelta(microseconds=1e6 - microsecond)

    return dt


def ExtractPatterns(myCP, myFs, lat, long):
    CTNum = 0
    ColNames = ['CTNum', 'Date', 'DayNight', 'Length', 'Species', 'Behav', 'Calf', 'Notes']
    myCTInfo = pd.DataFrame(data=None, columns=ColNames)
    myCP = myCP.drop(myCP[myCP.duration > 450].index)
    myCP.reset_index(inplace=True, drop=True)
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

    # myCP = myCP.drop(myCP[(myCP.CPS.diff() > 80.0)].index)
    myCP.reset_index(inplace=True, drop=True)
    myCP = NewICI(myCP, myFs)
    ColNames = list(myCP.columns)
    CTrains = pd.DataFrame(data=None, columns=ColNames)
    # # Find click trains
    TimeGaps = myCP[(myCP.ICI > 700.0) | (myCP.ICI < 0.0)].index.to_series()
    TimeGaps.reset_index(inplace=True, drop=True)
    DiffTimeGaps = TimeGaps.diff()

    # Find very long CT and reduce them to CT with fewer than 1000 clicks
    LongCTs = DiffTimeGaps[DiffTimeGaps > 1000].index
    if len(LongCTs) > 0:
        for m in range(0, len(LongCTs) - 1):
            Length = int(DiffTimeGaps[LongCTs[m]])
            CTsInside = math.floor(Length / 1000) + 1  # integer less than
            # Add Positions to TimeGaps
            PosInCP = int(TimeGaps[LongCTs[m]])
            NextPos = int(TimeGaps[LongCTs[m] + 1])
            NewPos = np.linspace(start=PosInCP, stop=NextPos, num=CTsInside + 2).astype(int)
            NewPos = pd.Series(NewPos)
            TimeGaps.append(NewPos[1:-1 - 1])
        DiffTimeGaps = TimeGaps.diff()
        CTs = DiffTimeGaps[DiffTimeGaps > 9].index
    else:
        CTs = DiffTimeGaps[DiffTimeGaps > 9].index

    for j in range(0, len(CTs)):  # j runs through all the CT c
        print('Processing click train', j, 'of', len(CTs))
        if j == 0:
            Start = 0
            End = TimeGaps[CTs[j]]
        else:
            Start = TimeGaps[CTs[j] - 1]
            End = TimeGaps[CTs[j]]
        CT = myCP[Start:End]
        CT.reset_index(inplace=True, drop=True)
        L1 = len(CT)
        L2 = 1
        while L2 != L1:
            L1 = len(CT)
            CT = NewICI(CT, myFs)
            CT = CT.drop(CT[(CT.start_sample.diff() < 500)
                                        & ((CT.amplitude.diff() < 0) | (CT.start_sample.diff() > 6))].index)
            CT.reset_index(inplace=True, drop=True)
            L2 = len(CT)
        CTNum = CTNum + 1
        CT = CT.assign(CT=CTNum)
        # Delete loose clicks
        L1 = len(CT)
        L2 = 1
        while L2 != L1:
            L1 = len(CT)
            CT = NewICI(CT, myFs)
            LooseClicks = CT[CT.ICI > 400].index.to_series()
            # Positions = find(diff(LooseClicks) == 1)
            RowsToDelete = LooseClicks[LooseClicks.diff() == 1]
            CT = CT.drop(RowsToDelete)
            CT.reset_index(inplace=True, drop=True)
            L2 = len(CT)
        # end
        # A large difference indicates echoes
        SortedCPS = CT.CPS.values.copy()
        SortedCPS.sort()
        DiffSorted = pd.DataFrame(SortedCPS).diff()
        DiffSorted = DiffSorted.drop([0])
        MaxDiffSorted = DiffSorted.max().values
        # ExploreCTTemp(CTTemp)

        if MaxDiffSorted <= 40:
            NewCT = CT.copy()
        else:
            # remove outliers (from https://stackoverflow.com/questions/62692771/
            # outlier-detection-based-on-the-moving-mean-in-python)
            # Calculate rolling median
            CT['rolling_CPS'] = CT['CPS'].rolling(window=9).median()
            # Calculate difference
            CT['diff_CPS'] = CT['CPS'] - CT['rolling_CPS']
            # Flag rows to be dropped as `1`
            # Set threshold for difference with rolling median (they have to be adaptable)
            upper_threshold = CT['CPS'].median() * 10
            lower_threshold = -upper_threshold
            CT['drop_flag'] = np.where(
                (CT['diff_CPS'] > upper_threshold) | (CT['diff_CPS'] < lower_threshold), 1, 0)
            # Drop flagged rows
            CT = CT[CT['drop_flag'] != 1]
            CT = CT.drop(['rolling_CPS', 'rolling_CPS', 'diff_CPS', 'drop_flag'], axis=1)
            CT.reset_index(inplace=True, drop=True)

            SortedCPS = CT.CPS.values.copy()
            SortedCPS.sort()
            DiffSorted = pd.DataFrame(SortedCPS).diff()
            DiffSorted = DiffSorted.drop([0])
            MaxDiffSorted = DiffSorted.max().values
            # ExploreCTTemp(CTTemp)
            if MaxDiffSorted <= 50:
                NewCT = CT.copy()
            elif len(CT) > 20:
                # Finding the stable areas
                CPSDiff = CT.CPS.diff()
                PercChangeS = (CPSDiff / CT.CPS[1:-1 - 1]) * 100
                PercChangeS = abs(PercChangeS[2:-1])
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
                    NewCT = CT.copy()
                    NewCT = NewICI(NewCT, myFs)
                else:  # go into the CT
                    RowN = StartRow[0]  # Low variability in CPS(for the next 4)
                    RowsToKeep = np.array(Here)  # # ok<FNDSB>
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
                        # end
                        RowsToKeep = np.append(RowsToKeep, RowN)
                    # end

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
                            # end
                            RowsToKeep = np.append(RowsToKeep, RowN)
                    # end
                    RowsToKeep = np.append(RowsToKeep, RowN)
                    # end

                    RowsToKeep = np.unique(RowsToKeep)
                    RowsToKeep = np.delete(RowsToKeep, np.where(RowsToKeep <= 0))
                    RowsToKeep = np.delete(RowsToKeep, np.where(RowsToKeep > len(CT) - 1))

                    if len(RowsToKeep) > 0:
                        NewCT = CT.iloc[RowsToKeep]
                        NewCT.reset_index(inplace=True, drop=True)
                        NewCT = NewICI(NewCT, myFs)
            else:
                NewCT = CT.copy()
                NewCT = NewICI(NewCT, myFs)
            # end
        # end
        NewCT.reset_index(inplace=True, drop=True)
        if len(NewCT) > 0:
            myCTInfo = CTInfoMaker(myCTInfo, NewCT, lat, long)
            # Put result into a final file
            if j == 0:
                CTrains = NewCT.copy()
            else:
                CTrains = CTrains.append(NewCT.copy())
            # end
    # end
    return CTrains, myCTInfo, myCP


def forceRange(v, maxi):
    # force v to be >= 0 and < max
    if v < 0:
        return v + maxi
    elif v >= maxi:
        return v - maxi
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
    cosH = (math.cos(TO_RAD * zenith) - (sinDec * math.sin(TO_RAD * latitude))) / (
            cosDec * math.cos(TO_RAD * latitude))

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


def UpdateWaterfall(XLimMin, YLimMin, ZLimMin, XLimMax, YLimMax, ZLimMax):
    pass
    # WaterfallAx.XLim = [XLimMin, XLimMax]
    # WaterfallAx.YLim = [YLimMin, YLimMax]
    # WaterfallAx.ZLim = [ZLimMin, ZLimMax]
    # WaterfallAx.ZTick = []


def CTInfoMaker(myCTInfo, myCTTemp, myLat, myLong):  # srise, sset
    # Store in CTInfo
    CTNum = myCTTemp.CT[0]
    # day/night
    TODAY = myCTTemp.datetime.iloc[0]
    day = int(TODAY[8:10])
    month = int(TODAY[5:7])
    year = int(TODAY[0:4])
    sriseH, sriseM = getSunriseTime(day, month, year, myLong, myLat)
    ssetH, ssetM = getSunsetTime(day, month, year, myLong, myLat)
    # I don't know which format time is returned here, need to correct when I do
    HH = TODAY[11:13]
    MM = TODAY[14:16]

    if (int(HH) > ssetH and int(MM) > ssetM) and (int(HH) < sriseH and int(MM) < sriseM):
        DayNight = 'Day'
    else:  # if DATE >= sset and DATE <= srise
        DayNight = 'Night'

    # Type
    Type = Species(myCTTemp)
    if Type == 'Non-NBHF':
        Behav = '-'
    else:
        Behav = Behaviour(myCTTemp)
    myCTInfo = myCTInfo.append({'CTNum': CTNum, 'Date': myCTTemp.datetime[0], 'Length': len(myCTTemp), 'Species': Type,
                                'DayNight': DayNight, 'Behav': Behav, 'Calf': '-', 'Notes': ' '}, ignore_index=True)
    return myCTInfo


def Species(CT):
    CFDiff = CT.CF.diff()
    PercChangeCF = (CFDiff / CT.CF[0:-1 - 1]) * 100
    MedianPercChangeCF = abs(PercChangeCF).median()
    SDCF = CT.CF.std()
    CPSDiff = CT.CPS.diff()
    CPSDiff = CPSDiff.drop([0])
    PercChange = (CPSDiff / CT.CPS[1:-1 - 1]) * 100
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


def Behaviour(CT):
    CFDiff = CT.CF.diff()
    PercChangeCF = (CFDiff / CT.CF[0:-1 - 1]) * 100
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


class WinTable(QtWidgets.QMainWindow):
    """
    CREATE pop up windows (menus)
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.title = ''
        self.top = 100
        self.left = 100
        self.width = 360
        self.height = 500
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)


class Ui_MainWindow(object):
    """
    MAIN WINDOW
    """

    def __init__(self):
        self.font = QtGui.QFont()

        self.Fig3D = WinTable()
        self.Pan3D = QtWidgets.QFrame(self.Fig3D)
        self.Plot3D = gl.GLViewWidget(self.Pan3D)

        self.font2 = QtGui.QFont()
        self.SummaryTable = None
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuDownloads = QtWidgets.QMenu(self.menubar)
        self.menuClick_Trains = QtWidgets.QMenu(self.menubar)
        # self.menuMetrics_Display = QtWidgets.QMenu(self.menubar)
        self.menuMain_Display = QtWidgets.QMenu(self.menubar)

        # self.MenuSetPorCC = QtWidgets.QAction(MainWindow)
        self.MenuSetDetector = QtWidgets.QAction(MainWindow)
        self.MenuOpenCT = QtWidgets.QAction(MainWindow)
        self.MenuNewCT = QtWidgets.QAction(MainWindow)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)

        self.MainDisplayTab = QtWidgets.QWidget()
        self.DisplaySettings = QtWidgets.QFrame(self.MainDisplayTab)
        self.AxesPan = QtWidgets.QFrame(self.MainDisplayTab)

        # Metrics
        self.MetricsDisplayTab = QtWidgets.QWidget()
        self.MetricsDisplayPan = QtWidgets.QFrame(self.MetricsDisplayTab)
        self.SummaryPan = QtWidgets.QFrame(self.MetricsDisplayTab)
        self.CPMAxesMetrics = pg.PlotWidget(self.SummaryPan)
        self.CPMAxesTot = pg.PlotWidget(self.SummaryPan)
        self.CPMAxesTimes = pg.PlotWidget(self.SummaryPan)

        self.MetricsTablePan = QtWidgets.QFrame(self.MetricsDisplayTab)
        self.MetricsUploadPan = QtWidgets.QFrame(self.MetricsDisplayTab)
        self.UploadData = QtWidgets.QPushButton(self.MetricsUploadPan)
        self.SelectFoldLab = QtWidgets.QLabel(self.MetricsUploadPan)
        self.SelectFolderMetricEdit = QtWidgets.QLineEdit(self.MetricsUploadPan)
        self.MetricBrowseButton = QtWidgets.QPushButton(self.MetricsUploadPan)

        self.ClearButton = QtWidgets.QPushButton(self.MetricsDisplayPan)
        self.DisplayMetricsButton = QtWidgets.QPushButton(self.MetricsDisplayPan)
        self.MetricsTable = QtWidgets.QTableWidget(self.MetricsTablePan)
        self.MetricPPMPan = QtWidgets.QFrame(self.MetricsDisplayPan)
        self.AllDayPan = QtWidgets.QFrame(self.MetricPPMPan)
        self.NightDayButton = QtWidgets.QRadioButton(self.AllDayPan)
        self.AllButtonMetrics = QtWidgets.QRadioButton(self.AllDayPan)
        self.SelectDatePan = QtWidgets.QFrame(self.MetricsDisplayPan)
        self.PPMPan = QtWidgets.QFrame(self.MetricPPMPan)
        self.TypeofCTButton = QtWidgets.QRadioButton(self.PPMPan)
        self.PPMButton = QtWidgets.QRadioButton(self.PPMPan)
        self.BehaviourButton = QtWidgets.QRadioButton(self.PPMPan)
        self.SelectDateLab = QtWidgets.QLabel(self.SelectDatePan)
        self.SelectadateDropDown = QtGui.QComboBox(self.SelectDatePan)
        # Validate pan
        self.ActionPan = QtWidgets.QFrame(self.MainDisplayTab)
        self.DisplayAllNBHF = QtWidgets.QFrame(self.ActionPan)
        self.ActionLabel = QtWidgets.QLabel(self.ActionPan)
        # Validate
        self.ValidatePan = QtWidgets.QFrame(self.ActionPan)
        self.ValButton = QtWidgets.QPushButton(self.ValidatePan)
        self.ValBehavLab = QtWidgets.QLabel(self.ValidatePan)
        self.BehaviourDropDown = QtGui.QComboBox(self.ValidatePan)
        self.ValTypeLab = QtWidgets.QLabel(self.ValidatePan)
        self.CTTypeDropDown = QtGui.QComboBox(self.ValidatePan)
        self.ValidateLab = QtWidgets.QLabel(self.ValidatePan)
        self.CT3DPan = QtWidgets.QFrame(self.ActionPan)
        self.SpectrogramButton = QtWidgets.QPushButton(self.CT3DPan)
        self.IndClicksAnd3D = QtWidgets.QPushButton(self.CT3DPan)
        self.IndCTLabel = QtWidgets.QLabel(self.CT3DPan)
        self.NBHFButton = QtWidgets.QRadioButton(self.DisplayAllNBHF)
        self.AllButton = QtWidgets.QRadioButton(self.DisplayAllNBHF)
        self.ActionTextLabel = QtWidgets.QLabel(self.DisplayAllNBHF)
        # Axes
        self.AmpAxesCT = pg.PlotWidget(self.AxesPan)
        self.ICIAxesCT = pg.PlotWidget(self.AxesPan)  # used to be QGraphicsView
        self.FreqAxesCT = pg.PlotWidget(self.AxesPan)
        self.AmpAxesCT.setBackground(background='w')
        self.ICIAxesCT.setBackground(background='w')
        self.FreqAxesCT.setBackground(background='w')
        self.AmpAxesCT.showAxis('bottom', show=False)
        self.ICIAxesCT.showAxis('bottom', show=False)
        self.AmplitudeDB = QtWidgets.QLabel(self.AxesPan)

        self.FreqPan = QtWidgets.QFrame(self.AxesPan)
        self.DirectionofarrivalButton = QtWidgets.QRadioButton(self.FreqPan)
        self.CentroidfrequencykHzButton = QtWidgets.QRadioButton(self.FreqPan)
        self.ICIPan = QtWidgets.QFrame(self.AxesPan)
        self.ClickspersecondButton = QtWidgets.QRadioButton(self.ICIPan)
        self.InterclickintervalmsButton = QtWidgets.QRadioButton(self.ICIPan)

        self.NotesPan = QtWidgets.QFrame(self.DisplaySettings)
        self.SaveupdatesButton = QtWidgets.QPushButton(self.DisplaySettings)
        self.SeeAddButton = QtWidgets.QPushButton(self.NotesPan)
        self.NoteLabel = QtWidgets.QLabel(self.NotesPan)
        self.CalfPan = QtWidgets.QFrame(self.DisplaySettings)
        self.CalfLabel = QtWidgets.QLabel(self.CalfPan)
        self.CalfTextLabel = QtWidgets.QLabel(self.CalfPan)
        self.BehavPan = QtWidgets.QFrame(self.DisplaySettings)
        self.BehavLabel = QtWidgets.QLabel(self.BehavPan)
        self.BehavTextLabel = QtWidgets.QLabel(self.BehavPan)
        self.CTPan = QtWidgets.QFrame(self.DisplaySettings)
        self.CTTypeLabel = QtWidgets.QLabel(self.CTPan)
        self.CTLabel = QtWidgets.QLabel(self.CTPan)
        self.LengthPan = QtWidgets.QFrame(self.DisplaySettings)
        self.LengthLabel = QtWidgets.QLabel(self.LengthPan)
        self.LengthText = QtWidgets.QLabel(self.LengthPan)
        self.CTInfoPan = QtWidgets.QFrame(self.DisplaySettings)
        self.TotalLabel = QtWidgets.QLabel(self.CTInfoPan)
        self.CTForw = QtWidgets.QPushButton(self.CTInfoPan)
        self.CTNumD = QtWidgets.QLineEdit(self.CTInfoPan)
        self.CTNumD.setStyleSheet("QLineEdit{background: white}")
        self.CTBack = QtWidgets.QPushButton(self.CTInfoPan)
        self.CTNumLabel = QtWidgets.QLabel(self.CTInfoPan)
        self.DatePan = QtWidgets.QFrame(self.DisplaySettings)
        self.DateandtimeofCTLabel = QtWidgets.QLabel(self.DatePan)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.MainTab = QtWidgets.QTabWidget(self.centralwidget)
        self.DayPan = QtWidgets.QFrame(self.DisplaySettings)
        self.DayLabel = QtWidgets.QLabel(self.DayPan)
        self.DayNightTextLabel = QtWidgets.QLabel(self.DayPan)
        self.DateLabel = QtWidgets.QLabel(self.DatePan)

        """
        MENU
        """
        # Detector Menu
        self.MenuDetSetFig = WinTable()
        self.centralwidgetdet = QtWidgets.QWidget(self.MenuDetSetFig)
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidgetdet)
        self.statusbartabs = QtWidgets.QStatusBar(self.MenuDetSetFig)
        self.actionExit = QtWidgets.QAction(self.MenuDetSetFig)
        self.actionNew = QtWidgets.QAction(self.MenuDetSetFig)

        self.ProjectTab = QtWidgets.QWidget()
        self.ProjectPan = QtWidgets.QFrame(self.ProjectTab)
        self.FolderPathDet = QtWidgets.QLineEdit(self.ProjectPan)
        self.DefaultBut = QtWidgets.QPushButton(self.ProjectPan)
        self.LengthFileEdit = QtWidgets.QLineEdit(self.ProjectPan)
        self.LengthFileSave = QtWidgets.QLabel(self.ProjectPan)
        self.SelectFolderTextSave = QtWidgets.QLabel(self.ProjectPan)
        self.BrowseDet = QtWidgets.QPushButton(self.ProjectPan)
        self.ZipMode = QtWidgets.QCheckBox(self.ProjectPan)
        self.CheckAllFolders = QtWidgets.QCheckBox(self.ProjectPan)

        self.HydTab = QtWidgets.QWidget()

        self.HydPan = QtWidgets.QFrame(self.HydTab)
        self.ChannelLab = QtWidgets.QLabel(self.HydPan)
        self.SerialNoEdit = QtWidgets.QLineEdit(self.HydPan)
        self.SerialNoLab = QtWidgets.QLabel(self.HydPan)
        self.DataTypeLbl = QtWidgets.QLabel(self.HydPan)
        self.DataTypeDD = QtGui.QComboBox(self.HydPan)
        self.HydParLbl = QtWidgets.QLabel(self.HydPan)
        self.DAQppEdit = QtWidgets.QLineEdit(self.HydPan)
        self.DAQLab = QtWidgets.QLabel(self.HydPan)
        self.GainEditDet = QtWidgets.QLineEdit(self.HydPan)
        self.GainLab1 = QtWidgets.QLabel(self.HydPan)
        self.SensEditDet = QtWidgets.QLineEdit(self.HydPan)
        self.SenLabel = QtWidgets.QLabel(self.HydPan)
        self.EditHydNDet = QtWidgets.QLineEdit(self.HydPan)

        self.PorCCTab = QtWidgets.QWidget()
        self.PorCCPantab = QtWidgets.QFrame(self.PorCCTab)
        self.HQThresDet = QtWidgets.QLineEdit(self.PorCCPantab)
        self.ProbHQLab = QtWidgets.QLabel(self.PorCCPantab)
        self.LQThresDet = QtWidgets.QLineEdit(self.PorCCPantab)
        self.LQThresLab = QtWidgets.QLabel(self.PorCCPantab)
        self.PorCCLbl = QtWidgets.QLabel(self.PorCCPantab)
        self.MinSepEd = QtWidgets.QLineEdit(self.PorCCPantab)
        self.MinSepLbl = QtWidgets.QLabel(self.PorCCPantab)
        self.MinLengthEd = QtWidgets.QLineEdit(self.PorCCPantab)
        self.MinLengthLbl = QtWidgets.QLabel(self.PorCCPantab)
        self.MaxLengthEd = QtWidgets.QLineEdit(self.PorCCPantab)
        self.MaxLengthLbl = QtWidgets.QLabel(self.PorCCPantab)
        self.PostSampEd = QtWidgets.QLineEdit(self.PorCCPantab)
        self.PostSampLbl = QtWidgets.QLabel(self.PorCCPantab)
        self.PreSampEd = QtWidgets.QLineEdit(self.PorCCPantab)
        self.ClicksLab = QtWidgets.QLabel(self.PorCCPantab)
        self.PreSampLbl = QtWidgets.QLabel(self.PorCCPantab)

        self.DetectorTab = QtWidgets.QWidget()
        self.DetPan = QtWidgets.QFrame(self.DetectorTab)
        self.PreFiltFreq = QtWidgets.QLineEdit(self.DetPan)
        self.PreFiltLbl2 = QtWidgets.QLabel(self.DetPan)
        self.PreFiltPole = QtWidgets.QLineEdit(self.DetPan)
        self.PreFiltPoleLbl = QtWidgets.QLabel(self.DetPan)
        self.PreFiltDD = QtWidgets.QComboBox(self.DetPan)
        self.PreFiltDDLbl = QtWidgets.QLabel(self.DetPan)
        self.PreFiltLbl = QtWidgets.QLabel(self.DetPan)
        self.DetThreshold = QtWidgets.QLineEdit(self.DetPan)
        self.DetThresLbl = QtWidgets.QLabel(self.DetPan)
        self.MaxFreqEd = QtWidgets.QLineEdit(self.DetPan)
        self.MaxFreqLbl = QtWidgets.QLabel(self.DetPan)
        self.MinFreqEd = QtWidgets.QLineEdit(self.DetPan)
        self.MinFreqLbl = QtWidgets.QLabel(self.DetPan)
        self.TriggerFilLbl = QtWidgets.QLabel(self.DetPan)
        self.LongFiltEdit2 = QtWidgets.QLineEdit(self.DetPan)
        self.LongFiltLbl2 = QtWidgets.QLabel(self.DetPan)
        self.LongFiltEdit = QtWidgets.QLineEdit(self.DetPan)
        self.LongFiltLbl = QtWidgets.QLabel(self.DetPan)
        self.ShortFiltEdit = QtWidgets.QLineEdit(self.DetPan)
        self.ShortFiltLbl = QtWidgets.QLabel(self.DetPan)
        self.DetectorLbl = QtWidgets.QLabel(self.DetPan)
        self.CheckClassify = QtWidgets.QCheckBox(self.DetPan)

        # PorCC menu
        self.MenuPorCCF = WinTable()

        self.PorCCSetPan = QtWidgets.QFrame(self.MenuPorCCF)
        self.BrowsePorCC = QtWidgets.QPushButton(self.PorCCSetPan)
        self.DataTypeLabel = QtWidgets.QLabel(self.PorCCSetPan)
        self.FileDD = QtGui.QComboBox(self.PorCCSetPan)
        self.FolderPathPorCC = QtWidgets.QLineEdit(self.PorCCSetPan)
        self.SelectFolderText = QtWidgets.QLabel(self.PorCCSetPan)
        self.ParPanPorCC = QtWidgets.QFrame(self.PorCCSetPan)
        self.ProbPan = QtWidgets.QFrame(self.PorCCSetPan)

        self.dBLab = QtWidgets.QLabel(self.ParPanPorCC)
        self.SensEdit = QtWidgets.QLineEdit(self.ParPanPorCC)
        self.SensLab = QtWidgets.QLabel(self.ParPanPorCC)
        self.dBLab2 = QtWidgets.QLabel(self.ParPanPorCC)
        self.GainEdit = QtWidgets.QLineEdit(self.ParPanPorCC)
        self.GainLab = QtWidgets.QLabel(self.ParPanPorCC)
        self.kHzLab = QtWidgets.QLabel(self.ParPanPorCC)
        self.FsEdit = QtWidgets.QLineEdit(self.ParPanPorCC)
        self.SFreqLab = QtWidgets.QLabel(self.ParPanPorCC)
        self.EditHydN = QtWidgets.QLineEdit(self.ParPanPorCC)
        self.HydLabel = QtWidgets.QLabel(self.ParPanPorCC)
        self.DAQEditPorCC = QtWidgets.QLineEdit(self.ParPanPorCC)
        self.DAQLabPorCC = QtWidgets.QLabel(self.ParPanPorCC)


        self.LQLab = QtWidgets.QLabel(self.ProbPan)
        self.LQThres = QtWidgets.QLineEdit(self.ProbPan)
        self.LQ = QtWidgets.QLabel(self.ProbPan)
        self.HQLab = QtWidgets.QLabel(self.ProbPan)
        self.HQThres = QtWidgets.QLineEdit(self.ProbPan)
        self.HQ = QtWidgets.QLabel(self.ProbPan)

        self.NewCTFig = WinTable()
        self.NewCTPan = QtWidgets.QFrame(self.NewCTFig)

        self.CancelButtonCT = QtWidgets.QPushButton(self.NewCTPan)
        self.OKButtonCT = QtWidgets.QPushButton(self.NewCTPan)

        self.LocationPan = QtWidgets.QFrame(self.NewCTPan)
        self.NewCTBrowseButton = QtWidgets.QPushButton(self.NewCTPan)
        self.InclSubFoldersNewCT = QtWidgets.QCheckBox(self.NewCTPan)
        self.FolderPathNewCT = QtWidgets.QLineEdit(self.NewCTPan)
        self.SelectFolderText2 = QtWidgets.QLabel(self.NewCTPan)
        self.ParPanCT = QtWidgets.QFrame(self.NewCTPan)
        self.MinLLabel = QtWidgets.QLabel(self.ParPanCT)
        self.LongEdit = QtWidgets.QLineEdit(self.LocationPan)
        self.LongLable = QtWidgets.QLabel(self.LocationPan)
        self.LatEdit = QtWidgets.QLineEdit(self.LocationPan)
        self.LatLabel = QtWidgets.QLabel(self.LocationPan)

        self.OpenCTFig = WinTable()
        self.OpenCTPan = QtWidgets.QFrame(self.OpenCTFig)
        self.OpenCTCancelButton = QtWidgets.QPushButton(self.OpenCTPan)
        self.OpenCTButton = QtWidgets.QPushButton(self.OpenCTPan)


        self.SelectFolderText1 = QtWidgets.QLabel(self.OpenCTPan)

        self.SpectWindow = WinTable()

        self.SpecPan = QtWidgets.QFrame(self.SpectWindow)
        self.UpdateSpec = QtWidgets.QPushButton(self.SpecPan)
        self.OverSpecLbl = QtWidgets.QLabel(self.SpecPan)
        self.OverSpec = QtWidgets.QLineEdit(self.SpecPan)
        self.FFTSpecLbl = QtWidgets.QLabel(self.SpecPan)
        self.FFTSpec = QtWidgets.QLineEdit(self.SpecPan)
        self.SpecLbl = QtWidgets.QLabel(self.SpecPan)
        self.SpectAxes = pg.ImageView(self.SpecPan)
        self.WaveformLbl = QtWidgets.QLabel(self.SpecPan)
        self.WaveformSpec = pg.PlotWidget(self.SpecPan)
        # self.SpectAxes = pg.PlotWidget(self.SpecPan)

        self.NotesFig = WinTable()
        self.PanNotes = QtWidgets.QFrame(self.NotesFig)
        self.CancelNotes = QtWidgets.QPushButton(self.PanNotes)
        self.SaveNotes = QtWidgets.QPushButton(self.PanNotes)
        self.NotesText = QtWidgets.QTextEdit(self.PanNotes)
        self.NotesCTLabel = QtWidgets.QLabel(self.PanNotes)
        """
        TABLES TAB
        """
        self.TableDisplayTab = QtWidgets.QWidget()
        self.TablesPan = QtWidgets.QFrame(self.TableDisplayTab)

    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1700, 950)
        MainWindow.setMinimumSize(QtCore.QSize(1600, 950))
        MainWindow.setMaximumSize(QtCore.QSize(1600, 950))

        ###########################################################
        # MAIN TABS
        ###########################################################
        self.centralwidget.setObjectName("centralwidget")
        self.MainTab.setGeometry(QtCore.QRect(0, 0, 1422, 888))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.MainTab.sizePolicy().hasHeightForWidth())
        self.MainTab.setSizePolicy(sizePolicy)
        self.MainTab.setMinimumSize(QtCore.QSize(1600, 950))
        self.MainTab.setMaximumSize(QtCore.QSize(1600, 950))
        self.MainTab.setStyleSheet("background-color: rgb(240, 240, 238)")
        self.MainTab.setObjectName("MainTab")
        self.MainDisplayTab.setObjectName("MainDisplayTab")


        ####################
        # DISPLAY PARAMETERS
        ####################
        self.font = QtGui.QFont()
        self.font.setBold(True)
        self.font.setWeight(75)
        self.DisplaySettings.setGeometry(QtCore.QRect(10, 15, 1570, 122))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.DisplaySettings.sizePolicy().hasHeightForWidth())
        self.DisplaySettings.setSizePolicy(sizePolicy)
        self.DisplaySettings.setFrameShape(QtWidgets.QFrame.Box)
        self.DisplaySettings.setFrameShadow(QtWidgets.QFrame.Raised)
        self.DisplaySettings.setObjectName("DisplaySettings")
        # Date area
        self.DatePan.setGeometry(QtCore.QRect(10, 8, 250, 106))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.DatePan.sizePolicy().hasHeightForWidth())
        self.DatePan.setSizePolicy(sizePolicy)
        self.DatePan.setFrameShape(QtWidgets.QFrame.Box)
        self.DatePan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.DatePan.setObjectName("DatePan")
        self.DateandtimeofCTLabel.setGeometry(QtCore.QRect(20, 60, 200, 30))
        self.DateandtimeofCTLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.DateandtimeofCTLabel.setObjectName("DateandtimeofCTLabel")
        self.DateLabel.setGeometry(QtCore.QRect(100, 20, 50, 30))
        self.DateLabel.setTextFormat(QtCore.Qt.RichText)
        self.DateLabel.setObjectName("DateLabel")
        # Day / Night area
        self.DayPan.setGeometry(QtCore.QRect(270, 8, 120, 106))
        self.DayPan.setFrameShape(QtWidgets.QFrame.Box)
        self.DayPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.DayPan.setObjectName("DayPan")
        self.DayNightTextLabel.setGeometry(QtCore.QRect(20, 20, 90, 30))
        self.DayNightTextLabel.setTextFormat(QtCore.Qt.RichText)
        self.DayNightTextLabel.setObjectName("DayNightTextLabel")
        self.DayLabel.setGeometry(QtCore.QRect(30, 60, 40, 30))
        self.DayLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.DayLabel.setObjectName("DayLabel")
        # CT Info area
        self.CTInfoPan.setGeometry(QtCore.QRect(400, 8, 340, 106))
        self.CTInfoPan.setFrameShape(QtWidgets.QFrame.Box)
        self.CTInfoPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.CTInfoPan.setObjectName("CTInfoPan")
        self.CTNumLabel.setGeometry(QtCore.QRect(60, 20, 250, 30))
        self.CTNumLabel.setTextFormat(QtCore.Qt.RichText)
        self.CTNumLabel.setObjectName("CTNumLabel")
        self.CTBack.setGeometry(QtCore.QRect(30, 60, 40, 30))
        self.CTBack.setObjectName("CTBack")
        self.CTBack.clicked.connect(self.CTBackCB)
        self.CTNumD.setGeometry(QtCore.QRect(73, 60, 80, 30))
        self.CTNumD.setObjectName("CTNumD")
        self.CTForw.setGeometry(QtCore.QRect(156, 60, 40, 30))
        self.CTForw.setObjectName("CTForw")
        self.CTForw.clicked.connect(self.CTForwCB)
        self.TotalLabel.setGeometry(QtCore.QRect(200, 60, 50, 30))
        self.TotalLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.TotalLabel.setObjectName("TotalLabel")
        # Length area
        self.LengthPan.setGeometry(QtCore.QRect(750, 8, 90, 106))
        self.LengthPan.setFrameShape(QtWidgets.QFrame.Box)
        self.LengthPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.LengthPan.setObjectName("LengthPan")
        self.LengthText.setGeometry(QtCore.QRect(10, 20, 60, 30))
        self.LengthText.setTextFormat(QtCore.Qt.RichText)
        self.LengthText.setAlignment(QtCore.Qt.AlignCenter)
        self.LengthText.setObjectName("LengthText")
        self.LengthLabel.setGeometry(QtCore.QRect(20, 60, 50, 30))
        self.LengthLabel.setTextFormat(QtCore.Qt.PlainText)
        self.LengthLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.LengthLabel.setObjectName("LengthLabel")
        # CT Type area
        self.CTPan.setGeometry(QtCore.QRect(850, 8, 150, 106))
        self.CTPan.setFrameShape(QtWidgets.QFrame.Box)
        self.CTPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.CTPan.setObjectName("CTPan")
        self.CTLabel.setGeometry(QtCore.QRect(40, 20, 100, 30))
        self.CTLabel.setTextFormat(QtCore.Qt.RichText)
        self.CTLabel.setObjectName("CTLabel")
        self.CTTypeLabel.setGeometry(QtCore.QRect(20, 60, 100, 30))
        self.CTTypeLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.CTTypeLabel.setObjectName("CTTypeLabel")
        # Behaviour area
        self.BehavPan.setGeometry(QtCore.QRect(1010, 8, 150, 106))
        self.BehavPan.setFrameShape(QtWidgets.QFrame.Box)
        self.BehavPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.BehavPan.setObjectName("BehavPan")
        self.BehavTextLabel.setGeometry(QtCore.QRect(35, 20, 100, 30))
        self.BehavTextLabel.setTextFormat(QtCore.Qt.RichText)
        self.BehavTextLabel.setObjectName("BehavTextLabel")
        self.BehavLabel.setGeometry(QtCore.QRect(20, 60, 100, 30))
        self.BehavLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.BehavLabel.setObjectName("BehavLabel")
        # Calf area
        self.CalfPan.setGeometry(QtCore.QRect(1170, 8, 90, 106))
        self.CalfPan.setFrameShape(QtWidgets.QFrame.Box)
        self.CalfPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.CalfPan.setObjectName("CalfPan")
        self.CalfTextLabel.setGeometry(QtCore.QRect(25, 20, 50, 30))
        self.CalfTextLabel.setTextFormat(QtCore.Qt.RichText)
        self.CalfTextLabel.setObjectName("CalfTextLabel")
        self.CalfLabel.setGeometry(QtCore.QRect(25, 60, 30, 30))
        self.CalfLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.CalfLabel.setObjectName("CalfLabel")
        # Notes area
        self.NotesPan.setGeometry(QtCore.QRect(1270, 8, 160, 106))
        self.NotesPan.setFrameShape(QtWidgets.QFrame.Box)
        self.NotesPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.NotesPan.setObjectName("NotesPan")
        self.NoteLabel.setGeometry(QtCore.QRect(60, 20, 50, 30))
        self.NoteLabel.setObjectName("NoteLabel")
        self.SeeAddButton.setGeometry(QtCore.QRect(22, 55, 124, 40))
        self.SeeAddButton.setObjectName("SeeAddButton")
        self.SeeAddButton.clicked.connect(self.NotesCT)
        # Save button
        self.SaveupdatesButton.setGeometry(QtCore.QRect(1440, 8, 120, 106))
        self.SaveupdatesButton.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.SaveupdatesButton.setAutoDefault(False)
        self.SaveupdatesButton.setDefault(False)
        self.SaveupdatesButton.setFlat(False)
        self.SaveupdatesButton.setObjectName("SaveupdatesButton")
        self.SaveupdatesButton.clicked.connect(self.SaveUpdates)

        """
        AXES AREA
        """
        pg.setConfigOption('background', 'w')
        self.AxesPan.setGeometry(QtCore.QRect(10, 150, 1200, 700))
        self.AxesPan.setFrameShape(QtWidgets.QFrame.Box)
        self.AxesPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.AxesPan.setObjectName("AxesPan")
        # AMPLITUDE
        self.AmplitudeDB.setGeometry(QtCore.QRect(550, 7, 150, 20))
        self.AmplitudeDB.setObjectName("AmplitudeDB")
        # Axes to plot STEM
        self.AmpAxesCT.setGeometry(QtCore.QRect(20, 35, 1160, 192))
        self.AmpAxesCT.setObjectName("AmpAxesCT")
        # ICI / CPS
        self.ICIAxesCT.setGeometry(QtCore.QRect(20, 239, 1160, 214))
        self.ICIAxesCT.setObjectName("ICIAxesCT")
        self.ICIPan.setGeometry(QtCore.QRect(300, 214, 650, 42))
        self.ICIPan.setFrameShape(QtWidgets.QFrame.Box)
        self.ICIPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.ICIPan.setObjectName("ICIPan")
        # CPS radio button
        self.ClickspersecondButton.setGeometry(QtCore.QRect(50, 10, 200, 20))
        self.ClickspersecondButton.setFont(self.font)
        self.ClickspersecondButton.setObjectName("ClickspersecondButton")
        self.ClickspersecondButton.setChecked(True)
        # ICI radio button
        self.InterclickintervalmsButton.setGeometry(QtCore.QRect(400, 10, 200, 20))
        self.InterclickintervalmsButton.setFont(self.font)
        self.InterclickintervalmsButton.setObjectName("InterclickintervalmsButton")
        # FREQUENCY / BEARING
        self.FreqAxesCT.setGeometry(QtCore.QRect(20, 464, 1160, 210))
        self.FreqAxesCT.setObjectName("FreqAxesCT")
        self.FreqPan.setGeometry(QtCore.QRect(300, 434, 650, 42))
        self.FreqPan.setFrameShape(QtWidgets.QFrame.Box)
        self.FreqPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.FreqPan.setObjectName("FreqPan")
        # Frequency radio button
        self.CentroidfrequencykHzButton.setGeometry(QtCore.QRect(50, 10, 200, 20))
        self.CentroidfrequencykHzButton.setFont(self.font)
        self.CentroidfrequencykHzButton.setObjectName("CentroidfrequencykHzButton")
        self.CentroidfrequencykHzButton.setChecked(True)
        # Bearing radio button
        self.DirectionofarrivalButton.setGeometry(QtCore.QRect(400, 10, 200, 20))
        self.DirectionofarrivalButton.setFont(self.font)
        self.DirectionofarrivalButton.setObjectName("DirectionofarrivalButton")
        """
        ACTION AREA
        """
        self.ActionPan.setGeometry(QtCore.QRect(1220, 150, 360, 700))
        self.ActionPan.setFrameShape(QtWidgets.QFrame.Box)
        self.ActionPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.ActionPan.setObjectName("ActionPan")
        self.ActionLabel.setGeometry(QtCore.QRect(50, 10, 200, 20))
        self.ActionLabel.setFont(self.font)
        self.ActionLabel.setObjectName("ActionLabel")
        #  Display CT type
        self.DisplayAllNBHF.setGeometry(QtCore.QRect(20, 38, 320, 74))
        self.DisplayAllNBHF.setFrameShape(QtWidgets.QFrame.Box)
        self.DisplayAllNBHF.setFrameShadow(QtWidgets.QFrame.Raised)
        self.DisplayAllNBHF.setObjectName("DisplayAllNBHF")
        self.ActionTextLabel.setGeometry(QtCore.QRect(120, 10, 80, 20))
        self.ActionTextLabel.setFont(self.font)
        self.ActionTextLabel.setObjectName("ActionTextLabel")
        self.AllButton.setGeometry(QtCore.QRect(40, 35, 100, 20))
        self.AllButton.setObjectName("AllButton")
        self.AllButton.setChecked(True)
        self.NBHFButton.setGeometry(QtCore.QRect(180, 35, 100, 20))
        self.NBHFButton.setObjectName("NBHFButton")
        # Individual click trains
        self.CT3DPan.setGeometry(QtCore.QRect(20, 140, 320, 242))
        self.CT3DPan.setFrameShape(QtWidgets.QFrame.Box)
        self.CT3DPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.CT3DPan.setObjectName("CT3DPan")
        self.IndCTLabel.setGeometry(QtCore.QRect(30, 10, 250, 20))
        self.IndCTLabel.setFont(self.font)
        self.IndCTLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.IndCTLabel.setObjectName("IndCTLabel")
        self.IndClicksAnd3D.setGeometry(QtCore.QRect(30, 60, 260, 62))
        self.font2.setBold(False)
        self.font2.setWeight(50)
        # Click train in 3D
        self.IndClicksAnd3D.setFont(self.font2)
        self.IndClicksAnd3D.setObjectName("IndClicksAnd3D")
        self.IndClicksAnd3D.clicked.connect(self.OpenIndClicksAnd3D)
        self.SpectrogramButton.setGeometry(QtCore.QRect(30, 140, 260, 62))
        # Spectrogram area
        self.SpectrogramButton.setFont(self.font2)
        self.SpectrogramButton.setObjectName("SpectrogramButton")
        self.SpectrogramButton.clicked.connect(self.CreateSpectrogram)
        # Validation
        self.ValidatePan.setGeometry(QtCore.QRect(20, 402, 320, 270))
        self.ValidatePan.setFrameShape(QtWidgets.QFrame.Box)
        self.ValidatePan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.ValidatePan.setObjectName("ValidatePan")
        self.ValidateLab.setGeometry(QtCore.QRect(50, 20, 200, 20))
        self.ValidateLab.setFont(self.font)
        self.ValidateLab.setAlignment(QtCore.Qt.AlignCenter)
        self.ValidateLab.setObjectName("ValidateLab")
        # CT Type
        self.CTTypeDropDown.setGeometry(QtCore.QRect(120, 60, 170, 40))
        self.CTTypeDropDown.setObjectName("CTTypeDropDown")
        self.CTTypeDropDown.addItem("Select")
        self.CTTypeDropDown.addItem("NBHF")
        self.CTTypeDropDown.addItem("LQ-NBHF")
        self.CTTypeDropDown.addItem("Non-NBHF")
        self.CTTypeDropDown.addItem("Sonar")
        self.ValTypeLab.setGeometry(QtCore.QRect(20, 70, 70, 20))
        self.ValTypeLab.setObjectName("ValTypeLab")
        # Behaviour
        self.BehaviourDropDown.setGeometry(QtCore.QRect(120, 130, 170, 40))
        self.BehaviourDropDown.setObjectName("BehaviourDropDown")
        self.BehaviourDropDown.addItem("Select")
        self.BehaviourDropDown.addItem("Orientation")
        self.BehaviourDropDown.addItem("Feeding")
        self.BehaviourDropDown.addItem("Socialising")
        self.BehaviourDropDown.addItem("Unknown")
        self.ValBehavLab.setGeometry(QtCore.QRect(20, 140, 70, 20))
        self.ValBehavLab.setObjectName("ValBehavLab")
        # Validate button
        self.ValButton.setGeometry(QtCore.QRect(30, 200, 260, 50))
        self.ValButton.setObjectName("ValButton")
        self.ValButton.clicked.connect(self.Validate)
        ########################################################
        # METRICS DISPLAY
        #######################################################
        self.MainTab.addTab(self.MainDisplayTab, "")
        self.MetricsDisplayTab.setObjectName("MetricsDisplayTab")
        self.MetricsUploadPan.setGeometry(QtCore.QRect(10, 8, 560, 122))
        self.MetricsUploadPan.setFrameShape(QtWidgets.QFrame.Box)
        self.MetricsUploadPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.MetricsUploadPan.setObjectName("MetricsUploadPan")
        self.UploadData.setGeometry(QtCore.QRect(440, 8, 112, 106))
        self.UploadData.setFont(self.font)
        self.UploadData.setObjectName("UploadData")
        self.UploadData.clicked.connect(self.UploadMetricData)
        self.SelectFolderMetricEdit.setGeometry(QtCore.QRect(20, 40, 400, 30))
        self.SelectFolderMetricEdit.setStyleSheet("QLineEdit{background: white}")
        self.SelectFolderMetricEdit.setObjectName("SelectFolderMetricEdit")
        self.SelectFolderMetricEdit.setText('C:/')
        self.SelectFoldLab.setGeometry(QtCore.QRect(20, 10, 180, 20))
        self.SelectFoldLab.setFont(self.font)
        self.SelectFoldLab.setObjectName("SelectFoldLab")
        # Browse area
        self.MetricBrowseButton.setGeometry(QtCore.QRect(300, 74, 120, 40))
        self.MetricBrowseButton.setFont(self.font)
        self.MetricBrowseButton.setObjectName("MetricBrowseButton")

        # TODO SHOULD CALL MetricBrowse.py here
        self.MetricBrowseButton.clicked.connect(self.MetricsBrowse)

        """
        Metrics Table
        """
        self.MetricsTablePan.setGeometry(QtCore.QRect(10, 140, 560, 700))
        self.MetricsTablePan.setFrameShape(QtWidgets.QFrame.Box)
        self.MetricsTablePan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.MetricsTablePan.setObjectName("MetricsTablePan")
        self.MetricsTable.setGeometry(QtCore.QRect(20, 20, 520, 660))
        self.MetricsTable.setObjectName("MetricsTable")
        self.MetricsTable.setRowCount(20)
        self.MetricsTable.setColumnCount(12)
        item = QtWidgets.QTableWidgetItem()
        self.MetricsTable.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.MetricsTable.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.MetricsTable.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.MetricsTable.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.MetricsTable.setHorizontalHeaderItem(4, item)
        item = QtWidgets.QTableWidgetItem()
        self.MetricsTable.setHorizontalHeaderItem(5, item)
        item = QtWidgets.QTableWidgetItem()
        self.MetricsTable.setHorizontalHeaderItem(6, item)
        item = QtWidgets.QTableWidgetItem()
        self.MetricsTable.setHorizontalHeaderItem(7, item)
        item = QtWidgets.QTableWidgetItem()
        self.MetricsTable.setHorizontalHeaderItem(8, item)
        item = QtWidgets.QTableWidgetItem()
        self.MetricsTable.setHorizontalHeaderItem(9, item)
        item = QtWidgets.QTableWidgetItem()
        self.MetricsTable.setHorizontalHeaderItem(10, item)
        item = QtWidgets.QTableWidgetItem()
        self.MetricsTable.setHorizontalHeaderItem(11, item)

        """
        Metrics summary
        """
        self.SummaryPan.setGeometry(QtCore.QRect(580, 140, 1000, 700))
        self.SummaryPan.setFrameShape(QtWidgets.QFrame.Box)
        self.SummaryPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.SummaryPan.setObjectName("SummaryPan")
        self.CPMAxesTimes.setGeometry(10, 10, 480, 680)
        self.CPMAxesTimes.setVisible(False)
        self.CPMAxesTot.setGeometry(500, 10, 480, 680)
        self.CPMAxesTot.setVisible(False)
        self.CPMAxesMetrics.setGeometry(10, 10, 980, 680)
        self.MetricsDisplayPan.setGeometry(QtCore.QRect(580, 8, 1000, 122))
        self.MetricsDisplayPan.setFrameShape(QtWidgets.QFrame.Box)
        self.MetricsDisplayPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.MetricsDisplayPan.setObjectName("MetricsDisplayPan")
        self.SelectDatePan.setGeometry(QtCore.QRect(10, 8, 270, 106))
        self.SelectDatePan.setFrameShape(QtWidgets.QFrame.Box)
        self.SelectDatePan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.SelectDatePan.setObjectName("SelectDatePan")
        self.SelectadateDropDown.setGeometry(QtCore.QRect(20, 54, 220, 40))
        self.SelectadateDropDown.setObjectName("SelectadateDropDown")
        self.SelectadateDropDown.addItem("All days")
        self.SelectDateLab.setGeometry(QtCore.QRect(20, 20, 100, 20))
        self.SelectDateLab.setObjectName("SelectDateLab")
        self.MetricPPMPan.setGeometry(QtCore.QRect(290, 8, 500, 106))
        self.MetricPPMPan.setFrameShape(QtWidgets.QFrame.Box)
        self.MetricPPMPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.MetricPPMPan.setObjectName("MetricPPMPan")
        """
        Positive porpoise minutes
        """
        self.PPMPan.setGeometry(QtCore.QRect(10, 8, 310, 90))
        self.PPMPan.setFrameShape(QtWidgets.QFrame.Box)
        self.PPMPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.PPMPan.setObjectName("PPMPan")
        self.PPMButton.setGeometry(QtCore.QRect(20, 55, 250, 16))
        self.PPMButton.setFont(self.font)
        self.PPMButton.setObjectName("PPMButton")
        self.TypeofCTButton.setGeometry(QtCore.QRect(20, 20, 150, 20))
        self.TypeofCTButton.setFont(self.font)
        self.TypeofCTButton.setObjectName("TypeofCTButton")
        self.TypeofCTButton.setChecked(True)
        self.BehaviourButton.setGeometry(QtCore.QRect(150, 20, 150, 20))
        self.BehaviourButton.setFont(self.font)
        self.BehaviourButton.setObjectName("BehaviourButton")
        # Day / Night
        self.AllDayPan.setGeometry(QtCore.QRect(330, 8, 160, 90))
        self.AllDayPan.setFrameShape(QtWidgets.QFrame.Box)
        self.AllDayPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.AllDayPan.setObjectName("AllDayPan")
        self.AllButtonMetrics.setGeometry(QtCore.QRect(20, 20, 50, 20))
        self.AllButtonMetrics.setObjectName("AllButtonMetrics")
        self.AllButtonMetrics.setFont(self.font)
        self.AllButtonMetrics.setChecked(True)
        self.NightDayButton.setGeometry(QtCore.QRect(20, 55, 130, 20))
        self.NightDayButton.setObjectName("NightDayButton")
        self.NightDayButton.setFont(self.font)
        # Update and clear buttons
        self.DisplayMetricsButton.setGeometry(QtCore.QRect(796, 8, 88, 106))
        self.DisplayMetricsButton.setObjectName("DisplayMetricsButton")
        self.DisplayMetricsButton.clicked.connect(self.DisplayMetrics)
        self.ClearButton.setGeometry(QtCore.QRect(890, 8, 100, 106))
        self.ClearButton.setObjectName("ClearButton")
        self.MainTab.addTab(self.MetricsDisplayTab, "")
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        # TABLES
        self.TablesPan.setGeometry(QtCore.QRect(10, 140, 560, 700))
        self.TablesPan.setFrameShape(QtWidgets.QFrame.Box)
        self.TablesPan.setFrameShadow(QtWidgets.QFrame.Raised)
        self.TablesPan.setObjectName("TablesPan")

        """
        MENU
        """
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1422, 18))
        self.menubar.setObjectName("menubar")
        self.menuMain_Display.setObjectName("menuMain_Display")
        # self.menuMetrics_Display.setObjectName("menuMetrics_Display")
        # Click trains menu
        self.menuClick_Trains.setObjectName("menuClick_Trains")
        self.MenuNewCT.setObjectName("MenuNewCT")
        self.MenuOpenCT.setObjectName("MenuOpenCT")
        # Downloads
        self.menuDownloads.setObjectName("menuDownloads")
        self.menuHelp.setObjectName("menuHelp")
        MainWindow.setMenuBar(self.menubar)
        self.MenuSetDetector.setObjectName("MenuSetDetector")
        # self.MenuSetPorCC.setObjectName("MenuSetPorCC")
        self.menuMain_Display.addSeparator()
        self.menuMain_Display.addAction(self.MenuSetDetector)
        # self.menuMetrics_Display.addAction(self.MenuSetPorCC)
        self.menuClick_Trains.addAction(self.MenuNewCT)
        self.menuClick_Trains.addAction(self.MenuOpenCT)
        self.menubar.addAction(self.menuMain_Display.menuAction())
        # self.menubar.addAction(self.menuMetrics_Display.menuAction())
        self.menubar.addAction(self.menuClick_Trains.menuAction())
        self.menubar.addAction(self.menuDownloads.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())
        self.retranslateUi(MainWindow)
        self.MainTab.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        self.MenuOpenCT.triggered.connect(lambda: self.OpenCTMenu())
        self.MenuNewCT.triggered.connect(lambda: self.NewCTMenu())
        # self.MenuSetPorCC.triggered.connect(lambda: self.OpenPorCCSetMenu())
        self.MenuSetDetector.triggered.connect(lambda: self.OpenDetSetMenu())

    """
    FUNCTIONS
    """
    def PushPorCCButton(self):
        pass
        # self.MenuPorCCF.close()
        # # TODO Message that it is working
        # # f = uifigure('Position', [300 300 300 150])
        # # UpdateMes = uilabel(f, 'Position', [10 10 280 130])
        # # Message = 'PorCC is at work! \n A message will pop up when ready.'
        # ##UpdateMes.Text = Message
        # # uialert(UIFigure, Message, 'PorCC Running', 'Icon', 'warning')
        # # TabGroup.Enable = 'off'
        #
        # ##Create a table with 0s to make it faster and more effective
        # if HydN > 1:
        #     ColunmNames = ['ID', 'Date', 'start_sample', 'Duration', 'Amp', 'PF', 'CF' \
        #                'BW', 'RMSBW', 'Bearing']
        # else:
        #     ColunmNames = ['ID', 'Date', 'start_sample', 'Duration', 'Amp', 'PF', 'CF' \
        #                'BW', 'RMSBW', 'ICI', 'CT', 'InFile', 'ClickN', 'Class']
        #
        #         # end
        # CP = pd.DataFrame(0, index=np.arange(2000000), columns=ColunmNames)
        #
        # # Get list of all subfolders. \
        # allSubFolders = os.listdir(SelectedFolder) # wasgenpath(SelectedFolder)
        # # Parse into a cell array. \
        # remain = allSubFolders
        # listOfFolderNames = []
        # WhileVar = 1
        # while WhileVar == 1:
        #     [singleSubFolder, remain] = [int(x) for x in remain.split('') if x.strip()]
        #     if not singleSubFolder:
        #         WhileVar = 2
        #     # end
        #     listOfFolderNames = [listOfFolderNames, singleSubFolder]
        # # end
        # numberOfFolders = len(listOfFolderNames)
        #
        # # Get this folder and print it out.
        # thisFolder = listOfFolderNames[1]
        #
        # # how many binary files in this folder
        # clickfiles = self.findBinaryFiles(thisFolder)
        # # Now we have a list of all files in this folder.
        # NumberBinaryFiles = len(clickfiles)
        #
        # # Running the classifier FFT = 512 RowN = 0 for p in NumberBinaryFiles: # each pgdf file click_id = 0  #
        # Click Id within the pgdf file # print out because lots of files and script takes a very long time. #
        # Message = strcat('Loading file pgdf', ' ', num2str(p), ' of a total of  ', ' ', num2str(length(
        # clickfiles))) # UpdateMes.Text = Message
        #
        #     # LOAD THE FILES
        #     fileName = SelectedFolder + '/' + clickfiles[p]
        #    # fid = fopen(fileName, 'r', 'ieee-be.l64')
        #
        #     ### TRY THIS
        #     fid = open(fileName, 'rb')
        #     # dim = np.fromfile(fid, dtype='>u4')
        #
        #     # initialize variables
        #     prevPos = -1
        #     Cols = []
        #     fileInfo = ['readModuleHeader', 'readModuleFooter']
        #     fileInfo.iloc[0,0] = self.readStdModuleHeader(fid) # @readStdModuleHeader
        #     fileInfo.iloc[0,1] = self.readStdModuleFooter(fid) # @readStdModuleFooter
        #
        #     # Create a file with information about each file
        #     # FilesInfo.FilNum(p, 1) = p
        #     # FilesInfo.FilName(p, 1) = filename
        #     # FilesInfo.Fs(p, 1) = Fs
        #     # FilesInfo.Lat(p, 1) = Lat
        #     # FilesInfo.Long(p, 1) = Long
        #     # Main loop = goes into.pgdf files to find clicks
        #     while 1:
        #         # if for some reason we're stuck at one byte, warn the user and abort
        #         pos = fid.tell() # ftell()
        #         if pos == prevPos:
        #             print('File stuck at byte %d', pos)
        #             break
        #         # end
        #         prevPos = pos
        #
        #         # read in the length of the object and the type. Move the menufile
        #         # pointer back to the beginning of the length variable, and switch
        #         # on the type. If we couldn't read nextLen or nextType, assume
        #         # that means we've hit the end of the menufile and break out of loop
        #         # [~, nL] = fread(fid, 1, 'int32')
        #         # [nextType, nT] = fread(fid, 1, 'int32')
        #         nLtemp = fid.read(1)
        #         nL = nLtemp[0]
        #         NTemp = fid.read(2)
        #         nT = NTemp[0]
        #         nextType = NTemp[1]
        #         if nL == 0 or nT == 0:
        #             break
        #
        #         fid.seek(-8, 0)
        #
        #         if nextType == -1:
        #             # Case 1: MenuFile Header.Read in the menufile header, and then set the function handles
        #             # depending on what type of data the binary menufile holds. The module
        #             # type is typically specified in the package class that extends PamControlledUnit, and is a
        #             # string unique to that module.
        #             fileInfo.fileHeader = self.readFileHeader(fid)
        #             Type = fileInfo.fileHeader.moduleType
        #             # Click Detector Module
        #                 if Type == 'Click Detector':
        #                     NewCase = fileInfo.fileHeader.streamName
        #                         if NewCase == 'Clicks':
        #                             fileInfo.objectType = 1000
        #                             fileInfo.readModuleData = self.readClickData()
        #                             fileInfo.readModuleFooter = self.readClickFooter()
        #                         elif NewCase == 'Trigger Background':
        #                             fileInfo.objectType = 0
        #                             fileInfo.readModuleData = self.readClickTriggerData()
        #                             fileInfo.readModuleHeader = self.readClickTriggerHeader()
        #                     # end
        #
        #                 else:
        #                     print(['Don''t recognize type ' fileInfo.fileHeader.moduleType '.  Aborting load'])
        #                     break
        #             # end
        #
        #         elif nextType == -2:
        #             # Case 2: MenuFile Footer.The menufile version should have been set
        #             # when we read the menufile header. If the menufile header is empty,
        #             # something has gone wrong so warn the user and exit
        #             if not fileInfo.fileHeader:
        #                 print('Error: found file footer before file header.  Aborting load')
        #                 break
        #             # end
        #             fileInfo.fileFooter = readFileFooterInfo(fid, fileInfo.fileHeader.fileFormat)
        #             # Case 3: Module Header. The correct function handle should have
        #             # been set when we read the menufile header. If the menufile header
        #             # is empty, something has gone wrong so warn the user and exit
        #
        #         elif nextType == -3:
        #             if not fileInfo.fileHeader:
        #                 print('Error: found module header before file header.  Aborting load')
        #                 break
        #
        #             fileInfo.moduleHeader = fileInfo.readModuleHeader(fid)
        #         elif nextType == -4:
        #             if not fileInfo.fileHeader:
        #                 print('Error: found module footer before file header.  Aborting load')
        #                 break
        #
        #             fileInfo.moduleFooter = fileInfo.readModuleFooter(fid)
        #         elif nextType == -5:
        #             # Case 5: Data. The correct function handle should have been set when we
        #             # read the menufile header. If the menufile header is empty something
        #             # has gone wrong so warn the user and exit
        #         else:
        #             if not fileInfo.fileHeader:
        #                 print('Error: found data before file header.  Aborting load')
        #                 break
        #             # end
        #             dataPoint = self.readPamData(fid, fileInfo)
        #             # check if the datapoint has a click in it
        #             a = isfield(dataPoint, 'wave')
        #             if a == 0:
        #                 continue
        #
        #           # If there are two hydrophones, use the hydrophone where the click
        #         # impinges first, because it is expected to be of the best quality
        #         if HydN > 1:
        #             if dataPoint.delays < 0:
        #                 click = dataPoint.wave[:, 1]
        #             else:
        #                 click = dataPoint.wave[:, 2]
        #             # end
        #         else:
        #             click = dataPoint.wave
        #         # end
        #         # Start counting clicks (impulsive sound), regardless of whether we use it
        #         click_id = click_id + 1
        #
        # # PorCC Classification # PF = peak frequency, CF = centroid frequency # PSD = power spectrum, f = frequency
        # [PSD, f, CF_Val, PF_Val] = self.CentroidAndPeakFrq1(click, FFT, Fs) if 160000 > PF_Val > 100000 < CF_Val <
        # 160000  # Møhl & Andersen, 1973 # RMS - BW = RMS bandwidth, ratio = PF / CF, BW = -3 dB bandwidth #
        # duration = estimated using 80 #energy of the clip # Q = CF / RMSBW, XC = peak cross - correlation
        # coefficient with a model click RMSBW_Val = math.sqrt(sum((f - CF_Val)**2. * PSD**2) / sum(PSD**2))) / 1000
        # durationV = self.clickDuration(click) Q_Val = (CF_Val / RMSBW_Val) / 1000 Ratio_Val = PF_Val/ CF_Val XCorr
        # = np.correlate(click, click_model) XCoeff = XCorr[0] XC_Val = max(XCoeff) BW_Val = powerbw(click,
        # Fs) / 1000 if Q_Val >= 4 and durationV < 450: # estimate probabilities for HQ ProbHQ = glmval(logitCoefHQ,
        # [Q_Val durationV], 'logit') # assign each clip to a category if ProbHQ >= ThresHQ: Porps_Val = 1  # HQ
        # click else: ProbLQ = glmval(logitCoefLQ, [Q_Val durationV Ratio_Val XC_Val CF_Val BW_Val], 'logit')  # Log
        # Model 2 if ProbLQ > ThresLQ: Porps_Val = 2  # LQ click else: Porps_Val = 3  # noise # end # end else:
        # Porps_Val = 3  # if Q < 4 and duration > 450 # end else: Porps_Val = 3  # if peak is outside porpoise range
        # end #if Porps_Val... # store in ClickParameters if Class HQ or LQ if Porps_Val == 1 or Porps_Val == 2 or
        # Porps_Val == 3: RowN = RowN + 1  # our ID ###### save HQ & LQ clicks #####% CP.ID[RowN] = RowN CP.Date[
        # RowN] = dataPoint.date CP.start_sample[RowN] = dataPoint.start_sample CP.Duration[RowN] = durationV
        # ##DURATION CP.Amp_Val = self.clickAmplitude(click, HydSens, Gain, DAQpp) CP.Amp[RowN] = Amp_Val  # MAXIMUM
        # AMPLITUDE CP.PF[RowN] = PF_Val  # PEAK FREQUENCY CP.CF[RowN] = CF_Val  # CENTROID FREQUENCY CP.BW[RowN] =
        # BW_Val  # BANDWIDTH - 3 dB CP.RMSBW[RowN] = RMSBW_Val  # RMS - BW if HydN > 1:  # and exist(
        # dataPoint.delays) CP.Bearing[RowN] = dataPoint.angles * 57.2958 # end if RowN == 1: ici = 0 else: ici = (
        # CP.start_sample[RowN] - CP.start_sample[RowN-1])/(Fs/1000) #in microseconds if ici < 0: ici = (CP.Date[RowN]
        # - CP.Date[RowN-]) * 120 * 24 #in milliseconds # end # end P.ICI[RowN] = ici # PORPOISE OR NOT? CT[RowN] = 0
        # b = fileName[end - 56:-1] InFile[RowN] = b ClickN[RowN] = click_id Class[RowN] = Porps_Val # end if
        # Porps_Val ~ = 3 # end switch # end # while (1) fclose(fid) # end # for Binary files
        #
        #     ###### Remove the ends
        #     ID[RowN+1:-1] = []
        #     Date[RowN+1:-1] = []
        #     start_sample[RowN+1:-1] = []
        #     Duration[RowN+1:-1] = []  ##DURATION
        #     Amp[RowN+1:-1] = []   # MAXIMUM AMPLITUDE
        #     PF[RowN+1:-1] = []  # PEAK FREQUENCY
        #     CF[RowN+1:-1] = []   # CENTROID FREQUENCY
        #     BW[RowN+1:-1] = []  # BANDWIDTH - 3 dB
        #     RMSBW[RowN+1:-1] = []
        #     if RowN < 2000000:
        #         ICI[RowN+1:-1] = []
        #     else:
        #         ICI[RowN+1:-1] = []              # end
        #     CT[RowN+1:-1] = []
        #     InFile[RowN+1:-1] = []
        #     ClickN[RowN+1:-1] = []
        #     Class[RowN+1:-1] = []
        #
        # if HydN == 1: ClickParameters = table(ID, Date, start_sample, Duration, Amp, PF, CF, BW, RMSBW, ICI, CT,
        # InFile, ClickN, Class) elif HydN == 2: Bearing[RowN+1:-1] = [] ClickParameters = table(ID, Date,
        # start_sample, Duration, Amp, PF, CF, BW, RMSBW, ICI, Bearing, CT, InFile, ClickN, Class) # end # S CP[
        # RowN+1:-1] = [] ###REMOVE CLICKS WITH THE SAME DATE # ZeroRep = find(diff(ClickParameters.Date) == 0) # if
        # size(ZeroRep, 1) > 0 # d = 1 # while size(d, 1) > 0 # a = find(diff(ClickParameters.Date) == 0) # b = find(
        # diff(ClickParameters.Amp) < -3) # c = arrayfun( @ (x) find(b == x, 1), a, 'UniformOutput', False) # c =
        # cell2mat(c) # d = b(c) # d = d + 1 ## remove the same date but amp < first # ClickParameters(d,
        # :) = [] # end # d = 1 # while size(d, 1) > 0 # a = find(diff(ClickParameters.Date) == 0) # b = find(diff(
        # ClickParameters.Amp) > 3) # c = arrayfun( @ (x) find(b == x, 1), a, 'UniformOutput', false) # c = cell2mat(
        # c) # d = b(c) ##remove the same date but amp < first # ClickParameters(d,:) = [] # end # ICI = diff(
        # ClickParameters.start_sample) / (Fs / 1000) # ClickParameters.ICI(1, 1) = 0 # ClickParameters.ICI(2: end,
        # 1) = ICI # else # ICI = diff(ClickParameters.start_sample) / (Fs / 1000) # ClickParameters.ICI(1,
        # 1) = 0 # ClickParameters.ICI(2: end, 1) = ICI # end
        #
        # # ICI CP = self.NewICI(self, CP, Fs) DateTime = strftime(CP.Date.loc[1]) yyyy = DateTime[8:11] mmm =
        # DateTime[4:6] dd = DateTime[1:2] hh = DateTime[13:14] mm = DateTime[16:17] ss = DateTime[19:20] DateTimeEnd
        # = strftime(CP.Date.loc[-1]) hhEnd = DateTimeEnd[13:14] mmEnd = DateTimeEnd[16:17] ssEnd = DateTimeEnd[
        # 19:20] FileDateTimePar = strcat(SelectedFolder, '\ClickParameters', yyyy, mmm, dd, '_', hh, mm, ss, '_',
        # hhEnd, mmEnd, ssEnd, '.mat') ClickParameters = self.CP save(FileDateTimePar, 'ClickParameters')
        #
        #         # message to let know it's finished
        #         # message = sprintf('The file was now generated, you can start seeing the data.')
        #         # uialert(UIFigure, message, 'Done!', 'Icon', 'success')
        #         # TabGroup.Enable = 'on'
        #         # end Of click OK Button for PorCC

    #
    # def MetricsOK(self):
    #     pass
    # #     # search dates within the selected folder
    # #     # Name of folders are 20150812 -> 12 Aug 2015
    # #     # Create a list of days and display in dropdown button
    # #     # Create table with dates
    # #     Metrics = pd.DataFrame()
    # #     FoldersIn = uigetdir(MetricSelectedFolder)
    # #     Folders = FoldersIn([FoldersIn.isdir] == 1)
    # #     if len(Folders) > 3:
    # #         for d in Folders:
    # #             NameDate = Folders.name[d]
    # #             yyyy = int(NameDate[1:4])
    # #             mmm = int(NameDate[5:6])
    # #             dd = int(NameDate([7:8])
    # #             Metrics.Date[d-2] = strftime(datetime(yyyy, mmm, dd))
    # #             # Continue later getting more metrics
    # #             # Metrics(d - 2, 2) =
    # #             # Metrics(d - 2, 3) =
    # #     else:
    # #         NameDate = MetricSelectedFolder[end-7:-1]
    # #         yyyy = int(NameDate[1:4])
    # #         mmm = int(NameDate[5:6])
    # #         dd = int(NameDate[7:8])
    # #         Metrics.Date[0,0] = strftime(datetime(yyyy, mmm, dd))
    # #
    # # # Upload data button pushed

    # #         ## SAVE POSITIVE PORPOISE MINUTES
    # #         AllPPMFile = [topLevelFolderMetrics + '\PosPorMin.mat']
    # #         save (AllPPMFile, 'PosPorMin')
    # #
    # #         RowNCT = 0
    # #         SummaryTable = pd.DataFrame()
    # #         for FolderN in numberOfFoldersMetrics:
    # #             # Get this folder and print it out.
    # #             thisFolder = listOfFolderNamesMetrics.iloc[FolderN]
    # #             ClickParNames = [thisFolder + '\ClickParameters*.mat']
    # #             ClickParFilesInFolder = uigetdir(ClickParNames)
    # #             ClickParFileName = ClickParFilesInFolder(1).name
    # #             ClickParFile = strcat(thisFolder, '\', ClickParFileName)
    # #             CP = pd.read_csv(ClickParFile)
    # #
    # #             LDate = len(CP)
    # #             DATE = CP.Date(round(LDate / 2), 1)
    # #             if LAT != 0 and LONG != 0
    # #                 [srise, sset, ~] = sunrise(LAT, LONG, 0, [], DATE)
    # #
    # #                 # Sunrise
    # #                 SunriseTimFe = strftime(srise)
    # #                 HourSunR = SunriseTime(13:14)
    # #                 HourSunR = int(HourSunR)
    # #                 MinSunR = SunriseTime(16:17)
    # #                 MinSunR = int(MinSunR)
    # #                 SunRSum = (HourSunR * 60) + MinSunR
    # #                 # Sunset
    # #                 SunsetTime = strftime(srise)
    # #                 HourSunS = SunsetTime(13:14)
    # #                 HourSunS = int(HourSunS)
    # #                 MinSunS = SunsetTime(16:17)
    # #                 MinSunS = int(MinSunS)
    # #                 SunSSum = (HourSunS * 60) + MinSunS
    # #             # end
    # #
    # #         # Open CTInfo and extract them
    # #         CTInfoName = [thisFolder + '\CTInfoNew.mat']
    # #         CTInfo = pd.read_csv(CTInfoName)
    # #
    # #         ## CTInfo NBHF decision
    # #         RowNCT = RowNCT + 2
    # #         # CT types
    # #         NBHF = CTInfo[CTInfo.Species == 'NBHF']
    # #         LQNBHF = CTInfo[CTInfo.Species == 'LQ-NBHF']
    # #         NonNBHF = CTInfo(strcmp({CTInfo.Species}, 'Non-NBHF'))
    # #         Sonar = CTInfo(strcmp({CTInfo.Species}, 'Sonar'))
    # #         Orient = CTInfo(strcmp({CTInfo.Behav}, 'Orientation'))
    # #         Forage = CTInfo(strcmp({CTInfo.Behav}, 'Foraging'))
    # #         Comm = CTInfo(strcmp({CTInfo.Behav}, 'Socialising'))
    # #         Unknown = CTInfo(strcmp({CTInfo.Behav}, 'Unknown'))
    # #
    # #         # Calves = CTInfo(strcmp({CTInfo.Calf}, 'Yes'))
    # #         # NoCalves = length(NBHF) - length(Calves) % CTInfo(strcmp({CTInfo.Calf}, 'No'))
    # #         cla(CPMAxesMetrics, 'reset')
    # #
    # #         # Day time
    # #         SummaryTable.Date[RowNCT] = DATE
    # #         SummaryTable.RealDate[RowNCT] = {strftime(SummaryTable.Date(RowNCT, 1), 'dd-mmm-yyyy')}
    # #         SummaryTable.NBHF[RowNCT] = sum(strcmp({NBHF.DayNight}, 'Day'))
    # #         SummaryTable.LQNBHF[RowNCT] = sum(strcmp({LQNBHF.DayNight}, 'Day'))
    # #         SummaryTable.NonNBHF[RowNCT] = sum(strcmp({NonNBHF.DayNight}, 'Day'))
    # #         SummaryTable.Sonar[RowNCT] = sum(strcmp({Sonar.DayNight}, 'Day'))
    # #         SummaryTable.Orient[RowNCT] = sum(strcmp({Orient.DayNight}, 'Day'))
    # #         SummaryTable.Forage[RowNCT] = sum(strcmp({Forage.DayNight}, 'Day'))
    # #         SummaryTable.Social[RowNCT] = sum(strcmp({Comm.DayNight}, 'Day'))
    # #         SummaryTable.Unknown[RowNCT] = sum(strcmp({Unknown.DayNight}, 'Day'))
    # #         SummaryTable.Day[RowNCT] = 1
    # #         SummaryTable.Night[RowNCT] = 0
    # #         # ClickTrains.Calves[RowNCT] = sum(strcmp({Calves.DayNight}, 'Day'))
    # #         # night time
    # #         SummaryTable.Date[RowNCT+1] = DATE
    # #         SummaryTable.RealDate[RowNCT+1] = {strftime(SummaryTable.Date(RowNCT + 1, 1), 'dd-mmm-yyyy')}
    # #         SummaryTable.NBHF[RowNCT+1] = [NBHF.DayNight == 'Night'].sum()
    # #         SummaryTable.LQNBHF[RowNCT+1] = LQNBHF.DayNight == 'Night'].sum()
    # #         SummaryTable.NonNBHF[RowNCT+1] = onNBHF.DayNight =='Night'].sum()
    # #         SummaryTable.Sonar[RowNCT+1] = Sonar.DayNight == 'Night'].sum()
    # #         SummaryTable.Orient[RowNCT+1] = Orient.DayNight == 'Night'].sum()
    # #         SummaryTable.Forage[RowNCT+1] = Forage.DayNight == 'Night'].sum()
    # #         SummaryTable.Social[RowNCT+1] = Comm.DayNight == 'Night'].sum()
    # #         SummaryTable.Unknown[RowNCT+1] = Unknown.DayNight == 'Night'].sum()
    # #         SummaryTable.Day[RowNCT+1] = 0
    # #         SummaryTable.Night[RowNCT+1] = 1
    # #         # ClickTrains.Calves(RowNCT + 1, 1) = sum(strcmp({Calves.DayNight}, 'Night'))
    # #
    # #         #end # For - going into each folder end
    # #         # If it already exists or not \
    # #         ClickTrainsName = [topLevelFolderMetrics + '\SummaryTable.mat']
    # #         save(ClickTrainsName, 'SummaryTable')
    # #
    # #         # DataDislplayLabel.Text = OriginalText
    # #         DaysList = SummaryTable.RealDate.to_list()
    # #         DaysList = np.unique(DaysList)
    # #         DaysList = sort(DaysList)
    # #
    # #         # list of days
    # #         DropDown = ['All days' + DaysList]
    # #         self.SelectadateDropDown.Items = DropDown
    # #         self.SelectadateDropDown.Enable = 'on'
    # #         DisplayMetricsButton.Enable = 'on'
    # #
    # #         ## Variables I need to display
    # #         NBHF = SummaryTable['NBHF'].sum()
    # #         LQNBHF = SummaryTable['LQNBHF'].sum()
    # #         NonNBHF = SummaryTable['NonNBHF'].sum()
    # #         Sonar = SummaryTable['Sonar'].sum()
    # #         Orient = SummaryTable['Orient'].sum()
    # #         Forage = SummaryTable['Forage'].sum()
    # #         Comm = SummaryTable['Social'].sum()
    # #         Unknown = SummaryTable['Unknown'].sum()
    # #         # Calves = sum(ClickTrains.Calves)
    # #         # NoCalves = length(NoCalves)
    # #
    # #         # Day vs Night
    # #         Night = SummaryTable(SummaryTable.Night == 1,:)
    # #         Day = SummaryTable(SummaryTable.Day == 1,:)
    # #         # Day
    # #         self.NBHFd = Day['NBHF'].sum()
    # #         self.LQNBHFd = sum(Day.LQNBHF)
    # #         self.NonNBHFd = sum(Day.NonNBHF)
    # #         self.NoiseD = sum(Day.Sonar)
    # #         self.TravelD = sum(Day.Orient)
    # #         self.ForageD = sum(Day.Forage)
    # #         self.CommD = sum(Day.Social)
    # #         self.UnknownD = sum(Day.Unknown)
    # #         self.CalvesD = sum(Day.Calves)
    # #         # Night
    # #         self.NBHFn = sum(Night.NBHF)
    # #         self.LQNBHFn = sum(Night.LQNBHF)
    # #         self.NonNBHFn = sum(Night.NonNBHF)
    # #         self.NoiseN = sum(Night.Sonar)
    # #         self.TravelN = sum(Night.Orient)
    # #         self.ForageN = sum(Night.Forage)
    # #         self.CommN = sum(Night.Social)
    # #         self.UnknownN = sum(Night.Unknown)
    # #         self.CalvesN = sum(Night.Calves)

    def DisplayMetrics(self):
        # clear axis (not possible)
        # What to display
        DaysToPlot = self.SelectadateDropDown.getText()
        if self.TypeofCTButton:
            if DaysToPlot == 'All Days':
                if self.AllButtonMetrics:
                    TypesCT = [sum(self.SummaryTable.NBHF), sum(self.SummaryTable.LQNBHF),
                               sum(self.SummaryTable.NonNBHF), sum(self.SummaryTable.Sonar)]
                    plt.bar(self.CPMAxesMetrics, TypesCT)
                else:
                    TypesCT = [sum(self.SummaryTable.NBHF[SummaryTable.Day == 1]),
                               sum(self.SummaryTable.LQNBHF[SummaryTable.Day == 1]),
                               sum(self.SummaryTable.NonNBHF[SummaryTable.Day == 1]),
                               sum(self.SummaryTable.Sonar[SummaryTable.Day == 1])]
                    plt.bar(self.CPMAxesMetrics, TypesCT)
        elif self.BehaviourButton:
            if DaysToPlot == 'All Days':
                if self.AllButtonMetrics:
                    BehavCT = [sum(self.SummaryTable.Orient), sum(self.SummaryTable.Forage),
                               sum(self.SummaryTable.Social)]
                    plt.bar(self.CPMAxesMetrics, BehavCT)
                else:
                    BehavCT = [sum(self.SummaryTable.Orient[SummaryTable.Day == 1]),
                               sum(self.SummaryTable.Forage[SummaryTable.Day == 1]),
                               sum(self.SummaryTable.Social[SummaryTable.Day == 1])]
                    plt.bar(self.CPMAxesMetrics, BehavCT)
            else:
                ThisDay = self.SummaryTable[SummaryTable.Date == DaysToPlot]
                if self.AllButtonMetrics:
                    BehavCT = [sum(ThisDay.Orient), sum(ThisDay.Forage), sum(ThisDay.Social)]
                    plt.bar(self.CPMAxesMetrics, BehavCT)
                else:
                    BehavCT = [sum(ThisDay.Orient[ThisDay.Day == 1]), sum(ThisDay.Forage[SummaryTable.Day == 1]),
                               sum(ThisDay.Social[ThisDay.Day == 1])]
                    plt.bar(self.CPMAxesMetrics, BehavCT)
        else:  # self.PosPorMinButton
            for j in range(0, numberOfFoldersMetrics - 1):
                for i in range(0, 1440):
                    if PosPorMin.Positive[i, j] == 1:
                        plt.plot(self.CPMAxesTimes, self.SummaryTable.Minutes[i, 0], j, '.', 'color', 'black')
                    # end
                # end
            # end
            for j in range(0, numberOfFoldersMetrics - 1):
                Folder = self.listOfFolderNamesMetrics[j + 1]
                DateName = Folder[-1 - 7:-1]
                FDay = DateName[7:8]
                Month = DateName[5:6]
                FDate = FDay + '/' + Month
                FDate = {FDate}
                self.CPMAxesTimes.setGraphXLabel(FDate)
            # end
            self.CPMAxesTimes.setGraphXLabel('Min of the day')

            for j in range(0, numberOfFoldersMetrics - 1):
                plt.bar(self.CPMAxesTot, j, sum(PosPorMin.Positive[:, j]))
            # end
            self.CPMAxesTot.setGraphYLabel('Positive minutes')

    # Translate
    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "D-PorCCA"))
        self.DateandtimeofCTLabel.setText(_translate("MainWindow", "12 Aug 2015, 12.35.32"))
        self.DayNightTextLabel.setText(_translate("MainWindow",
                                                  "<html><head/><body><p><span style=\" "
                                                  "font-weight:600\">Day/Night</span></p></body></html>"))
        self.DayLabel.setText(_translate("MainWindow", "Day"))
        self.CTNumLabel.setText(_translate("MainWindow",
                                           "<html><head/><body><p><span style=\" font-weight:600\">Click train number "
                                           "(Total)</span></p></body></html>"))
        self.CTForw.setText(_translate("MainWindow", ">"))
        self.CTBack.setText(_translate("MainWindow", "<"))
        self.TotalLabel.setText(_translate("MainWindow", "()"))
        self.LengthText.setText(_translate("MainWindow",
                                           "<html><head/><body><p><span style=\" "
                                           "font-weight:600\">Length</span></p></body></html>"))
        self.LengthLabel.setText(_translate("MainWindow", "20"))
        self.SaveupdatesButton.setText(_translate("MainWindow", " Save \n"
                                                                "updates"))
        self.CTLabel.setText(_translate("MainWindow",
                                        "<html><head/><body><p><span style=\" font-weight:600\">CT "
                                        "Type</span></p></body></html>"))
        self.CTTypeLabel.setText(_translate("MainWindow", "NBHF"))
        self.BehavTextLabel.setText(_translate("MainWindow",
                                               "<html><head/><body><p><span style=\" "
                                               "font-weight:600\">Behaviour</span></p></body></html>"))
        self.BehavLabel.setText(_translate("MainWindow", "Socialising"))
        self.CalfTextLabel.setText(_translate("MainWindow",
                                              "<html><head/><body><p><span style=\" "
                                              "font-weight:600\">Calf</span></p></body></html>"))
        self.CalfLabel.setText(_translate("MainWindow", "No"))
        self.NoteLabel.setText(_translate("MainWindow",
                                          "<html><head/><body><p><span style=\" "
                                          "font-weight:600\">Notes</span></p></body></html>"))
        self.SeeAddButton.setText(_translate("MainWindow", "Add/Modify"))
        self.DateLabel.setText(_translate("MainWindow",
                                          "<html><head/><body><p><span style=\" "
                                          "font-weight:600\">Date</span></p></body></html>"))
        self.ActionLabel.setText(_translate("MainWindow", "ACTIONS - CLICK TRAINS"))
        self.ActionTextLabel.setText(_translate("MainWindow", "Display"))
        self.AllButton.setText(_translate("MainWindow", "All"))
        self.NBHFButton.setText(_translate("MainWindow", "NBHF"))
        self.IndCTLabel.setText(_translate("MainWindow", "Individual click trains - 3D"))
        self.IndClicksAnd3D.setText(_translate("MainWindow", "Click train in 3D"))
        self.SpectrogramButton.setText(_translate("MainWindow", "Spectrogram"))
        self.ValidateLab.setText(_translate("MainWindow", "Validation/changes"))
        self.ValButton.setText(_translate("MainWindow", "Validate"))
        # self.BehaviourDropDown.setText(_translate("MainWindow", "Select"))
        # self.CTTypeDropDown.setText(_translate("MainWindow", "Select"))
        self.ValTypeLab.setText(_translate("MainWindow", "CT Type"))
        self.ValBehavLab.setText(_translate("MainWindow", "Behaviour"))
        self.AmplitudeDB.setText(_translate("MainWindow",
                                            "<html><head/><body><p><span style=\" font-weight:600\">Amplitude ("
                                            "dB)</span></p></body></html>"))
        self.InterclickintervalmsButton.setText(_translate("MainWindow", "Inter-click interval"))
        self.ClickspersecondButton.setText(_translate("MainWindow", "Clicks per second"))
        self.CentroidfrequencykHzButton.setText(_translate("MainWindow", "Centroid Frequency"))
        self.DirectionofarrivalButton.setText(_translate("MainWindow", "Direction of arrival"))
        self.MainTab.setTabText(self.MainTab.indexOf(self.MainDisplayTab), _translate("MainWindow", "Main Display"))
        self.UploadData.setText(_translate("MainWindow", "UPLOAD \n"
                                                         " DATA"))
        self.SelectFoldLab.setText(_translate("MainWindow", "Select project folder"))
        self.MetricBrowseButton.setText(_translate("MainWindow", "Browse"))
        #  self.SelectadateDropDown.setText(_translate("MainWindow", "All Days"))
        self.SelectDateLab.setText(_translate("MainWindow", "Select dates"))
        self.AllButtonMetrics.setText(_translate("MainWindow", "All"))
        self.NightDayButton.setText(_translate("MainWindow", "Day/Night"))
        self.PPMButton.setText(_translate("MainWindow", "Positive Porpoise Minute"))
        self.TypeofCTButton.setText(_translate("MainWindow", "Type of CT"))
        self.BehaviourButton.setText(_translate("MainWindow", "Behaviour"))
        self.DisplayMetricsButton.setText(_translate("MainWindow", "UPDATE"))
        self.ClearButton.setText(_translate("MainWindow", "CLEAR ALL"))
        self.MainTab.setTabText(self.MainTab.indexOf(self.MetricsDisplayTab), _translate("MainWindow", "Metrics"))
        self.MainTab.setTabText(self.MainTab.indexOf(self.TableDisplayTab), _translate("MainWindow", "Tables"))
        self.menuMain_Display.setTitle(_translate("MainWindow", "Detector"))
        # self.menuMetrics_Display.setTitle(_translate("MainWindow", "PorCC"))
        self.menuClick_Trains.setTitle(_translate("MainWindow", "Click Trains"))
        self.menuDownloads.setTitle(_translate("MainWindow", "Downloads"))
        self.menuHelp.setTitle(_translate("MainWindow", "Help"))
        self.MenuSetDetector.setText(_translate("MainWindow", "Settings Detector"))
        # self.MenuSetPorCC.setText(_translate("MainWindow", "Settings PorCC"))
        self.MenuNewCT.setText(_translate("MainWindow", "New CT Project"))
        self.MenuOpenCT.setText(_translate("MainWindow", "Open existing project"))
        item = self.MetricsTable.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "Date"))
        item = self.MetricsTable.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "NBHF"))
        item = self.MetricsTable.horizontalHeaderItem(2)
        item.setText(_translate("MainWindow", "LQNBHF"))
        item = self.MetricsTable.horizontalHeaderItem(3)
        item.setText(_translate("MainWindow", "NonNBHF"))
        item = self.MetricsTable.horizontalHeaderItem(4)
        item.setText(_translate("MainWindow", "Sonar"))
        item = self.MetricsTable.horizontalHeaderItem(5)
        item.setText(_translate("MainWindow", "Orient"))
        item = self.MetricsTable.horizontalHeaderItem(6)
        item.setText(_translate("MainWindow", "Forage"))
        item = self.MetricsTable.horizontalHeaderItem(7)
        item.setText(_translate("MainWindow", "Social"))
        item = self.MetricsTable.horizontalHeaderItem(8)
        item.setText(_translate("MainWindow", "Unknown"))
        item = self.MetricsTable.horizontalHeaderItem(9)
        item.setText(_translate("MainWindow", "Day"))
        item = self.MetricsTable.horizontalHeaderItem(10)
        item.setText(_translate("MainWindow", "Night"))
        item = self.MetricsTable.horizontalHeaderItem(11)
        item.setText(_translate("MainWindow", "Calf"))
        # self.Update = UpdateCT.UpdatePlots

    """
    Display buttons
    """

    def CTBackCB(self):
        NumCT = int(self.CTNumD.text())
        First = CTInfo['CTNum'].iloc[0]
        if NumCT > First:
            self.AmpAxesCT.clear()
            self.ICIAxesCT.clear()
            self.FreqAxesCT.clear()
            RowCT = CTInfo[CTInfo.CTNum == NumCT].index[0]
            NumCT = int(CTInfo.CTNum[RowCT - 1])
            self.UpdateCT(NumCT, CP, CTInfo)

    def CTForwCB(self):
        NumCT = int(self.CTNumD.text())
        Tot = CTInfo['CTNum'].iloc[-1]
        if NumCT == Tot:
            a = 1  # do nothing
        elif NumCT < Tot:
            self.AmpAxesCT.clear()
            self.ICIAxesCT.clear()
            self.FreqAxesCT.clear()
            RowCT = CTInfo[CTInfo.CTNum == NumCT].index[0]
            # print(RowCT)
            NumCT = int(CTInfo.CTNum[RowCT + 1])
            # print(NumCT)
            self.UpdateCT(NumCT, CP, CTInfo)

    def UpdateCT(self, NumCT, Cp, myCTInfo):
        global CTTemp
        CTTemp = Cp[Cp.CT == NumCT]
        # print(CTTemp)
        if NumCT != 1:
            CTTemp.reset_index(inplace=True, drop=True)
        if len(CTTemp) > 9:
            fs = 1 / (CTTemp.iloc[3]['ICI'] / (
                    1000 * (CTTemp.iloc[3]["start_sample"] - CTTemp.iloc[2]["start_sample"])))
            CTTemp = NewICI(CTTemp, fs)
            CTTemp.loc[:, 'SumMs'] = int(0)
            # print(CTTemp)
            for i in range(1, len(CTTemp)):
                CTTemp.SumMs[i] = int(CTTemp.SumMs[i - 1]) + int(CTTemp.ICI[i])
            CTTemp.SumMs = CTTemp.SumMs / 1000
            CT1HQ = CTTemp[CTTemp['pyPorCC'] == 1]
            CT1LQ = CTTemp[CTTemp['pyPorCC'] == 2]
            self.CTNumD.setText(str(NumCT))
            self.CTTypeLabel.setText(str(myCTInfo.Species[myCTInfo.CTNum == NumCT].values[0]))
            self.DateandtimeofCTLabel.setText(str(myCTInfo.Date[myCTInfo.CTNum == NumCT].values[0]))
            self.LengthLabel.setText(str(len(CTTemp)))
            self.BehavLabel.setText(str(myCTInfo.Behav[myCTInfo.CTNum == NumCT].values[0]))
            self.CalfLabel.setText(str(myCTInfo.Calf[myCTInfo.CTNum == NumCT].values[0]))
            self.TotalLabel.setText('(' + str(myCTInfo['CTNum'].iloc[-1]) + ')')
            self.DayLabel.setText(str(myCTInfo.DayNight[myCTInfo.CTNum == NumCT].values[0]))
            # if CTInfo["Saved"][[CTInfo.CTNum] == NumCT] == 1:
            #    SelectedCTCheckBox.Value = 1

            self.AmpAxesCT.clear()
            # TODO set the max and min in all axis (think of SoundTrap data)
            WidthBar = max(CTTemp.SumMs) / 500
            AmpLinesLQ = pg.BarGraphItem(x=CT1LQ.SumMs, height=CT1LQ.amplitude, brush='b', width=WidthBar)
            AmpLinesHQ = pg.BarGraphItem(x=CT1HQ.SumMs, height=CT1HQ.amplitude, brush='r', width=WidthBar)
            AmpDotsLQ = pg.ScatterPlotItem(x=CT1LQ.SumMs, y=CT1LQ.amplitude, symbol='o', brush='b', width=0.2)
            AmpDotsHQ = pg.ScatterPlotItem(x=CT1HQ.SumMs, y=CT1HQ.amplitude, symbol='o', brush='r', width=0.2)
            self.AmpAxesCT.addItem(AmpLinesLQ)
            self.AmpAxesCT.addItem(AmpDotsLQ)
            self.AmpAxesCT.addItem(AmpLinesHQ)
            self.AmpAxesCT.addItem(AmpDotsHQ)
            self.AmpAxesCT.setXRange(0, max(CTTemp.SumMs) + 0.1)
            self.AmpAxesCT.setYRange(80, 170)

            # plot click per second (default) or ICI
            ICIorCPS = self.InterclickintervalmsButton.isChecked()

            if ICIorCPS == 1:
                ICILQ = CT1LQ.ICI.to_list()
                ICIHQ = CT1HQ.ICI.to_list()
                ICIDotsLQ = pg.ScatterPlotItem(x=CT1LQ.SumMs, y=ICILQ, symbol='o', brush='b', width=2)
                ICIDotsHQ = pg.ScatterPlotItem(x=CT1HQ.SumMs, y=ICIHQ, symbol='o', brush='r', width=2)
                self.ICIAxesCT.addItem(ICIDotsLQ)
                self.ICIAxesCT.addItem(ICIDotsHQ)
                self.ICIAxesCT.setXRange(0, max(CTTemp.SumMs) + 0.1)
                self.ICIAxesCT.setYRange(0, max(CTTemp.ICI[2:-1]) + 10)
            else:
                # plot clicks per second
                CPSLQ = CT1LQ.CPS.to_list()
                CPSHQ = CT1HQ.CPS.to_list()
                CPSDotsLQ = pg.ScatterPlotItem(x=CT1LQ.SumMs, y=CPSLQ, symbol='o', brush='b', width=2)
                CPSDotsHQ = pg.ScatterPlotItem(x=CT1HQ.SumMs, y=CPSHQ, symbol='o', brush='r', width=2)
                self.ICIAxesCT.addItem(CPSDotsLQ)
                self.ICIAxesCT.addItem(CPSDotsHQ)
                self.ICIAxesCT.setXRange(0, max(CTTemp.SumMs) + 0.1)
                self.ICIAxesCT.setYRange(0, max(CTTemp.CPS[2:-1]) + 30)
            a = 1
            if a == 0:  # If bearing exist
                # TODO identifying if bearing exists and plot
                CTTemp.Bearing = CTTemp.Bearing.to_list()
                self.FreqAxesCT.plot(CTTemp.SumMs, CTTemp.Bearing, pen=None, symbol='o', color='b')
                self.FreqAxesCT.setXRange(0, max(CTTemp.SumMs) + 0.1)
                self.FreqAxesCT.setYRange(0, 180)
            #     FreqOrBearing = DirectionofarrivalButton.Value
            #     BearExist = strcmpi('Bearing', CP.Properties.VariableNames)
            #     if FreqOrBearing == 1 and len(BearExist) > 0:
            else:
                FreqLQ = CT1LQ.CF / 1000
                FreqHQ = CT1HQ.CF / 1000
                FreqLQ = FreqLQ.to_list()
                FreqHQ = FreqHQ.to_list()
                FreqLinesLQ = pg.BarGraphItem(x=CT1LQ.SumMs, height=FreqLQ, brush='b', width=WidthBar)
                FreqLinesHQ = pg.BarGraphItem(x=CT1HQ.SumMs, height=FreqHQ, brush='r', width=WidthBar)
                FreqDotsLQ = pg.ScatterPlotItem(x=CT1LQ.SumMs, y=FreqLQ, symbol='o', brush='b', width=2)
                FreqDotsHQ = pg.ScatterPlotItem(x=CT1HQ.SumMs, y=FreqHQ, symbol='o', brush='r', width=2)
                self.FreqAxesCT.addItem(FreqLinesLQ)
                self.FreqAxesCT.addItem(FreqDotsLQ)
                self.FreqAxesCT.addItem(FreqLinesHQ)
                self.FreqAxesCT.addItem(FreqDotsHQ)
                self.FreqAxesCT.setXRange(0, max(CTTemp.SumMs) + 0.1)
                self.FreqAxesCT.setYRange(50, 180)

    """
    Buttons
    """

    def NotesCT(self):
        global CTInfo
        CTNum = int(self.CTNumD.text())
        Notes = CTInfo.Notes[CTInfo.CTNum == CTNum].values[0]
        self.NotesFig.setGeometry(200, 200, 400, 500)
        self.NotesFig.setWindowTitle('Notes')
        self.NotesFig.setFixedSize(400, 500)
        # Create panel within figure
        self.PanNotes.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.PanNotes.setGeometry(10, 10, 380, 480)
        # Label with information about the click train itself
        self.NotesCTLabel = QtWidgets.QLabel(self.PanNotes)
        self.NotesCTLabel2 = QtWidgets.QLabel(self.PanNotes)
        self.NotesCTLabel.setText('Date: ' + str(CTInfo.Date[CTInfo.CTNum == CTNum].values[0]))
        self.NotesCTLabel.setGeometry(20, 15, 250, 30)
        self.NotesCTLabel2.setText('Click train: ' + str(CTNum))
        self.NotesCTLabel2.setGeometry(20, 50, 150, 30)
        # Create panel within figure
        self.NotesText.setGeometry(10, 100, 360, 330)
        self.NotesText.setText(Notes)
        self.SaveNotes.setGeometry(40, 440, 100, 30)
        self.SaveNotes.setText('Save')
        self.SaveNotes.clicked.connect(self.PushSaveNButton)
        self.CancelNotes.setGeometry(230, 440, 100, 30)
        self.CancelNotes.setText('Cancel')
        self.CancelNotes.clicked.connect(self.PushCanNButton)
        self.NotesFig.show()

    def PushCanNButton(self):
        self.NotesFig.close()

    def PushSaveNButton(self):
        CTNum = int(self.CTNumD.text())
        NotesUpdate = self.NotesText.toPlainText()
        CTInfo.Notes[CTInfo.CTNum == CTNum] = NotesUpdate
        self.NotesFig.close()

    def OpenIndClicksAnd3D(self):
        global CTTemp
        x1 = []
        y1 = []
        z1 = []
        self.Fig3D.setGeometry(50, 50, 1200, 800)
        self.Fig3D.setWindowTitle('Click train in 3D')
        self.Pan3D.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.Pan3D.setGeometry(10, 10, 1180, 780)
        self.Plot3D.setGeometry(600, 300, 570, 460)

        ## WATERFALL PLOT
        FFT = 512
        CTTemp['iciSum'] = 0
        for i in range(1, len(CTTemp)):
            CTTemp.iciSum[i] = CTTemp.ICI[i] + CTTemp.iciSum[i - 1]
        WavFileToOpen = CTTemp.filename[0]

        # CREATE empty x, y, and z for the waterfall plot
        z1 = np.zeros((257, len(CTTemp)), dtype=float)  # power of each click in each row

        # FILL the variables
        # X = frequency
        click1, fs = soundfile.read(WavFileToOpen, start=1, stop=150)
        freqs, psd = signal.welch(click1, fs=fs, window='hann', nfft=512)
        x1 = freqs / 1000
        # x1.to_numpy()
        # Y = time( in secs)
        y1 = CTTemp.iciSum.to_numpy()
        y1 = y1 / 1000
        # y1 = y1.T

        # Normalised Amplitude
        for i in range(0, len(CTTemp)):
            Start = CTTemp.start_sample[i]
            End = Start + CTTemp.duration_samples[i]
            click, Fs = soundfile.read(WavFileToOpen, start=int(Start), stop=int(End))
            freqs, psd = signal.welch(click, fs=Fs, window='hann', nfft=512)
            z1[:, i] = psd
        a = z1.max()
        z1 = z1 / a
        # ax = plt.axes(projection='3d')

        # Data for a three-dimensional line
        # ax.contour3D(x1, y1, z1)
        # ax.set_xlabel('Frequency (kHz)')
        # ax.set_ylabel('Time (ms)')
        # ax.set_zlabel('Amplitude')

        # Data for three-dimensional scattered points
        # zdata = 15 * np.random.random(100)
        # xdata = np.sin(zdata) + 0.1 * np.random.randn(100)
        # ydata = np.cos(zdata) + 0.1 * np.random.randn(100)
        # ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens')

        # self.Plot3D.add_subplot(111, projection='3d')

        # n = 0
        print(len(y1), len(x1), len(z1))
        p2 = gl.GLSurfacePlotItem(y=y1, x=x1, z=z1, shader='shaded', color=(0.5, 0.5, 1, 1))
        p2.translate(-10, -30, 10)
        p2.scale(1.0, 1.0, 0.5)
        self.Plot3D.addItem(p2)

        # self.Plot3D.contour3D(x1, y1, z1)
        # self.Plot3D.set_xlabel('Frequency (kHz)')
        # self.Plot3D.set_ylabel('Time (ms)')
        # self.Plot3D.set_zlabel('Amplitude')
        # plt.scatter(self.Plot3D, x1, y1, z1)
        self.Fig3D.show()
        self.Plot3D.show()

    def CreateSpectrogram(self):
        global CTTemp, Fs
        self.SpectWindow.setGeometry(200, 200, 1000, 600)
        self.SpectWindow.setWindowTitle('Spectrogram')
        # Create panel within figure
        self.SpecPan.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.SpecPan.setGeometry(10, 10, 980, 580)
        # where the plot will be
        self.WaveformSpec.setGeometry(160, 20, 800, 250)
        self.WaveformLbl.setGeometry(200, 0, 150, 20)
        self.WaveformLbl.setText('Waveform')
        self.SpectAxes.setGeometry(QtCore.QRect(160, 300, 800, 250))
        self.SpecLbl.setGeometry(200, 270, 150, 20)
        self.SpecLbl.setText('Spectrogram')
        self.FFTSpec.setGeometry(90, 300, 50, 30)
        self.FFTSpec.setText('512')
        self.FFTSpecLbl.setGeometry(20, 300, 60, 30)
        self.FFTSpecLbl.setText('FFT')
        self.OverSpec.setGeometry(90, 340, 50, 30)
        self.OverSpec.setText('128')
        self.OverSpecLbl.setGeometry(20, 340, 60, 30)
        self.OverSpecLbl.setText('Overlap')
        self.UpdateSpec.setGeometry(20, 390, 120, 40)
        self.UpdateSpec.setText('Update')
        self.UpdateSpec.clicked.connect(self.UpdateSpect)

        FileToOpen = CTTemp.filename[0]
        Start = CTTemp.start_sample.iloc[0] - 50000
        if Start < 0:
            Start = 0
        End = CTTemp.start_sample.iloc[-1] + 50000
        Signal, Fs = soundfile.read(FileToOpen, start=int(Start), stop=int(End))
        MeanSig = sum(Signal) / len(Signal)
        Signal = Signal - MeanSig
        # Signal = Signal/max(Signal)
        sos = signal.butter(10, 2000, 'hp', fs=Fs, output='sos')
        self.FiltSig = signal.sosfilt(sos, Signal)
        # the signal
        Duration = len(Signal) / Fs
        t = np.arange(0.0, Duration, 1 / Fs)
        NFFT = int(self.FFTSpec.text())  # length of the windowing segments
        Overlap = 128  # int(self.OverSpec.text())
        self.WaveformSpec.plot(t, self.FiltSig)
        # window = signal.get_window('hann', NFFT)

        fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
        ax1.plot(t, self.FiltSig)
        Pxx, freqs, bins, im = ax2.specgram(self.FiltSig, NFFT=NFFT, Fs=Fs, noverlap=Overlap)  # , cmap='jet')
        ax1.grid(False)
        ax2.grid(False)
        # ax1.ylabel('Amplitude')
        # ax1.title('Waveform')
        # ax2.ylabel('Frequency (Hz)')
        # ax2.xlabel('Time (s)')
        # ax2.title('Spectrogram')
        plt.show()

        # freq, t, Pxx = signal.spectrogram(self.FiltSig, fs=Fs, nfft=NFFT, window=window, scaling='density',
        # noverlap=Overlap) Pxx = 10*np.log10(Pxx**2) self.SpectAxes.setImage(Pxx.T, autoRange=False, scale=(100,
        # 600)) ColorMap = pg.ColorMap(pos=[0, 0.25, 0.5, 0.75, 1], color=[(0, 0, 255), (200, 255, 255), (100, 255,
        # 150), (255, 255, 0), (255, 0, 100)]) self.SpectAxes.setColorMap(ColorMap) self.SpectWindow.show()

    def UpdateSpect(self):
        NFFT = int(self.FFTSpec.text())  # length of the windowing segments
        Overlap = int(self.OverSpec.text())
        window = signal.get_window('hann', NFFT)
        freq, t, Pxx = signal.spectrogram(self.FiltSig, fs=Fs, nfft=NFFT, window=window, scaling='density',
                                          noverlap=Overlap)
        Pxx = 10 * np.log10(Pxx ** 2)
        Ratio = NFFT / Overlap
        x = (Ratio ** 2) * 7
        self.SpectAxes.setImage(Pxx.T, autoRange=False, scale=(x, 600))
        ColorMap = pg.ColorMap(pos=[0, 0.25, 0.5, 0.75, 1], color=[(0, 0, 255), (200, 255, 255), (100, 255, 150),
                                                                   (255, 255, 0), (255, 0, 100)])
        self.SpectAxes.setColorMap(ColorMap)

    def Validate(self):
        global CTInfo
        CTNum = int(self.CTNumD.text())
        Species = self.CTTypeDropDown.currentText()
        Behaviour = self.BehaviourDropDown.currentText()
        # if Species == 'Select' and Behaviour == 'Select':
        #     a = 1
        # el
        if not Species == 'Select' and Behaviour == 'Select':
            CTInfo.Species[CTInfo.CTNum == CTNum] = Species
            self.CTTypeLabel.setText(str(Species))
        elif Species == 'Select' and not Behaviour == 'Select':
            CTInfo.Behav[CTInfo.CTNum == CTNum] = Behaviour
            self.BehavLabel.setText(str(Behaviour))
        elif not Species == 'Select' and not Behaviour == 'Select':
            CTInfo.Species[CTInfo.CTNum == CTNum] = Species
            CTInfo.Behav[CTInfo.CTNum == CTNum] = Behaviour
            self.BehavLabel.setText(str(Behaviour))
            self.CTTypeLabel.setText(str(Species))

    """
    METRICS STUFF
    """

    def MetricsBrowse(self):
        global BrowseSelectedFolder
        root = tk.Tk()
        root.withdraw()
        BrowseSelectedFolder = filedialog.askdirectory()
        self.SelectFolderMetricEdit.setText(BrowseSelectedFolder)

    """
    MENUS
    """

    def OpenPorCCSetMenu(self):
        pass
        # self.MenuPorCCF.setGeometry(200, 200, 360, 550)
        # self.MenuPorCCF.setWindowTitle('PorCC Settings')
        # # Create panel within figure
        # self.PorCCSetPan.setFrameShape(QtWidgets.QFrame.StyledPanel)
        # self.PorCCSetPan.setGeometry(10, 10, 340, 530)
        # self.SelectFolderText.setGeometry(20, 10, 100, 30)
        # self.SelectFolderText.setText('Select Folder')
        # # Edit path
        # self.FolderPathPorCC.setGeometry(20, 40, 300, 30)
        # self.FolderPathPorCC.setText('C://')
        # # Dropdown
        # self.FileDD.setGeometry(150, 75, 170, 30)
        # self.FileDD.addItem("PAMGuard output")
        # self.FileDD.addItem("SoundTrap output")
        # #   self.FileDD.activated[str].connect(self.style_choice)
        # self.DataTypeLabel.setGeometry(20, 75, 100, 30)
        # self.DataTypeLabel.setText('Data type')
        # # browse
        # self.BrowsePorCC.setGeometry(220, 115, 100, 30)
        # self.BrowsePorCC.setText('Browse')
        # self.BrowsePorCC.clicked.connect(self.BrowseButtonPorCC)
        #
        # # PARAMETERS PANEL
        # self.ParPanPorCC.setGeometry(10, 160, 320, 170)
        # self.ParPanPorCC.setFrameShape(QtWidgets.QFrame.StyledPanel)
        # # HydrophoneS section
        # self.HydLabel.setText('Number of Hydrophones')
        # self.HydLabel.setGeometry(10, 10, 180, 25)
        # self.EditHydN.setGeometry(210, 10, 50, 25)
        # self.EditHydN.setText("1")
        # # Sampling Frequency
        # self.SFreqLab.setText('Sampling Frequency')
        # self.SFreqLab.setGeometry(10, 40, 180, 25)
        # self.FsEdit.setText('576')
        # self.FsEdit.setGeometry(210, 40, 50, 25)
        # self.kHzLab.setText('kHz')
        # self.kHzLab.setGeometry(270, 40, 30, 25)
        # # Sensitivity
        # self.SensLab = QtWidgets.QLabel(self.ParPanPorCC)
        # self.SensLab.setText('Hyd Sensitivity')
        # self.SensLab.setGeometry(10, 70, 180, 25)
        # # self.SensEdit = QtWidgets.QLineEdit(self.ParPanPorCC)
        # self.SensEdit.setText('-172')
        # self.SensEdit.setGeometry(210, 70, 50, 25)
        # self.dBLab = QtWidgets.QLabel(self.ParPanPorCC)
        # self.dBLab.setText('dB')
        # self.dBLab.setGeometry(270, 70, 30, 25)
        # # Gain
        # self.GainLab.setText('Gain')
        # self.GainLab.setGeometry(10, 100, 180, 25)
        # self.GainEdit.setText('0')
        # self.GainEdit.setGeometry(210, 100, 50, 25)
        # self.dBLab2.setText('dB')
        # self.dBLab2.setGeometry(270, 100, 30, 25)
        # # DAQ Peak to Peak (Clipping level)
        # self.DAQLabPorCC.setText('DAQ pp Clipping level')
        # self.DAQLabPorCC.setGeometry(10, 130, 180, 25)
        # self.DAQEditPorCC.setText('2')
        # self.DAQEditPorCC.setGeometry(210, 130, 50, 25)
        # self.dBLab.setText('Volts')
        # self.dBLab.setGeometry(270, 130, 30, 25)
        # # Classifier panel
        # self.ProbPan.setFrameShape(QtWidgets.QFrame.StyledPanel)
        # self.ProbPan.setGeometry(10, 340, 320, 140)
        # # HQ
        # self.HQ.setText('Prob threshold')
        # self.HQ.setGeometry(10, 10, 100, 25)
        # self.HQThres.setText('0.999999')
        # self.HQThres.setGeometry(170, 10, 100, 25)
        # self.HQLab.setText('HQ')
        # self.HQLab.setGeometry(280, 10, 50, 25)
        # # LQ
        # self.LQ.setText('Prob threshold')
        # self.LQ.setGeometry(10, 50, 100, 25)
        # self.LQThres.setText('0.6')
        # self.LQThres.setGeometry(170, 50, 100, 25)
        # self.LQLab.setText('LQ')
        # self.LQLab.setGeometry(280, 50, 50, 25)
        # # Buttons
        # # OK button
        # self.OpenCTButton = QtWidgets.QPushButton(self.PorCCSetPan)
        # self.OpenCTButton.setText('OK')
        # self.OpenCTButton.setGeometry(40, 490, 100, 30)
        # self.OpenCTButton.clicked.connect(self.OKButtonPorCC)
        # # Cancel button
        # self.OpenCTCancelButton = QtWidgets.QPushButton(self.PorCCSetPan)
        # self.OpenCTCancelButton.setText('Cancel')
        # self.OpenCTCancelButton.setGeometry(210, 490, 100, 30)
        # self.OpenCTCancelButton.clicked.connect(self.CancelButtonPorCC)
        # # Disable resize and show menu
        # self.MenuPorCCF.setFixedSize(360, 550)
        # self.MenuPorCCF.show()

    def OKButtonPorCC(self):
        pass
        # DataType = self.FileDD.text()
        #
        # if DataType == 'SoundTrap output':
        #     name = 'SoundTrap'
        #     model = 1
        #     serial_number = 67416073
        #     Vpp = 2
        #     hydrophone = pyhy.SoundTrapHF(name, model, serial_number, Vpp)
        #     # clicks_df = hydrophone.read_HFclicks(sound_folder_path)
        #
        # # TODO make the OK PorCC function

    def CancelButtonPorCC(self):
        self.MenuPorCCF.close()

    def ApplyButtonPorCC(self):
        pass
        # TODO make the apply PorCC function

    def OpenDetSetMenu(self):
        # TODO finish the menu
        self.font.setBold(True)
        # CREATE MENU FIGURE
        self.centralwidgetdet.setObjectName("centralwidgetdet")
        self.tabWidget.setGeometry(QtCore.QRect(0, 0, 550, 490))
        self.tabWidget.setObjectName("tabWidget")
        self.ProjectTab.setObjectName("ProjectTab")
        self.tabWidget.addTab(self.ProjectTab, "General settings")
        self.HydTab.setObjectName("HydTab")
        self.tabWidget.addTab(self.HydTab, "Hydrophone")
        self.DetectorTab.setObjectName("DetectorTab")
        self.tabWidget.addTab(self.DetectorTab, "Detector")
        self.PorCCTab.setObjectName("PorCCTab")
        self.tabWidget.addTab(self.PorCCTab, "Clicks")
        self.MenuDetSetFig.setCentralWidget(self.centralwidgetdet)

        self.retranslateUi(self.MenuDetSetFig)
        QtCore.QMetaObject.connectSlotsByName(self.MenuDetSetFig)

        # Create panel within figure
        self.ProjectPan.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.ProjectPan.setGeometry(5, 5, 535, 460)
        self.SelectFolderText = QtWidgets.QLabel(self.ProjectPan)
        self.SelectFolderText.setGeometry(10, 10, 120, 30)
        self.SelectFolderText.setText('Select Folder')
        self.SelectFolderText.setFont(self.font)
        # Browse button
        self.FolderPathDet.setGeometry(10, 50, 400, 30)
        self.FolderPathDet.setText("C:/")
        # Checkbox
        self.CheckAllFolders.setGeometry(10, 85, 180, 30)
        self.CheckAllFolders.setText('Include subfolders')
        self.CheckAllFolders.setChecked(True)
        # Checkbox zipfile
        self.ZipMode.setGeometry(10, 115, 180, 30)
        self.ZipMode.setText('Zipped files')
        self.ZipMode.setChecked(False)
        # Browse button
        self.BrowseDet.setGeometry(300, 85, 110, 30)
        self.BrowseDet.setText("Browse")
        self.BrowseDet.clicked.connect(self.BrowseButtonDet)

        # saving
        self.SelectFolderTextSave.setGeometry(10, 140, 300, 30)
        self.SelectFolderTextSave.setText('(Files are saved in the same (sub)folder)')

        # length files
        self.LengthFileSave.setGeometry(10, 190, 220, 30)
        self.LengthFileSave.setText('Length limit - No of clicks')
        self.LengthFileEdit.setGeometry(300, 190, 110, 30)
        self.LengthFileEdit.setText('100000')
        self.LengthFileEdit.setAlignment(QtCore.Qt.AlignRight)
        self.LengthFileEdit.HorizontalAlignment = 'right'
        # Return to default values
        self.DefaultBut.setGeometry(20, 260, 150, 30)
        self.DefaultBut.setText("Default values")
        self.DefaultBut.clicked.connect(self.SetDefaults)

        # Detector
        self.DetPan.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.DetPan.setGeometry(5, 5, 535, 460)
        self.CheckClassify.setGeometry(300, 330, 180, 30)
        self.CheckClassify.setText('Classify with PorCC')
        self.CheckClassify.setChecked(True)

        # Trigger Filter
        self.DetectorLbl.setGeometry(20, 10, 150, 30)
        self.DetectorLbl.setText('Detector settings')
        self.DetectorLbl.setFont(self.font)
        self.ShortFiltLbl.setGeometry(10, 40, 120, 30)
        self.ShortFiltLbl.setText('Short filter')
        self.ShortFiltEdit.setGeometry(160, 40, 80, 30)
        self.ShortFiltEdit.setText('0.1')
        self.ShortFiltEdit.setAlignment(QtCore.Qt.AlignRight)
        self.LongFiltLbl.setGeometry(10, 80, 120, 30)
        self.LongFiltLbl.setText('Long filter')
        self.LongFiltEdit.setGeometry(160, 80, 80, 30)
        self.LongFiltEdit.setText('0.00001')
        self.LongFiltEdit.setAlignment(QtCore.Qt.AlignRight)
        self.LongFiltLbl2.setGeometry(10, 120, 120, 30)
        self.LongFiltLbl2.setText('Long filter 2')
        self.LongFiltEdit2.setGeometry(160, 120, 80, 30)
        self.LongFiltEdit2.setText('0.000001')
        self.LongFiltEdit2.setAlignment(QtCore.Qt.AlignRight)

        self.TriggerFilLbl.setGeometry(320, 10, 150, 30)
        self.TriggerFilLbl.setText('Trigger filter')
        self.TriggerFilLbl.setFont(self.font)
        self.MinFreqLbl.setGeometry(300, 50, 150, 30)
        self.MinFreqLbl.setText('Min Freq (kHz)')
        self.MinFreqEd.setGeometry(430, 50, 80, 30)
        self.MinFreqEd.setText('100')
        self.MinFreqEd.setAlignment(QtCore.Qt.AlignRight)
        self.MaxFreqLbl.setGeometry(300, 90, 150, 30)
        self.MaxFreqLbl.setText('Max Freq (kHz)')
        self.MaxFreqEd.setGeometry(430, 90, 80, 30)
        self.MaxFreqEd.setText('150')
        self.MaxFreqEd.setAlignment(QtCore.Qt.AlignRight)
        self.DetThresLbl.setGeometry(300, 130, 120, 30)
        self.DetThresLbl.setText('Threshold (dB)')
        self.DetThreshold.setGeometry(430, 130, 80, 30)
        self.DetThreshold.setText('10')
        self.DetThreshold.setAlignment(QtCore.Qt.AlignRight)

        # Prefilter
        self.PreFiltLbl.setGeometry(20, 160, 150, 30)
        self.PreFiltLbl.setText('Pre-filter')
        self.PreFiltLbl.setFont(self.font)
        self.PreFiltDDLbl.setGeometry(10, 200, 120, 30)
        self.PreFiltDDLbl.setText('Filter type')
        self.PreFiltDD.setGeometry(100, 200, 180, 30)
        self.PreFiltDD.addItem('Butterworth')
        self.PreFiltDD.addItem('Chebyshev I')
        self.PreFiltDD.addItem('Chebyshev II')
        self.PreFiltDD.addItem('Bessel')
        self.PreFiltPoleLbl.setGeometry(10, 240, 120, 30)
        self.PreFiltPoleLbl.setText('Filter order')
        self.PreFiltPole.setGeometry(160, 240, 80, 30)
        self.PreFiltPole.setText('4')
        self.PreFiltLbl2.setGeometry(10, 280, 120, 30)
        self.PreFiltLbl2.setText('Cut off freq')
        self.PreFiltFreq.setGeometry(160, 280, 80, 30)
        self.PreFiltFreq.setText('20000')

        # Clicks tab
        self.PorCCPantab.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.PorCCPantab.setGeometry(5, 5, 535, 460)
        self.ClicksLab.setGeometry(30, 10, 150, 30)
        self.ClicksLab.setText('Clicks (samples)')
        self.ClicksLab.setFont(self.font)
        self.PreSampLbl.setGeometry(10, 40, 150, 30)
        self.PreSampLbl.setText('Pre-Samples')
        self.PreSampEd.setGeometry(200, 40, 80, 30)
        self.PreSampEd.setText('40')
        self.PreSampEd.setAlignment(QtCore.Qt.AlignRight)
        self.PostSampLbl.setGeometry(10, 80, 150, 30)
        self.PostSampLbl.setText('Post-Samples')
        self.PostSampEd.setGeometry(200, 80, 80, 30)
        self.PostSampEd.setText('40')
        self.PostSampEd.setAlignment(QtCore.Qt.AlignRight)
        self.MaxLengthLbl.setGeometry(10, 120, 150, 30)
        self.MaxLengthLbl.setText('Max length')
        self.MaxLengthEd.setGeometry(200, 120, 80, 30)
        self.MaxLengthEd.setText('1024')
        self.MaxLengthEd.setAlignment(QtCore.Qt.AlignRight)
        self.MinLengthLbl.setGeometry(10, 160, 150, 30)
        self.MinLengthLbl.setText('Min Length')
        self.MinLengthEd.setGeometry(200, 160, 80, 30)
        self.MinLengthEd.setText('90')
        self.MinLengthEd.setAlignment(QtCore.Qt.AlignRight)
        self.MinSepLbl.setGeometry(10, 200, 150, 30)
        self.MinSepLbl.setText('Min separation')
        self.MinSepEd.setGeometry(200, 200, 80, 30)
        self.MinSepEd.setText('100')
        self.MinSepEd.setAlignment(QtCore.Qt.AlignRight)

        # label
        self.PorCCLbl.setGeometry(20, 240, 150, 30)
        self.PorCCLbl.setText('PorCC settings')
        self.PorCCLbl.setFont(self.font)
        # #LQ
        self.LQThresLab.setText('Prob threshold - LQ')
        self.LQThresLab.setGeometry(10, 280, 160, 30)
        self.LQThresDet.setText('0.6')
        self.LQThresDet.setGeometry(200, 280, 80, 30)
        self.LQThresDet.setAlignment(QtCore.Qt.AlignRight)
        # #HQ
        self.ProbHQLab.setText('Prob threshold - HQ')
        self.ProbHQLab.setGeometry(10, 320, 160, 20)
        self.HQThresDet.setText('0.999999')
        self.HQThresDet.setGeometry(200, 320, 80, 30)
        self.HQThresDet.setAlignment(QtCore.Qt.AlignRight)

        # Hydrophones section
        self.HydPan.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.HydPan.setGeometry(5, 5, 535, 460)
        self.HydParLbl.setGeometry(20, 10, 250, 30)
        self.HydParLbl.setText('Hydrophone settings')
        self.HydParLbl.setFont(self.font)
        self.DataTypeDD.setGeometry(200, 40, 150, 30)
        self.DataTypeDD.addItem("ST300HF")
        self.DataTypeDD.addItem("ST500HF")
        self.DataTypeDD.addItem('ST Click Detector')
        self.DataTypeDD.addItem("PAMGuard")
        #   self.FileDD.activated[str].connect(self.style_choice)
        self.DataTypeLbl.setGeometry(10, 40, 200, 30)
        self.DataTypeLbl.setText('Recorder/Hydrophone')

        self.SerialNoLab.setText('Serial number')
        self.SerialNoLab.setGeometry(10, 80, 150, 30)
        self.SerialNoEdit.setGeometry(200, 80, 150, 30)
        self.SerialNoEdit.setText('0')
        self.SerialNoEdit.setAlignment(QtCore.Qt.AlignRight)

        self.ChannelLab.setText('No of Channels')
        self.ChannelLab.setGeometry(10, 150, 150, 30)
        self.EditHydNDet.setGeometry(200, 150, 150, 30)
        self.EditHydNDet.setText('1')
        self.EditHydNDet.setAlignment(QtCore.Qt.AlignRight)
        # Sensitivity
        self.SenLabel.setText('Hyd Sensitivity (dB)')
        self.SenLabel.setGeometry(10, 190, 150, 30)
        self.SensEditDet.setText('-172')
        self.SensEditDet.setGeometry(200, 190, 150, 30)
        self.SensEditDet.setAlignment(QtCore.Qt.AlignRight)
        # #Gain
        self.GainLab1.setText('Gain (dB)')
        self.GainLab1.setGeometry(10, 230, 150, 30)
        self.GainEditDet.setText('0')
        self.GainEditDet.setGeometry(200, 230, 150, 30)
        self.GainEditDet.setAlignment(QtCore.Qt.AlignRight)
        self.DAQLab.setText('Clip level (p-p, Volts)')
        self.DAQLab.setGeometry(10, 270, 150, 30)
        self.DAQppEdit.setText('2')
        self.DAQppEdit.setGeometry(200, 270, 150, 30)
        self.DAQppEdit.setAlignment(QtCore.Qt.AlignRight)

        #         ## if soundtrap
        # TODO the soundtrap needs serial number
        # read_HF()

        # Check button to save everything or just LQ and HQ
        # self.LQHQSaveLbl = QtWidgets.QLabel(self.PorCCPantab)
        # self.LQHQSaveLbl.setGeometry(340, 220, 150, 30)
        # self.LQHQSaveLbl.setText('Save:')
        # self.LQHQSaveLbl.setFont(self.font)
        # self.LQHQCheck = QtWidgets.QCheckBox(self.PorCCPantab)
        # self.LQHQCheck.setText('LQ, HQ')
        # self.LQHQCheck.setGeometry(340, 260, 160, 30)
        # self.LQHQCheck.setChecked(True)
        # self.LQHQNoise = QtWidgets.QCheckBox(self.PorCCPantab)
        # self.LQHQNoise.setText('LQ, HQ, Noise')
        # self.LQHQNoise.setGeometry(340, 300, 160, 30)

        # OK and Cancel buttons
        # OK button in Project
        self.OKButtonDet = QtWidgets.QPushButton(self.ProjectPan)
        self.OKButtonDet.setText('Run detector')
        self.OKButtonDet.setGeometry(60, 380, 150, 40)
        self.CancelButtonDet = QtWidgets.QPushButton(self.ProjectPan)
        self.CancelButtonDet.setText('Cancel')
        self.CancelButtonDet.setGeometry(320, 380, 150, 40)
        self.OKButtonDet.clicked.connect(self.RunDetector)
        self.CancelButtonDet.clicked.connect(self.CancelDetector)
        # OK button in Hydrophone
        self.OKButtonDet = QtWidgets.QPushButton(self.HydPan)
        self.OKButtonDet.setText('Run detector')
        self.OKButtonDet.setGeometry(60, 380, 150, 40)
        self.CancelButtonDet = QtWidgets.QPushButton(self.HydPan)
        self.CancelButtonDet.setText('Cancel')
        self.CancelButtonDet.setGeometry(320, 380, 150, 40)
        self.OKButtonDet.clicked.connect(self.RunDetector)
        self.CancelButtonDet.clicked.connect(self.CancelDetector)
        # in Detector
        self.OKButtonDet = QtWidgets.QPushButton(self.DetPan)
        self.OKButtonDet.setText('Run detector')
        self.OKButtonDet.setGeometry(60, 380, 150, 40)
        self.CancelButtonDet = QtWidgets.QPushButton(self.DetPan)
        self.CancelButtonDet.setText('Cancel')
        self.CancelButtonDet.setGeometry(320, 380, 150, 40)
        self.OKButtonDet.clicked.connect(self.RunDetector)
        self.CancelButtonDet.clicked.connect(self.CancelDetector)
        # in PorCC
        self.OKButtonDet = QtWidgets.QPushButton(self.PorCCPantab)
        self.CancelButtonDet = QtWidgets.QPushButton(self.PorCCPantab)
        self.OKButtonDet.setText('Run detector')
        self.OKButtonDet.setGeometry(60, 380, 150, 40)
        self.CancelButtonDet.setText('Cancel')
        self.CancelButtonDet.setGeometry(320, 380, 150, 40)
        self.OKButtonDet.clicked.connect(self.RunDetector)
        self.CancelButtonDet.clicked.connect(self.CancelDetector)

        self.MenuDetSetFig.setFixedSize(550, 490)
        self.MenuDetSetFig.show()

    def SetDefaults(self):
        self.LengthFileEdit.setText('1000000')
        self.EditHydN.setText('1')
        self.SensEditDet.setText('-172')
        self.GainEditDet.setText('0')
        self.DAQppEdit.setText('2')
        self.ShortFiltEdit.setText('0.1')
        self.LongFiltEdit.setText('0.00001')
        self.LongFiltEdit2.setText('0.000001')
        self.PreFiltPole.setText('4')
        self.PreFiltFreq.setText('20000')
        self.MinFreqEd.setText('100')
        self.MaxFreqEd.setText('150')
        self.DetThreshold.setText('10')
        self.CheckClassify.setChecked(True)
        self.PreSampEd.setText('40')
        self.PostSampEd.setText('40')
        self.MaxLengthEd.setText('1024')
        self.MinLengthEd.setText('90')
        self.MinSepEd.setText('100')
        self.HQThresDet.setText('0.999999')
        self.LQThresDet.setText('0.6')

    def RunDetector(self):
        self.MenuDetSetFig.close()
        MainFolder = self.FolderPathDet.text()
        ModelSel = self.DataTypeDD.currentText()
        classifier = porcc.PorCC(load_type='manual', config_file='default')
        LQ = float(self.LQThresDet.text())
        HQ = float(self.HQThresDet.text())
        # update the thresholds
        classifier.th1 = HQ
        classifier.th2 = LQ
        if ModelSel == 'ST300HF' or ModelSel == 'ST500HF':
            name = 'SoundTrap'
            serial_number = int(self.SerialNoEdit.text())  # 0
            if serial_number == 0:
                Sens = int(self.SensEditDet.text())
                hydrop = pyhy.soundtrap.SoundTrap(name=name, model=ModelSel, sensitivity=Sens,
                                                  serial_number=serial_number)
            else:
                hydrop = pyhy.soundtrap.SoundTrap(name=name, model=ModelSel, serial_number=serial_number)

            MaxLenFile = int(self.LengthFileEdit.text())
            LongFilt = float(self.LongFiltEdit.text())
            LongFilt2 = float(self.LongFiltEdit2.text())
            ShortFilt = float(self.ShortFiltEdit.text())
            DetThres = float(self.DetThreshold.text())
            PreSam = int(self.PreSampEd.text())
            PostSam = int(self.PostSampEd.text())
            MaxLenClick = int(self.MaxLengthEd.text())
            MinLenClick = int(self.MinLengthEd.text())
            MinSep = int(self.MinSepEd.text())
            MinFrq = float(self.MinFreqEd.text()) * 1000
            MaxFrq = float(self.MaxFreqEd.text()) * 1000

            pfilter = click_detector.Filter(filter_name='butter', filter_type='bandpass', order=4,
                                            frequencies=[MinFrq, MaxFrq])
            dfilter = click_detector.Filter(filter_name='butter', filter_type='high', order=4, frequencies=20000)
            detector = click_detector.ClickDetector(hydrophone=hydrop, long_filt=LongFilt, long_filt2=LongFilt2,
                                                    short_filt=ShortFilt, threshold=DetThres, min_separation=MinSep,
                                                    max_length=MaxLenClick, min_length=MinLenClick, pre_samples=PreSam,
                                                    post_samples=PostSam, prefilter=pfilter, dfilter=dfilter,
                                                    save_max=MaxLenFile, save_folder=MainFolder, convert=True,
                                                    click_model_path=None, classifier=classifier, save_noise=False)
            blocksize = 3456000
        elif ModelSel == 'ST Click Detector':
            name = 'SoundTrap'
            serial_number = int(self.SerialNoEdit.text())  # 738496579
            if serial_number == 0:
                Sens = int(self.SensEditDet.text())
                hydrop = pyhy.soundtrap.SoundTrapHF(name=name, model=ModelSel, sensitivity=Sens,
                                                    serial_number=serial_number)
            else:
                hydrop = pyhy.soundtrap.SoundTrapHF(name=name, model=ModelSel, serial_number=serial_number)
            MinFrq = float(self.MinFreqEd.text()) * 1000
            MaxFrq = float(self.MaxFreqEd.text()) * 1000

            pfilter = click_detector.Filter(filter_name='butter', filter_type='bandpass', order=4,
                                            frequencies=[MinFrq, MaxFrq])
            detector = click_detector.ClickDetectorSoundTrapHF(hydrophone=hydrop, prefilter=pfilter,
                                                               save_folder=MainFolder, convert=True,
                                                               click_model_path=None, classifier=classifier,
                                                               save_noise=False)
            blocksize = None

        else:
            blocksize = 3456000
            detector = None

        # for loop to go into subfolders
        zip_mode = self.ZipMode.isChecked()
        if self.CheckAllFolders.isChecked():
            FilesAndFolders = os.listdir(MainFolder)
            if zip_mode:
                Folders = [s for s in FilesAndFolders if 'zip' in s]
                detector.save_folder = MainFolder
            else:
                Folders = [s for s in FilesAndFolders if os.path.isdir(os.path.join(MainFolder, s))]
            for myFolder in Folders:
                ThisFolderSave = MainFolder + '/' + myFolder
                print('Detecting clicks in', ThisFolderSave)
                if not zip_mode:
                    detector.save_folder = ThisFolderSave
                detector.detect_click_clips_folder(ThisFolderSave, blocksize=blocksize, zip_mode=zip_mode)
        else:
            detector.detect_click_clips_folder(MainFolder, blocksize=blocksize, zip_mode=zip_mode)

    def CancelDetector(self):
        self.MenuDetSetFig.close()

    def OpenCTMenu(self):
        self.OpenCTFig.setGeometry(350, 350, 380, 210)
        self.OpenCTFig.setWindowTitle("Open Project")
        # create a panel inside the figure
        self.OpenCTPan.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.OpenCTPan.setGeometry(QtCore.QRect(10, 10, 360, 190))
        # label "Select folder"
        self.SelectFolderText1.setGeometry(20, 10, 200, 30)
        self.SelectFolderText1.setText('Select Folder')
        # edit - selected folder to be read
        self.FolderPathOpenCT = QtWidgets.QLineEdit(self.OpenCTPan)
        self.FolderPathOpenCT.setGeometry(20, 50, 310, 30)
        # browse button
        self.BrowseButtonOpen = QtWidgets.QPushButton(self.OpenCTPan)
        self.BrowseButtonOpen.setGeometry(230, 85, 100, 30)
        self.BrowseButtonOpen.setText('Browse')
        self.BrowseButtonOpen.clicked.connect(self.PushBrowseButtonOpenCT)
        # OK button
        self.OpenCTButton.setText('OK')
        self.OpenCTButton.setGeometry(40, 140, 100, 30)
        self.OpenCTButton.clicked.connect(self.OpenCT)
        # Cancel button
        self.OpenCTCancelButton.setText('Cancel')
        self.OpenCTCancelButton.setGeometry(210, 140, 100, 30)
        self.OpenCTCancelButton.clicked.connect(self.OpenCTCancel)
        # Disable resize and show menu
        self.OpenCTFig.setFixedSize(380, 210)
        self.OpenCTFig.show()

    def BrowseButtonDet(self):
        root = tk.Tk()
        root.withdraw()
        self.SelectedFolder = filedialog.askdirectory()
        self.FolderPathDet.setText(self.SelectedFolder)

    def BrowseButtonPorCC(self):
        root = tk.Tk()
        root.withdraw()
        SelectedFolderPorCC = filedialog.askdirectory()
        self.FolderPathPorCC.setText(SelectedFolderPorCC)

    def NewCTMenu(self):
        self.NewCTFig.setGeometry(100, 100, 360, 320)
        self.NewCTFig.setWindowTitle('New Project')
        self.NewCTFig.sizePolicy = 'Fixed'
        self.NewCTFig.setFixedSize(360, 320)
        # create a panel inside the figure
        self.NewCTPan.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.NewCTPan.setGeometry(10, 10, 340, 300)
        self.SelectFolderText2.setGeometry(20, 10, 100, 30)
        self.SelectFolderText2.setText('Select Folder')
        self.FolderPathNewCT.setGeometry(10, 40, 300, 30)
        # self.FolderPathNewCT.setText('C:\Users')
        self.InclSubFoldersNewCT.setGeometry(10, 80, 180, 20)
        self.InclSubFoldersNewCT.setText("Include subfolders")
        # browse button
        self.NewCTBrowseButton.setText("Browse")
        self.NewCTBrowseButton.setGeometry(QtCore.QRect(210, 90, 100, 30))
        self.NewCTBrowseButton.clicked.connect(self.NewCTBrowse)
        # LATITUDE AND LONGITUDE PANEL
        self.LocationPan.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.LocationPan.setGeometry(10, 140, 320, 100)
        # Latitude
        self.LatLabel.setGeometry(10, 10, 80, 30)
        self.LatLabel.setText("Latitude")
        self.LatEdit.setGeometry(200, 10, 100, 30)
        self.LatEdit.setText("55.576197")
        # Longitude
        self.LongLable.setText('Longitude')
        self.LongLable.setGeometry(10, 50, 100, 30)
        self.LongEdit.setText('9.848555')
        self.LongEdit.setGeometry(200, 50, 100, 30)
        # Buttons
        self.OKButtonCT.setGeometry(QtCore.QRect(10, 250, 100, 30))
        self.OKButtonCT.setText("OK")
        self.OKButtonCT.clicked.connect(self.IdentifyCT)
        self.CancelButtonCT.setGeometry(QtCore.QRect(230, 250, 100, 30))
        self.CancelButtonCT.setText("Cancel")
        self.CancelButtonCT.clicked.connect(self.PushCancelButtonNewCT)

        self.NewCTFig.show()

    def PushCancelButtonNewCT(self):
        self.NewCTFig.close()
        # TODO make only one cancel function

    def NewCTBrowse(self):
        root = tk.Tk()
        root.withdraw()
        SelectedFolderNewCT = filedialog.askdirectory()
        self.FolderPathNewCT.setText(SelectedFolderNewCT)

    def IdentifyCT(self):
        global CP  # sset, srise
        self.NewCTFig.close()
        MainFolder = self.FolderPathNewCT.text()
        latitude = float(self.LatEdit.text())
        longitude = float(self.LongEdit.text())
        if self.InclSubFoldersNewCT.isChecked():
            FilesAndFolders = os.listdir(MainFolder)
            myFolders = [s for s in FilesAndFolders if "." not in s]
        else:
            myFolders = MainFolder
            print(myFolders)

        for myFolder in myFolders:
            if myFolders == MainFolder:
                FilesInFolder = os.listdir(MainFolder)
            else:
                FilesInFolder = os.listdir(MainFolder + '/' + myFolder)
            if any("." in s for s in FilesInFolder):
                CPFile = [s for s in FilesInFolder if "CP.csv" in s]
                print(len(CPFile))
                if len(CPFile) > 0:
                    if myFolders == MainFolder:
                        CPFileName = MainFolder + '/CP.csv'
                    else:
                        CPFileName = MainFolder + '/' + thisFolder + '/CP.csv'
                    CP = pd.read_csv(CPFileName)
                else:
                    DetClicks = [s for s in FilesInFolder if "clips.csv" in s]
                    CP = pd.DataFrame()
                    for thisFile in DetClicks:
                        if myFolders == MainFolder:
                            CPFile = MainFolder + '/' + thisFile
                        else:
                            CPFile = MainFolder + '/' + thisFolder + '/' + thisFile
                        TempCP = pd.read_csv(CPFile, index_col=0)
                        TempCP['id'] = TempCP.index.values
                        CP = CP.append(TempCP, ignore_index=True)
                    if myFolders == MainFolder:
                        CPFileName = MainFolder + '/CP.csv'
                    else:
                        CPFileName = MainFolder + '/' + thisFolder + '/CP.csv'
                    CP.to_csv(CPFileName, index=False)

            fs = 576000  # I will need a settings file
            CTrains, thisCTInfo, CP = ExtractPatterns(CP, fs, latitude, longitude)
            if myFolders == MainFolder:
                CTInfoFileName = MainFolder + '/CTInfo.csv'
                CTrainsFileName = MainFolder + '/CTrains.csv'
                thisCTInfo.to_csv(CTInfoFileName, index=False)
                CTrains.to_csv(CTrainsFileName, index=False)
                break
            else:
                CTInfoFileName = MainFolder + '/' + thisFolder + '/CTInfo.csv'
                CTrainsFileName = MainFolder + '/' + thisFolder + '/CTrains.csv'
                thisCTInfo.to_csv(CTInfoFileName, index=False)
                CTrains.to_csv(CTrainsFileName, index=False)

    def OpenCTCancel(self):
        self.OpenCTFig.close()

    # OPEN CT
    def PushBrowseButtonOpenCT(self):
        root = tk.Tk()
        root.withdraw()
        self.SelectedFolderCT = filedialog.askdirectory()
        self.FolderPathOpenCT.setText(self.SelectedFolderCT)

    def OpenCT(self):
        global CP, CTInfo
        self.OpenCTFig.close()
        CPFileName = self.SelectedFolder + "/CTrains.csv"
        CP = pd.read_csv(CPFileName)
        CTInfoFileName = self.SelectedFolder + "/CTInfo.csv"
        CTInfo = pd.read_csv(CTInfoFileName)
        CTInfo = CTInfo[CTInfo.Length > 9]
        CTInfo.reset_index(inplace=True, drop=True)
        NumCT = int(CTInfo.CTNum.iloc[0])
        self.UpdateCT(NumCT, CP, CTInfo)

    def UploadMetricData(self):
        global topLevelFolderMetrics, CTInfo, CP, thisFolder
        topLevelFolderMetrics = self.SelectFolderMetricEdit.text()
        FilesAndFolders = os.listdir(topLevelFolderMetrics)
        PPMFile = [s for s in FilesAndFolders if "PosPorMin" in s]
        SummaryFile = [s for s in FilesAndFolders if "SummaryTable" in s]
        if PPMFile and SummaryFile:
            FileToOpenPPM = topLevelFolderMetrics + '/PosPorMin.csv'
            PPM = pd.read_csv(FileToOpenPPM)
            FileToOpenCT = topLevelFolderMetrics + '/SummaryTable.csv'
            SumTable = pd.read_csv(FileToOpenCT)
        else:  # If the data needs to be prepared
            if any("." not in s for s in FilesAndFolders):
                listOfFolderNamesMetrics = [s for s in FilesAndFolders if "." not in s]
                numberOfFoldersMetrics = len(listOfFolderNamesMetrics)
                ColName = ['Minute']
                Minutes = pd.DataFrame(0, index=np.arange(1440), columns=ColName)
                Minutes.iloc[0, 0] = '00:00'
                Minutes.iloc[1:1440, 0] = range(1, 1440)
                for i in range(1, len(Minutes)):
                    if int(i) < 60:
                        a = time(hour=0, minute=int(i))
                        Minutes.iloc[i, 0] = a.strftime('%H:%M')
                    else:
                        h = int(int(i) / 60)
                        m = int(i) % 60
                        a = time(hour=int(h), minute=int(m))
                        Minutes.iloc[i, 0] = a.strftime('%H:%M')
                # Put in Clicks Per Minute
                #    self.TimesRec = pd.DataFrame()
                numberOfFoldersMetrics = len(listOfFolderNamesMetrics)
                PPM = pd.DataFrame(0, index=np.arange(1440), columns=['Minute', listOfFolderNamesMetrics])
                for FolderN in range(0, numberOfFoldersMetrics - 1):
                    # load file
                    thisFolder = listOfFolderNamesMetrics[FolderN]
                    CTInfoFile = thisFolder + '/CTInfo.csv'
                    CTInfo = pd.read_csv(CTInfoFile)
                    # Load ClickParameters
                    ClickParFile = thisFolder + '/CPar.csv'
                    CP = pd.read_csv(ClickParFile)
                    CTInfo.CT = int(CTInfo.CT)
                    # check the positive min
                    AllNBHF = CTInfo[CTInfo.Species != 'Non-NBHF']
                    for i in range(0, len(AllNBHF)):
                        Time = AllNBHF.Date.iloc[i]
                        Time = str(Time)
                        Hour = int(Time[11:13]) * 60
                        Min = int(Time[14:16])
                        RowN = Hour + Min + 1
                        PPM.Positive[RowN, FolderN - 1] = 1
                        PPM.NumCT[RowN, FolderN - 1] = PPM.NumCT[RowN, FolderN - 1] + 1
                        PPM.NumClicks[RowN, FolderN - 1] = PPM.NumClicks[RowN, FolderN - 1] + AllNBHF.Length[i]
                    FullName = thisFolder + 'PosPorMin.csv'
                    PPM.to_csv(FullName)

    def SaveUpdates(self):
        FullName = self.SelectedFolder + '/CTInfo.csv'
        CTInfo.to_csv(FullName)


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("fusion")
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
