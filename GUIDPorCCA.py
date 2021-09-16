"""
Module : GUIDPorCCA.py
Authors : Mel Cosentino
Institution : University of Strathclyde (UK), Aarhus University (Denmark)
Last Accessed : 27/03/2021
"""

__author__ = "Mel Cosentino"
__version__ = "0"
__credits__ = "Mel Cosentino"
__email__ = "orcinus.orca.1758@gmail.com"
__status__ = "Development"

import os
import tkinter as tk
import warnings
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
import zipfile
import pathlib

import click_trains

pd.options.mode.chained_assignment = None

global BrowseSelectedFolder, CTInfo, CTTemp, NHyd, CP, topLevelFolderMetrics
global thisFolder, sset, srise, Fs, PosPorMin, SummaryTable, numberOfFoldersMetrics


def UpdateWaterfall(XLimMin, YLimMin, ZLimMin, XLimMax, YLimMax, ZLimMax):
    pass


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
        self.CheckAllFoldersMetr = QtWidgets.QCheckBox(self.MetricsUploadPan)

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
        self.open_ct_cancel_b = QtWidgets.QPushButton(self.OpenCTPan)
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

        # MAIN TABS
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

        # DISPLAY PARAMETERS
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
        self.SaveupdatesButton.clicked.connect(self.save_updates)

        # AXES AREA
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

        # ACTION AREA
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
        self.CTTypeDropDown.addItem("Noise")
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

        # METRICS DISPLAY
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
        # CheckBox
        self.CheckAllFoldersMetr.setGeometry(20, 74, 180, 30)
        self.CheckAllFoldersMetr.setText('Include subfolders')
        self.CheckAllFoldersMetr.setChecked(True)

        self.MetricBrowseButton.clicked.connect(self.MetricsBrowse)

        # Metrics Table
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

        # Metrics summary
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

        # Positive porpoise minutes
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
        self.DisplayMetricsButton.clicked.connect(self.display_metrics)
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

        # MAIN MENU
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

    # FUNCTIONS ###
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

    def metrics(self):
        pass

    def display_metrics(self):
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

    def update_ct(self, NumCT, Cp, myCTInfo):
        """
            Displays the click train number 'NumCT'

            :param:
                NumCT: click train number
                Cp:  pandas dataframe containing information about all clicks detected in the recordings
                myCTInfo: contains summary data of all click trains
        """
        global CTTemp
        CTTemp = Cp[Cp.CT == NumCT]
        if NumCT != 1:
            CTTemp.reset_index(inplace=True, drop=True)
        if len(CTTemp) > 9:
            CTTemp = click_trains.new_ici(CTTemp)
            CTTemp.loc[:, 'SumMs'] = int(0)
            for i in range(1, len(CTTemp)):
                CTTemp.SumMs[i] = int(CTTemp.SumMs[i - 1]) + int(CTTemp.ICI[i])
            CTTemp.SumMs = CTTemp.SumMs / 1000
            CT1HQ = CTTemp[CTTemp['pyPorCC'] == 1]
            CT1LQ = CTTemp[CTTemp['pyPorCC'] == 2]
            self.CTNumD.setText(str(NumCT))
            self.CTTypeLabel.setText(str(myCTInfo.CTType[myCTInfo.CTNum == NumCT].values[0]))
            self.DateandtimeofCTLabel.setText(str(myCTInfo.Date[myCTInfo.CTNum == NumCT].values[0]))
            self.LengthLabel.setText(str(len(CTTemp)))
            # self.BehavLabel.setText(str(myCTInfo.Behav[myCTInfo.CTNum == NumCT].values[0]))
            # self.CalfLabel.setText(str(myCTInfo.Calf[myCTInfo.CTNum == NumCT].values[0]))
            self.BehavLabel.setText('-')
            self.CalfLabel.setText('-')
            self.TotalLabel.setText('(' + str(myCTInfo['CTNum'].iloc[-1]) + ')')
            self.DayLabel.setText(str(myCTInfo.DayNight[myCTInfo.CTNum == NumCT].values[0]))
            # if CTInfo["Saved"][[CTInfo.CTNum] == NumCT] == 1:
            #    SelectedCTCheckBox.Value = 1
            self.AmpAxesCT.clear()

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

    # Display buttons
    def CTBackCB(self):
        """
            Moves to the previous click train
        """
        NumCT = int(self.CTNumD.text())
        First = CTInfo['CTNum'].iloc[0]
        if NumCT > First:
            self.AmpAxesCT.clear()
            self.ICIAxesCT.clear()
            self.FreqAxesCT.clear()
            RowCT = CTInfo[CTInfo.CTNum == NumCT].index[0]
            NumCT = int(CTInfo.CTNum[RowCT - 1])
            self.update_ct(NumCT, CP, CTInfo)

    def CTForwCB(self):
        """
            Moves to the next click train
        """
        NumCT = int(self.CTNumD.text())
        Tot = CTInfo['CTNum'].iloc[-1]
        if NumCT == Tot:
            a = 1  # do nothing
        elif NumCT < Tot:
            self.AmpAxesCT.clear()
            self.ICIAxesCT.clear()
            self.FreqAxesCT.clear()
            RowCT = CTInfo[CTInfo.CTNum == NumCT].index[0]
            NumCT = int(CTInfo.CTNum[RowCT + 1])
            self.update_ct(NumCT, CP, CTInfo)

    def NotesCT(self):
        """
            Pops up a window where the user can write comments about the click train, which is then saved into the
            CTInfo dataframe, which contains summary information of all click trains.
        """
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
        pass
        # global CTTemp
        # x1 = []
        # y1 = []
        # z1 = []
        # self.Fig3D.setGeometry(50, 50, 1200, 800)
        # self.Fig3D.setWindowTitle('Click train in 3D')
        # self.Pan3D.setFrameShape(QtWidgets.QFrame.StyledPanel)
        # self.Pan3D.setGeometry(10, 10, 1180, 780)
        # self.Plot3D.setGeometry(600, 300, 570, 460)
        #
        # ## WATERFALL PLOT
        # FFT = 512
        # CTTemp['iciSum'] = 0
        # for i in range(1, len(CTTemp)):
        #     CTTemp.iciSum[i] = CTTemp.ICI[i] + CTTemp.iciSum[i - 1]
        # WavFileToOpen = CTTemp.filename[0]
        #
        # # CREATE empty x, y, and z for the waterfall plot
        # z1 = np.zeros((257, len(CTTemp)), dtype=float)  # power of each click in each row
        #
        # # FILL the variables
        # # X = frequency
        # click1, fs = soundfile.read(WavFileToOpen, start=1, stop=150)
        # freqs, psd = signal.welch(click1, fs=fs, window='hann', nfft=512)
        # x1 = freqs / 1000
        # # x1.to_numpy()
        # # Y = time( in secs)
        # y1 = CTTemp.iciSum.to_numpy()
        # y1 = y1 / 1000
        # # y1 = y1.T
        #
        # # Normalised Amplitude
        # for i in range(0, len(CTTemp)):
        #     Start = CTTemp.start_sample[i]
        #     End = Start + CTTemp.duration_samples[i]
        #     click, Fs = soundfile.read(WavFileToOpen, start=int(Start), stop=int(End))
        #     freqs, psd = signal.welch(click, fs=Fs, window='hann', nfft=512)
        #     z1[:, i] = psd
        # a = z1.max()
        # z1 = z1 / a
        # # ax = plt.axes(projection='3d')
        #
        # # Data for a three-dimensional line
        # # ax.contour3D(x1, y1, z1)
        # # ax.set_xlabel('Frequency (kHz)')
        # # ax.set_ylabel('Time (ms)')
        # # ax.set_zlabel('Amplitude')
        #
        # # Data for three-dimensional scattered points
        # # zdata = 15 * np.random.random(100)
        # # xdata = np.sin(zdata) + 0.1 * np.random.randn(100)
        # # ydata = np.cos(zdata) + 0.1 * np.random.randn(100)
        # # ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens')
        #
        # # self.Plot3D.add_subplot(111, projection='3d')
        #
        # # n = 0
        # print(len(y1), len(x1), len(z1))
        # p2 = gl.GLSurfacePlotItem(y=y1, x=x1, z=z1, shader='shaded', color=(0.5, 0.5, 1, 1))
        # p2.translate(-10, -30, 10)
        # p2.scale(1.0, 1.0, 0.5)
        # self.Plot3D.addItem(p2)
        #
        # # self.Plot3D.contour3D(x1, y1, z1)
        # # self.Plot3D.set_xlabel('Frequency (kHz)')
        # # self.Plot3D.set_ylabel('Time (ms)')
        # # self.Plot3D.set_zlabel('Amplitude')
        # # plt.scatter(self.Plot3D, x1, y1, z1)
        # self.Fig3D.show()
        # self.Plot3D.show()

    def CreateSpectrogram(self):
        """
            Creates a pop up figure and plots the waveform and spectrogram of the click train
        """
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
        if '.zip' in str(FileToOpen):
            if isinstance(FileToOpen, str):
                FileToOpen = pathlib.Path(FileToOpen)
            zip_file = zipfile.ZipFile(FileToOpen.parent)
            FileToOpen = zip_file.open(FileToOpen.name)
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
        """
            Updates the spectrogram based on a series of parameters provided by the user.
        """
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
        """
            Allows the user to change the classification of the click train manually.
        """
        global CTInfo
        CTNum = int(self.CTNumD.text())
        CTType = self.CTTypeDropDown.currentText()
        Behaviour = self.BehaviourDropDown.currentText()
        # if CTType == 'Select' and Behaviour == 'Select':
        #     a = 1
        # el
        if not CTType == 'Select' and Behaviour == 'Select':
            CTInfo.CTType[CTInfo.CTNum == CTNum] = CTType
            self.CTTypeLabel.setText(str(CTType))
        elif CTType == 'Select' and not Behaviour == 'Select':
            CTInfo.Behav[CTInfo.CTNum == CTNum] = Behaviour
            self.BehavLabel.setText(str(Behaviour))
        elif not CTType == 'Select' and not Behaviour == 'Select':
            CTInfo.CTType[CTInfo.CTNum == CTNum] = CTType
            CTInfo.Behav[CTInfo.CTNum == CTNum] = Behaviour
            self.BehavLabel.setText(str(Behaviour))
            self.CTTypeLabel.setText(str(CTType))

    # METRICS STUFF

    def MetricsBrowse(self):
        global BrowseSelectedFolder
        self.root_metric_browse_b = tk.Tk()
        self.root_metric_browse_b.withdraw()
        BrowseSelectedFolder = filedialog.askdirectory()
        self.SelectFolderMetricEdit.setText(BrowseSelectedFolder)
        self.root_metric_browse_b.mainloop(0)

    def UploadMetricData(self):
        """
            Generates summary data of the project and uploads it for visualisation purposes

            PosPorMin: positive porpoise minute - contains at least one NBHF click train (strict criteria) or a NBHF
                or a LQ-NBHF click train (relaxed criteria).
        """
        global topLevelFolderMetrics, CTInfo, CP, thisFolder
        topLevelFolderMetrics = self.SelectFolderMetricEdit.text()
        # Check whether the PosPorMin and SummaryTable files have been created
        FilesAndFolders = os.listdir(topLevelFolderMetrics)
        PPMFile = [s for s in FilesAndFolders if "PosPorMin.csv" in s]
        SummaryFile = [s for s in FilesAndFolders if "SummaryTable.csv" in s]
        # Open if ready
        if PPMFile and SummaryFile:
            FileToOpenPPM = os.path.join(topLevelFolderMetrics, 'PosPorMin.csv')
            PPM = pd.read_csv(FileToOpenPPM)
            FileToOpenCT = os.path.join(topLevelFolderMetrics, 'SummaryTable.csv')
            SumTable = pd.read_csv(FileToOpenCT)
        else:  # If the data needs to be prepared
            PPM = pd.DataFrame()
            if self.CheckAllFoldersMetr.isChecked():
                # List subfolders
                Folders = [item for item in os.listdir(topLevelFolderMetrics) if
                           os.path.isdir(os.path.join(topLevelFolderMetrics, item))]
                numberOfFoldersMetrics = len(Folders)
                for myFolder in Folders:
                    Files = os.listdir(os.path.join(topLevelFolderMetrics, myFolder))
                    CTInfoFile = [s for s in Files if "CTInfo.csv" in s]
                    if len(CTInfoFile) > 0:
                        CTInfo = pd.read_csv(os.path.join(topLevelFolderMetrics, myFolder, "CTInfo.csv"))
                    # Clicks = [s for s in Files if "Clicks.csv" in s]
                    Date = CTInfo.Date[0]
                    NewDate = int(Date[0:4] + Date[5:7] + Date[8:10])

                    # Metrics.Date[d-2] = strftime(datetime(yyyy, mmm, dd))

                # if any("." not in s for s in FilesAndFolders):
                #     listOfFolderNamesMetrics = [s for s in FilesAndFolders if "." not in s]
                #     numberOfFoldersMetrics = len(listOfFolderNamesMetrics)
                #     ColName = ['Minute']
                #     Minutes = pd.DataFrame(0, index=np.arange(1440), columns=ColName)
                #     Minutes.iloc[0, 0] = '00:00'
                #     Minutes.iloc[1:1440, 0] = range(1, 1440)
                #     for i in range(1, len(Minutes)):
                #         if int(i) < 60:
                #             a = time(hour=0, minute=int(i))
                #             Minutes.iloc[i, 0] = a.strftime('%H:%M')
                #         else:
                #             h = int(int(i) / 60)
                #             m = int(i) % 60
                #             a = time(hour=int(h), minute=int(m))
                #             Minutes.iloc[i, 0] = a.strftime('%H:%M')
                #     # Put in Clicks Per Minute
                #     #    self.TimesRec = pd.DataFrame()
                #     numberOfFoldersMetrics = len(listOfFolderNamesMetrics)
                #     PPM = pd.DataFrame(0, index=np.arange(1440), columns=['Minute', listOfFolderNamesMetrics])
                #     for FolderN in range(0, numberOfFoldersMetrics - 1):
                #         # load file
                #         thisFolder = listOfFolderNamesMetrics[FolderN]
                #         CTInfoFile = thisFolder + '/CTInfo.csv'
                #         CTInfo = pd.read_csv(CTInfoFile)
                #         # Load ClickParameters
                #         ClickParFile = thisFolder + '/CPar.csv'
                #         CP = pd.read_csv(ClickParFile)
                #         CTInfo.CT = int(CTInfo.CT)
                #         # check the positive min
                #         AllNBHF = CTInfo[CTInfo.CTType != 'Noise']
                #         for i in range(0, len(AllNBHF)):
                #             Time = AllNBHF.Date.iloc[i]
                #             Time = str(Time)
                #             Hour = int(Time[11:13]) * 60
                #             Min = int(Time[14:16])
                #             RowN = Hour + Min + 1
                #             PPM.Positive[RowN, FolderN - 1] = 1
                #             PPM.NumCT[RowN, FolderN - 1] = PPM.NumCT[RowN, FolderN - 1] + 1
                #             PPM.NumClicks[RowN, FolderN - 1] = PPM.NumClicks[RowN, FolderN - 1] + AllNBHF.Length[i]
                #         FullName = thisFolder + 'PosPorMin.csv'
                #         PPM.to_csv(FullName)

    """
    MENUS
    """

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
        """
            Detects and saves individual impulsive sounds from audio files (wav). It is an adaptation of the click
                detector module in PAMGuard.
            We recommend using the default values.

            Additional parameters and functions
                PorCC: Porpoise Click Classifier (Cosentino et al, 2019). Classifies the detected sounds into either of
                    3 categories: high-quality porpoise click (HQ), low-quality porpoise click (LQ), and high-frequency
                    noise (N).
        """
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
            serial_number = int(self.SerialNoEdit.text())
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
                ThisFolderSave = os.path.join(MainFolder, myFolder)
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
        self.open_ct_cancel_b.setText('Cancel')
        self.open_ct_cancel_b.setGeometry(210, 140, 100, 30)
        self.open_ct_cancel_b.clicked.connect(self.open_ct_cancel)
        # Disable resize and show menu
        self.OpenCTFig.setFixedSize(380, 210)
        self.OpenCTFig.show()

    def BrowseButtonDet(self):
        self.root_browse_button_det = tk.Tk()
        self.root_browse_button_det.withdraw()
        SelectedFolderDet = filedialog.askdirectory()
        self.FolderPathDet.setText(SelectedFolderDet)
        self.root_browse_button_det.mainloop(0)

    def BrowseButtonPorCC(self):
        self.root_browse_button_porcc = tk.Tk()
        self.root_browse_button_porcc.withdraw()
        SelectedFolderPorCC = filedialog.askdirectory()
        self.FolderPathPorCC.setText(SelectedFolderPorCC)
        self.root_browse_button_porcc.mainloop(0)

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
        self.FolderPathNewCT.setText('C:/')
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
        self.root_new_ct_browse = tk.Tk()
        self.root_new_ct_browse.withdraw()
        SelectedFolderNewCT = tk.filedialog.askdirectory()
        self.FolderPathNewCT.setText(SelectedFolderNewCT)
        self.root_new_ct_browse.mainloop(0)

    def IdentifyCT(self):
        """
            Prepares the data to locate and classify click trains, which is done via extract_patterns() and ct_type()
        """
        global CP  # sset, srise
        self.NewCTFig.close()
        MainFolder = self.FolderPathNewCT.text()
        latitude = float(self.LatEdit.text())
        longitude = float(self.LongEdit.text())
        if self.InclSubFoldersNewCT.isChecked():
            myFolders = [item for item in os.listdir(MainFolder) if os.path.isdir(os.path.join(MainFolder, item))]
            if len(myFolders) == 0:
                warnings.warn("You selected 'Include subfolders' but there are no subfolders - please try again")
                myFolders = "None"
            else:
                myFolders = [item for item in os.listdir(MainFolder) if os.path.isdir(os.path.join(MainFolder, item))]
        else:
            FilesInFolder = os.listdir(MainFolder)
            ClipsFiles = [s for s in FilesInFolder if "clips.csv" in s]
            if len(ClipsFiles) == 0:
                CPFile = [s for s in FilesInFolder if "CP.csv" in s]
                if len(CPFile) == 0:
                    warnings.warn("You did not select 'Include subfolders' and there are no files to process in this "
                                  "folder - please try again")
                    myFolders = 'None'
                else:
                    myFolders = [MainFolder]
            else:
                myFolders = [MainFolder]

        for myFolder in myFolders:
            a = 1
            if myFolders == 'None':
                warnings.warn("There were no files to analyse. Select another folder")
                FilesInFolder = []
            elif myFolders == MainFolder:
                FilesInFolder = os.listdir(MainFolder)
            else:
                FilesInFolder = os.listdir(os.path.join(MainFolder, myFolder))
            if any("." in s for s in FilesInFolder):
                CPFile = [s for s in FilesInFolder if "CP.csv" in s]
                if len(CPFile) > 0:
                    if myFolders == MainFolder:
                        print('Processing folder', MainFolder)
                        CPFileName = os.path.join(MainFolder, 'CP.csv')
                    else:
                        print('Processing folder', myFolders)
                        CPFileName = os.path.join(MainFolder, myFolder, 'CP.csv')
                    CP = pd.read_csv(CPFileName, sep=',')
                else:
                    DetClicks = [s for s in FilesInFolder if "clips.csv" in s]
                    if len(DetClicks) > 0:
                        CP = pd.DataFrame()
                        for thisFile in DetClicks:
                            if myFolders == MainFolder:
                                CPFile = os.path.join(MainFolder, thisFile)
                            else:
                                CPFile = os.path.join(MainFolder, myFolder, thisFile)
                            TempCP = pd.read_csv(CPFile, index_col=0)
                            TempCP['id'] = TempCP.index.values
                            CP = CP.append(TempCP, ignore_index=True)
                        if myFolders == MainFolder:
                            CPFileName = os.path.join(MainFolder, 'CP.csv')
                        else:
                            CPFileName = os.path.join(MainFolder, myFolder, 'CP.csv')
                        CP.to_csv(CPFileName, index=False)
                        CP = pd.read_csv(CPFileName, sep=',')

                Clicks, thisCTInfo, CP = click_trains.extract_patterns(CP, latitude, longitude)
                # TODO FREQUENCY IS HARDCODED!
                if myFolders == MainFolder:
                    CTInfoFileName = os.path.join(MainFolder, 'CTInfo.csv')
                    ClicksFileName = os.path.join(MainFolder, 'Clicks.csv')
                    thisCTInfo.to_csv(CTInfoFileName, index=False)
                    Clicks.to_csv(ClicksFileName, index=False)
                    CPFileName = os.path.join(MainFolder, 'CP.csv')
                    CP.to_csv(CPFileName, index=False)
                    self.root_new_ct_browse.mainloop(0)
                    a = a + 1
                    if a == 2:
                        break
                else:
                    CTInfoFileName = os.path.join(MainFolder, myFolder, 'CTInfo.csv')
                    ClicksFileName = os.path.join(MainFolder, myFolder, 'Clicks.csv')
                    thisCTInfo.to_csv(CTInfoFileName, index=False)
                    Clicks.to_csv(ClicksFileName, index=False)
                    CPFileName = os.path.join(MainFolder, myFolder, 'CP.csv')
                    CP.to_csv(CPFileName, index=False)
        self.root_new_ct_browse.mainloop(0)

    def open_ct_cancel(self):
        self.OpenCTFig.close()
        self.root_new_ct_browse.mainloop(0)

    # OPEN CT
    def PushBrowseButtonOpenCT(self):
        self.root_open_ct_browse_b = tk.Tk()
        self.root_open_ct_browse_b.destroy()
        self.SelectedFolderCT = filedialog.askdirectory()
        self.FolderPathOpenCT.setText(self.SelectedFolderCT)
        self.root_open_ct_browse_b.mainloop(0)

    def OpenCT(self):
        global CP, CTInfo
        self.OpenCTFig.close()
        CPFileName = os.path.join(self.SelectedFolderCT, "Clicks.csv")
        CP = pd.read_csv(CPFileName)
        CTInfoFileName = os.path.join(self.SelectedFolderCT, "CTInfo.csv")
        CTInfo = pd.read_csv(CTInfoFileName)
        CTInfo = CTInfo[CTInfo.Length > 9]
        CTInfo.reset_index(inplace=True, drop=True)
        NumCT = int(CTInfo.CTNum.iloc[0])
        self.update_ct(NumCT, CP, CTInfo)
        self.root_open_ct_browse_b.mainloop(0)

    def save_updates(self):
        FullName = os.path.join(self.SelectedFolderCT, 'CTInfo.csv')
        CTInfo.to_csv(FullName)

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
        # self.open_ct_cancel_b = QtWidgets.QPushButton(self.PorCCSetPan)
        # self.open_ct_cancel_b.setText('Cancel')
        # self.open_ct_cancel_b.setGeometry(210, 490, 100, 30)
        # self.open_ct_cancel_b.clicked.connect(self.CancelButtonPorCC)
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


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("fusion")
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
