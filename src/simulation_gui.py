import sys, os
import subprocess
import datetime

import numpy as np
import pandas as pd
from matplotlib import cm
from PIL import Image
from PIL.ImageQt import ImageQt

from PyQt5.QtWidgets import QMainWindow
from PyQt5 import QtCore, QtGui, QtWidgets

import xml.etree.cElementTree as ET



class Ui_MainWindow( QMainWindow ):
    
    def __init__( self ):
        super().__init__()
        self.setObjectName("MainWindow")
        self.resize(771, 584)
        self.centralwidget = QtWidgets.QWidget(self)
        self.centralwidget.setObjectName("centralwidget")
        self.labelForPixmap = QtWidgets.QLabel(self.centralwidget)
        self.labelForPixmap.setGeometry(QtCore.QRect(9, 6, 550, 550))
        self.labelForPixmap.setText("")
        self.labelForPixmap.setObjectName("labelForPixmap")
        self.dimensionalityLabel = QtWidgets.QLabel(self.centralwidget)
        self.dimensionalityLabel.setGeometry(QtCore.QRect(580, 10, 161, 21))
        self.dimensionalityLabel.setObjectName("dimensionalityLabel")
        self.dimensionalityComboBox = QtWidgets.QComboBox(self.centralwidget)
        self.dimensionalityComboBox.setGeometry(QtCore.QRect(580, 30, 161, 31))
        self.dimensionalityComboBox.setObjectName("dimensionalityComboBox")
        self.sizeLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.sizeLineEdit.setGeometry(QtCore.QRect(580, 60, 151, 21))
        self.sizeLineEdit.setObjectName("sizeLineEdit")
        self.neighborhoodLabel = QtWidgets.QLabel(self.centralwidget)
        self.neighborhoodLabel.setGeometry(QtCore.QRect(580, 90, 161, 21))
        self.neighborhoodLabel.setObjectName("neighborhoodLabel")
        self.neighborhoodComboBox = QtWidgets.QComboBox(self.centralwidget)
        self.neighborhoodComboBox.setGeometry(QtCore.QRect(580, 110, 161, 31))
        self.neighborhoodComboBox.setObjectName("neighborhoodComboBox")
        self.boundaryConditionsLabel = QtWidgets.QLabel(self.centralwidget)
        self.boundaryConditionsLabel.setGeometry(QtCore.QRect(580, 140, 161, 21))
        self.boundaryConditionsLabel.setObjectName("boundaryConditionsLabel")
        self.boundaryConditionsComboBox = QtWidgets.QComboBox(self.centralwidget)
        self.boundaryConditionsComboBox.setGeometry(QtCore.QRect(580, 160, 161, 31))
        self.boundaryConditionsComboBox.setObjectName("boundaryConditionsComboBox")
        self.nucleusLabel = QtWidgets.QLabel(self.centralwidget)
        self.nucleusLabel.setGeometry(QtCore.QRect(580, 190, 161, 21))
        self.nucleusLabel.setObjectName("nucleusLabel")
        self.simulationComboBox = QtWidgets.QComboBox(self.centralwidget)
        self.simulationComboBox.setGeometry(QtCore.QRect(580, 260, 161, 31))
        self.simulationComboBox.setObjectName("simulationComboBox")
        self.nucleusLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.nucleusLineEdit.setGeometry(QtCore.QRect(580, 210, 151, 21))
        self.nucleusLineEdit.setObjectName("nucleusLineEdit")
        self.simulationLabel = QtWidgets.QLabel(self.centralwidget)
        self.simulationLabel.setGeometry(QtCore.QRect(580, 240, 161, 21))
        self.simulationLabel.setObjectName("simulationLabel")
        self.iterationsNumberLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.iterationsNumberLineEdit.setGeometry(QtCore.QRect(580, 290, 151, 21))
        self.iterationsNumberLineEdit.setObjectName("iterationsNumberLineEdit")
        self.runPushButton = QtWidgets.QPushButton(self.centralwidget)
        self.runPushButton.setGeometry(QtCore.QRect(660, 320, 70, 32))
        self.runPushButton.setObjectName("runPushButton")
        self.createPushButton = QtWidgets.QPushButton(self.centralwidget)
        self.createPushButton.setGeometry(QtCore.QRect(610, 345, 100, 32))
        self.createPushButton.setObjectName("createPushButton")
        self.dispPushButton = QtWidgets.QPushButton(self.centralwidget)
        self.dispPushButton.setGeometry(QtCore.QRect(590, 320, 80, 32))
        self.dispPushButton.setObjectName("dispPushButton")
        self.displayTextEdit = QtWidgets.QTextEdit(self)
        self.displayTextEdit.setGeometry(QtCore.QRect(560, 371, 210, 188))
        self.displayTextEdit.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.displayTextEdit.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.displayTextEdit.setReadOnly(True)
        self.displayTextEdit.setObjectName("displayTextEdit")
        self.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(self)
        self.statusbar.setObjectName("statusbar")
        self.setStatusBar(self.statusbar)

        QtCore.QMetaObject.connectSlotsByName(self)

        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.dimensionalityLabel.setText(_translate("MainWindow", "Dimensionality type:"))
        self.sizeLineEdit.setText(_translate("MainWindow", ""))
        self.sizeLineEdit.setPlaceholderText("Dimensionality size") 
        self.neighborhoodLabel.setText(_translate("MainWindow", "Neighborhood type:"))
        self.boundaryConditionsLabel.setText(_translate("MainWindow", "Boundary conditions type:"))
        self.nucleusLabel.setText(_translate("MainWindow", "Number of nucleus:"))
        self.nucleusLineEdit.setText(_translate("MainWindow", ""))
        self.nucleusLineEdit.setPlaceholderText("Nucleus")
        self.simulationLabel.setText(_translate("MainWindow", "Simulation:"))
        self.iterationsNumberLineEdit.setText(_translate("MainWindow", ""))
        self.iterationsNumberLineEdit.setPlaceholderText("Iterations MC")
        self.runPushButton.setText(_translate("MainWindow", "RUN"))
        self.createPushButton.setText(_translate("MainWindow", "Create sim."))
        self.dispPushButton.setText(_translate("MainWindow", "Disp Img"))
        
        self.dimensionalityComboBox.addItems(["2D", "3D"])
        self.boundaryConditionsComboBox.addItems(["absorbing", "periodic"])
        self.neighborhoodComboBox.addItems(["von_neumann", "moore"])
        self.simulationComboBox.addItems(["none", "monte_carlo"])

        self.runPushButton.clicked.connect( self.runPushButtonClicked )
        self.createPushButton.clicked.connect( self.createPushButtonClicked )
        self.dispPushButton.clicked.connect( self.loadFunction )
    
    
    def createPushButtonClicked(self):
        file_name = self.save_data_to_xml_file()
        
        
    def save_data_to_xml_file(self) -> str:
        dimensionality_type = self.dimensionalityComboBox.currentText()
        dim_size = self.sizeLineEdit.text() if self.sizeLineEdit.text() else "2"
        neighborhood_type = self.neighborhoodComboBox.currentText()
        boundary_conditions = self.boundaryConditionsComboBox.currentText()
        nucleus_num = self.nucleusLineEdit.text() if self.nucleusLineEdit.text() else "2"
        simulation_type = self.simulationComboBox.currentText()
        iterations_MC = self.iterationsNumberLineEdit.text() if self.iterationsNumberLineEdit.text() else "2"
        
        dir_path = "data/"
        now = datetime.datetime.now()
        file_name = str( now )[:10] +'-'+ str( now )[11:-7] +':'+ str( now )[-6:-4] +'_'+ neighborhood_type +'_' \
                    + boundary_conditions +'_'+ simulation_type +'_'+ dimensionality_type +'_'+ dim_size
        #file_name = "data" ### tmp ###
        xml_file = dir_path + file_name + ".xml"
        
        GG_config = ET.Element("GG_config")

        ET.SubElement(GG_config, "dimensionalityType").text = dimensionality_type
        ET.SubElement(GG_config, "dimSize").text = dim_size
        ET.SubElement(GG_config, "CA_NeighborhoodType").text = neighborhood_type
        ET.SubElement(GG_config, "boundaryConditionsType").text = boundary_conditions
        ET.SubElement(GG_config, "numberOfNucleus").text = nucleus_num
        ET.SubElement(GG_config, "simulationType").text = simulation_type
        ET.SubElement(GG_config, "iterationsMC").text = iterations_MC
        ET.SubElement(GG_config, "outputFile").text = dir_path + file_name + ".csv"
        ET.SubElement(GG_config, "measurementsFile").text = dir_path + file_name + ".txt"

        tree = ET.ElementTree(GG_config)
        tree.write( xml_file, xml_declaration=True )
        
        return xml_file
    
    
    def runPushButtonClicked(self):
        files_list = self.getFilesList()
        self.runSimulations(files_list)
        self.setTextOnDisplayEdit(files_list)
        self.removeFiles(files_list)
    
    
    def getFilesList(self):
        dir_path = "data/"
        files_list = [ dir_path + file for file in os.listdir(dir_path) if file.endswith(".xml") ]
        return files_list
    
    
    def runSimulations(self, files_list):
        for file_path in files_list:
            subprocess.call(["./main", file_path])
    
    
    def setTextOnDisplayEdit(self, files_list):
        text = ""
        
        for file_path in files_list:
            
            text += file_path[5:-4].replace("_", " ") + "\n"
            
            with open(file_path[:-4] + ".txt") as time_file:
                for line in time_file:
                    text += " " + line
            text += "\n\n"
        
        self.displayTextEdit.clear()
        self.displayTextEdit.setPlainText(text)
    
    
    def removeFiles(self, files_list):
        for file_path in files_list:
            os.remove(file_path)
    
    
    def loadFunction(self):
        try:
            filePath, _ = QtWidgets.QFileDialog.getOpenFileName()
            df = pd.read_csv(filePath, sep=',')
            df = df[['x', 'y', 'z', 'grainId']]
        except Exception as e:
            return
        
        if (df['z'] == 0).all():
            self.disp2Dimage(df)
        else:
            self.disp3Dimage(df)
    
    
    def loadFile(self):
        try:
            filePath, _ = QtWidgets.QFileDialog.getOpenFileName()
            if not filePath or filePath[-4:] != '.csv':
                return
            df = pd.read_csv(filePath, sep=',')
            df = df[['x', 'y', 'z', 'grainId']]
        except OSError as err:
            print("OS error: {0}".format(err))
            raise
        
        return df
    
    
    def disp2Dimage(self, df):
        image_matrix = self.makePartitionMxt( df, 'x', 'y', df['x'].max()+1 )
        image_matrix = image_matrix / df['grainId'].max() # normalization <0,1>
        image  = Image.fromarray( np.uint8( cm.gist_earth(image_matrix) * 255 ) )
        pixmap = QtGui.QPixmap.fromImage( ImageQt(image) )
        self.labelForPixmap.setPixmap( pixmap.scaled( self.labelForPixmap.size(), QtCore.Qt.KeepAspectRatio ) )
    
    
    def disp3Dimage(self, df):
        size = df['x'].max() + 1
        zeros_mtx = np.zeros( (size, size), dtype=int )
        
        
        tmp_mtx_0 = np.flipud( self.makePartitionMxt(df[ df["x"]==0 ], "z", "y", size) ) # buttom
        image_matrix = np.concatenate( (zeros_mtx, tmp_mtx_0, zeros_mtx), axis=1 )
        
        tmp_mtx_0 = np.fliplr( self.makePartitionMxt(df[ df["y"]==0 ], "x", "z", size) )
        tmp_mtx_1 = self.makePartitionMxt(df[ df["z"]==0 ], "x", "y", size) # front
        tmp_mtx_2 = self.makePartitionMxt(df[ df["y"]==size-1 ], "x", "z", size)
        tmp_mtx_row = np.concatenate( (tmp_mtx_0, tmp_mtx_1, tmp_mtx_2), axis=1 )
        image_matrix = np.concatenate( (image_matrix, tmp_mtx_row), axis=0 )
        
        tmp_mtx_0 = self.makePartitionMxt(df[ df["x"]==size-1 ], "z", "y", size) # down
        tmp_mtx_row = np.concatenate( (zeros_mtx, tmp_mtx_0, zeros_mtx), axis=1 )
        image_matrix = np.concatenate( (image_matrix, tmp_mtx_row), axis=0 )
        
        tmp_mtx_0 = np.flipud( self.makePartitionMxt(df[ df["z"]==size-1 ], "x", "y", size) ) # back
        tmp_mtx_row = np.concatenate( (zeros_mtx, tmp_mtx_0, zeros_mtx), axis=1 )
        image_matrix = np.concatenate( (image_matrix, tmp_mtx_row), axis=0 )
        
        
        max_fraction = df['grainId'].max()
        image_matrix[ image_matrix==0 ] = max_fraction+1
        image_matrix = image_matrix / (max_fraction+1) # normalization <0,1>
        
        image  = Image.fromarray( np.uint8( cm.gist_earth(image_matrix) * 255 ) )
        pixmap = QtGui.QPixmap.fromImage( ImageQt(image) )
        self.labelForPixmap.setPixmap( pixmap.scaled( self.labelForPixmap.size(), QtCore.Qt.KeepAspectRatio ) )
    
    
    def makePartitionMxt(self, df, i, j, size):
        matrix = np.zeros( (size, size), dtype=int )
        
        for index, row in df.iterrows():
            matrix[ row[i] ][ row[j] ] = row['grainId']
        
        return matrix



if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    ui = Ui_MainWindow()
    ui.show()
    sys.exit(app.exec_())
