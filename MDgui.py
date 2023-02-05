# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MDgui.ui'
#
# Created by: PyQt5 UI code generator 5.14.1
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(1116, 703)
        self.frame = QtWidgets.QFrame(Dialog)
        self.frame.setGeometry(QtCore.QRect(10, 140, 221, 151))
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.label = QtWidgets.QLabel(self.frame)
        self.label.setGeometry(QtCore.QRect(10, 10, 181, 20))
        self.label.setObjectName("label")
        self.textEdit = QtWidgets.QTextEdit(self.frame)
        self.textEdit.setGeometry(QtCore.QRect(10, 40, 191, 75))
        self.textEdit.setObjectName("textEdit")
        self.textEdit.setText("1ubq.pdb")
        self.label_14 = QtWidgets.QLabel(self.frame)
        self.label_14.setGeometry(QtCore.QRect(10, 120, 191, 21))
        self.label_14.setObjectName("label_14")
        self.frame_2 = QtWidgets.QFrame(Dialog)
        self.frame_2.setGeometry(QtCore.QRect(10, 300, 221, 161))
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.textEdit_2 = QtWidgets.QTextEdit(self.frame_2)
        self.textEdit_2.setGeometry(QtCore.QRect(10, 50, 191, 75))
        self.textEdit_2.setObjectName("textEdit_2")
        self.textEdit_2.setText("leaprc.protein.ff14SB")
        self.label_15 = QtWidgets.QLabel(self.frame_2)
        self.label_15.setGeometry(QtCore.QRect(10, 130, 201, 21))
        self.label_15.setObjectName("label_15")
        self.label_2 = QtWidgets.QLabel(self.frame_2)
        self.label_2.setGeometry(QtCore.QRect(10, 20, 201, 20))
        self.label_2.setObjectName("label_2")
        self.frame_3 = QtWidgets.QFrame(Dialog)
        self.frame_3.setGeometry(QtCore.QRect(260, 120, 311, 361))
        self.frame_3.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3.setObjectName("frame_3")
        self.lineEdit_3 = QtWidgets.QLineEdit(self.frame_3)
        self.lineEdit_3.setGeometry(QtCore.QRect(20, 70, 113, 28))
        self.lineEdit_3.setObjectName("lineEdit_3")
        self.lineEdit_3.setText("12")
        self.label_8 = QtWidgets.QLabel(self.frame_3)
        self.label_8.setGeometry(QtCore.QRect(20, 50, 271, 21))
        self.label_8.setObjectName("label_8")
        self.lineEdit_4 = QtWidgets.QLineEdit(self.frame_3)
        self.lineEdit_4.setGeometry(QtCore.QRect(20, 140, 113, 28))
        self.lineEdit_4.setObjectName("lineEdit_4")
        self.lineEdit_4.setText("1")
        self.label_9 = QtWidgets.QLabel(self.frame_3)
        self.label_9.setGeometry(QtCore.QRect(20, 110, 261, 21))
        self.label_9.setObjectName("label_9")
        self.label_10 = QtWidgets.QLabel(self.frame_3)
        self.label_10.setGeometry(QtCore.QRect(20, 170, 271, 21))
        self.label_10.setObjectName("label_10")
        self.lineEdit_5 = QtWidgets.QLineEdit(self.frame_3)
        self.lineEdit_5.setGeometry(QtCore.QRect(20, 190, 113, 28))
        self.lineEdit_5.setObjectName("lineEdit_5")
        self.lineEdit_5.setText("50")
        self.label_11 = QtWidgets.QLabel(self.frame_3)
        self.label_11.setGeometry(QtCore.QRect(20, 230, 281, 21))
        self.label_11.setObjectName("label_11")
        self.lineEdit_6 = QtWidgets.QLineEdit(self.frame_3)
        self.lineEdit_6.setGeometry(QtCore.QRect(20, 260, 113, 28))
        self.lineEdit_6.setObjectName("lineEdit_6")
        self.lineEdit_6.setText("10")
        self.lineEdit_7 = QtWidgets.QLineEdit(self.frame_3)
        self.lineEdit_7.setGeometry(QtCore.QRect(20, 320, 261, 28))
        self.lineEdit_7.setObjectName("lineEdit_7")
        self.lineEdit_7.setText("/home/gabsbi/Desktop/amber18/dat/leap/cmd/")
        self.label_12 = QtWidgets.QLabel(self.frame_3)
        self.label_12.setGeometry(QtCore.QRect(20, 290, 281, 21))
        self.label_12.setObjectName("label_12")
        self.label_19 = QtWidgets.QLabel(self.frame_3)
        self.label_19.setGeometry(QtCore.QRect(20, 10, 271, 31))
        self.label_19.setObjectName("label_19")
        self.frame_5 = QtWidgets.QFrame(Dialog)
        self.frame_5.setGeometry(QtCore.QRect(10, 500, 311, 181))
        self.frame_5.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_5.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_5.setObjectName("frame_5")
        self.pushButton_2 = QtWidgets.QPushButton(self.frame_5)
        self.pushButton_2.setGeometry(QtCore.QRect(110, 130, 81, 28))
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton = QtWidgets.QPushButton(self.frame_5)
        self.pushButton.setGeometry(QtCore.QRect(20, 50, 281, 28))
        self.pushButton.setObjectName("pushButton")
        self.pushButton_3 = QtWidgets.QPushButton(self.frame_5)
        self.pushButton_3.setGeometry(QtCore.QRect(20, 90, 281, 28))
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton.clicked.connect(self.run_ambertools)       
        self.pushButton_2.clicked.connect(self.closeIt)
        self.pushButton_3.clicked.connect(self.run_openmm)
        self.label_18 = QtWidgets.QLabel(self.frame_5)
        self.label_18.setGeometry(QtCore.QRect(80, 10, 181, 31))
        self.label_18.setObjectName("label_18")
        self.frame_6 = QtWidgets.QFrame(Dialog)
        self.frame_6.setGeometry(QtCore.QRect(560, 280, 511, 201))
        self.frame_6.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_6.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_6.setObjectName("frame_6")
        self.label_13 = QtWidgets.QLabel(self.frame_6)
        self.label_13.setGeometry(QtCore.QRect(20, 10, 261, 21))
        self.label_13.setObjectName("label_13")
        self.textBrowser_2 = QtWidgets.QTextBrowser(self.frame_6)
        self.textBrowser_2.setGeometry(QtCore.QRect(20, 40, 471, 151))
        self.textBrowser_2.setObjectName("textBrowser_2")
        self.frame_7 = QtWidgets.QFrame(Dialog)
        self.frame_7.setGeometry(QtCore.QRect(580, 130, 121, 141))
        self.frame_7.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_7.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_7.setObjectName("frame_7")
        self.label_16 = QtWidgets.QLabel(self.frame_7)
        self.label_16.setGeometry(QtCore.QRect(10, 10, 101, 31))
        self.label_16.setObjectName("label_16")
        self.radioButton = QtWidgets.QRadioButton(self.frame_7)
        self.radioButton.setGeometry(QtCore.QRect(10, 60, 104, 21))
        self.radioButton.setObjectName("radioButton")
        self.radioButton.setChecked(False)
        self.radioButton_2 = QtWidgets.QRadioButton(self.frame_7)
        self.radioButton_2.setGeometry(QtCore.QRect(10, 100, 104, 21))
        self.radioButton_2.setObjectName("radioButton_2")
        self.radioButton_2.setChecked(True)
        self.frame_8 = QtWidgets.QFrame(Dialog)
        self.frame_8.setGeometry(QtCore.QRect(360, 500, 711, 181))
        self.frame_8.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_8.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_8.setObjectName("frame_8")
        self.checkBox = QtWidgets.QCheckBox(self.frame_8)
        self.checkBox.setGeometry(QtCore.QRect(20, 40, 801, 31))
        self.checkBox.setObjectName("checkBox")
        self.checkBox.setChecked(True)
        self.checkBox_2 = QtWidgets.QCheckBox(self.frame_8)
        self.checkBox_2.setGeometry(QtCore.QRect(20, 70, 791, 21))
        self.checkBox_2.setObjectName("checkBox_2")
        self.checkBox_2.setChecked(False)
        self.checkBox_3 = QtWidgets.QCheckBox(self.frame_8)
        self.checkBox_3.setGeometry(QtCore.QRect(20, 130, 791, 21))
        self.checkBox_3.setObjectName("checkBox_3")
        self.checkBox_3.setChecked(True)
        self.checkBox_4 = QtWidgets.QCheckBox(self.frame_8)
        self.checkBox_4.setGeometry(QtCore.QRect(20, 100, 781, 21))
        self.checkBox_4.setObjectName("checkBox_4")
        self.checkBox_4.setChecked(False)
        self.label_17 = QtWidgets.QLabel(self.frame_8)
        self.label_17.setGeometry(QtCore.QRect(30, 10, 371, 31))
        self.label_17.setObjectName("label_17")
        self.label_7 = QtWidgets.QLabel(Dialog)
        self.label_7.setGeometry(QtCore.QRect(140, 10, 791, 91))
        self.label_7.setStyleSheet("background-color: rgb(0, 0, 0);")
        self.label_7.setObjectName("label_7")
        self.label_4 = QtWidgets.QLabel(Dialog)
        self.label_4.setGeometry(QtCore.QRect(20, 10, 211, 91))
        self.label_4.setText("")
        self.label_4.setPixmap(QtGui.QPixmap("droid.png"))
        self.label_4.setScaledContents(True)
        self.label_4.setObjectName("label_4")
        self.label_3 = QtWidgets.QLabel(Dialog)
        self.label_3.setGeometry(QtCore.QRect(860, 10, 211, 91))
        self.label_3.setText("")
        self.label_3.setPixmap(QtGui.QPixmap("maxdemon.png"))
        self.label_3.setScaledContents(True)
        self.label_3.setObjectName("label_3")
        self.frame_4 = QtWidgets.QFrame(Dialog)
        self.frame_4.setGeometry(QtCore.QRect(720, 130, 351, 141))
        self.frame_4.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_4.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_4.setObjectName("frame_4")
        self.textBrowser = QtWidgets.QTextBrowser(self.frame_4)
        self.textBrowser.setGeometry(QtCore.QRect(20, 50, 261, 71))
        self.textBrowser.setObjectName("textBrowser")
        self.label_20 = QtWidgets.QLabel(self.frame_4)
        self.label_20.setGeometry(QtCore.QRect(20, 20, 261, 21))
        self.label_20.setObjectName("label_20")

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.label.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-weight:600;\">PDB file list (up to 5)</span></p></body></html>"))
        self.label_14.setText(_translate("Dialog", "<html><head/><body><p>e.g. 1ubq.pdb</p></body></html>"))
        
        self.label_15.setText(_translate("Dialog", "<html><head/><body><p>e.g. leaprc.protein.ff14SB</p></body></html>"))
        self.label_2.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-weight:600;\">force field list (up to 5)</span></p></body></html>"))
        self.label_8.setText(_translate("Dialog", "<html><head/><body><p>size of water box (nm-octrahedral)</p></body></html>"))
        self.label_9.setText(_translate("Dialog", "<html><head/><body><p>length of MD heating (ns)</p></body></html>"))
        self.label_10.setText(_translate("Dialog", "<html><head/><body><p>length of MD equilibration (ns)</p></body></html>"))
        self.label_11.setText(_translate("Dialog", "<html><head/><body><p>length of MD production run (ns)</p></body></html>"))
        self.label_12.setText(_translate("Dialog", "<html><head/><body><p>path to force field folder</p></body></html>"))
        self.label_19.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600; color:#0000ca;\">MD run parameters</span></p></body></html>"))
        self.pushButton_2.setText(_translate("Dialog", "exit"))
        self.pushButton.setText(_translate("Dialog", "pre-processing (AmberTools)"))
        self.pushButton_3.setText(_translate("Dialog", "run MD simulation (openMM)"))
        self.label_18.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600; color:#d00000;\">program control</span></p></body></html>"))
        self.label_13.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-weight:600;\">system dependencies</span></p></body></html>"))
        self.textBrowser_2.setHtml(_translate("Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">BabbittLab at RIT</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">https://people.rit.edu/gabsbi/</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:14pt; font-weight:600;\">dependencies</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; line-height:100%; background-color:transparent;\"><span style=\" font-weight:600; background-color:transparent;\">CUDA graphics toolkit library</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; line-height:100%; background-color:transparent;\"><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\">https://developer.nvidia.com/cuda-downloads</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; line-height:100%; background-color:transparent;\"><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt; font-weight:600;\">AmberTools</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; line-height:100%; background-color:transparent;\"><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\">https://ambermd.org/</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; line-height:100%; background-color:transparent;\"><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\">install via conda (preferred)</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; line-height:100%; background-color:transparent;\"><span style=\" font-family:\'monospace\'; font-size:large; color:#000000;\">conda create --name AmberTools22   conda activate AmberTools22</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; line-height:100%; background-color:transparent;\"><span style=\" font-family:\'monospace\'; font-size:large; color:#000000;\">conda install -c conda-forge ambertools=22 compilers</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; line-height:100%; background-color:transparent;\"><span style=\" font-family:\'Open Sans, serif\'; font-size:12pt; font-weight:600;\">openMM </span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\">(also conda installed)</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; line-height:100%; background-color:transparent;\"><span style=\" font-family:\'Consolas,Menlo,DejaVu Sans Mono,Bitstream Vera Sans Mono,monospace\'; color:#3e4349;\">conda install </span><span style=\" font-family:\'Consolas,Menlo,DejaVu Sans Mono,Bitstream Vera Sans Mono,monospace\'; color:#666666;\">-</span><span style=\" font-family:\'Consolas,Menlo,DejaVu Sans Mono,Bitstream Vera Sans Mono,monospace\'; color:#3e4349;\">c conda</span><span style=\" font-family:\'Consolas,Menlo,DejaVu Sans Mono,Bitstream Vera Sans Mono,monospace\'; color:#666666;\">-</span><span style=\" font-family:\'Consolas,Menlo,DejaVu Sans Mono,Bitstream Vera Sans Mono,monospace\'; color:#3e4349;\">forge openmm cudatoolkit</span><span style=\" font-family:\'Consolas,Menlo,DejaVu Sans Mono,Bitstream Vera Sans Mono,monospace\'; color:#666666;\">=</span><span style=\" font-family:\'Consolas,Menlo,DejaVu Sans Mono,Bitstream Vera Sans Mono,monospace\'; color:#208050;\">10.0</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">for anaconda go here</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">https://www.anaconda.com/products/distribution</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.label_16.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600; color:#121212;\">solvation</span></p></body></html>"))
        self.radioButton.setText(_translate("Dialog", "implicit"))
        self.radioButton_2.setText(_translate("Dialog", "explicit"))
        self.checkBox.setText(_translate("Dialog", "reduce PDB structure (add H) and remove waters (pdb4amber)"))
        self.checkBox_2.setText(_translate("Dialog", "run force field modifications for small molecule via sqm (antechamber)"))
        self.checkBox_3.setText(_translate("Dialog", "create topology and input coordinates for explicit solvent system (tleap)"))
        self.checkBox_4.setText(_translate("Dialog", "create topology and input coordinates for implicit solvent system (tleap)"))
        self.label_17.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600; color:#222dc8;\">MD pre-processing options</span></p></body></html>"))
        self.label_7.setText(_translate("Dialog", "<html><head/><body><p align=\"center\"><span style=\" font-size:16pt; color:#ffffff;\">AmberTools/openMM</span></p><p align=\"center\"><span style=\" font-size:16pt; color:#ffffff;\">Comparative Molecular Dynamics Simulator</span></p></body></html>"))
        self.textBrowser.setHtml(_translate("Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">/path2file/1cdw_bound.pdb</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">/path2file/1cdw_unbound.pdb</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">/path2file/1cdw_ortholog.pdb</p></body></html>"))
        self.label_20.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-weight:600;\">file list example </span></p></body></html>"))





##########################################################################################
#################   button actions    ####################################################
##########################################################################################
# create run program action  (make DROIDS.ctl and launch pipeline)
    def run_ambertools(self):
        # check and set graph style
        print("writing control file")
        box_size = self.lineEdit_3.text()
        box_size = int(box_size)
        #print(box_size)
        heat_len = self.lineEdit_4.text()
        heat_len = int(heat_len)
        heat_len = heat_len*1000000
        #print(heat_len)
        eq_len = self.lineEdit_5.text()
        eq_len = int(eq_len)
        eq_len = eq_len*1000000
        #print(eq_len)
        prod_len = self.lineEdit_6.text()
        prod_len = int(prod_len)
        prod_len = prod_len*1000000
        #print(prod_len)
            
        tleap_path = self.lineEdit_7.text()
        #print(tleap_path)
        pdb_list = self.textEdit.toPlainText()
        #print(pdb_list)
        pdb_list = str.split(pdb_list, "\n")
        #print(pdb_list)
        n_pdb = len(pdb_list)
        if(n_pdb == 1 or n_pdb == 2 or n_pdb == 3 or n_pdb == 4 or n_pdb == 5):
            first_pdb = pdb_list[0]
        if(n_pdb == 2 or n_pdb == 3 or n_pdb == 4 or n_pdb == 5):
            second_pdb = pdb_list[1]
        if(n_pdb == 3 or n_pdb == 4 or n_pdb == 5):
            third_pdb = pdb_list[2]
        if(n_pdb == 4 or n_pdb == 5):
            fourth_pdb = pdb_list[3]
        if(n_pdb == 5):
            fifth_pdb = pdb_list[4]
        ff_list = self.textEdit_2.toPlainText()
        #print(ff_list)
        ff_list = str.split(ff_list, "\n")
        #print(ff_list)
        n_ff = len(ff_list)
        if(n_ff == 1 or n_ff == 2 or n_ff == 3 or n_ff == 4 or n_ff == 5):
            first_ff = ff_list[0]
        if(n_ff == 2 or n_ff == 3 or n_ff == 4 or n_ff == 5):
            second_ff = ff_list[1]
        if(n_ff == 3 or n_ff == 4 or n_ff == 5):
            third_ff = ff_list[2]
        if(n_ff == 4 or n_ff == 5):
            fourth_ff = ff_list[3]
        if(n_ff == 5):
            fifth_ff = ff_list[4]
        if self.radioButton_2.isChecked() == True:
            solvation = "explicit"
        else:
            solvation = "implicit"
        
        # analyses to run
        if self.checkBox.isChecked() == True:
            pdb4amber = "yes"
        else:
            pdb4amber = "no"
        if self.checkBox_2.isChecked() == True:
            antechamber = "yes"
        else:
            antechamber = "no"
        if self.checkBox_3.isChecked() == True:
            tleap_explicit = "yes"
        else:
            tleap_explicit = "no"
        if self.checkBox_4.isChecked() == True:
            tleap_implicit = "yes"
        else:
            tleap_implicit = "no"
        
        # write file
        out = open("./MDr.ctl", "w") 
        if(n_pdb == 1 or n_pdb == 2 or n_pdb == 3 or n_pdb == 4 or n_pdb == 5):
            out.write("firstID,%s,#pdb id for first structure\n" % first_pdb)
        if(n_pdb == 2 or n_pdb == 3 or n_pdb == 4 or n_pdb == 5):
            out.write("secondID,%s,#pdb id for second structure\n" % second_pdb)
        if(n_pdb == 3 or n_pdb == 4 or n_pdb == 5):
            out.write("thirdID,%s,#pdb id for third structure\n" % third_pdb)
        if(n_pdb == 4 or n_pdb == 5):
            out.write("fourthID,%s,#pdb id for fourth structure\n" % fourth_pdb)
        if(n_pdb == 5):
            out.write("fifthID,%s,#pdb id for fifth structure\n" % fifth_pdb)
        if(n_ff == 1 or n_ff == 2 or n_ff == 3 or n_ff == 4 or n_ff == 5):
            out.write("firstFF,%s,#pdb id for first force field\n" % first_ff)
        if(n_ff == 2 or n_ff == 3 or n_ff == 4 or n_ff == 5):
            out.write("secondFF,%s,#pdb id for second force field\n" % second_ff)
        if(n_ff == 3 or n_ff == 4 or n_ff == 5):
            out.write("thirdFF,%s,#pdb id for third force field\n" % third_ff)
        if(n_ff == 4 or n_ff == 5):
            out.write("fourthFF,%s,#pdb id for fourth force field\n" % fourth_ff)
        if(n_ff == 5):
            out.write("fifthFF,%s,#pdb id for fifth force field\n" % fifth_ff)
        out.write("box_size,%s,#water box size\n" % box_size)
        out.write("heat_len,%s,#length of heating\n" % heat_len)
        out.write("eq_len,%s,#length of equilibration\n" % eq_len)
        out.write("prod_len,%s,#length of production\n" % prod_len)
        out.write("pdb4amber,%s,#run pdb4amber\n" % pdb4amber)
        out.write("antechamber,%s,#run antechamber\n" % antechamber)
        out.write("tleap_implicit,%s,#run tleap\n" % tleap_implicit)        
        out.write("tleap_explicit,%s,#run tleap\n" % tleap_explicit)
        out.write("tleap_path,%s,#path to force field folder\n" % tleap_path)
        out.close()    
                
        print("\nIMPORTANT: CHECK THAT RESIDUES IN ALL CHAINS ARE NUMBERED SEQUENTIALLY STARTING FROM 1 to END OF LAST CHAIN\n")
        print("\nALSO: SET STARTING CHAIN LABEL TO A, put additional things (ligands, DNA etc) at end of PDB file and SAVE AS SINGLE MODEL (i.e. may need to remove model X and endmdl lines manually\n")
               
        user_input = input("\nIs AmberToolsXX installed in the AmberToolsXX conda environment? (yes/no)\n\n")
        
        if(user_input == "yes" or user_input == "y"):
            #what_version = input("\n Enter AmberTools version number (e.g. 22)")
            #os.system("conda activate AmberTools%s\n" % what_version)
            print("\nRUN THE FOLLOWING CMD'S/SCRIPTS IN THE NEW TERMINAL\n")
            print("conda activate OR conda activate AmberTools22 (as set on your system)")
            print("python3 MD_protein_pdb4amber.py")
            print("python3 MD_protein_antechamber.py (if needed)")
            print("python3 MD_protein_tleap_explicit.py")
            print("conda deactivate")
            print("\n\n")
            #os.system("conda deactivate\n")
            print("CLOSE TERMINAL WHEN MD SIMULATION PREPARATIONS ARE COMPLETED\n\n")
            os.system("x-terminal-emulator\n")
            
        if(user_input == "no" or user_input == "n"):
            if(pdb4amber == "yes"):
                print("reducing/drying PDB structure")
                os.system("python3 MD_protein_pdb4amber.py")
            if(antechamber == "yes"):
                print("modifying GAFF2 force field")
                os.system("python3 MD_protein_antechamber.py")
            if(tleap_explicit == "yes" or tleap_implicit == "yes"):
                print("building final files for MD simulation") 
                os.system("python3 MD_protein_tleap.py")
                
        print("\nPLEASE CHECK TERMINAL OUTPUT FOR WARNINGS\n")
        print("\nMD simulation prep stages are completed\n")
                
    def run_openmm(self):
        
        print("running molecular dynamics simulation")
        user_input = input("\nIs openMM installed in the base conda environment? (yes/no)\n\n")
        
        if(user_input == "yes" or user_input == "y"):
            #os.system("conda config --set auto_activate_base true\n")
            print("\nRUN THE FOLLOWING CMD's/SCRIPTS IN THE NEW TERMINAL\n")
            print("conda activate OR conda activate openmm (as set on your system)")
            print("python3 MD_protein_openMM.py")
            print("conda deactivate\n")
            print("optional- to monitor system during simulation, open additional terminals and run 'top' or 'htop' and 'nvidia-smi -l 10'")
            print("\n\n")
            #os.system("conda config --set auto_activate_base false\n")
            print("CLOSE TERMINAL WHEN MD SIMULATION IS COMPLETED\n\n")
            os.system("x-terminal-emulator\n")
        
        if(user_input == "no" or user_input == "n"):
            os.system("x-terminal-emulator -e top\n")
            os.system("x-terminal-emulator -e nvidia-smi -l 10\n")
            cmd = "python3 MD_protein_openMM.py"
            os.system(cmd)
        
        print("\nPLEASE CHECK TERMINAL OUTPUT FOR WARNINGS\n")
        print("\nMD simulation run stages are completed\n")         
            
    def closeIt(self):
        print("MDgui.py program closed")
        sys.exit(app.exec_())

#########################################################################################
#########################################################################################


if __name__ == "__main__":
    import sys
    import os
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())
