# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'droidsTEMPLATE.ui'
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
        #self.textEdit.setText("1cdw_bound.pdb\n1cdw_bound.prmtop\n1cdw_bound.nc")
        self.label_14 = QtWidgets.QLabel(self.frame)
        self.label_14.setGeometry(QtCore.QRect(10, 120, 191, 16))
        self.label_14.setObjectName("label_14")
        self.frame_2 = QtWidgets.QFrame(Dialog)
        self.frame_2.setGeometry(QtCore.QRect(10, 300, 221, 161))
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.textEdit_2 = QtWidgets.QTextEdit(self.frame_2)
        self.textEdit_2.setGeometry(QtCore.QRect(10, 50, 191, 75))
        self.textEdit_2.setObjectName("textEdit_2")
        #self.textEdit_2.setText("1cdw_unbound.pdb\n1cdw_unbound.prmtop\n1cdw_unbound.nc")
        self.label_15 = QtWidgets.QLabel(self.frame_2)
        self.label_15.setGeometry(QtCore.QRect(10, 130, 201, 16))
        self.label_15.setObjectName("label_15")
        self.label_2 = QtWidgets.QLabel(self.frame_2)
        self.label_2.setGeometry(QtCore.QRect(10, 20, 201, 20))
        self.label_2.setObjectName("label_2")
        self.frame_3 = QtWidgets.QFrame(Dialog)
        self.frame_3.setGeometry(QtCore.QRect(240, 100, 311, 381))
        self.frame_3.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3.setObjectName("frame_3")
        self.lineEdit = QtWidgets.QLineEdit(self.frame_3)
        self.lineEdit.setGeometry(QtCore.QRect(20, 30, 113, 28))
        self.lineEdit.setObjectName("lineEdit")
        #self.lineEdit.setText("25")
        self.label_5 = QtWidgets.QLabel(self.frame_3)
        self.label_5.setGeometry(QtCore.QRect(20, 10, 261, 16))
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(self.frame_3)
        self.label_6.setGeometry(QtCore.QRect(20, 60, 271, 16))
        self.label_6.setObjectName("label_6")
        self.lineEdit_2 = QtWidgets.QLineEdit(self.frame_3)
        self.lineEdit_2.setGeometry(QtCore.QRect(20, 80, 113, 28))
        self.lineEdit_2.setObjectName("lineEdit_2")
        #self.lineEdit_2.setText("100")
        self.lineEdit_3 = QtWidgets.QLineEdit(self.frame_3)
        self.lineEdit_3.setGeometry(QtCore.QRect(20, 130, 113, 28))
        self.lineEdit_3.setObjectName("lineEdit_3")
        #self.lineEdit_3.setText("5000")
        self.label_8 = QtWidgets.QLabel(self.frame_3)
        self.label_8.setGeometry(QtCore.QRect(20, 110, 271, 21))
        self.label_8.setObjectName("label_8")
        self.lineEdit_4 = QtWidgets.QLineEdit(self.frame_3)
        self.lineEdit_4.setGeometry(QtCore.QRect(20, 180, 113, 28))
        self.lineEdit_4.setObjectName("lineEdit_4")
        #self.lineEdit_4.setText("1")
        self.label_9 = QtWidgets.QLabel(self.frame_3)
        self.label_9.setGeometry(QtCore.QRect(20, 160, 261, 16))
        self.label_9.setObjectName("label_9")
        self.label_10 = QtWidgets.QLabel(self.frame_3)
        self.label_10.setGeometry(QtCore.QRect(20, 210, 271, 16))
        self.label_10.setObjectName("label_10")
        self.lineEdit_5 = QtWidgets.QLineEdit(self.frame_3)
        self.lineEdit_5.setGeometry(QtCore.QRect(20, 230, 113, 28))
        self.lineEdit_5.setObjectName("lineEdit_5")
        #self.lineEdit_5.setText("179")
        self.label_11 = QtWidgets.QLabel(self.frame_3)
        self.label_11.setGeometry(QtCore.QRect(20, 260, 281, 16))
        self.label_11.setObjectName("label_11")
        self.lineEdit_6 = QtWidgets.QLineEdit(self.frame_3)
        self.lineEdit_6.setGeometry(QtCore.QRect(20, 280, 113, 28))
        self.lineEdit_6.setObjectName("lineEdit_6")
        #self.lineEdit_6.setText("1")
        self.lineEdit_7 = QtWidgets.QLineEdit(self.frame_3)
        self.lineEdit_7.setGeometry(QtCore.QRect(20, 330, 113, 28))
        self.lineEdit_7.setObjectName("lineEdit_7")
        self.lineEdit_7.setText("/usr/lib/ucsf-chimerax/bin/")
        self.label_12 = QtWidgets.QLabel(self.frame_3)
        self.label_12.setGeometry(QtCore.QRect(20, 310, 281, 16))
        self.label_12.setObjectName("label_12")
        self.label_19 = QtWidgets.QLabel(self.frame_3)
        self.label_19.setGeometry(QtCore.QRect(230, 10, 111, 31))
        self.label_19.setObjectName("label_19")
        self.frame_5 = QtWidgets.QFrame(Dialog)
        self.frame_5.setGeometry(QtCore.QRect(10, 500, 211, 181))
        self.frame_5.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_5.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_5.setObjectName("frame_5")
        self.pushButton_2 = QtWidgets.QPushButton(self.frame_5)
        self.pushButton_2.setGeometry(QtCore.QRect(70, 140, 81, 28))
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.clicked.connect(self.closeIt)
        self.pushButton = QtWidgets.QPushButton(self.frame_5)
        self.pushButton.setGeometry(QtCore.QRect(20, 50, 171, 28))
        self.pushButton.setObjectName("pushButton")
        self.pushButton.clicked.connect(self.run_program)
        self.pushButton_3 = QtWidgets.QPushButton(self.frame_5)
        self.pushButton_3.setGeometry(QtCore.QRect(40, 90, 131, 28))
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_3.clicked.connect(self.run_analyses)
        self.label_18 = QtWidgets.QLabel(self.frame_5)
        self.label_18.setGeometry(QtCore.QRect(20, 10, 181, 31))
        self.label_18.setObjectName("label_18")
        self.frame_6 = QtWidgets.QFrame(Dialog)
        self.frame_6.setGeometry(QtCore.QRect(560, 280, 511, 201))
        self.frame_6.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_6.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_6.setObjectName("frame_6")
        self.label_13 = QtWidgets.QLabel(self.frame_6)
        self.label_13.setGeometry(QtCore.QRect(20, 10, 261, 16))
        self.label_13.setObjectName("label_13")
        self.textBrowser_2 = QtWidgets.QTextBrowser(self.frame_6)
        self.textBrowser_2.setGeometry(QtCore.QRect(20, 40, 471, 141))
        self.textBrowser_2.setObjectName("textBrowser_2")
        self.frame_7 = QtWidgets.QFrame(Dialog)
        self.frame_7.setGeometry(QtCore.QRect(580, 130, 121, 141))
        self.frame_7.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_7.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_7.setObjectName("frame_7")
        self.label_16 = QtWidgets.QLabel(self.frame_7)
        self.label_16.setGeometry(QtCore.QRect(10, 10, 81, 31))
        self.label_16.setObjectName("label_16")
        self.radioButton = QtWidgets.QRadioButton(self.frame_7)
        self.radioButton.setGeometry(QtCore.QRect(10, 60, 104, 21))
        self.radioButton.setObjectName("radioButton")
        self.radioButton.setChecked(True)
        self.radioButton_2 = QtWidgets.QRadioButton(self.frame_7)
        self.radioButton_2.setGeometry(QtCore.QRect(10, 100, 104, 21))
        self.radioButton_2.setObjectName("radioButton_2")
        self.radioButton_2.setChecked(False)
        self.frame_8 = QtWidgets.QFrame(Dialog)
        self.frame_8.setGeometry(QtCore.QRect(250, 500, 821, 181))
        self.frame_8.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_8.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_8.setObjectName("frame_8")
        self.checkBox = QtWidgets.QCheckBox(self.frame_8)
        self.checkBox.setGeometry(QtCore.QRect(20, 40, 801, 31))
        self.checkBox.setObjectName("checkBox")
        self.checkBox.setChecked(True)
        self.checkBox_2 = QtWidgets.QCheckBox(self.frame_8)
        self.checkBox_2.setGeometry(QtCore.QRect(20, 80, 791, 21))
        self.checkBox_2.setObjectName("checkBox_2")
        self.checkBox_2.setChecked(False)
        self.checkBox_3 = QtWidgets.QCheckBox(self.frame_8)
        self.checkBox_3.setGeometry(QtCore.QRect(20, 140, 791, 21))
        self.checkBox_3.setObjectName("checkBox_3")
        self.checkBox_3.setChecked(False)
        self.checkBox_4 = QtWidgets.QCheckBox(self.frame_8)
        self.checkBox_4.setGeometry(QtCore.QRect(20, 110, 781, 21))
        self.checkBox_4.setObjectName("checkBox_4")
        self.checkBox_4.setChecked(False)
        self.label_17 = QtWidgets.QLabel(self.frame_8)
        self.label_17.setGeometry(QtCore.QRect(160, 0, 371, 31))
        self.label_17.setObjectName("label_17")
        #self.checkBox_5 = QtWidgets.QCheckBox(self.frame_8)
        #self.checkBox_5.setGeometry(QtCore.QRect(20, 130, 791, 21))
        #self.checkBox_5.setObjectName("checkBox_5")
        #self.checkBox_5.setChecked(False)
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
        self.label_20.setGeometry(QtCore.QRect(20, 20, 261, 16))
        self.label_20.setObjectName("label_20")

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.label.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:11pt; font-weight:600;\">query protein file list</span></p></body></html>"))
        self.label_14.setText(_translate("Dialog", "<html><head/><body><p>i.e. bound or mutant</p></body></html>"))
        self.label_15.setText(_translate("Dialog", "<html><head/><body><p>i.e.unbound or wildtype</p></body></html>"))
        self.label_2.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-weight:600;\">reference protein file list</span></p></body></html>"))
        self.label_5.setText(_translate("Dialog", "<html><head/><body><p>number of subsamples</p></body></html>"))
        self.label_6.setText(_translate("Dialog", "<html><head/><body><p>frames iper subsample (e.g. 100)</p></body></html>"))
        self.label_8.setText(_translate("Dialog", "<html><head/><body><p>total number of frames</p></body></html>"))
        self.label_9.setText(_translate("Dialog", "<html><head/><body><p>number of protein chains</p></body></html>"))
        self.label_10.setText(_translate("Dialog", "<html><head/><body><p>number of AA sites in protein</p></body></html>"))
        self.label_11.setText(_translate("Dialog", "<html><head/><body><p>start number on N terminus (e.g. 1)</p></body></html>"))
        self.label_12.setText(_translate("Dialog", "<html><head/><body><p>chimerax path /usr/lib/ucsf-chimerax/bin/</p></body></html>"))
        self.label_19.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600;\">inputs</span></p></body></html>"))
        self.pushButton_2.setText(_translate("Dialog", "exit"))
        self.pushButton.setText(_translate("Dialog", "run MD sampling"))
        self.pushButton_3.setText(_translate("Dialog", "run analyses"))
        self.label_18.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600;\">program control</span></p></body></html>"))
        self.label_13.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-weight:600;\">more information</span></p></body></html>"))
        self.textBrowser_2.setHtml(_translate("Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">BabbittLab at RIT</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">https://people.rit.edu/gabsbi/</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">citations</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; line-height:100%; background-color:transparent;\"><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\">Babbitt G.A.</span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt; font-style:italic;\"> Coppola E.E</span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\">. </span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt; font-style:italic;\">Mortensen J.S.</span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\"> </span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt; font-style:italic;\">Adams L.E. Liao J. K.</span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\"> 2018. DROIDS 1.2 – a GUI-based pipeline for GPU-accelerated comparative protein dynamics. BIOPHYSICAL JOURNAL 114: 1009-1017. CELL Press.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; line-height:100%; background-color:transparent;\"><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\">Babbitt G.A.</span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt; font-style:italic;\"> </span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\">Fokoue E.</span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt; font-style:italic;\"> Evans J.R. Diller K.I. Adams L.E. </span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\">2020. DROIDS 3.0 - Detection of genetic and drug class variant impact on conserved protein binding dynamics. BIOPHYSICAL JOURNAL 118: 541-551 CELL Press.</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; line-height:100%; background-color:transparent;\"><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\">Babbitt G.A. Fokoue E.P. </span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt; font-style:italic;\">Srivastava H.R. Callahan B. Rajendran M</span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\">. 2022.</span><span style=\" font-family:\'ArialUnicodeMS, serif\'; font-size:14pt; color:#333666;\"> </span><span style=\" font-family:\'Open Sans, serif\'; font-size:11pt;\">Statistical machine learning for comparative protein dynamics with the DROIDS/maxDemon software pipeline. STAR PROTOCOLS 3(1):101194</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.label_16.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600;\">graphs</span></p></body></html>"))
        self.radioButton.setText(_translate("Dialog", "light"))
        self.radioButton_2.setText(_translate("Dialog", "dark"))
        self.checkBox.setText(_translate("Dialog", "site-wise comparison of atom fluctuations (signed symmetric KL divergence)"))
        self.checkBox_2.setText(_translate("Dialog", "site-wise comparison of atom correlations across the protein (MMD on learned features)"))
        self.checkBox_3.setText(_translate("Dialog", "site-wise identification of conserved molecular dynamics (ortholog vs neutral learning profiles)"))
        self.checkBox_4.setText(_translate("Dialog", "site-wise identification of coordinated dynamics (mutual information on learned classifications)"))
        #self.checkBox_5.setText(_translate("Dialog", "site-wise comparison of genetic/drug class variants (divergence metrics on conserved dynamics)"))
        self.label_17.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600;\">machine learning analyses to be run</span></p></body></html>"))
        self.label_7.setText(_translate("Dialog", "<html><head/><body><p align=\"center\"><span style=\" font-size:16pt; color:#ffffff;\">DROIDS/maxDemon 5.0 </span></p><p align=\"center\"><span style=\" font-size:16pt; color:#ffffff;\">AI-assisted Comparative Molecular Dynamics</span></p></body></html>"))
        self.textBrowser.setHtml(_translate("Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">1cdw_bound.pdb</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">1cdw_bound.prmtop</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">1cdw_bound.nc</p></body></html>"))
        self.label_20.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-weight:600;\">file list example </span></p></body></html>"))


##########################################################################################
#################   button actions    ####################################################
##########################################################################################
# create run program action  (make DROIDS.ctl and launch pipeline)
    def run_program(self):
        
        print("writing control file")
        n_subsamples = self.lineEdit.text()
        #print(n_subsamples)
        frs_subsample = self.lineEdit_2.text()
        #print(frs_subsample)
        frs_total = self.lineEdit_3.text()
        #print(frs_total)
        n_chains = self.lineEdit_4.text()
        #print(n_chains)
        n_sites = self.lineEdit_5.text()
        #print(n_sites)
        start_site = self.lineEdit_6.text()
        #print(start_site)
        chx_path = self.lineEdit_7.text()
        #print(chx_path)
        query_list = self.textEdit.toPlainText()
        #print(query_list)
        query_list = str.split(query_list, "\n")
        #print(query_list)
        query_pdb = query_list[0]
        query_id = query_pdb[:-4]
        query_top = query_list[1]
        query_traj = query_list[2]
        reference_list = self.textEdit_2.toPlainText()
        #print(reference_list)
        reference_list = str.split(reference_list, "\n")
        #print(reference_list)
        reference_pdb = reference_list[0]
        reference_id = reference_pdb[:-4]
        reference_top = reference_list[1]
        reference_traj = reference_list[2]
        
        # check and set graph style
        if self.radioButton.isChecked() == True:
            graph_color = "light"
        else:
            graph_color = "dark"
        
        # analyses to run
        if self.checkBox.isChecked() == True:
            divergence = "yes"
        else:
            divergence = "no"
        if self.checkBox_2.isChecked() == True:
            discrepancy = "yes"
        else:
            discrepancy = "no"
        if self.checkBox_3.isChecked() == True:
            conservation = "yes"
        else:
            conservation = "no"
        if self.checkBox_4.isChecked() == True:
            coordination = "yes"
        else:
            coordination = "no"
        #if self.checkBox_5.isChecked() == True:
        #    variants = "yes"
        #else:
        #    variants = "no"
        
        # write file
        f = open("./DROIDS.ctl", "w") 
        f.write("queryID,%s,#pdb id for query structure\n" % query_id)
        f.write("referenceID,%s,#pdb id for ref structure\n" % reference_id)
        f.write("queryPDB,%s,#pdb file for query structure\n" % query_pdb)
        f.write("referencePDB,%s,#pdb file for ref structure\n" % reference_pdb)
        f.write("queryTOP,%s,#topology for query structure\n" % query_top)
        f.write("referenceTOP,%s,#topology for ref structure\n" % reference_top)
        f.write("queryTRAJ,%s,#trajectory for query structure\n" % query_traj)
        f.write("referenceTRAJ,%s,#trajectory for ref structure\n" % reference_traj)
        f.write("subsamples,%s,#number of subsamples\n" % n_subsamples)
        f.write("frame_size,%s,#number of frames per subsample\n" % frs_subsample)
        f.write("n_frames,%s,#total number of frames in simulation (5000 per ns)\n" % frs_total)
        f.write("num_chains,%s,#number of protein chains\n" % n_chains)
        f.write("length,%s,#total length of protein\n" % n_sites)
        f.write("start,%s,#Nterminal AA starts at position...\n" % start_site)
        f.write("chimerax,%s,#path to chimerax binary\n" % chx_path)
        f.write("bgcolor,%s,#background color graphics\n" % graph_color)
        f.write("divergence,%s,#run KL divergence analysis\n" % divergence)
        f.write("discrepancy,%s,#run MMD analysis\n" % discrepancy)
        f.write("conservation,%s,#run conserved dynamics analysis\n" % conservation)
        f.write("coordination,%s,#run coordinated dynamics analysis\n" % coordination)
        #f.write("variants,%s,#run coordinated dynamics analysis\n" % variants)
        f.close()
        
        
        print("running DROIDS/maxDemon 5.0 parsers")
        cmd1 = "python3 cpptraj_sampler.py"
        os.system(cmd1)
        
    def run_analyses(self):
        print("running DROIDS/maxDemon 5.0 analyses")
        cmd2 = "python3 chimerax_analyzer.py"
        os.system(cmd2)
    
    def closeIt(self):
        print("DROIDS/maxDemon program closed")
        sys.exit(app.exec_())


#########################################################################################
##########################################################################################
if __name__ == "__main__":
    import sys
    import os
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())
