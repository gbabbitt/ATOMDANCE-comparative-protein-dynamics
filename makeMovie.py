# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'makeMovie.ui'
#
# Created by: PyQt5 UI code generator 5.14.1
#
# WARNING! All changes made in this file will be lost!

import numpy as np
import random as rd
from sklearn.preprocessing import MinMaxScaler
from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(761, 287)
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(110, 20, 141, 21))
        self.label.setObjectName("label")
        self.textEdit = QtWidgets.QTextEdit(Dialog)
        self.textEdit.setGeometry(QtCore.QRect(10, 60, 361, 161))
        self.textEdit.setObjectName("textEdit")
        self.label_2 = QtWidgets.QLabel(Dialog)
        self.label_2.setGeometry(QtCore.QRect(460, 20, 201, 21))
        self.label_2.setObjectName("label_2")
        self.textBrowser = QtWidgets.QTextBrowser(Dialog)
        self.textBrowser.setGeometry(QtCore.QRect(390, 60, 351, 161))
        self.textBrowser.setObjectName("textBrowser")
        self.pushButton = QtWidgets.QPushButton(Dialog)
        self.pushButton.setGeometry(QtCore.QRect(430, 230, 251, 41))
        self.pushButton.setObjectName("pushButton")
        self.pushButton.clicked.connect(self.morphIt)
        self.pushButton_2 = QtWidgets.QPushButton(Dialog)
        self.pushButton_2.setGeometry(QtCore.QRect(690, 230, 61, 41))
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.clicked.connect(self.closeIt)
        self.pushButton_3 = QtWidgets.QPushButton(Dialog)
        self.pushButton_3.setGeometry(QtCore.QRect(10, 230, 401, 41))
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_3.clicked.connect(self.convertIt)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "generate MMD weighted traj files"))
        self.label.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600;\">input file list</span></p></body></html>"))
        self.label_2.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600;\">file list example</span></p></body></html>"))
        self.textBrowser.setHtml(_translate("Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">1ubq.pdb</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">1ubq.prmtop</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">1ubq.nc</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">maxMeanDiscrepancy_flux.txt</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.pushButton.setText(_translate("Dialog", "rewrite pdb weighted by MMD"))
        self.pushButton_2.setText(_translate("Dialog", "exit"))
        self.pushButton_3.setText(_translate("Dialog", "create multi-frame pdb from .nc trajectory file"))

#####################################################################

    def closeIt(self):
        print("makeMovie program closed")
        sys.exit(app.exec_())

    def convertIt(self):
        print("reading input")
        reference_list = self.textEdit.toPlainText()
        print(reference_list)
        reference_list = str.split(reference_list, "\n")
        print(reference_list)
        pdbFILE = reference_list[0]
        pdbID = pdbFILE[:-4]
        topFILE = reference_list[1]
        trajFILE = reference_list[2]
        mmdFILE = reference_list[3]
        
        user_input = input("\nIs AmberToolsXX or cpptraj installed in a conda environment? (yes/no)\n\n")
        
        if(user_input == "yes" or user_input == "y"):
            #what_version = input("\n Enter AmberTools version number (e.g. 22)")
            #os.system("conda activate AmberTools%s\n" % what_version)
            print("\nRUN THE FOLLOWING CMD'S/SCRIPTS IN THE NEW TERMINAL\n")
            print("conda activate OR conda activate AmberTools22 (as set on your system)")
            print("python3 MD_protein_pdb4amber.py")
            print("python3 MD_protein_antechamber.py (if needed)")
            print("python3 MD_protein_tleap.py")
            print("conda deactivate")
            print("\n\n")
            #os.system("conda deactivate\n")
            print("CLOSE TERMINAL WHEN MD SIMULATION PREPARATIONS ARE COMPLETED\n\n")
            os.system("x-terminal-emulator\n")
        
        if(user_input != "yes" or user_input != "y"):
            print("converting .nc file to multiframe .pdb file")
            cmd1 = "cpptraj -p %s -y %s -x mfMMD_%s" % (topFILE, trajFILE, pdbFILE)
            os.system(cmd1)
    
    
    def morphIt(self):
        print("reading input")
        reference_list = self.textEdit.toPlainText()
        print(reference_list)
        reference_list = str.split(reference_list, "\n")
        print(reference_list)
        pdbFILE = reference_list[0]
        pdbID = pdbFILE[:-4]
        topFILE = reference_list[1]
        trajFILE = reference_list[2]
        mmdFILE = reference_list[3]
        
        print("deleting waters in multiframe .pdb")
        outfile = open("./mfMMD_noWAT_%s" % pdbFILE, "w") 
        infile = open("./mfMMD_%s" % pdbFILE, "r") 
        infile_lines = infile.readlines()
        skip = "no"
        for x in range(len(infile_lines)):
            if(x==0):
                continue
            skip = "no"
            infile_line = infile_lines[x]
            infile_line_array = str.split(infile_line)
            #print(infile_line_array)
            header = infile_line_array[0]
            #print(header)
            if(header == "MODEL"):
                modelN = infile_line_array[1]
                print("removing waters from model %s" % modelN)
            if(header == "ATOM"):
                residue = infile_line_array[3]
                #print(residue)
                if(residue == "WAT" or residue == "Na+" or residue == "Cl-"):
                    skip = "yes"
            if(header == "TER"):
                residue = infile_line_array[2]
                #print(residue)
                if(residue == "WAT" or residue == "Na+" or residue == "Cl-"):
                    skip = "yes"
            if(skip == "no"):
                outfile.write(infile_line)
        outfile.close
        infile.close
        
        print("weighting MMD to multiframe .pdb")
        outfile = open("./mfMMD_adjust_%s" % pdbFILE, "w") 
        infile = open("./mfMMD_noWAT_%s" % pdbFILE, "r")
        infile_lines = infile.readlines()
        infile2 = open("./%s" % mmdFILE, "r")
        infile2_lines = infile2.readlines()
        for x in range(len(infile_lines)):
            infile_line = infile_lines[x]
            infile_line_array = str.split(infile_line)
            #print(infile_line_array)
            header = infile_line_array[0]
            #print(header)
            if(header == "MODEL"):
                modelN = infile_line_array[1]
                print("adjusting trajectories according to MMD in model %s" % modelN)
                
            if(header == "ATOM"):
                atomN1 = infile_line_array[1]
                atomType1 = infile_line_array[2]
                res1 = infile_line_array[3]
                #print(residue)
                pos1 = infile_line_array[4]
                mmd="NA"
                mmd_wt=1.000 # does not alter positions in DNA or ligand etc
                myMMD = []
                myPOS = []
                myRES = []
                for x in range(len(infile2_lines)): # normalize MMD from 0 to 1
                    if(x==0):
                        continue
                    infile2_line = infile2_lines[x]
                    infile2_line_array = str.split(infile2_line)
                    #print(infile2_line_array)
                    pos2 = infile2_line_array[0]
                    res2 = infile2_line_array[1]
                    mmd2 = (infile2_line_array[2])
                    #print(mmd2)
                    mmd2 = float(mmd2)
                    myPOS.append(pos2)
                    myRES.append(res2)
                    myMMD.append(mmd2)
                myMMD = np.array(myMMD)
                myMMD = myMMD.reshape(-1, 1) # reshape because it is 1D feature
                scaler = MinMaxScaler()
                scaler.fit(myMMD)
                myMMD_norm = scaler.transform(myMMD)
                myMMD_norm = myMMD_norm*2 # range from 0,2
                #print(myMMD) # untransformed
                #print(myMMD_norm) # transformed
                for x in range(len(myMMD_norm)):     
                    pos_test = myPOS[x]
                    #print(pos_test)
                    #print(pos1)
                    if(pos1==pos_test):
                        mmd_wt=myMMD_norm[x]
                        mmd_wt = float(mmd_wt)
                        #print("yes")
                        #print(mmd_wt)
                atomX = infile_line_array[5]
                atomY = infile_line_array[6]
                atomZ = infile_line_array[7]
                atomX = float(atomX)
                atomY = float(atomY)
                atomZ = float(atomZ)
                #print("values XYZ %s %s %s and mmd %s" % (atomX, atomY, atomZ, mmd_wt))
                ###### recalculate XYX #######
                noise = rd.gauss(0, 0.125)
                atomX_adj = atomX+(noise*mmd_wt**4)
                atomY_adj = atomY+(noise*mmd_wt**4)
                atomZ_adj = atomZ+(noise*mmd_wt**4)
                atomX_adj = round(atomX_adj, 3)
                atomY_adj = round(atomY_adj, 3)
                atomZ_adj = round(atomZ_adj, 3)
                #print("adjust XYZ %s %s %s and mmd %s" % (atomX_adj, atomY_adj, atomZ_adj, mmd_wt))
                atomX_adj = '%.3f' % atomX_adj # convert to fixed length string
                atomY_adj = '%.3f' % atomY_adj
                atomZ_adj = '%.3f' % atomZ_adj
                infile_line_segment1 = infile_line[:32]
                infile_line_segment2 = infile_line[32:55]
                infile_line_segment3 = infile_line[55:79]
                #print(infile_line_segment1)
                #print(infile_line_segment2)
                #print(infile_line_segment3)
                infile_line_segment2_adj = "%s  %s  %s " %(atomX_adj, atomY_adj,atomZ_adj )
                adj_infile_line = "%s%s%s\n" % (infile_line_segment1,infile_line_segment2_adj,infile_line_segment3)
                ##############################
            if(header == "ATOM"):
                outfile.write(adj_infile_line)
            if(header != "ATOM"):
                outfile.write(infile_line)    
        outfile.close
        infile.close
        
        print("to create movie from multiframe .pdb")
        print("close this program then...")
        print("open Tools/General/MolecularDynamicsViewer and load multi-frame .pdb and attribute file for MMD")
        #cmd2 = "chimerax"
        #os.system(cmd2)
#####################################################################

if __name__ == "__main__":
    import sys
    import os
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())
