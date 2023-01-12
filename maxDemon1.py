# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'maxDemon1.ui'
#
# Created by: PyQt5 UI code generator 5.14.1
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(526, 326)
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(10, 20, 231, 21))
        self.label.setObjectName("label")
        self.textEdit = QtWidgets.QTextEdit(Dialog)
        self.textEdit.setGeometry(QtCore.QRect(10, 50, 211, 211))
        self.textEdit.setObjectName("textEdit")
        self.label_2 = QtWidgets.QLabel(Dialog)
        self.label_2.setGeometry(QtCore.QRect(300, 20, 201, 21))
        self.label_2.setObjectName("label_2")
        self.textBrowser = QtWidgets.QTextBrowser(Dialog)
        self.textBrowser.setGeometry(QtCore.QRect(250, 50, 251, 111))
        self.textBrowser.setObjectName("textBrowser")
        self.pushButton = QtWidgets.QPushButton(Dialog)
        self.pushButton.setGeometry(QtCore.QRect(300, 170, 141, 41))
        self.pushButton.setObjectName("pushButton")
        self.pushButton.clicked.connect(self.run_sampler)
        self.pushButton_2 = QtWidgets.QPushButton(Dialog)
        self.pushButton_2.setGeometry(QtCore.QRect(300, 220, 141, 41))
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.clicked.connect(self.analyzeIt)
        self.pushButton_3 = QtWidgets.QPushButton(Dialog)
        self.pushButton_3.setGeometry(QtCore.QRect(300, 270, 141, 41))
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_3.clicked.connect(self.closeIt)
        
        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Enter 5 files for ortholog/genetic variant comparison to query"))
        self.label.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600;\">ortholog/mutant files</span></p></body></html>"))
        self.label_2.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600;\">file list example</span></p></body></html>"))
        self.textBrowser.setHtml(_translate("Dialog", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Ubuntu\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">1ubq_ortholog.pdb</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">1ubq_ortholog.prmtop</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">1ubq_ortholog.nc</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">1ubq.fasta</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:12pt;\">1ubq_ortholog.fasta</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.pushButton.setText(_translate("Dialog", "run sampling"))
        self.pushButton_2.setText(_translate("Dialog", "analyze"))
        self.pushButton_3.setText(_translate("Dialog", "close afterwards"))

#########################################################################################
       
    
    def run_sampler(self):
        print("writing control file")
        ortho_list = self.textEdit.toPlainText()
        print(ortho_list)
        ortho_list = str.split(ortho_list, "\n")
        print(ortho_list)
        ortho_pdb = ortho_list[0]
        ortho_id = ortho_pdb[:-4]
        ortho_top = ortho_list[1]
        ortho_traj = ortho_list[2]
        query_fas = ortho_list[3]
        ortho_fas = ortho_list[4]
                
        # write file
        f = open("./maxDemon.ctl", "w") 
        f.write("orthoID,%s,#pdb id for ortholog structure\n" % ortho_id)
        f.write("orthoPDB,%s,#pdb file for ortholog structure\n" % ortho_pdb)
        f.write("orthoTOP,%s,#topology for ortholog structure\n" % ortho_top)
        f.write("orthoTRAJ,%s,#trajectory for ortholog structure\n" % ortho_traj)
        f.write("queryFAS,%s,#alignment sequence for query structure\n" % query_fas)
        f.write("orthoFAS,%s,#alignment sequence for ortholog structure\n" % ortho_fas)
        f.close()
        
        ############# run ortholog sampling program ###################
        print("running DROIDS/maxDemon 5.0 ortholog/mutant trajectory sampler")
        cmd1 = "python3 cpptraj_ortholog_sampler.py"
        os.system(cmd1)
        print("close user interface when completed")
    
    def analyzeIt(self):
        print("running DROIDS/maxDemon 5.0 ortholog/mutant trajectory sampler")
        cmd1 = "python3 chimerax_consdyn.py"
        os.system(cmd1)
    
        
    def closeIt(self):
        print("maxDemon sampler program closed")
        sys.exit(app.exec_())



if __name__ == "__main__":
    import sys
    import os
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())
