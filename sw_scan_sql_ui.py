# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'sw_scan_sql.ui'
#
# Created by: PyQt5 UI code generator 5.15.4
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1518, 603)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.splitter = QtWidgets.QSplitter(self.centralwidget)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.frame = QtWidgets.QFrame(self.splitter)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.frame)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.pb_save_tsv = QtWidgets.QPushButton(self.frame)
        self.pb_save_tsv.setObjectName("pb_save_tsv")
        self.gridLayout.addWidget(self.pb_save_tsv, 10, 0, 1, 1)
        self.label = QtWidgets.QLabel(self.frame)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.le_description2 = QtWidgets.QLineEdit(self.frame)
        self.le_description2.setObjectName("le_description2")
        self.gridLayout.addWidget(self.le_description2, 3, 1, 1, 2)
        self.pte_sql = QtWidgets.QPlainTextEdit(self.frame)
        self.pte_sql.setObjectName("pte_sql")
        self.gridLayout.addWidget(self.pte_sql, 8, 0, 1, 3)
        self.pb_filter = QtWidgets.QPushButton(self.frame)
        self.pb_filter.setObjectName("pb_filter")
        self.gridLayout.addWidget(self.pb_filter, 6, 0, 1, 1)
        self.pb_save_fbs = QtWidgets.QPushButton(self.frame)
        self.pb_save_fbs.setObjectName("pb_save_fbs")
        self.gridLayout.addWidget(self.pb_save_fbs, 10, 1, 1, 1)
        self.le_align_length = QtWidgets.QLineEdit(self.frame)
        self.le_align_length.setObjectName("le_align_length")
        self.gridLayout.addWidget(self.le_align_length, 5, 1, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.frame)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 3, 0, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.frame)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 5, 0, 1, 1)
        self.cb_identity_relation = QtWidgets.QComboBox(self.frame)
        self.cb_identity_relation.setObjectName("cb_identity_relation")
        self.cb_identity_relation.addItem("")
        self.cb_identity_relation.addItem("")
        self.cb_identity_relation.addItem("")
        self.gridLayout.addWidget(self.cb_identity_relation, 4, 1, 1, 1)
        self.le_id = QtWidgets.QLineEdit(self.frame)
        self.le_id.setObjectName("le_id")
        self.gridLayout.addWidget(self.le_id, 0, 1, 1, 2)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem, 6, 1, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.frame)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.le_description1 = QtWidgets.QLineEdit(self.frame)
        self.le_description1.setObjectName("le_description1")
        self.gridLayout.addWidget(self.le_description1, 1, 1, 1, 2)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem1, 7, 0, 1, 1)
        self.pb_run_query = QtWidgets.QPushButton(self.frame)
        self.pb_run_query.setObjectName("pb_run_query")
        self.gridLayout.addWidget(self.pb_run_query, 9, 0, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.frame)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 4, 0, 1, 1)
        self.le_identity = QtWidgets.QLineEdit(self.frame)
        self.le_identity.setObjectName("le_identity")
        self.gridLayout.addWidget(self.le_identity, 4, 2, 1, 1)
        self.pb_clear = QtWidgets.QPushButton(self.frame)
        self.pb_clear.setObjectName("pb_clear")
        self.gridLayout.addWidget(self.pb_clear, 6, 2, 1, 1)
        self.pb_save_fasta = QtWidgets.QPushButton(self.frame)
        self.pb_save_fasta.setObjectName("pb_save_fasta")
        self.gridLayout.addWidget(self.pb_save_fasta, 10, 2, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.tw = QtWidgets.QTableView(self.splitter)
        self.tw.setObjectName("tw")
        self.horizontalLayout.addWidget(self.splitter)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1518, 22))
        self.menubar.setObjectName("menubar")
        self.menuQuit = QtWidgets.QMenu(self.menubar)
        self.menuQuit.setObjectName("menuQuit")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.action_quit = QtWidgets.QAction(MainWindow)
        self.action_quit.setObjectName("action_quit")
        self.action_load_file = QtWidgets.QAction(MainWindow)
        self.action_load_file.setObjectName("action_load_file")
        self.actionAbout = QtWidgets.QAction(MainWindow)
        self.actionAbout.setObjectName("actionAbout")
        self.menuQuit.addAction(self.action_load_file)
        self.menuQuit.addAction(self.action_quit)
        self.menuHelp.addAction(self.actionAbout)
        self.menubar.addAction(self.menuQuit.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.pb_save_tsv.setText(_translate("MainWindow", "Save TSV"))
        self.label.setText(_translate("MainWindow", "ID contains"))
        self.pb_filter.setText(_translate("MainWindow", "Filter"))
        self.pb_save_fbs.setText(_translate("MainWindow", "Save \"FBS\""))
        self.label_5.setText(_translate("MainWindow", "Description does not contain"))
        self.label_4.setText(_translate("MainWindow", "Alignment length"))
        self.cb_identity_relation.setItemText(0, _translate("MainWindow", "="))
        self.cb_identity_relation.setItemText(1, _translate("MainWindow", ">="))
        self.cb_identity_relation.setItemText(2, _translate("MainWindow", "<="))
        self.label_2.setText(_translate("MainWindow", "Description contains"))
        self.pb_run_query.setText(_translate("MainWindow", "Run query"))
        self.label_3.setText(_translate("MainWindow", "identity"))
        self.pb_clear.setText(_translate("MainWindow", "Clear"))
        self.pb_save_fasta.setText(_translate("MainWindow", "Save in FASTA format"))
        self.menuQuit.setTitle(_translate("MainWindow", "SW Scan"))
        self.menuHelp.setTitle(_translate("MainWindow", "Help"))
        self.action_quit.setText(_translate("MainWindow", "Quit"))
        self.action_load_file.setText(_translate("MainWindow", "Load file"))
        self.actionAbout.setText(_translate("MainWindow", "About"))