#!/usr/bin/env python3
"""

"""


from PyQt5.QtWidgets import (QMainWindow, QApplication, QLabel, QWidget,
                             QTableWidgetItem,
                             QFileDialog)
# from PyQt5.QtGui import (QFileDialog)


from sw_scan_ui import Ui_MainWindow
import sys
import pandas as pd

class SW_Scan(QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None):
        super(SW_Scan, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("SW Scan")

        self.initialize_var()

        self.connections()

        self.tw.setRowCount(3)
        self.tw.setColumnCount(3)
 
 
        self.tw.setItem(0, 0, QTableWidgetItem("Name"))
        self.tw.setItem(0, 1, QTableWidgetItem("Email"))
        self.tw.setItem(0, 2, QTableWidgetItem("Phone"))
 
        #adding items
 
        self.tw.setItem(1,0, QTableWidgetItem("Parwiz"))
        self.tw.setItem(1, 1, QTableWidgetItem("parwiz@gmail.com"))
        self.tw.setItem(1, 2, QTableWidgetItem("565656"))
 
        self.tw.setItem(2, 0, QTableWidgetItem("John"))
        self.tw.setItem(2, 1, QTableWidgetItem("john@gmail.com"))
        self.tw.setItem(2, 2, QTableWidgetItem("3232323"))



    def initialize_var(self):
        self.file_name = ""
        self.align = None

    def connections(self):

        self.action_load_file.triggered.connect(self.load_file)
        self.pb_filter.clicked.connect(self.filter)
        self.action_quit.triggered.connect(self.close)


    def load_file(self):
        self.file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "")
        if not self.file_name:
            return
        self.align = pd.read_table(self.file_name)
        self.statusBar().showMessage(f"File {self.file_name} loaded")
        print(self.align['accession'])


    def filter(self):

        FIELDS_NAME = ["accession", "description", "frame", "identity", "score"]

        if self.le_id.text():
            results = self.align[self.align.accession.str.contains(self.le_id.text())]
        if self.le_description.text():
            results = results[self.align.accession.str.contains(self.le_id.text())]

        self.tw.setRowCount(len(results))
        self.tw.setColumnCount(len(FIELDS_NAME))

        for idx in range(len(results)):
            for i, name in enumerate(FIELDS_NAME):
                self.tw.setItem(idx, i, QTableWidgetItem(results.iloc[idx][name]))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    program = SW_Scan()
    program.show()
    program.raise_()
    sys.exit(app.exec_())
    

