#!/usr/bin/env python3
"""

"""

FIELDS_NAME = ("accession", "description", "frame", "identity",
"score", "align_length", "target_length", "aligned_query_sequence",
"aligned_target_sequence",
"query_begin", "query_end", "target_begin", "target_end_optimal")


from PyQt5.QtWidgets import (QMainWindow, QApplication,
                             QTableWidgetItem,
                             QFileDialog)


from sw_scan_ui import Ui_MainWindow
import sys
import pandas as pd
import csv

class SW_Scan(QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None):
        super(SW_Scan, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("SW Scan")

        self.initialize_var()

        self.connections()


    def initialize_var(self):
        self.file_name = ""
        self.align = None

    def connections(self):

        self.action_load_file.triggered.connect(self.load_file)
        self.pb_filter.clicked.connect(self.filter2)
        self.pb_save_tsv.clicked.connect(self.save_tsv)
        self.pb_save_fbs.clicked.connect(self.save_fbs)
        self.action_quit.triggered.connect(self.close)


    def load_file(self):

        #self.align = pd.read_table("/data/tmp/orf8_vrl.sorted.sw.1e5")
        #return

        self.file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "")
        if not self.file_name:
            return

        with open(self.file_name, "r") as f:
            reader = csv.reader(f, delimiter = "\t")
            row_count = len(list(reader))

        #self.align = pd.read_table(self.file_name)
        self.statusBar().showMessage(f"File {self.file_name} loaded. {row_count} sequences found")


    def filter(self):

        self.tw.clear()

        results = None
        if self.le_id.text():
            results = self.align[self.align.accession.str.contains(self.le_id.text())]
        if self.le_description.text():
            if results is not None:
                results = results[self.align.description.str.contains(self.le_id.text())]
            else:
                results = self.align[self.align.description.str.contains(self.le_id.text())]

        self.tw.setRowCount(len(results))
        self.tw.setColumnCount(len(FIELDS_NAME))
        self.tw.setHorizontalHeaderLabels(FIELDS_NAME)

        for idx in range(len(results)):
            for i, name in enumerate(FIELDS_NAME):
                self.tw.setItem(idx, i, QTableWidgetItem(str(results.iloc[int(idx)][name])))


    def filter2(self):

        def identity(value, relation, data):
            if relation == '=':
                return value == data
            if relation == '>=':
                return float(data) >= value
            if relation == '<=':
                return float(data) <= value

        self.tw.clear()
        self.statusBar().showMessage("Filtering...")
        app.processEvents()
        if self.file_name == "":
            return
        infile = open(self.file_name, "r")
        read = csv.reader(infile, delimiter='\t')
        headers = next(read)

        self.tw.setRowCount(0)
        self.tw.setColumnCount(len(FIELDS_NAME))
        self.tw.setHorizontalHeaderLabels(FIELDS_NAME)

        id = self.le_id.text().upper() if self.le_id.text() else ""
        description = self.le_description.text().upper() if self.le_description.text() else ""
        identity_pc = float(self.le_identity.text()) if self.le_identity.text() else 0
        min_align_length = int(self.le_align_length.text()) if self.le_align_length.text() else 0

        for row in read:

            if (id in row[0].upper()) and \
               (description in row[1].upper()) and \
               ((identity_pc == 0) or (identity_pc and identity(identity_pc, self.cb_identity_relation.currentText(), row[3]))) and \
               ((min_align_length == 0) or (min_align_length and int(row[5]) >= min_align_length)):

                self.tw.setRowCount(self.tw.rowCount() + 1)
                for i, _ in enumerate(FIELDS_NAME):
                    self.tw.setItem(self.tw.rowCount() - 1, i, QTableWidgetItem(row[i]))

        self.statusBar().showMessage(f"{self.tw.rowCount()} sequences found")

    def save_tsv(self):
        with open("results.tsv", "w") as f_out:
            for row in range(self.tw.rowCount()):
                for col, _ in enumerate(FIELDS_NAME):
                    print(self.tw.item(row, col).text(), file=f_out, end="\t")
                print(file=f_out, end="\n")


    def save_fbs(self):

        # all uniq aligned query
        '''
        aligned_query_list = []
        for row in range(self.tw.rowCount()):
            for col in [7]:
                if self.tw.item(row, col).text() not in aligned_query_list:
                   aligned_query_list.append(self.tw.item(row, col).text()) 
        '''

        aligned_query_list = []
        for row in range(self.tw.rowCount()):
            if aligned_query_list == []:
                aligned_query_list.append(self.tw.item(row, 7).text())
            else:
                flag_in = False
                for x in aligned_query_list:
                    if self.tw.item(row, 7).text() in x:
                        flag_in = True
                        break
                if not flag_in:
                    aligned_query_list.append(self.tw.item(row, 7).text()) 

        print(aligned_query_list)

        with open("results_fbs.tsv", "w") as f_out:
            out = ""
            for aligned_query in aligned_query_list:
                out += "query" + (" " *  25) + aligned_query + "\n"
                for row in range(self.tw.rowCount()):

                    # check aligned query
                    if self.tw.item(row, 7).text() not in aligned_query:
                        continue

                    for col, _ in enumerate(FIELDS_NAME):

                        if col == 0:
                            id = self.tw.item(row, col).text()
                            while len(id) < 30:
                                id += ' '
                            out += id

                        if col == 8: # align target sequence
                            ats = ""
                            aligned_target_sequence = self.tw.item(row, col).text()
                            aligned_target_sequence = (" " * aligned_query.index(self.tw.item(row, 7).text())) + aligned_target_sequence
                            for nq, nt in zip(aligned_query, aligned_target_sequence):
                                if nq == nt:
                                    ats += "."
                                elif nt == " ":
                                    ats += " "
                                else:
                                    ats += nt

                            out += ats + "\n"
                out += "\n\n"

            f_out.write(out) 



if __name__ == "__main__":
    app = QApplication(sys.argv)
    program = SW_Scan()
    program.show()
    program.raise_()
    sys.exit(app.exec_())
    

