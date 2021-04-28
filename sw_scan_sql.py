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
# import csv
import sqlite3

class SW_Scan(QMainWindow, Ui_MainWindow):

    def __init__(self, input_file_name: str, parent=None):
        super(SW_Scan, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("SW Scan")

        self.initialize_var()

        self.connections()
        
        app.processEvents()

        if input_file_name:
            self.load_file(input_file_name)


    def initialize_var(self):
        self.file_name = ""
        self.align = None
        self.connection = None


    def connections(self):
        self.action_load_file.triggered.connect(self.load_file)
        self.pb_filter.clicked.connect(self.filter)
        self.pb_save_tsv.clicked.connect(self.save_tsv)
        self.pb_save_fbs.clicked.connect(self.save_fbs)
        self.action_quit.triggered.connect(self.close)


    def load_file(self, file_name):

        if not file_name:
            self.file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "")
            if not self.file_name:
                return
        else:
            self.file_name = file_name

        self.statusBar().showMessage("Loading file...")
        app.processEvents()

        self.connection = sqlite3.connect(':memory:')
        cur = self.connection.cursor()
        cur.execute('create table sequences (accession text, description text, frame text, identity float, score float, align_length int, target_length int, aligned_query_sequence text, aligned_target_sequence text, query_begin int, query_end int, target_begin int, target_end_optimal int)')
        self.connection.commit()

        df = pd.read_csv(self.file_name, sep='\t')
        df.to_sql('sequences', self.connection, if_exists='replace', index=False)
        cur.execute('SELECT count(*) FROM sequences')
        row = cur.fetchone()
        row_count = row[0]

        cur = self.connection.cursor()
        cur.execute("SELECT * FROM sequences ORDER by identity DESC LIMIT 1000")

        self.tw.setRowCount(0)
        self.tw.setColumnCount(len(FIELDS_NAME))
        self.tw.setHorizontalHeaderLabels(FIELDS_NAME)

        for row in cur.fetchall():
            idx = 0
            self.tw.setRowCount(self.tw.rowCount() + 1)
            for r in row:
                self.tw.setItem(self.tw.rowCount() - 1, idx, QTableWidgetItem(str(r)))
                idx += 1

        self.statusBar().showMessage(f"File {self.file_name} loaded. {row_count} sequences found ({self.tw.rowCount()} displayed)")


    def filter(self):

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

        self.tw.setRowCount(0)
        self.tw.setColumnCount(len(FIELDS_NAME))
        self.tw.setHorizontalHeaderLabels(FIELDS_NAME)

        id = self.le_id.text().upper() if self.le_id.text() else ""
        description1 = self.le_description.text().upper() if self.le_description.text() else ""
        identity_pc = float(self.le_identity.text()) if self.le_identity.text() else 0
        min_align_length = int(self.le_align_length.text()) if self.le_align_length.text() else 0

        sql = "SELECT * FROM sequences WHERE "

        sql2 = ""        
        if id:
            sql2 += f" description LIKE '%{id}%' "

        if description1:
            if sql2: sql2 += " AND "
            sql2 += f" description LIKE '%{description1}%' "

        if identity_pc:
            if sql2: sql2 += " AND "
            sql2 += f" identity {self.cb_identity_relation.currentText()} {identity_pc} "

        if min_align_length:
            if sql2: sql2 += " AND "
            sql2 += f" align_length >= {min_align_length} "

        cur = self.connection.cursor()
        cur.execute(sql + sql2 + " ORDER by identity")

        for row in cur.fetchall():
            idx = 0
            self.tw.setRowCount(self.tw.rowCount() + 1)
            for r in row:
                self.tw.setItem(self.tw.rowCount() - 1, idx, QTableWidgetItem(str(r)))
                idx += 1

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
    if len(sys.argv) > 1:
        file_name = sys.argv[1]
    else:
        file_name = ""
    program = SW_Scan(file_name)
    program.show()
    program.raise_()
    sys.exit(app.exec_())
    

