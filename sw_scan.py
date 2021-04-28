#!/usr/bin/env python3
"""
Display and filter results of sw.py
(c) Olivier Friard 2021

"""

FIELDS_NAME = ("accession", "description", "frame", "identity",
"score", "align_length", "target_length", "aligned_query_sequence",
"aligned_target_sequence",
"query_begin", "query_end", "target_begin", "target_end_optimal")

SEQ_LIMIT_NB = 1000
SEQ_ORDER = " ORDER BY identity DESC "

from PyQt5.QtWidgets import (QMainWindow, QApplication,
                             QTableWidgetItem,
                             QFileDialog, QMessageBox)

from sw_scan_ui import Ui_MainWindow
import sys
import pandas as pd
import sqlite3

__version__ = '2'
__version_date__ = "2021-04-28"

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
        self.n_seq = 0


    def connections(self):
        self.action_load_file.triggered.connect(lambda: self.load_file(""))
        self.pb_filter.clicked.connect(self.filter)
        self.pb_clear.clicked.connect(self.clear)
        self.pb_run_query.clicked.connect(self.run_query)
        self.pb_save_tsv.clicked.connect(self.save_tsv)
        self.pb_save_fbs.clicked.connect(self.save_fbs)
        self.actionAbout.triggered.connect(self.about)
        self.action_quit.triggered.connect(self.close)

    def about(self):

        about_dialog = QMessageBox()
        #about_dialog.setIconPixmap(QPixmap(":/small_logo"))

        about_dialog.setWindowTitle("About SW Scan")
        about_dialog.setStandardButtons(QMessageBox.Ok)
        about_dialog.setDefaultButton(QMessageBox.Ok)
        about_dialog.setEscapeButton(QMessageBox.Ok)

        about_dialog.setInformativeText((
            f"<b>SW Scan</b> v. {__version__} - {__version_date__}"
            "<p>Copyright &copy; 2021 Olivier Friard<br>"
            "Department of Life Sciences and Systems Biology<br>"
            "University of Torino - Italy<br>"
        ))

        about_dialog.setDetailedText("")

        _ = about_dialog.exec_()



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
        self.n_seq = row[0]

        cur = self.connection.cursor()
        sql = f"SELECT * FROM sequences {SEQ_ORDER} LIMIT {SEQ_LIMIT_NB}"
        self.pte_sql.setPlainText(sql)
        cur.execute(sql)

        self.tw.setRowCount(0)
        self.tw.setColumnCount(len(FIELDS_NAME))
        self.tw.setHorizontalHeaderLabels(FIELDS_NAME)

        for row in cur.fetchall():
            idx = 0
            self.tw.setRowCount(self.tw.rowCount() + 1)
            for r in row:
                self.tw.setItem(self.tw.rowCount() - 1, idx, QTableWidgetItem(str(r)))
                idx += 1

        self.statusBar().showMessage(f"File {self.file_name} loaded. {self.n_seq} sequences found ({self.tw.rowCount()} displayed)")


    def clear(self):
        for w in (self.le_id, self.le_description1, self.le_description2,
                  self.le_identity, self.le_align_length):
            w.clear()
        self.filter()


    def filter(self):

        self.tw.clear()
        self.statusBar().showMessage("Filtering...")
        app.processEvents()
        if self.file_name == "":
            return

        id = self.le_id.text().upper() if self.le_id.text() else ""
        description1 = self.le_description1.text().upper() if self.le_description1.text() else ""
        description2 = self.le_description2.text().upper() if self.le_description2.text() else ""
        identity_pc = float(self.le_identity.text()) if self.le_identity.text() else 0
        min_align_length = int(self.le_align_length.text()) if self.le_align_length.text() else 0

        sql1 = "SELECT * FROM sequences WHERE "

        sql2 = ""        
        if id:
            sql2 += f" description LIKE '%{id}%' "

        if description1:
            for term in description1.split(","):
                if sql2: sql2 += " AND "
                sql2 += f" description LIKE '%{term}%' "

        if description2:
            for term in description2.split(","):
                if sql2: sql2 += " AND "
                sql2 += f" description NOT LIKE '%{term}%' "

        if identity_pc:
            if sql2: sql2 += " AND "
            sql2 += f" identity {self.cb_identity_relation.currentText()} {identity_pc} "

        if min_align_length:
            if sql2: sql2 += " AND "
            sql2 += f" align_length >= {min_align_length} "

        if sql2:
            sql = sql1 + sql2 + SEQ_ORDER
            flag_all = False
        else:
            sql = f"SELECT * FROM sequences {SEQ_ORDER} LIMIT {SEQ_LIMIT_NB}" 
            flag_all = True

        self.pte_sql.setPlainText(sql)

        self.tw.setRowCount(0)
        self.tw.setColumnCount(len(FIELDS_NAME))
        self.tw.setHorizontalHeaderLabels(FIELDS_NAME)

        cur = self.connection.cursor()
        cur.execute(sql)
        for row in cur.fetchall():
            idx = 0
            self.tw.setRowCount(self.tw.rowCount() + 1)
            for r in row:
                self.tw.setItem(self.tw.rowCount() - 1, idx, QTableWidgetItem(str(r)))
                idx += 1

        if flag_all:
            self.statusBar().showMessage(f"{self.n_seq} sequences found ({self.tw.rowCount()} displayed)")
        else:
            self.statusBar().showMessage(f"{self.tw.rowCount()} sequences found")


    def run_query(self):
        if self.pte_sql.toPlainText():

            self.tw.setRowCount(0)
            self.tw.setColumnCount(len(FIELDS_NAME))
            self.tw.setHorizontalHeaderLabels(FIELDS_NAME)

            cur = self.connection.cursor()
            try:
                cur.execute(self.pte_sql.toPlainText())
            except sqlite3.OperationalError:
                self.statusBar().showMessage("Query error!")
                return

            for row in cur.fetchall():
                idx = 0
                self.tw.setRowCount(self.tw.rowCount() + 1)
                for r in row:
                    self.tw.setItem(self.tw.rowCount() - 1, idx, QTableWidgetItem(str(r)))
                    idx += 1

            self.statusBar().showMessage(f"{self.tw.rowCount()} sequences found")


    def save_tsv(self):
        file_name, _ = QFileDialog.getSaveFileName(self, "Save file TSV", "")
        if not file_name:
            return
        with open(file_name, "w") as f_out:
            for row in range(self.tw.rowCount()):
                for col, _ in enumerate(FIELDS_NAME):
                    print(self.tw.item(row, col).text(), file=f_out, end="\t")
                print(file=f_out, end="\n")


    def save_fbs(self):

        # all uniq aligned query
        file_name, _ = QFileDialog.getSaveFileName(self, "Save file FBS", "")
        if not file_name:
            return

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

        with open(file_name, "w") as f_out:
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
    

