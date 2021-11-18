#!/usr/bin/env python3
"""
Display and filter results of sw (SQLite version)
(c) Olivier Friard 2021

"""

FIELDS_NAME = ("accession", "description", "frame", "identity",
"score", "align_length", "target_length", "aligned_query_sequence",
"aligned_target_sequence",
"query_begin", "query_end", "target_begin", "target_end_optimal")


SEQ_ORDER = " ORDER BY identity DESC "

from PyQt5.QtWidgets import (QMainWindow, QApplication,
                             QTableWidgetItem,
                             QFileDialog, QMessageBox)
from PyQt5 import QtSql

from PyQt5.QtCore import QObject, QThread, pyqtSignal

from PyQt5 import QtCore

from sw_scan_ui import Ui_MainWindow
import sys

__version__ = '4'
__version_date__ = "2021-11-17"

'''
class Worker(QObject):
    finished = pyqtSignal()

    def run(self):
        print("running")
        q = QtSql.QSqlQuery(self.query)
        if q.exec_():
            q.first()
            print(q.value('n'))
            #self.statusBar().showMessage(f"{q.value('n'):,} sequence(s) filtered")

        self.finished.emit()
'''

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

        self.sql2 = ""


    def connections(self):
        self.action_load_file.triggered.connect(lambda: self.load_file(""))
        self.pb_filter.clicked.connect(self.filter)
        self.pb_clear.clicked.connect(self.clear)
        self.pb_run_query.clicked.connect(self.run_query)
        self.pb_save_fasta.clicked.connect(self.save_fasta)
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


        self.db = QtSql.QSqlDatabase.addDatabase("QSQLITE")
        self.db.setDatabaseName(self.file_name)
        if not self.db.open():
            sys.exit(-1)

        self.model = QtSql.QSqlTableModel()
        self.tw.setModel(self.model)
        self.model.setTable("sequences")
        self.model.select()


        q = QtSql.QSqlQuery("SELECT count(*) as n FROM sequences")
        if q.exec_():
            q.first()
            self.statusBar().showMessage(f"SW results loaded: {q.value('n'):,} sequence(s)")



    def clear(self):
        for w in (self.le_id, self.le_description1, self.le_description2,
                  self.le_identity, self.le_align_length):
            w.clear()

        self.filter()


    def filter(self):

        self.statusBar().showMessage("")
        if self.file_name == "":
            return

        self.statusBar().showMessage(f"Filtering sequences")
        app.processEvents()

        id = self.le_id.text().upper() if self.le_id.text() else ""
        description1 = self.le_description1.text().upper() if self.le_description1.text() else ""
        description2 = self.le_description2.text().upper() if self.le_description2.text() else ""
        identity_pc = float(self.le_identity.text()) if self.le_identity.text() else 0
        min_align_length = int(self.le_align_length.text()) if self.le_align_length.text() else 0

        #sql1 = "SELECT * FROM sequences WHERE "

        self.sql2 = ""
        if id:
            self.sql2 += f" `id` LIKE '%{id}%' "

        if description1:
            for term in description1.split(","):
                if self.sql2: self.sql2 += " AND "
                self.sql2 += f" description LIKE '%{term}%' "

        if description2:
            for term in description2.split(","):
                if self.sql2: self.sql2 += " AND "
                self.sql2 += f" description NOT LIKE '%{term}%' "

        if identity_pc:
            if self.sql2: self.sql2 += " AND "
            self.sql2 += f" identity {self.cb_identity_relation.currentText()} {identity_pc} "

        if min_align_length:
            if self.sql2: self.sql2 += " AND "
            self.sql2 += f" align_length >= {min_align_length} "

        if self.sql2:
            self.sql = self.sql2 + SEQ_ORDER
        else:
            self.sql = ""

        self.pte_sql.setPlainText(self.sql)

        self.model.setFilter(self.sql)
        self.model.select()

        q = QtSql.QSqlQuery(f"SELECT count(*) as n FROM sequences WHERE {self.sql2}")
        if q.exec_():
            q.first()
            self.statusBar().showMessage(f"{q.value('n'):,} sequence(s) filtered")

        '''
        self.thread = QThread()
        self.worker = Worker()
        self.worker.db = self.db
        self.worker.query = f"SELECT count(*) as n FROM sequences WHERE {sql2}"
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)

        #self.worker.finished.connect(self.thread.quit)

        self.worker.finished.connect(lambda: self.statusBar().showMessage("OK"))
        self.thread.start()
        '''



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


    def save_fasta(self):
        """
        save sequences from results in FASTA format
        """
        file_name, _ = QFileDialog.getSaveFileName(self, "Save file in FASTA format", "")
        if not file_name:
            return
        with open(file_name, "w") as f_out:

            if self.sql2:
                conditions = f"WHERE {self.sql2}"
            else:
                conditions = ""
            q = QtSql.QSqlQuery(f"SELECT id, frame, aligned_target_sequence FROM sequences {conditions} ORDER BY id, frame")
            if q.exec_():
                while q.next():
                    print(f">{q.value('id')}_{q.value('frame')}\n{q.value('aligned_target_sequence').replace('-', '')}", file=f_out)



    def save_tsv(self):
        """
        save sequences in TSV format (same as input)
        """

        file_name, _ = QFileDialog.getSaveFileName(self, "Save file TSV", "")
        if not file_name:
            return

        self.statusBar().showMessage(f"Saving TSV file")
        app.processEvents()
        with open(file_name, "w") as f_out:

            if self.sql2:
                conditions = f"WHERE {self.sql2}"
            else:
                conditions = ""


            q = QtSql.QSqlQuery(f"SELECT * FROM sequences {conditions} ORDER BY id, frame")
            if q.exec_():
                while q.next():
                    app.processEvents()
                    print("\t".join([str(q.value(self.model.headerData(i, QtCore.Qt.Horizontal))) for i in range(self.model.columnCount())]), file=f_out)

        self.statusBar().showMessage(f"Saving TSV file done")


    def message_dialog(self, title, text, buttons):
        message = QMessageBox()
        message.setWindowTitle(title)
        message.setText(text)
        message.setIcon(QMessageBox.Question)
        for button in buttons:
            message.addButton(button, QMessageBox.YesRole)

        message.exec_()
        return message.clickedButton().text()


    def save_fbs(self):
        """
        Export sequences in query anchored format
        """

        flag_group = (self.message_dialog("SW Scan", "Group identical sequences?", ["Yes", "No"]) == "Yes")

        # all uniq aligned query
        file_name, _ = QFileDialog.getSaveFileName(self, "Save file FBS", "")
        if not file_name:
            return

        self.statusBar().showMessage(f"Saving FBS file")
        app.processEvents()

        if self.sql2:
            conditions = f"WHERE {self.sql2}"
        else:
            conditions = ""


        q = QtSql.QSqlQuery(f"SELECT MAX(length(id) + length(description)) AS max_id_descr_len FROM sequences {conditions}")
        if not q.exec_():
            self.statusBar().showMessage(f"SQL error")

        q.first()
        max_id_len = q.value('max_id_descr_len')
        app.processEvents()


        q = QtSql.QSqlQuery(f"SELECT count(distinct aligned_query_sequence) as n FROM sequences {conditions}")
        if  q.exec_():
            q.first()
            n_distinct_aligned_query = q.value('n')


        with open(file_name, "w") as f_out:

            q = QtSql.QSqlQuery(f"SELECT distinct aligned_query_sequence FROM sequences {conditions} ORDER BY score DESC, LENGTH(aligned_query_sequence) DESC")
            if not q.exec_():
                self.statusBar().showMessage(f"SQL error")

            count = 0
            out = ""
            count_group = 1
            while q.next():

                self.statusBar().showMessage(f"Saving FBS file {count} / {n_distinct_aligned_query}")
                app.processEvents()

                out += "query" + (" " *  (max_id_len + 5 - 5)) + q.value('aligned_query_sequence') + "\n"
                count += 1
                if conditions:
                    conditions2 = " AND " + conditions.replace("WHERE", "")
                else:
                    conditions2 = ""

                '''
                query2 = (f"SELECT id, description, frame, aligned_query_sequence, aligned_target_sequence FROM sequences "
                                      f"where aligned_query_sequence like '%{q.value('aligned_query_sequence')}%' {conditions2} "
                                      "ORDER BY id")
                '''

                query2 = (f"SELECT id, description, frame, aligned_query_sequence, aligned_target_sequence FROM sequences "
                                      f"where aligned_query_sequence = '{q.value('aligned_query_sequence')}' {conditions2} "
                                      "ORDER BY score DESC, id")
                                    

                q2 = QtSql.QSqlQuery(query2)

                if not q2.exec_():
                    self.statusBar().showMessage(f"SQL error in q2")

                id_dict = {}
                cleaned_seq = {}

                while q2.next():

                    id_descr = f'{q2.value("id")} {q2.value("description")} {q2.value("frame")}'
                    aligned_target_sequence = q2.value("aligned_target_sequence")

                    if aligned_target_sequence not in id_dict:
                        id_dict[aligned_target_sequence] = [id_descr]
                    else:
                        id_dict[aligned_target_sequence].append(id_descr)


                    formated_aligned_target_sequence = (" " * q2.value("aligned_query_sequence").index(q.value("aligned_query_sequence"))) + aligned_target_sequence

                    ats = ""
                    for nq, nt in zip(q.value("aligned_query_sequence"), formated_aligned_target_sequence):
                        if nq == nt:
                            ats += "."
                        elif nt == " ":
                            ats += " "
                        else:
                            ats += nt

                    cleaned_seq[aligned_target_sequence] = ats

                for seq in id_dict:
                    if flag_group:
                        if len(id_dict[seq]) > 1:
                            id_descr = f"group #{count_group} ({len(id_dict[seq]):,} seq)"
                            count_group += 1
                        else:
                            id_descr = id_dict[seq][0]

                        while len(id_descr) < max_id_len + 5:
                            id_descr += " "

                        out += id_descr

                        out += cleaned_seq[seq] + "\n"

                    else:
                        for id in id_dict[seq]:
                            id_descr = id
                            while len(id_descr) < max_id_len + 5:
                                id_descr += " "
                            out += id_descr
                            out += cleaned_seq[seq] + "\n"

                out += "\n\n"

            f_out.write(out)

        self.statusBar().showMessage(f"Saving FBS file done")



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


