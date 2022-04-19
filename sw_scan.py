#!/usr/bin/env python3

"""
Display and filter results of the sw (SQLite version) and the blast.py programs
See https://github.com/olivierfriard/sw_scan

(c) Olivier Friard 2021-2022


v.8
-------
Added check for duplicate sequences in save in FASTA format operation
Added check for error during the 'save file in FASTA format' and 'save file in TSV format' operation

v.7
-------
Added file name on title bar

v.6
-------
Added blast.py output


"""

FIELDS_NAME = (
    "accession",
    "description",
    "frame",
    "identity",
    "score",
    "align_length",
    "target_length",
    "aligned_query_sequence",
    "aligned_target_sequence",
    "query_begin",
    "query_end",
    "target_begin",
    "target_end_optimal",
)

SEQ_ORDER = " ORDER BY identity DESC "

from PyQt5.QtWidgets import QMainWindow, QApplication, QFileDialog, QMessageBox
from PyQt5 import QtSql
from PyQt5 import QtCore

from sw_scan_ui import Ui_MainWindow
import sys
import time
import traceback

__version__ = "8"
__version_date__ = "2022-04-19"


class SW_Scan(QMainWindow, Ui_MainWindow):
    def __init__(self, input_file_name: str, parent=None):

        super(SW_Scan, self).__init__(parent)

        self.setupUi(self)

        sys.excepthook = self.excepthook

        self.setWindowTitle(f"SW Scan v.{__version__}")
        self.lb_copyright.setText("(c) 2021-2022 Olivier Friard")
        self.initialize_var()
        self.connections()
        app.processEvents()

        if input_file_name:
            self.load_file(input_file_name)

    def excepthook(self, exception_type, exception_value, traceback_object):
        """
        error management
        """

        exception_text = "".join(traceback.format_exception(exception_type, exception_value, traceback_object))

        error_text = exception_text.replace("\r\n", "\n").replace("\n", "<br>")
        text = f"SW Scan version: {__version__}<br><br>" f"<b>An error has occured</b>:<br>" f"{error_text}<br><br>"

        errorbox = QMessageBox()
        errorbox.setWindowTitle("SW Scan error occured")
        errorbox.setText(text)
        errorbox.setTextFormat(Qt.RichText)
        errorbox.setStandardButtons(QMessageBox.Abort)

        continueButton = errorbox.addButton("Ignore and try to continue", QMessageBox.RejectRole)

        ret = errorbox.exec_()

        if ret == QMessageBox.Abort:
            sys.exit(1)

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
        # about_dialog.setIconPixmap(QPixmap(":/small_logo"))

        about_dialog.setWindowTitle("About SW Scan")
        about_dialog.setStandardButtons(QMessageBox.Ok)
        about_dialog.setDefaultButton(QMessageBox.Ok)
        about_dialog.setEscapeButton(QMessageBox.Ok)
        about_dialog.setInformativeText(
            (
                f"<b>SW Scan</b> v. {__version__} - {__version_date__}"
                "<p>&copy; 2021-2022 Olivier Friard<br>"
                "Department of Life Sciences and Systems Biology<br>"
                "University of Torino - Italy<br>"
            )
        )

        about_dialog.setDetailedText("")

        _ = about_dialog.exec_()

    def load_file(self, file_name):
        """
        Load alignment file (sw or blast.py)
        """

        if not file_name:
            self.file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "")
            if not self.file_name:
                return
        else:
            self.file_name = file_name

        t1 = time.time()

        self.db = QtSql.QSqlDatabase.addDatabase("QSQLITE")
        self.db.setDatabaseName(self.file_name)
        if not self.db.open():
            self.statusBar().showMessage(f"Error opening the alignments file")
            sys.exit(-1)

        self.model = QtSql.QSqlTableModel()
        self.tw.setModel(self.model)
        self.model.setTable("sequences")
        self.model.select()

        t2 = time.time()
        print("file loading time: ", t2 - t1)

        q = QtSql.QSqlQuery("SELECT count(*) as n FROM sequences")
        if q.exec_():
            q.first()
            self.statusBar().showMessage(
                f"Alignments file loaded: {q.value('n'):,} sequence{'s' if q.value('n') > 1 else ''}"
            )
            self.setWindowTitle(f"{self.file_name} - SW Scan v.{__version__}")
            self.lb_results_file.setText(
                f"Alignments file: <b>{self.file_name}</b> ({q.value('n'):,} sequence{'s' if q.value('n') > 1 else ''})"
            )
        else:
            self.statusBar().showMessage(f"Error opening the alignments file")

    def clear(self):
        for w in (self.le_id, self.le_description1, self.le_description2, self.le_identity, self.le_align_length):
            w.clear()

        self.filter()

    def filter(self):
        """
        filter sequences
        """

        self.statusBar().showMessage("")
        if self.file_name == "":
            return

        self.statusBar().showMessage(f"Filtering sequences")

        id = self.le_id.text().upper() if self.le_id.text() else ""
        description1 = self.le_description1.text().upper() if self.le_description1.text() else ""
        description2 = self.le_description2.text().upper() if self.le_description2.text() else ""
        identity_pc = float(self.le_identity.text()) if self.le_identity.text() else 0
        min_align_length = int(self.le_align_length.text()) if self.le_align_length.text() else 0

        self.sql2 = ""
        if id:
            self.sql2 += f" `id` LIKE '%{id}%' "

        if description1:
            for term in description1.split(","):
                if self.sql2:
                    self.sql2 += " AND "
                self.sql2 += f" description LIKE '%{term}%' "

        if description2:
            for term in description2.split(","):
                if self.sql2:
                    self.sql2 += " AND "
                self.sql2 += f" description NOT LIKE '%{term}%' "

        if identity_pc:
            if self.sql2:
                self.sql2 += " AND "
            self.sql2 += f" identity {self.cb_identity_relation.currentText()} {identity_pc} "

        if min_align_length:
            if self.sql2:
                self.sql2 += " AND "
            self.sql2 += f" align_length >= {min_align_length} "

        if self.sql2:
            self.sql = self.sql2 + SEQ_ORDER
        else:
            self.sql = ""

        self.frame.setEnabled(False)
        self.pte_sql.setPlainText(self.sql2)
        self.le_order.setText(SEQ_ORDER)

        app.processEvents()

        t1 = time.time()
        self.model.setFilter(self.sql)
        self.model.select()
        t2 = time.time()
        print("filter:", t2 - t1)

        q = QtSql.QSqlQuery(f"SELECT count(*) as n FROM sequences WHERE {self.sql2}")
        if q.exec_():
            q.first()
            self.statusBar().showMessage(f"{q.value('n'):,} sequence{'s' if q.value('n') > 1 else ''} filtered")

        self.frame.setEnabled(True)

    def run_query(self):
        """
        Run SQL query defined by user
        """
        if self.pte_sql.toPlainText():

            self.statusBar().showMessage(f"Running query")
            self.frame.setEnabled(False)
            app.processEvents()
            t1 = time.time()
            self.model.setFilter(self.pte_sql.toPlainText())
            self.model.select()
            t2 = time.time()
            print("run query", t2 - t1)

            q = QtSql.QSqlQuery(f"SELECT count(*) as n FROM sequences WHERE {self.pte_sql.toPlainText()}")
            if q.exec_():
                q.first()
                self.statusBar().showMessage(f"{q.value('n'):,} sequence{'s' if q.value('n') > 1 else ''} filtered")

            self.sql2 = self.pte_sql.toPlainText()
            self.frame.setEnabled(True)

        else:
            self.sql2 = ""
            self.statusBar().showMessage(f"No query to run")

    def save_fasta(self):
        """
        save filtered sequences in FASTA format
        modify the sequence id in case of duplicates
        """

        file_name, _ = QFileDialog.getSaveFileName(self, "Save sequences in FASTA format", "")
        if not file_name:
            return
        try:

            # save sequences
            with open(file_name, "w") as f_out:

                if self.sql2:
                    conditions = f"WHERE {self.sql2}"
                else:
                    conditions = ""

                # check for duplicates sequence id
                duplicate_seq_id = []
                q = QtSql.QSqlQuery(f"SELECT id, frame FROM sequences {conditions} GROUP BY id HAVING COUNT(id) > 1")
                if q.exec_():
                    while q.next():
                        duplicate_seq_id.append(f"{q.value('id')}_{q.value('frame')}")

                print(f"saving sequences in FASTA format {conditions}")
                self.statusBar().showMessage(f"Saving sequences in FASTA format. Please wait...")
                self.frame.setEnabled(False)
                t1 = time.time()
                while time.time() - t1 < 1:
                    app.processEvents()

                q = QtSql.QSqlQuery(
                    f"SELECT id, frame, aligned_target_sequence FROM sequences {conditions} {self.le_order.text()}"
                )
                if q.exec_():
                    while q.next():
                        if f"{q.value('id')}_{q.value('frame')}" not in duplicate_seq_id:
                            print(
                                f">{q.value('id')}_{q.value('frame')}\n{q.value('aligned_target_sequence').replace('-', '')}",
                                file=f_out,
                            )
                        else:
                            # duplicate sequence
                            id = 1
                            while f"{q.value('id')}_{q.value('frame')}#{id}" in duplicate_seq_id:
                                id += 1
                            print(
                                f">{q.value('id')}_{q.value('frame')}#{id}\n{q.value('aligned_target_sequence').replace('-', '')}",
                                file=f_out,
                            )
                            duplicate_seq_id.append(f"{q.value('id')}_{q.value('frame')}#{id}")

                print(f"Saving sequences in FASTA format done")
                self.statusBar().showMessage(f"Sequences saved in FASTA format in {file_name}")

        except:
            QMessageBox.critical(
                self,
                "SW Scan",
                ("An error occured during the file saving operation."),
            )

        self.frame.setEnabled(True)

    def save_tsv(self):
        """
        save sequences in TSV format (same as input)
        """

        file_name, _ = QFileDialog.getSaveFileName(self, "Save file TSV", "")
        if not file_name:
            return
        try:
            with open(file_name, "w") as f_out:

                if self.sql2:
                    conditions = f"WHERE {self.sql2}"
                else:
                    conditions = ""

                print(f"saving sequences in FASTA format {conditions}")
                self.statusBar().showMessage(f"Saving TSV file. Please wait")
                self.frame.setEnabled(False)
                t1 = time.time()
                while time.time() - t1 < 1:
                    app.processEvents()

                q = QtSql.QSqlQuery(f"SELECT * FROM sequences {conditions} {self.le_order.text()}")
                if q.exec_():
                    while q.next():
                        print(
                            "\t".join(
                                [
                                    str(q.value(self.model.headerData(i, QtCore.Qt.Horizontal)))
                                    for i in range(self.model.columnCount())
                                ]
                            ),
                            file=f_out,
                        )

            self.statusBar().showMessage(f"Saving TSV file done")
        except Exception:
            QMessageBox.critical(
                self,
                "SW Scan",
                ("An error occured during the file saving operation."),
            )
        self.frame.setEnabled(True)

    def message_dialog(self, title, text, buttons):
        """
        display message box and returns the clicked button
        """
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

        flag_group = self.message_dialog("SW Scan", "Group identical sequences?", ["Yes", "No"]) == "Yes"

        # all uniq aligned query
        file_name, _ = QFileDialog.getSaveFileName(self, "Save file FBS", "")
        if not file_name:
            return

        print("Saving sequences in FBS format")
        self.statusBar().showMessage(f"Saving FBS file")
        self.frame.setEnabled(False)

        t1 = time.time()
        while time.time() - t1 < 1:
            app.processEvents()

        t1 = time.time()
        if self.sql2:
            conditions = f"WHERE {self.sql2}"
        else:
            conditions = ""

        # check for duplicates sequence id
        duplicate_seq_id = []
        q = QtSql.QSqlQuery(f"SELECT id, frame FROM sequences {conditions} GROUP BY id HAVING COUNT(id) > 1")
        if q.exec_():
            while q.next():
                duplicate_seq_id.append(f"{q.value('id')}")

        print(f"{duplicate_seq_id=}")

        # get max length of id and description
        q = QtSql.QSqlQuery(
            f"SELECT MAX(length(id) + length(description)) AS max_id_descr_len FROM sequences {conditions}"
        )
        if not q.exec_():
            self.statusBar().showMessage(f"SQL error")

        q.first()
        max_id_len = q.value("max_id_descr_len")

        q = QtSql.QSqlQuery(f"SELECT count(distinct aligned_query_sequence) as n FROM sequences {conditions}")
        if q.exec_():
            q.first()
            n_distinct_aligned_query = q.value("n")

        with open(file_name, "w") as f_out:

            q = QtSql.QSqlQuery(
                f"SELECT DISTINCT aligned_query_sequence FROM sequences {conditions} ORDER BY score DESC, LENGTH(aligned_query_sequence) DESC"
            )
            if not q.exec_():
                self.statusBar().showMessage(f"SQL error")

            count = 0
            out = ""
            count_group = 1
            while q.next():

                self.statusBar().showMessage(f"Saving FBS file {count} / {n_distinct_aligned_query}")
                app.processEvents()

                out += "query" + (" " * (max_id_len + 5 - 5)) + q.value("aligned_query_sequence") + "\n"
                count += 1
                if conditions:
                    conditions2 = " AND " + conditions.replace("WHERE", "")
                else:
                    conditions2 = ""

                query2 = (
                    f"SELECT id, description, frame, aligned_query_sequence, aligned_target_sequence FROM sequences "
                    f"WHERE aligned_query_sequence = '{q.value('aligned_query_sequence')}' {conditions2} "
                    "ORDER BY score DESC, id"
                )

                q2 = QtSql.QSqlQuery(query2)

                if not q2.exec_():
                    self.statusBar().showMessage(f"SQL error in q2")

                id_dict = {}
                cleaned_seq = {}

                while q2.next():

                    # remove id if contained in description (for BLAST output)
                    descr = q2.value("description").replace(q2.value("id") + " ", "")

                    if q2.value("id") in duplicate_seq_id:
                        # duplicate sequence
                        print(q2.value("id"), "duplicate")
                        idx = 1
                        while f"{q2.value('id')}#{idx}" in duplicate_seq_id:
                            idx += 1
                        id_ = f"{q2.value('id')}#{idx}"
                        print("replaced by", id_)
                        duplicate_seq_id.append(id_)
                    else:
                        id_ = q2.value("id")

                    id_descr = f'{id_} {descr} {q2.value("frame")}'

                    aligned_target_sequence = q2.value("aligned_target_sequence")

                    if aligned_target_sequence not in id_dict:
                        id_dict[aligned_target_sequence] = [id_descr]
                    else:
                        id_dict[aligned_target_sequence].append(id_descr)

                    formated_aligned_target_sequence = (
                        " " * q2.value("aligned_query_sequence").index(q.value("aligned_query_sequence"))
                    ) + aligned_target_sequence

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

        print(f"Sequences saved in FBS format in {round(time.time() - t1)} s")
        self.statusBar().showMessage(f"Saving FBS file done")
        self.frame.setEnabled(True)


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
