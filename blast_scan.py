#!/usr/bin/env python3
"""
Display and filter results of sw (SQLite version)
(c) Olivier Friard 2021

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

from PyQt5.QtWidgets import QMainWindow, QApplication, QTableWidgetItem, QFileDialog, QMessageBox
from PyQt5 import QtSql

from PyQt5.QtCore import QObject, QThread, pyqtSignal

from PyQt5 import QtCore

from blast_scan_ui import Ui_MainWindow
import sys
import pandas as pd

__version__ = "1"
__version_date__ = "2022-03-14"


class SW_Scan(QMainWindow, Ui_MainWindow):
    def __init__(self, input_file_name: str, parent=None):

        super(SW_Scan, self).__init__(parent)

        self.setupUi(self)

        self.setWindowTitle(f"BLAST Scan v.{__version__}")

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
        # about_dialog.setIconPixmap(QPixmap(":/small_logo"))

        about_dialog.setWindowTitle("About BLAST Scan")
        about_dialog.setStandardButtons(QMessageBox.Ok)
        about_dialog.setDefaultButton(QMessageBox.Ok)
        about_dialog.setEscapeButton(QMessageBox.Ok)

        about_dialog.setInformativeText(
            (
                f"<b>BLAST Scan</b> v. {__version__} - {__version_date__}"
                "<p>Copyright &copy; 2021 Olivier Friard<br>"
                "Department of Life Sciences and Systems Biology<br>"
                "University of Torino - Italy<br>"
            )
        )

        about_dialog.setDetailedText("")

        _ = about_dialog.exec_()

    def load_file(self, file_name):

        if not file_name:
            self.file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "")
            if not self.file_name:
                return
        else:
            self.file_name = file_name

        """
        self.connection = sqlite3.connect(':memory:')
        cur = self.connection.cursor()
        cur.execute(("CREATE TABLE sequences ("
        "accession text, "
        "description text, "
        "frame text, "
        "identity float, "
        "score float, "
        "align_length int, "
        "target_length int, "
        "aligned_query_sequence text, "
        aligned_target_sequence text, query_begin int, query_end int, target_begin int, target_end_optimal int)")

        )
        self.connection.commit()
        """

        df = pd.read_csv(
            self.file_name, sep="\t"
        )  # https://towardsdatascience.com/%EF%B8%8F-load-the-same-csv-file-10x-times-faster-and-with-10x-less-memory-%EF%B8%8F-e93b485086c7
        print(df.columns)

        """
        sseqid stitle length pident evalue bitscore score qframe sframe qstart qend sstart send slen qseq sseq
        """

        """
        df.to_sql('sequences', self.connection, if_exists='replace', index=False)
        cur.execute('SELECT count(*) FROM sequences')
        row = cur.fetchone()
        self.n_seq = row[0]
        """

        self.db = QtSql.QSqlDatabase.addDatabase("QSQLITE")
        self.db.setDatabaseName(":memory:")
        if not self.db.open():
            print("Error opening the file")
            sys.exit(-1)

        createTableQuery = QtSql.QSqlQuery()

        createTableQuery.exec(
            (
                "CREATE TABLE sequences ("
                "id INTEGER PRIMARY KEY AUTOINCREMENT UNIQUE NOT NULL,"
                "accession text, "  # sseqid
                "description text, "  # stitle
                "frame text, "  # sframe
                "identity float, "  # pident
                "score float, "  # score
                "align_length int, "  # length
                "target_length int, "  # slen
                "aligned_query_sequence text, "  # qseq
                "aligned_target_sequence text, "  # sseq
                "query_begin int, "  # qstart
                "query_end int, "  # qend
                "target_begin int, "  # sstart
                "target_end_optimal int)"  # send
            )
        )

        print(self.db.tables())

        query = QtSql.QSqlQuery()

        return

        self.model = QtSql.QSqlTableModel()
        self.tw.setModel(self.model)
        self.model.setTable("sequences")
        self.model.select()

        q = QtSql.QSqlQuery("SELECT count(*) as n FROM sequences")
        if q.exec_():
            q.first()
            self.statusBar().showMessage(f"BLAST results loaded: {q.value('n'):,} sequence(s)")

    def clear(self):
        for w in (self.le_id, self.le_description1, self.le_description2, self.le_identity, self.le_align_length):
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

        self.pte_sql.setPlainText(self.sql2)
        self.le_order.setText(SEQ_ORDER)

        self.model.setFilter(self.sql)
        self.model.select()

        q = QtSql.QSqlQuery(f"SELECT count(*) as n FROM sequences WHERE {self.sql2}")
        if q.exec_():
            q.first()
            self.statusBar().showMessage(f"{q.value('n'):,} sequence(s) filtered")

    def run_query(self):
        """
        Run SQL query defined by user
        """
        if self.pte_sql.toPlainText():

            self.model.setFilter(self.pte_sql.toPlainText())
            self.model.select()

            q = QtSql.QSqlQuery(f"SELECT count(*) as n FROM sequences WHERE {self.pte_sql.toPlainText()}")
            if q.exec_():
                q.first()
                self.statusBar().showMessage(f"{q.value('n'):,} sequence(s) filtered")

            self.sql2 = self.pte_sql.toPlainText()

        else:

            self.sql2 = ""

            self.statusBar().showMessage(f"No query to run")

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
            q = QtSql.QSqlQuery(
                f"SELECT id, frame, aligned_target_sequence FROM sequences {conditions}  {self.le_order.text()}"
            )
            if q.exec_():
                while q.next():
                    print(
                        f">{q.value('id')}_{q.value('frame')}\n{q.value('aligned_target_sequence').replace('-', '')}",
                        file=f_out,
                    )

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

            q = QtSql.QSqlQuery(f"SELECT * FROM sequences {conditions} {self.le_order.text()}")
            if q.exec_():
                while q.next():
                    app.processEvents()
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

        flag_group = self.message_dialog("SW Scan", "Group identical sequences?", ["Yes", "No"]) == "Yes"

        # all uniq aligned query
        file_name, _ = QFileDialog.getSaveFileName(self, "Save file FBS", "")
        if not file_name:
            return

        self.statusBar().showMessage(f"Saving FBS file")
        app.processEvents()

        print(f"{self.sql2=}")

        if self.sql2:
            conditions = f"WHERE {self.sql2}"
        else:
            conditions = ""

        # get max of id and description
        q = QtSql.QSqlQuery(
            f"SELECT MAX(length(id) + length(description)) AS max_id_descr_len FROM sequences {conditions}"
        )
        if not q.exec_():
            self.statusBar().showMessage(f"SQL error")

        q.first()
        max_id_len = q.value("max_id_descr_len")
        app.processEvents()

        q = QtSql.QSqlQuery(f"SELECT count(distinct aligned_query_sequence) as n FROM sequences {conditions}")
        if q.exec_():
            q.first()
            n_distinct_aligned_query = q.value("n")

        with open(file_name, "w") as f_out:

            q = QtSql.QSqlQuery(
                f"SELECT distinct aligned_query_sequence FROM sequences {conditions} ORDER BY score DESC, LENGTH(aligned_query_sequence) DESC"
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

                """
                query2 = (f"SELECT id, description, frame, aligned_query_sequence, aligned_target_sequence FROM sequences "
                                      f"where aligned_query_sequence like '%{q.value('aligned_query_sequence')}%' {conditions2} "
                                      "ORDER BY id")
                """

                query2 = (
                    f"SELECT id, description, frame, aligned_query_sequence, aligned_target_sequence FROM sequences "
                    f"where aligned_query_sequence = '{q.value('aligned_query_sequence')}' {conditions2} "
                    "ORDER BY score DESC, id"
                )

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